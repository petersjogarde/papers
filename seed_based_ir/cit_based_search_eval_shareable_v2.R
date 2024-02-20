#------------------------------------------------------------------------------
# The script compares different approaches to retrieve related publications to a set of seeds. 
# 
# Data from the NIH OCC collection is used by this script: https://nih.figshare.com/collections/iCite_Database_Snapshots_NIH_Open_Citation_Collection_/4586573/42
# Data should be located in a Postgres database. Names of used tables are specified below. Column names originate from the downloaded data. 
# 
#------------------------------------------------------------------------------
#Set main working directory
folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(folder)

#Packages
library(dplyr)
library(httr)
library(jsonlite)
library(scales)
library(RPostgreSQL)
library(ggplot2)
library(DBI)
library(getPass)
library(xlsx)

source("getDCfromSeeds.R")
source("getBCfromSeeds.R")
source("getCCfromSeeds.R")
source("getRAfromSeeds.R")
source("getPrecisionRecall.R")

#------------------------------------------------------------------------------
#Variables
srYear <- 2022 #Publication year of the systematic reviews
refFirstYear <- 2010 #First year of references
refLastYear <- 2022 #Last year of references
minNumRefs <- 30 #Minimum number of references for a systematic review to be included
numSystematicReviews <- 3000 #Number of systematic reviews to include
numSeedRefs <- 5 #Number of seeds to be used for each systematic review
maxRank <- 2000 #Max rank to be shown in precision and recall graphs
minCCscore <- 2 #Minimum number of CC-relations to be used (for efficiency)  
minBCscore <- 2 #Minimum number of BC-relations to be used (for efficiency)
runQueries = T #Get new data from PubMed related article (data can be stored and imported instead of rerunning)
outputFolder <- paste0(Sys.Date(), "_output") 
readFolder <- "2023-10-13_output" #If runQueries is FALSE - this is the location to retrieve data from
dir.create(outputFolder, showWarnings = FALSE, recursive = FALSE)

#Names of the database tables including the NIHOCC data from figshare 
nihOccMetadataTable = "rawdata.nihocc_occ_metadata_pmids"
nihOccCitationsTable = "rawdata.nihocc_citations_pmids"

#------------------------------------------------------------------------------

#Connect to db
if(runQueries){
 db <- dbConnect(drv = dbDriver(""),
                 dbname = "",
                 host = '127.0.0.1',
                 port = 55432,
                 user = "",
                 pass = getPass(msg = "PASSWORD: ", noblank = FALSE, forcemask = FALSE))
}

#Get list of systematic reviews, their refs and make seed selection
#==============================
if(runQueries){
 sql <- paste0("
 SELECT 
   a.pmid, 
   a.title, 
   a.journal, 
   a.authors, 
   a.year, 
   a.doi, 
   count(DISTINCT cp.cited) as num_refs
 FROM ", nihOccMetadataTable, " a 
 INNER JOIN ", nihOccCitationsTable, " cp ON (cp.citing = a.pmid)
 INNER JOIN ", nihOccMetadataTable, " ac ON (ac.pmid = cp.cited)
 WHERE a.title ~* '\\msystematic review\\M'
 AND a.year = ", srYear,"
 AND ac.year BETWEEN ", refFirstYear," AND ", refLastYear,"
 GROUP BY 1 
 HAVING count(*) > ", minNumRefs,"
 ORDER BY random()
 LIMIT ", numSystematicReviews)
 
 sr <- dbGetQuery(db, sql)
 
 #Save systematic reviews to file
 write.csv(sr,paste0(outputFolder, "/systematic_reviews.csv"),row.names=F)
 saveRDS(sr,paste0(outputFolder, "/sr.rds"))
 
 print(paste0("Number of systematic reviews: ", nrow(sr)))
 #==============================
 
 #Get refs and select proxy seed articles
 #==============================
 # Create string of systematic review pmids
 srPmids <- paste0(sr$pmid, collapse=",")
 
 sql <- paste0("
 SELECT 
  cp.citing AS pmid, 
  cp.cited,
  ac.citation_count AS c
 FROM ", nihOccCitationsTable, " cp
 	INNER JOIN ", nihOccMetadataTable, " ac ON (ac.pmid = cp.cited)
 WHERE 
  cp.citing IN (", srPmids, ")
  AND ac.year >= ", refFirstYear,"
 ")
 
 srRefs <- dbGetQuery(db, sql)
 print(paste0("Total number of refs: ", nrow(srRefs)))
 
 #Save refs to file
 write.csv(srRefs, paste0(outputFolder, "/systematic_reviews_refs.csv"),row.names=F)
 saveRDS(srRefs,paste0(outputFolder, "/srRefs.rds"))
 
 #Gets sample
 set.seed(38)
 #srRefs$weight <- sqrt(srRefs$rcr) #Probability to be drawn
 srRefs$weight <- 1 #Probability to be drawn
 srSeeds <- srRefs %>% group_by(pmid) %>% sample_n(numSeedRefs, weight=weight)
 srSeeds <- srSeeds[,c("pmid","cited","c")]
 names(srSeeds) <- c("pmid","seed","c")
 srSeeds$is_seed <- T
 
 #Remove seeds from refs
 srRefsNoSeeds <- merge(srRefs, srSeeds, by.x=c("pmid","cited"), by.y=c("pmid","seed"), all.x=T)
 srRefsNoSeeds <- srRefsNoSeeds[is.na(srRefsNoSeeds$is_seed),]
 srRefsNoSeeds <- srRefsNoSeeds[,c("pmid","cited")]
 srRefsNoSeeds$isRef <- 1
 
 #Count refs for each sr
 srRefsCount <- srRefsNoSeeds %>% group_by(pmid) %>% summarise(num_refs = n())
 
 #Save seeds to file
 write.csv(srSeeds, paste0(outputFolder, "/systematic_reviews_seeds.csv"),row.names=F)
 saveRDS(srSeeds,paste0(outputFolder, "/srSeeds.rds"))
 saveRDS(srRefsNoSeeds,paste0(outputFolder, "/srRefsNoSeeds.rds"))
 saveRDS(srRefsCount,paste0(outputFolder, "/srRefsCount.rds"))
} else {
 srSeeds <- readRDS(paste0(readFolder, "/srSeeds.rds"))
 srRefs <- readRDS(paste0(readFolder, "/srRefs.rds"))
 sr <- readRDS(paste0(readFolder, "/sr.rds"))
 srRefsNoSeeds <- readRDS(paste0(readFolder, "/srRefsNoSeeds.rds"))
 srRefsCount <- readRDS(paste0(readFolder, "/srRefsCount.rds"))
} 
#==============================

#==============================
if(runQueries){
 # Results objects 
 
 # DC
 r_dc_pr <- vector("list", length = numSystematicReviews)
 dc_precisions <- double(numSystematicReviews)
 dc_recalls <- double(numSystematicReviews)
 dc_counts <- double(numSystematicReviews)
 dc_rnk_dfs <- vector("list", length = numSystematicReviews)
 
 # BC 
 r_bc_pr <- vector("list", length = numSystematicReviews)
 bc_precisions <- double(numSystematicReviews)
 bc_recalls <- double(numSystematicReviews)
 bc_counts <- double(numSystematicReviews)
 bc_rnk_dfs <- vector("list", length = numSystematicReviews)
 
 # CC
 r_cc_pr <- vector("list", length = numSystematicReviews)
 cc_precisions <- double(numSystematicReviews)
 cc_recalls <- double(numSystematicReviews)
 cc_counts <- double(numSystematicReviews)
 cc_rnk_dfs <- vector("list", length = numSystematicReviews)
 
 # RA
 r_ra_pr <- vector("list", length = numSystematicReviews)
 ra_precisions <- double(numSystematicReviews)
 ra_recalls <- double(numSystematicReviews)
 ra_counts <- double(numSystematicReviews)
 ra_rnk_dfs <- vector("list", length = numSystematicReviews)
 
 # DC_BC_CC
 r_dc_bc_cc_pr <- vector("list", length = numSystematicReviews)
 dc_bc_cc_precisions <- double(numSystematicReviews)
 dc_bc_cc_recalls <- double(numSystematicReviews)
 dc_bc_cc_counts <- double(numSystematicReviews)
 dc_bc_cc_rnk_dfs <- vector("list", length = numSystematicReviews)
 
 # DC_BC_CC_RA
 r_dc_bc_cc_ra_pr <- vector("list", length = numSystematicReviews)
 dc_bc_cc_ra_precisions <- double(numSystematicReviews)
 dc_bc_cc_ra_recalls <- double(numSystematicReviews)
 dc_bc_cc_ra_counts <- double(numSystematicReviews)
 dc_bc_cc_ra_rnk_dfs <- vector("list", length = numSystematicReviews)
 
 for(i in 1:length(unique(srSeeds$pmid))){
  pmid <- unique(srSeeds$pmid)[i]
  if(i %% 1 == 0){
   print(paste0("Row: ", i))
   print(paste0("PMID: ", pmid))
  }
  
  #Put seed pmids into string
  seeds <- srSeeds$seed[srSeeds$pmid==pmid]
  seedPmids_str <- paste0(seeds, collapse=",")
  pmidsToFind <- srRefsNoSeeds$cited[srRefsNoSeeds$pmid==pmid]
  
  #Direct citations 
  #==============================
  #Get related pmids by DC 
  dc_score <- getDCfromSeeds(seedPmids_str, c(pmid, seeds)) #Do not include the systematic review nor the seeds
  dc_score <- dc_score[dc_score$related_pmid != pmid, ]
  
  #Get precision and recall, also by rank
  r_dc_pr[[i]] <- getPrecisionRecall(dc_score, pmidsToFind, maxRank)
  dc_precisions[i] <- r_dc_pr[[i]]$total_precision
  dc_recalls[i] <- r_dc_pr[[i]]$total_recall
  dc_counts[i] <- r_dc_pr[[i]]$total_count
  r_dc_pr[[i]]$df$sr_pmid <- pmid
  dc_rnk_dfs[[i]] <- r_dc_pr[[i]]$df
  #==============================
  
  # Bibliographic coupling
  #==============================
  #Get related pmids by BC
  bc_score <- getBCfromSeeds(seedPmids_str, minBCscore) #Do not include the systematic review nor the seeds
  bc_score <- bc_score[bc_score$related_pmid != pmid, ]
  
  #Get precision and recall, also by rank
  r_bc_pr[[i]] <- getPrecisionRecall(bc_score, pmidsToFind, maxRank)
  bc_precisions[i] <- r_bc_pr[[i]]$total_precision
  bc_recalls[i] <- r_bc_pr[[i]]$total_recall
  bc_counts[i] <- r_bc_pr[[i]]$total_count
  r_bc_pr[[i]]$df$sr_pmid <- pmid
  bc_rnk_dfs[[i]] <- r_bc_pr[[i]]$df
  #==============================
  
  # Co-citations 
  #============================== 
  #Get related pmids by DC 
  cc_score <- getCCfromSeeds(seedPmids_str, removePmids=pmid, minCCscore=minCCscore) #Do not include the systematic review nor the seeds
  cc_score <- cc_score[cc_score$related_pmid != pmid, ]
  
  #Get precision and recall, also by rank
  r_cc_pr[[i]] <- getPrecisionRecall(cc_score, pmidsToFind, maxRank)
  
  cc_precisions[i] <- r_cc_pr[[i]]$total_precision
  cc_recalls[i] <- r_cc_pr[[i]]$total_recall
  cc_counts[i] <- r_cc_pr[[i]]$total_count
  r_cc_pr[[i]]$df$sr_pmid <- pmid
  cc_rnk_dfs[[i]] <- r_cc_pr[[i]]$df
  #============================== 
  
  # PubMed related articles
  #==============================
  #Get related pmids by RA
  ra_score <- getRAfromSeeds(seeds, sysSleepTime=0) #Do not include the systematic review nor the seeds
  ra_score <- ra_score[ra_score$related_pmid != pmid, ]
  
  #Get precision and recall, also by rank
  r_ra_pr[[i]] <- getPrecisionRecall(ra_score, pmidsToFind, maxRank)
  
  ra_precisions[i] <- r_ra_pr[[i]]$total_precision
  ra_recalls[i] <- r_ra_pr[[i]]$total_recall
  ra_counts[i] <- r_ra_pr[[i]]$total_count
  r_ra_pr[[i]]$df$sr_pmid <- pmid
  ra_rnk_dfs[[i]] <- r_ra_pr[[i]]$df
  #==============================
  
  # DC_BC_CC related articles (delimit to maxRank)
  #==============================
  #Normalize scores before combining
  dc_score$score_norm <- dc_score$score / 1
  bc_score$score_norm <- bc_score$score / 10
  cc_score$score_norm <- cc_score$score / 10
  
  #Sum dc bc and cc 
  dc_bc_cc_prep <- bind_rows(dc_score,bc_score,cc_score)
  dc_bc_cc_score <- dc_bc_cc_prep %>% group_by(related_pmid) %>% summarise(score=sum(score_norm))
  
  #Get precision and recall, also by rank
  r_dc_bc_cc_pr[[i]] <- getPrecisionRecall(dc_bc_cc_score, pmidsToFind, maxRank)
  
  dc_bc_cc_precisions[i] <- r_dc_bc_cc_pr[[i]]$total_precision
  dc_bc_cc_recalls[i] <- r_dc_bc_cc_pr[[i]]$total_recall
  dc_bc_cc_counts[i] <- r_dc_bc_cc_pr[[i]]$total_count
  r_dc_bc_cc_pr[[i]]$df$sr_pmid <- pmid
  dc_bc_cc_rnk_dfs[[i]] <- r_dc_bc_cc_pr[[i]]$df
  #==============================
  
  # DC_BC_CC_RA related articles (delimit to rank 100)
  #==============================
  #Normalize scores before 
  ra_score$score_norm <- rescale(ra_score$score, c(min(dc_bc_cc_score$score),max(dc_bc_cc_score$score)))
  
  #Sum dc bc and cc 
  dc_bc_cc_ra_prep <- bind_rows(dc_score,bc_score,cc_score,ra_score)
  dc_bc_cc_ra_score <- dc_bc_cc_ra_prep %>% group_by(related_pmid) %>% summarise(score=sum(score_norm))
  
  #Get precision and recall, also by rank
  r_dc_bc_cc_ra_pr[[i]] <- getPrecisionRecall(dc_bc_cc_ra_score, pmidsToFind, maxRank)
  
  dc_bc_cc_ra_precisions[i] <- r_dc_bc_cc_ra_pr[[i]]$total_precision
  dc_bc_cc_ra_recalls[i] <- r_dc_bc_cc_ra_pr[[i]]$total_recall
  dc_bc_cc_ra_counts[i] <- r_dc_bc_cc_ra_pr[[i]]$total_count
  r_dc_bc_cc_ra_pr[[i]]$df$sr_pmid <- pmid
  dc_bc_cc_ra_rnk_dfs[[i]] <- r_dc_bc_cc_ra_pr[[i]]$df
  #==============================
 }
 #Merge resulting data frames into one
 dc_rnk_df <- do.call("rbind", dc_rnk_dfs)
 bc_rnk_df <- do.call("rbind", bc_rnk_dfs)
 cc_rnk_df <- do.call("rbind", cc_rnk_dfs)
 ra_rnk_df <- do.call("rbind", ra_rnk_dfs)
 dc_bc_cc_rnk_df <- do.call("rbind", dc_bc_cc_rnk_dfs)
 dc_bc_cc_ra_rnk_df <- do.call("rbind", dc_bc_cc_ra_rnk_dfs)
 
 saveRDS(dc_rnk_df, paste0(outputFolder, "/dc_rnk_df.rds"))
 saveRDS(bc_rnk_df, paste0(outputFolder, "/bc_rnk_df.rds"))
 saveRDS(cc_rnk_df, paste0(outputFolder, "/cc_rnk_df.rds"))
 saveRDS(ra_rnk_df, paste0(outputFolder, "/ra_rnk_df.rds"))
 saveRDS(dc_bc_cc_rnk_df, paste0(outputFolder, "/dc_bc_cc_rnk_df.rds"))
 saveRDS(dc_bc_cc_ra_rnk_df, paste0(outputFolder, "/dc_bc_cc_ra_rnk_df.rds"))
 
 saveRDS(dc_precisions, paste0(outputFolder, "/dc_precisions.rds"))
 saveRDS(bc_precisions, paste0(outputFolder, "/bc_precisions.rds"))
 saveRDS(cc_precisions, paste0(outputFolder, "/cc_precisions.rds"))
 saveRDS(ra_precisions, paste0(outputFolder, "/ra_precisions.rds"))
 saveRDS(dc_bc_cc_precisions, paste0(outputFolder, "/dc_bc_cc_precisions.rds"))
 saveRDS(dc_bc_cc_ra_precisions, paste0(outputFolder, "/dc_bc_cc_ra_precisions.rds"))
 
 saveRDS(dc_recalls, paste0(outputFolder, "/dc_recalls.rds"))
 saveRDS(bc_recalls, paste0(outputFolder, "/bc_recalls.rds"))
 saveRDS(cc_recalls, paste0(outputFolder, "/cc_recalls.rds"))
 saveRDS(ra_recalls, paste0(outputFolder, "/ra_recalls.rds"))
 saveRDS(dc_bc_cc_recalls, paste0(outputFolder, "/dc_bc_cc_recalls.rds"))
 saveRDS(dc_bc_cc_ra_recalls, paste0(outputFolder, "/dc_bc_cc_ra_recalls.rds")) 
 
 saveRDS(dc_counts, paste0(outputFolder, "/dc_counts.rds"))
 saveRDS(bc_counts, paste0(outputFolder, "/bc_counts.rds"))
 saveRDS(cc_counts, paste0(outputFolder, "/cc_counts.rds"))
 saveRDS(ra_counts, paste0(outputFolder, "/ra_counts.rds"))
 saveRDS(dc_bc_cc_counts, paste0(outputFolder, "/dc_bc_cc_counts.rds"))
 saveRDS(dc_bc_cc_ra_counts, paste0(outputFolder, "/dc_bc_cc_ra_counts.rds")) 
} else {
 dc_rnk_df <- readRDS(paste0(readFolder, "/dc_rnk_df.rds"))
 bc_rnk_df <- readRDS(paste0(readFolder, "/bc_rnk_df.rds"))
 cc_rnk_df <- readRDS(paste0(readFolder, "/cc_rnk_df.rds"))
 ra_rnk_df <- readRDS(paste0(readFolder, "/ra_rnk_df.rds"))
 dc_bc_cc_rnk_df <- readRDS(paste0(readFolder, "/dc_bc_cc_rnk_df.rds"))
 dc_bc_cc_ra_rnk_df <- readRDS(paste0(readFolder, "/dc_bc_cc_ra_rnk_df.rds"))
 
 dc_precisions <- readRDS(paste0(readFolder, "/dc_precisions.rds"))
 bc_precisions <- readRDS(paste0(readFolder, "/bc_precisions.rds"))
 cc_precisions <- readRDS(paste0(readFolder, "/cc_precisions.rds"))
 ra_precisions <- readRDS(paste0(readFolder, "/ra_precisions.rds"))
 dc_bc_cc_precisions <- readRDS(paste0(readFolder, "/dc_bc_cc_precisions.rds"))
 dc_bc_cc_ra_precisions <- readRDS(paste0(readFolder, "/dc_bc_cc_ra_precisions.rds"))
 
 dc_recalls <- readRDS(paste0(readFolder, "/dc_recalls.rds"))
 bc_recalls <- readRDS(paste0(readFolder, "/bc_recalls.rds"))
 cc_recalls <- readRDS(paste0(readFolder, "/cc_recalls.rds"))
 ra_recalls <- readRDS(paste0(readFolder, "/ra_recalls.rds"))
 dc_bc_cc_recalls <- readRDS(paste0(readFolder, "/dc_bc_cc_recalls.rds"))
 dc_bc_cc_ra_recalls <- readRDS(paste0(readFolder, "/dc_bc_cc_ra_recalls.rds")) 
 
 dc_counts <- readRDS(paste0(readFolder, "/dc_counts.rds"))
 bc_counts <- readRDS(paste0(readFolder, "/bc_counts.rds"))
 cc_counts <- readRDS(paste0(readFolder, "/cc_counts.rds"))
 ra_counts <- readRDS(paste0(readFolder, "/ra_counts.rds"))
 dc_bc_cc_counts <- readRDS(paste0(readFolder, "/dc_bc_cc_counts.rds"))
 dc_bc_cc_ra_counts <- readRDS(paste0(readFolder, "/dc_bc_cc_ra_counts.rds")) 
}

#=========================================
#Plots 
#=========================================
colors <- c("#e6194B","#f58231","#3cb44b","#4363d8","#911eb4","#f032e6")
colors2 <- sort(rep(colors,2))
linetypes <- 1:6
linetypes2 <- rep(c(1,2), 6)
lineWidth <- 1
Approaches <- c("DC", "BC", "CC", "DC_BC_CC", "RA", "DC_BC_CC_RA")
#=========================================

# Number of retrieved publications - violine plot 
#=========================================
counts_df <- rbind(
 data.frame(Approach="DC", count=dc_counts),
 data.frame(Approach="BC", count=bc_counts),
 data.frame(Approach="CC", count=cc_counts),
 data.frame(Approach="RA", count=ra_counts),
 data.frame(Approach="DC_BC_CC", count=dc_bc_cc_counts),
 data.frame(Approach="DC_BC_CC_RA", count=dc_bc_cc_ra_counts)
)

#Control order
counts_df$Approach <- factor(counts_df$Approach, levels = Approaches)

png(paste0(outputFolder, "/violine_plot_count.png"), width=22, height=12, unit="cm", res=600)

#Plot (restrict to 100 relations or less)
counts_df %>% ggplot(aes(x=Approach, y=count, fill=Approach)) +
 geom_violin(width=1.2) +
 geom_boxplot(width=0.1, color="black", alpha=0) +
 scale_fill_manual(values=colors) +
 theme_classic() +
 theme(
  legend.position="none",
  plot.title = element_text(size=16),
  axis.text = element_text(size=14),
  axis.title = element_text(size=14)
 ) +
 ggtitle("") +
 xlab("") + ylab("# retrieved publ.") +
 scale_y_continuous(trans='log10')
# ylim(c(0,10))

dev.off()
#=========================================

# Recall violine plot 
#=========================================

# Bind data 
recall_df <- rbind(
 data.frame(Approach="DC", recall=dc_recalls),
 data.frame(Approach="BC", recall=bc_recalls),
 data.frame(Approach="CC", recall=cc_recalls),
 data.frame(Approach="RA", recall=ra_recalls),
 data.frame(Approach="DC_BC_CC", recall=dc_bc_cc_recalls),
 data.frame(Approach="DC_BC_CC_RA", recall=dc_bc_cc_ra_recalls)
)

#Control order
recall_df$Approach <- factor(recall_df$Approach, levels = Approaches)

png(paste0(outputFolder, "/violine_plot_recall.png"), width=22, height=12, unit="cm", res=600)

#Plot (restrict to 100 relations or less)
recall_df %>% ggplot(aes(x=Approach, y=recall, fill=Approach)) +
 geom_violin(width=1.2) +
 geom_boxplot(width=0.1, color="black", alpha=0) +
 scale_fill_manual(values=colors) +
 theme_classic() +
 theme(
  legend.position="none",
  plot.title = element_text(size=16),
  axis.text = element_text(size=14),
  axis.title = element_text(size=14)
 ) +
 ggtitle("") +
 xlab("") + ylab("Recall")

dev.off()
#=========================================

# Precision violine plot 
#=========================================
precision_df <- rbind(
 data.frame(Approach="DC", precision=dc_precisions),
 data.frame(Approach="BC", precision=bc_precisions),
 data.frame(Approach="CC", precision=cc_precisions),
 data.frame(Approach="RA", precision=ra_precisions),
 data.frame(Approach="DC_BC_CC", precision=dc_bc_cc_precisions),
 data.frame(Approach="DC_BC_CC_RA", precision=dc_bc_cc_ra_precisions)
)

#Control order
precision_df$Approach <- factor(precision_df$Approach, levels = Approaches)

png(paste0(outputFolder, "/violine_plot_precision.png"), width=22, height=12, unit="cm", res=600)

#Plot (restrict to 100 relations or less)
precision_df %>% ggplot(aes(x=Approach, y=precision, fill=Approach)) +
 geom_violin(width=1.2) +
 geom_boxplot(width=0.1, color="black", alpha=0) +
 scale_fill_manual(values=colors) +
 theme_classic() +
 theme(
  legend.position="none",
  plot.title = element_text(size=16),
  axis.text = element_text(size=14),
  axis.title = element_text(size=14)
 ) +
 ggtitle("") +
 xlab("") + ylab("Precision") +
 ylim(c(0,10))

dev.off()
#=========================================

# Total Recall Density plot
#=========================================
# Do not include RA and DC_BC_CC_RA because RA is restricted to top 100
recall_df <- rbind(
 data.frame(Approach="DC", recall=dc_recalls),
 data.frame(Approach="BC", recall=bc_recalls),
 data.frame(Approach="CC", recall=cc_recalls),
 data.frame(Approach="RA", recall=ra_recalls),
 data.frame(Approach="DC_BC_CC", recall=dc_bc_cc_recalls),
 data.frame(Approach="DC_BC_CC_RA", recall=dc_bc_cc_ra_recalls)
)

levels(recall_df$Approach) <- Approaches

#Control order
recall_df$Approach <- factor(recall_df$Approach, levels = Approaches)

#Density curve
recall_density_plot <- ggplot(recall_df, aes(x=recall, colour=Approach)) +
 geom_density(aes(linetype=Approach), size=lineWidth, show.legend=TRUE, alpha = 0.5) +
 theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
 theme(plot.margin = margin(0.3, 0.1, 0.2, 0.1, "cm")) +
 #stat_density(aes(x=recall, colour=Approach), size=lineWidth, geom="line",position="identity") +
 scale_linetype_manual(values = linetypes)

#Print
fileName <- paste0("density_recall.png")
png(paste0(outputFolder, "/", fileName), width=15, height=15, unit="cm", res=600)
print(recall_density_plot)
dev.off()
#=========================================

# Total Precision Density plot
#=========================================
precision_df <- rbind(
 data.frame(Approach="DC", precision=dc_precisions),
 data.frame(Approach="BC", precision=bc_precisions),
 data.frame(Approach="CC", precision=cc_precisions),
 data.frame(Approach="RA", precision=ra_precisions),
 data.frame(Approach="DC_BC_CC", precision=dc_bc_cc_precisions),
 data.frame(Approach="DC_BC_CC_RA", precision=dc_bc_cc_ra_precisions)
)

#Control order
precision_df$Approach <- factor(precision_df$Approach, levels = Approaches)

# Density curve
precision_density_plot <- ggplot(precision_df, aes(x=precision, colour=Approach)) +
 geom_density(aes(linetype=Approach), size=lineWidth, show.legend=TRUE) +
 scale_color_manual(values = c(colors)) +
 scale_linetype_manual(values = linetypes) +
 theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
 theme(plot.margin = margin(0.3, 0.1, 0.2, 0.1, "cm"))
#stat_density(aes(x=precision, colour=Approach), geom="line",position="identity")

#Print
fileName <- paste0("density_precision.png")
png(paste0(outputFolder, "/", fileName), width=15, height=15, unit="cm", res=600)
print(precision_density_plot)
dev.off()

# Density curve LIMITED X-axis
precision_density_plot_xlim <- ggplot(precision_df, aes(x=precision, colour=Approach)) +
 geom_density(aes(linetype=Approach), size=lineWidth, show.legend=TRUE) +
 scale_color_manual(values = c(colors)) +
 scale_linetype_manual(values = linetypes) +
 theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
 theme(plot.margin = margin(0.3, 0.1, 0.2, 0.1, "cm")) +
 xlim(c(0,5))
#stat_density(aes(x=precision, colour=Approach), geom="line",position="identity")

#Print
fileName <- paste0("density_precision_xlim0p1.png")
png(paste0(outputFolder, "/", fileName), width=15, height=15, unit="cm", res=600)
print(precision_density_plot_xlim)
dev.off()

# Recall/precision per rank - up to rank
#======================================
plotRank <- function(plotMaxRank, includeApproaches, measure, measure_name){
 
 dc_plot <- dc_rnk_df[dc_rnk_df$rnk <= plotMaxRank,]
 bc_plot <- bc_rnk_df[bc_rnk_df$rnk <= plotMaxRank,]
 cc_plot <- cc_rnk_df[cc_rnk_df$rnk <= plotMaxRank,]
 dc_bc_cc_plot <- dc_bc_cc_rnk_df[dc_bc_cc_rnk_df$rnk <= plotMaxRank,]
 ra_plot <- ra_rnk_df[ra_rnk_df$rnk <= plotMaxRank,]
 dc_bc_cc_ra_plot <- dc_bc_cc_ra_rnk_df[dc_bc_cc_ra_rnk_df$rnk <= plotMaxRank,]
 #dc_bc_cc_cluster_plot <- dc_bc_cc_cluster_rnk_df[dc_bc_cc_cluster_rnk_df$rnk <= plotMaxRank,]
 
 #If we do not have a value on a row we regard it as no hit because no related publication has been retrieved by the approach.
 pmids <- unique(srSeeds$pmid)
 ranks <- 1:plotMaxRank
 
 df_all_ranks <- expand.grid(pmids,ranks)
 names(df_all_ranks) <- c("sr_pmid","rnk")

 fillRanks <- function(df){
  
  df <- merge(df_all_ranks, df, by=c("sr_pmid","rnk"), all.x=T)
  df <- df[order(df$sr_pmid, df$rnk), ]
  
  if(measure == "rnk_recall"){
   #Fill with fill from tidyr. default decending. 
   df <- df %>% tidyr::fill(measure)
  }
  if(measure == "rnk_precision"){
   df <- df %>% tidyr::fill(hit_sum)
   df$rnk_precision <- df$hit_sum / df$rnk * 100
  }
  
  return(df)
 }
 dc_plot <- fillRanks(dc_plot)
 bc_plot <- fillRanks(bc_plot)
 cc_plot <- fillRanks(cc_plot)
 dc_bc_cc_plot <- fillRanks(dc_bc_cc_plot)
 ra_plot <- fillRanks(ra_plot)
 dc_bc_cc_ra_plot <- fillRanks(dc_bc_cc_ra_plot)
 #dc_bc_cc_cluster_plot <- fillRanks(dc_bc_cc_cluster_plot)
 
 #Add names and colors 
 dc_plot$Approach <- "DC"
 dc_plot$color <- colors[1]
 bc_plot$Approach <- "BC"
 bc_plot$color <- colors[2]
 cc_plot$Approach <- "CC"
 cc_plot$color <- colors[3]
 dc_bc_cc_plot$Approach <- "DC_BC_CC"
 dc_bc_cc_plot$color <- colors[4]
 #dc_bc_cc_cluster_plot$Approach <- "DC_BC_CC_CLUSTER"
 #dc_bc_cc_cluster_plot$color <- colors[5]
 ra_plot$Approach <- "RA"
 ra_plot$color <- colors[5]
 dc_bc_cc_ra_plot$Approach <- "DC_BC_CC_RA"
 dc_bc_cc_ra_plot$color <- colors[6]
 
 df_plot_prep <- rbind(dc_plot, bc_plot, cc_plot, dc_bc_cc_plot, ra_plot, dc_bc_cc_ra_plot)
 
 #Calculate mean per rank
 df_plot <- df_plot_prep %>% group_by(rnk, Approach) %>% summarise(mean=mean(!!sym(measure)), 
                                                                   sd=sd(!!sym(measure)),
                                                                   n=n())
 
 #Calculate error
 df_plot$error <- (1.96 * df_plot$sd/sqrt(df_plot$n)) / 2
 
 #Upper and lower error margins
 df_plot$ymin <- df_plot$mean - df_plot$error
 df_plot$ymax <- df_plot$mean + df_plot$error
 
 #Control order of approaches
 df_plot$Approach <- factor(df_plot$Approach, levels = Approaches)
 
 #Delimit approaches
 df_plot <- df_plot[df_plot$Approach %in% includeApproaches,]
 
 #Plot
 plot <- df_plot %>% ggplot(aes(x=rnk, y=mean, colour=Approach)) +
  geom_line(aes(linetype=Approach), size=lineWidth) + 
  scale_color_manual(values = c(colors)) +
  scale_linetype_manual(values = linetypes) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.margin = margin(0.3, 0.1, 0.2, 0.1, "cm")) +
  labs(title = "", y = paste0(measure_name, " (%)"), x = "Rank")
 plot <- plot + geom_ribbon(aes(ymin=df_plot$ymin, ymax=df_plot$ymax), linetype=2, alpha=0.1)
 
 #Print
 fileName <- paste0(measure, "_plot_", plotMaxRank, ".png")
 png(paste0(outputFolder, "/", fileName), width=15, height=15, unit="cm", res=600)
 print(plot)
 dev.off()
}
#======================================
plotRank(plotMaxRank=100, includeApproaches=Approaches, measure="rnk_recall",measure_name="Recall")
plotRank(plotMaxRank=maxRank, includeApproaches=c("DC","BC","CC","DC_BC_CC"), measure="rnk_recall",measure_name="Recall")
plotRank(plotMaxRank=100, includeApproaches=Approaches, measure="rnk_precision",measure_name="Precision")
plotRank(plotMaxRank=maxRank, includeApproaches=c("DC","BC","CC","DC_BC_CC"), measure="rnk_precision",measure_name="Precision")

dbDisconnect(db)
