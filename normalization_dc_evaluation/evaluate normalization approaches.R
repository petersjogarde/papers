#------------------------------------------------------------------------------
# Evaluate normalization approaches 
# Written by Peter Sjögårde
#------------------------------------------------------------------------------

#Set main working directory
folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(folder)

# Import libraries
#------------------------------------------------------------------------------
library(igraph)
library(ggplot2)
library(dplyr)
library(moments) #For skewness calculation
library(signal) #For interpolation
library(pdfCluster) #For Adjusted Rand Index (ARI)
library(hrbrthemes) #For boxplot
#------------------------------------------------------------------------------

# Variables
#------------------------------------------------------------------------------
getNewResults <- T #Choose if new results are to be calculated or if plots should be created from previous results
dataFolder <- "data"
outputFolder = paste0(Sys.Date(), "_plot_output") #The outputfolder to be created
outputFolder = "2023-12-14_data_output"
readFromFolder = "2023-12-14_data_output" #If getNewResults is false then read previous results from this folder. Otherwise this variable is not used. 

#Color scale used in plots
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

options(scipen=10) #Show values as normal numbers
min_cluster_size=1 #Used in silhouette (set to 1 if no restriction)
resolutions=c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05) #Resolutions to be tested
leiden_iterations = 100 #Number of iteration of the leiden algorithm
ariYear=2019 #Reviews from this year onward are used
ariMinTotalRefs=70 #Delimit to reviews with at least n total relations
ariMinRefsClusterShare=0.4 #Minimum proportion of references within fields for reviews
ariOverlapSameTopic=0.3 #This value determines when two reviews should be considered to address the same topic
#------------------------------------------------------------------------------

#Initialize
#------------------------------------------------------------------------------
#Start time
startTimeStamp <- Sys.time() # get start time

#Create output folder
dir.create(outputFolder, showWarnings = FALSE, recursive = FALSE) #Change showWarnings to TRUE if not to overwrite folder
if(getNewResults){
 dir.create(paste0(outputFolder, "/size distribution histograms"), showWarnings = FALSE, recursive = FALSE) #Change showWarnings to TRUE if not to overwrite folder
}

options(dplyr.summarise.inform = FALSE) #Do not show group by warning that have no effect on results
#------------------------------------------------------------------------------

#Functions
#===================

#Function to load and cluster data and get some of the results
#------------------------------------------------------------------------------
createOutput <- function(pmids_file, relations_file, field_name, consoleFile){
 
 #Print console to file
 con <- file(consoleFile)
 write(' ', file=con, append=F)
 sink(con, append=TRUE)
 sink(con, append=TRUE, type="message")
 
 print("----------------------------------------------------")
 print(paste0("Starting analysis for ", field_name))
 print("----------------------------------------------------")
 
 #Read files
 pubs <- read.table(pmids_file, header=T)
 relations <- read.table(relations_file, header=T)
 
 #Remove citations to the same document (some of these exists in the data)
 relations <- unique(relations[relations$citing != relations$cited,])
 
 #Remove pubs with no relations
 pubs <- pubs[pubs$cnt_relations > 0,]
 pubs <- pubs[pubs$pmid %in% unique(c(relations$citing, relations$cited)),]
 
 #Normalize citation relations 
 #============================
 relations$weigth_geo <- 1 / sqrt(relations$cnt_relations_citing * relations$cnt_relations_cited)
 
 relations$weight_frac <- (1 / relations$cnt_relations_citing + 1 / relations$cnt_relations_cited) / 2
 
 relations$weigth_geo_limit <- 1 / (
  sqrt(ifelse(relations$cnt_relations_citing < 5, 5, relations$cnt_relations_citing)) 
  * sqrt(ifelse(relations$cnt_relations_cited < 5, 5, relations$cnt_relations_cited))
 )
 
 relations$weight_unnormalized <- 1
 
 relations$weight_directed_frac <- (1 / relations$cnt_citing_refs + 1 / relations$cnt_cited_citations) / 2
 
 relations$weight_directed_geo  <- 1 / sqrt(relations$cnt_citing_refs * relations$cnt_cited_citations)
 
 #We now disregard the direction and put the lowest id first
 relations$tempciting <- ifelse(relations$citing < relations$cited, relations$citing, relations$cited)
 relations$tempcited <- ifelse(relations$citing < relations$cited, relations$cited, relations$citing)
 relations$citing <- relations$tempciting
 relations$cited <- relations$tempcited
 
 #If i cites j and j cites i use the average of the citation relation
 relations <- relations %>% group_by(citing, cited) %>% summarise(
  weight_unnormalized = mean(weight_unnormalized),
  weight_frac = mean(weight_frac),
  weigth_geo = mean(weigth_geo), 
  weigth_geo_limit = mean(weigth_geo_limit), 
  weight_directed_frac = mean(weight_directed_frac), 
  weight_directed_geo = mean(weight_directed_geo)
 )
 
 #Normalize to create approximately the same scale
 relations$weight_unnormalized = relations$weight_unnormalized / mean(relations$weight_unnormalized)
 relations$weigth_geo = relations$weigth_geo / mean(relations$weigth_geo)
 relations$weight_frac = relations$weight_frac / mean(relations$weight_frac)
 relations$weigth_geo_limit = relations$weigth_geo_limit / mean(relations$weigth_geo_limit)
 relations$weight_directed_frac = relations$weight_directed_frac / mean(relations$weight_directed_frac)
 relations$weight_directed_geo  = relations$weight_directed_geo  / mean(relations$weight_directed_geo )
 
 #Create graph
 g <- graph_from_data_frame(relations, directed=F, vertices=pubs)
 
 print("===================================================")
 print("Unnormalized approach")
 print("===================================================")
 unnormalized_silhouttes_means <- c()
 unnormalized_internal_citations <- c()
 unnormalized_granularities <- c()
 unnormalized_item_rel_dfs <- list()
 
 for(r in resolutions){
  #Get cluster solution
  r_leiden_unnormalized <- runLeiden(g, w=relations$weight_unnormalized, r=r, F, name="Unnormalized", field_name)
  
  #Add granularity
  unnormalized_granularities <- c(unnormalized_granularities, r_leiden_unnormalized$granularity)
  
  #Get silhoutte values and other metric
  unnormalized_item_rel_df <- getSilhouette(r_leiden_unnormalized$df, relations, min_cluster_size) # At the moment min_cluster_size must be 1, otherwise the calculation of silhoutte will not be correct
  
  #Add df to list
  unnormalized_item_rel_dfs <- append(unnormalized_item_rel_dfs, list(unnormalized_item_rel_df))
  
  #Calculate means and put in lists 
  unnormalized_silhouttes_means <- c(unnormalized_silhouttes_means, mean(unnormalized_item_rel_df$item_rel$silhouette, na.rm=TRUE))
  unnormalized_internal_citations <- c(unnormalized_internal_citations, mean(unnormalized_item_rel_df$item_rel$shr_rels_in_cluster, na.rm=TRUE))
  
  print("---")
 }
 print("===================================================")
 print("Fractional approach")
 print("===================================================")
 frac_silhouttes_means <- c()
 frac_internal_citations <- c()
 frac_granularities <- c()
 frac_item_rel_dfs <- list()
 
 for(r in resolutions){
  #Get cluster solution
  r_leiden_frac <- runLeiden(g, w=relations$weight_frac, r=r, F, name="Fractional", field_name)
  
  #Add granularity
  frac_granularities <- c(frac_granularities, r_leiden_frac$granularity)
  
  #Get silhoutte values and other metric
  frac_item_rel_df <- getSilhouette(r_leiden_frac$df, relations, min_cluster_size) # At the moment min_cluster_size must be 1, otherwise the calculation of silhoutte will not be correct
  
  #Add df to list
  frac_item_rel_dfs <- append(frac_item_rel_dfs, list(frac_item_rel_df))
  
  #Calculate means and put in lists 
  frac_silhouttes_means <- c(frac_silhouttes_means, mean(frac_item_rel_df$item_rel$silhouette, na.rm=TRUE))
  frac_internal_citations <- c(frac_internal_citations, mean(frac_item_rel_df$item_rel$shr_rels_in_cluster, na.rm=TRUE))
  
  print("---")
 }
 print("===================================================")
 print("Geometric approach")
 print("===================================================")
 geometric_silhouttes_means <- c()
 geometric_internal_citations <- c()
 geometric_granularities <- c()
 geometric_item_rel_dfs <- list() 
 
 for(r in resolutions){
  #Get cluster solution
  r_leiden_geometric <- runLeiden(g, w=relations$weigth_geo, r=r, F, name="Geometric", field_name)
  
  #Add granularity
  geometric_granularities <- c(geometric_granularities, r_leiden_geometric$granularity)
  
  #Get silhoutte values and other metric
  geometric_item_rel_df <- getSilhouette(r_leiden_geometric$df, relations, min_cluster_size) # At the moment min_cluster_size must be 1, otherwise the calculation of silhoutte will not be correct
  
  #Add df to list
  geometric_item_rel_dfs <- append(geometric_item_rel_dfs, list(geometric_item_rel_df))
  
  #Calculate means and put in lists 
  geometric_silhouttes_means <- c(geometric_silhouttes_means, mean(geometric_item_rel_df$item_rel$silhouette, na.rm=TRUE))
  geometric_internal_citations <- c(geometric_internal_citations, mean(geometric_item_rel_df$item_rel$shr_rels_in_cluster, na.rm=TRUE))
  
  print("---")
 }
 print("===================================================")
 print("Geometric approach with denominator minimum limit of 5")
 print("===================================================")
 geometric_lim5_silhouttes_means <- c()
 geometric_lim5_internal_citations <- c()
 geometric_lim5_granularities <- c()
 geometric_lim5_item_rel_dfs <- list() 
 
 for(r in resolutions){
  #Get cluster solution
  r_leiden_geometric_lim5 <- runLeiden(g, w=relations$weigth_geo_limit, r=r, F, name="Geometric mean with limit=5", field_name)
  
  #Add granularity
  geometric_lim5_granularities <- c(geometric_lim5_granularities, r_leiden_geometric_lim5$granularity)
  
  #Get silhoutte values and other metric
  geometric_lim5_item_rel_df <- getSilhouette(r_leiden_geometric_lim5$df, relations, min_cluster_size) # At the moment min_cluster_size must be 1, otherwise the calculation of silhoutte will not be correct
  
  #Add df to list
  geometric_lim5_item_rel_dfs <- append(geometric_lim5_item_rel_dfs, list(geometric_lim5_item_rel_df))
  
  #Calculate means and put in lists 
  geometric_lim5_silhouttes_means <- c(geometric_lim5_silhouttes_means, mean(geometric_lim5_item_rel_df$item_rel$silhouette, na.rm=TRUE))
  geometric_lim5_internal_citations <- c(geometric_lim5_internal_citations, mean(geometric_lim5_item_rel_df$item_rel$shr_rels_in_cluster, na.rm=TRUE))
  
  print("---")
 }
 
 print("===================================================")
 print("Directed approach - frac")
 print("===================================================")
 directed_frac_silhouttes_means <- c()
 directed_frac_internal_citations <- c()
 directed_frac_granularities <- c()
 directed_frac_item_rel_dfs <- list() 
 
 for(r in resolutions){
  #Get cluster solution
  r_leiden_directed_frac <- runLeiden(g, w=relations$weight_directed_frac, r=r, F, name="Directional-fractional", field_name)
  
  #Add granularity
  directed_frac_granularities <- c(directed_frac_granularities, r_leiden_directed_frac$granularity)
  
  #Get silhoutte values and other metric
  directed_frac_item_rel_df <- getSilhouette(r_leiden_directed_frac$df, relations, min_cluster_size) # At the moment min_cluster_size must be 1, otherwise the calculation of silhoutte will not be correct
  
  #Add df to list
  directed_frac_item_rel_dfs <- append(directed_frac_item_rel_dfs, list(directed_frac_item_rel_df))
  
  #Calculate means and put in lists 
  directed_frac_silhouttes_means <- c(directed_frac_silhouttes_means, mean(directed_frac_item_rel_df$item_rel$silhouette, na.rm=TRUE))
  directed_frac_internal_citations <- c(directed_frac_internal_citations, mean(directed_frac_item_rel_df$item_rel$shr_rels_in_cluster, na.rm=TRUE))
  
  print("---")
 }
 
 
 print("===================================================")
 print("Directed approach - geo")
 print("===================================================")
 directed_geo_silhouttes_means <- c()
 directed_geo_internal_citations <- c()
 directed_geo_granularities <- c()
 directed_geo_item_rel_dfs <- list() 
 
 for(r in resolutions){
  #Get cluster solution
  r_leiden_directed_geo  <- runLeiden(g, w=relations$weight_directed_geo , r=r, F, name="Directional-geometric", field_name)
  
  #Add granularity
  directed_geo_granularities <- c(directed_geo_granularities, r_leiden_directed_geo $granularity)
  
  #Get silhoutte values and other metric
  directed_geo_item_rel_df <- getSilhouette(r_leiden_directed_geo $df, relations, min_cluster_size) # At the moment min_cluster_size must be 1, otherwise the calculation of silhoutte will not be correct
  
  #Add df to list
  directed_geo_item_rel_dfs <- append(directed_geo_item_rel_dfs, list(directed_geo_item_rel_df))
  
  #Calculate means and put in lists 
  directed_geo_silhouttes_means <- c(directed_geo_silhouttes_means, mean(directed_geo_item_rel_df$item_rel$silhouette, na.rm=TRUE))
  directed_geo_internal_citations <- c(directed_geo_internal_citations, mean(directed_geo_item_rel_df$item_rel$shr_rels_in_cluster, na.rm=TRUE))
  
  print("---")
 }
 
 r=list("unnormalized_item_rel_dfs"=unnormalized_item_rel_dfs,
        "unnormalized_silhouttes_means"=unnormalized_silhouttes_means, 
        "unnormalized_internal_citations"=unnormalized_internal_citations, 
        "unnormalized_granularities"=unnormalized_granularities, 
        "frac_item_rel_dfs"=frac_item_rel_dfs,
        "frac_silhouttes_means"=frac_silhouttes_means, 
        "frac_internal_citations"=frac_internal_citations, 
        "frac_granularities"=frac_granularities, 
        "geometric_item_rel_dfs"=geometric_item_rel_dfs,
        "geometric_silhouttes_means"=geometric_silhouttes_means, 
        "geometric_internal_citations"=geometric_internal_citations, 
        "geometric_granularities"=geometric_granularities, 
        "geometric_lim5_item_rel_dfs"=geometric_lim5_item_rel_dfs,
        "geometric_lim5_silhouttes_means"=geometric_lim5_silhouttes_means, 
        "geometric_lim5_internal_citations"=geometric_lim5_internal_citations, 
        "geometric_lim5_granularities"=geometric_lim5_granularities,
        "directed_frac_item_rel_dfs"=directed_frac_item_rel_dfs,
        "directed_frac_silhouttes_means"=directed_frac_silhouttes_means, 
        "directed_frac_internal_citations"=directed_frac_internal_citations, 
        "directed_frac_granularities"=directed_frac_granularities, 
        "directed_geo_item_rel_dfs"=directed_geo_item_rel_dfs,
        "directed_geo_silhouttes_means"=directed_geo_silhouttes_means, 
        "directed_geo_internal_citations"=directed_geo_internal_citations, 
        "directed_geo_granularities"=directed_geo_granularities, 
        "pubs"=pubs,
        "relations"=relations
 )
 
 # Restore output to console
 sink() 
 sink(type="message")
 
 return=(r)
}
#------------------------------------------------------------------------------

#Function to run the Leiden algorithm and get distribution
#------------------------------------------------------------------------------
runLeiden <- function(g, w, r, printNetworks=F, name, field_name){
 
 leiden <- cluster_leiden(g, objective_function = "CPM", weights=w, resolution_parameter=r, n_iterations=leiden_iterations)
 
 t <- table(leiden$membership)
 
 t <- data.frame(t)
 
 #Create histogram
 h <- hist(t$Freq[t$Freq >= 10], breaks=100, plot=F)
 mainname=paste0(field_name, " - Histogram of cluster sizes\n", name, " (resolution=",r,")")
 
 #Output to file
 outputfilename=tolower(paste0(field_name, "_hist_sizes_", name, "_r",r))
 png(paste0(outputFolder, "/size distribution histograms/",outputfilename,".png"), width=22, height=12, unit="cm", res=600)
 plot(h, xlab="Cluster size", main=mainname)
 dev.off()
 
 nodes_clusters <- data.frame("id"=leiden$names, "cluster"=leiden$membership)
 
 cluster_count <- nodes_clusters %>% group_by(cluster) %>% summarize(cluster_size=n())
 
 print(paste0("Resolution: ", r))
 print(paste0("Number of clusters: ", nrow(cluster_count)))
 print(paste0("Average cluster size: ", mean(cluster_count$cluster_size)))
 print(paste0("Number of clusters with less than 10 publications: ", nrow(cluster_count[cluster_count$cluster_size < 10,])))
 print(paste0("Number of clusters with at least 10 publications: ", nrow(cluster_count[cluster_count$cluster_size >= 10,])))
 print(paste0("Largest cluster size: ", max(cluster_count$cluster_size)))
 granularity = sum(cluster_count$cluster_size)  / sum(cluster_count$cluster_size^2)
 print(paste0("Granularity: ", granularity))
 skewness = round(skewness(cluster_count$cluster_size),2)
 print(paste0("Skewness: ", skewness))
 
 r=list("df"=nodes_clusters, "skewness"=skewness, "granularity"=granularity)
 return(r)
}
#------------------------------------------------------------------------------

#Function to calculate silhouette value of clustering solution
#------------------------------------------------------------------------------
getSilhouette <- function(clusters, relations, min_cluster_size){
 
 clusters <- clusters[!is.na(clusters$cluster),]
 
 r1 <- relations[,c("cited","citing")]
 names(r1) <- c("id1","id2")
 r2 <- relations[,c("citing","cited")]
 names(r2) <- c("id1","id2")
 
 r <- rbind(r1,r2)
 rm(r1)
 rm(r2)
 
 r <- unique(r)
 
 r_cluster <- merge(r, clusters, by.x="id2", by.y="id")
 rm(r)
 
 cluster_count <- clusters %>% group_by(cluster) %>% summarize(Freq=n())
 names(cluster_count) <- c("cluster","cluster_size")
 
 #Count number of relations for each item
 item_relations_count <- r_cluster %>% group_by(id1, cluster) %>% summarize(Freq=n())
 names(item_relations_count) <- c("id","cluster","relations_count")
 rm(r_cluster)
 
 item_relations_count <- merge(item_relations_count, cluster_count, by="cluster")
 
 #Get cluster for the items 
 item_relations_count2 <- merge(item_relations_count, clusters, by="id")
 item_relations_count2 <- item_relations_count2 %>% rename(cited_cluster = cluster.x, item_cluster = cluster.y)
 
 #1 For each observation i, calculate the average dissimilarity ai 
 #  between i and all other points of the cluster to which i belongs.
 item_rel_in_cluster <- item_relations_count2[item_relations_count2$item_cluster == item_relations_count2$cited_cluster,]
 item_rel_in_cluster <- item_rel_in_cluster[item_rel_in_cluster$cluster_size >= min_cluster_size,]
 item_rel_in_cluster$ai_a <- 1 - (item_rel_in_cluster$relations_count / (item_rel_in_cluster$cluster_size - 1))
 
 #2 For all other clusters C, to which i does not belong, calculate the average dissimilarity d(i,C) of i 
 #  to all observations of C. The smallest of these d(i,C) is defined as bi=minCd(i,C). 
 #  The value of bi can be seen as the dissimilarity between i and its “neighbor” cluster, i.e., the nearest one to which it does not belong.
 item_rel_outside <- item_relations_count2[item_relations_count2$item_cluster != item_relations_count2$cited_cluster,]
 item_rel_outside <- item_rel_outside[item_rel_outside$cluster_size >= min_cluster_size,]
 item_rel_outside$d_ic <- 1 - item_rel_outside$relations_count / item_rel_outside$cluster_size
 item_rel_outside_min <- item_rel_outside %>% group_by(id) %>% summarize(b_i=min(d_ic)) 
 
 #3 Finally the silhouette width of the observation i is defined by the formula: Si=(bi−ai)/max(ai,bi).
 item_rel <- merge(item_rel_in_cluster, item_rel_outside_min, by="id", all.x=T, all.y=T)
 
 #Fill na with 1 = max dissimilarity
 item_rel$ai_a[is.na(item_rel$ai_a)] <- 1
 item_rel$b_i[is.na(item_rel$b_i)] <- 1
 
 #silhouette 
 item_rel$silhouette <- (item_rel$b_i - item_rel$ai_a) / pmax(item_rel$b_i, item_rel$ai_a)
 
 #Avoid division with 0. If both b_i and ai_a is 0 then set silhouett to 0 (same similarity) 
 item_rel$silhouette[item_rel$b_i == 0 & item_rel$ai_a == 0] <- 0
 
 #Compile data from items
 item_total_rels <- item_relations_count2 %>% group_by(id) %>% summarize(total_rels = sum(relations_count))
 item_rel <- merge(item_rel, item_total_rels, by="id")
 item_rel$shr_rels_in_cluster <- item_rel$relations_count / item_rel$total_rels
 item_rel$shr_rels_in_cluster[is.na(item_rel$shr_rels_in_cluster)] <- 0
 
 print(paste0("Average silhouette value: ", round(mean(item_rel$silhouette, na.rm=TRUE),3)))
 print(paste0("Average % relations within cluster: ", round(mean(item_rel$shr_rels_in_cluster, na.rm=TRUE)*100,1)))
 rm(item_relations_count)
 rm(item_relations_count2)
 
 return(list(item_rel=item_rel))
}
#------------------------------------------------------------------------------

#Function to get the number of inaccurate assignments in clustering solution
#------------------------------------------------------------------------------
getInaccurate <- function(dfs,
                          max_shr_rels_in_cluster=0.1, #If bad assignment, a low proportion of relations are within cluster
                          min_total_rels=20, #Only count as bad assignment if the publications have enough relations
                          max_silhouette=0){ #The silhouette width should be lower than 0, i.e another cluster assignment would be better
 r_PIA <- c()
 for(k in 1:length(resolutions)){
  t <- dfs[[k]]$item_rel
  bad <- t[t$shr_rels_in_cluster < max_shr_rels_in_cluster
           & t$total_rels >= min_total_rels 
           & t$silhouette <= max_silhouette,]
  #print(paste0("Bad: ", nrow(bad)))
  r_PIA <- c(r_PIA,nrow(bad))
 }
 return(r_PIA)
}
#------------------------------------------------------------------------------

#Function to get ARI value of clustering solutions
#------------------------------------------------------------------------------
getARI <- function(r, year=ariYear){
 
 #Get pubs from ariYear onwards
 baseline_pmids <- r$pubs[r$pubs$publication_year > year,]
 
 #Delimit to reviews with at least n total relations
 baseline_pmids <- baseline_pmids[baseline_pmids$cnt_total_refs >= ariMinTotalRefs,]
 baseline_pmids <- baseline_pmids[baseline_pmids$cnt_refs / baseline_pmids$cnt_total_refs >= ariMinRefsClusterShare,]
 baseline_pmids <- unique(baseline_pmids$pmid)
 
 print(paste0("Initial number of baseline publications:", length(baseline_pmids)))
 
 if(length(baseline_pmids)>0){
  #Get baseline refs 
  baseline <- r$relations[r$relations$citing %in% baseline_pmids, ]
  baseline <- baseline[,c("citing","cited")]
  names(baseline) <- c("review_pmid", "cited_pmid")
  
  #There should be only 1 review per topic. If overlap chose only 1 
  overlap <- merge(baseline, baseline, by="cited_pmid")
  overlap <- overlap[overlap$review_pmid.x > overlap$review_pmid.y, ]
  overlap2 <- overlap %>% group_by(review_pmid.x,review_pmid.y) %>% summarise("bc" = n())
  
  #Count refs 
  baseline_no_refs <- baseline %>% group_by(review_pmid) %>% summarise("no_refs" = n())
  
  overlap2 <- merge(overlap2, baseline_no_refs, by.x="review_pmid.x", by.y="review_pmid")
  overlap2 <- merge(overlap2, baseline_no_refs, by.x="review_pmid.y", by.y="review_pmid")
  
  overlap2$overlap <- 0.5 * (overlap2$bc / overlap2$no_refs.x + overlap2$bc / overlap2$no_refs.y)
  
  #Get groups of same topic
  same_topic_pairs <- overlap2[overlap2$overlap >= ariOverlapSameTopic, ]
  if(nrow(same_topic_pairs) > 0){
   g <- graph_from_data_frame(same_topic_pairs)
   
   c <- components(g)
   
   c_df <- as.data.frame(c$membership)
   names(c_df) <- c("m")
   c_df$pmid <- row.names(c_df)
   
   rand <- c_df %>% group_by(m) %>% slice_sample(n=1)
   
   #Delete baseline pmids not to be used
   baseline <- baseline[!baseline$review_pmid %in% c_df$pmid[!c_df$pmid %in% rand$pmid],]
  }
  
  #Now make sure we do not have overlapping references 
  overlapping_refs <- baseline %>% group_by(cited_pmid) %>% summarise("count" = n())
  overlapping_refs <- overlapping_refs[overlapping_refs$count > 1,]
  
  if(nrow(overlapping_refs) > 0){
   
   #Refer overlapping refs to the review to which it has most connections
   ref_rels1 = r$relations[r$relations$citing %in% overlapping_refs$cited_pmid, ]
   ref_rels2 = r$relations[r$relations$cited %in% overlapping_refs$cited_pmid, ]
   ref_rels1 <- ref_rels1[,c("citing", "cited")]
   ref_rels2 <- ref_rels2[,c("cited", "citing")]
   
   names(ref_rels1) <- c("overlapping_pmid", "related_pmid")
   names(ref_rels2) <- c("overlapping_pmid", "related_pmid")
   
   ref_rels <- rbind(ref_rels1, ref_rels2)
   
   bl_overlapping_refs <- baseline[baseline$cited_pmid %in% overlapping_refs$cited_pmid,]
   names(bl_overlapping_refs) <- c("review_pmid", "overlapping_pmid")
   refs_relations <- merge(bl_overlapping_refs, baseline, by="review_pmid")
   names(refs_relations) <- c("review_pmid", "overlapping_pmid", "related_pmid")
   refs_relations <- merge(refs_relations, ref_rels, by=c("overlapping_pmid", "related_pmid"))
   
   count_overlap_rels <- refs_relations %>% group_by(review_pmid, overlapping_pmid) %>% summarise("sum_rels" = n())
   
   #Refer to review_pmid with highest connection
   rank_overlapping_over_reviews <- count_overlap_rels %>% arrange(overlapping_pmid, -sum_rels) %>% group_by(overlapping_pmid) %>% mutate("rank" = row_number(-sum_rels))
   refer_overlapping_to_review <- rank_overlapping_over_reviews[rank_overlapping_over_reviews$rank == 1, ]
   refer_overlapping_to_review <- refer_overlapping_to_review[,1:2]
   names(refer_overlapping_to_review) <- names(baseline)
   
   baseline <- baseline[!baseline$cited_pmid %in% overlapping_refs$cited_pmid,]
   
   baseline <- rbind(baseline, refer_overlapping_to_review)
   
   names(baseline) <- c("class", "item")
  }
  print(paste0("Number of baseline publications after restrictions:", length(unique(baseline$class))))
  
  #Calculate ARI 
  r_unnormalized_ARI <- calculateARI(r$unnormalized_item_rel_dfs, baseline)
  r_frac_ARI <- calculateARI(r$frac_item_rel_dfs, baseline)
  r_geometric_ARI <- calculateARI(r$geometric_item_rel_dfs, baseline)
  r_geometric_lim5_ARI <- calculateARI(r$geometric_lim5_item_rel_dfs, baseline)
  r_directed_frac_ARI <- calculateARI(r$directed_frac_item_rel_dfs, baseline)
  r_directed_geo_ARI <- calculateARI(r$directed_geo_item_rel_dfs, baseline)
  
  return(list("r_unnormalized_ARI"=r_unnormalized_ARI, 
              "r_frac_ARI"=r_frac_ARI, 
              "r_geometric_ARI"=r_geometric_ARI, 
              "r_geometric_lim5_ARI"=r_geometric_lim5_ARI,
              "r_directed_frac_ARI"=r_directed_frac_ARI,
              "r_directed_geo_ARI"=r_directed_geo_ARI))
 }
}

#------------------------------------------------------------------------------

#Function to get ARI from different resolutions
#------------------------------------------------------------------------------
calculateARI <- function(dfs, baseline){
 
 #Vector for the result 
 r_ARI <- c()
 
 #For each resolution 
 for(i in 1:length(resolutions)){
  
  #Delimit data to baseline
  compare_set <- dfs[[i]]$item_rel
  compare_set <- compare_set[compare_set$id %in% baseline$item, c("id","item_cluster")]
  
  if(length(compare_set$id) != length(baseline$item)){
   stop()
   print("Baseline set is not the same length as in df")
  }
  
  #Sort sets
  baseline <- baseline %>% arrange(item)
  compare_set <- compare_set %>% arrange(id)
  
  r_ARI_i <- adj.rand.index(baseline$class, compare_set$item_cluster)
  
  r_ARI <- c(r_ARI, r_ARI_i)
 }
 
 return(r_ARI)
}

#Get skewness
getSkewness <- function(dfs){
 r_skew <- c()
 for(k in 1:length(resolutions)){
  t <- dfs[[k]]$item_rel
  
  t_cnt <- t %>% group_by(item_cluster) %>% summarise("cnt" = n())
  t_cnt <- t_cnt[!is.na(t_cnt$item_cluster),]
  
  skewvalue = round(skewness(t_cnt$cnt),2)
  
  r_skew = c(r_skew, skewvalue)
 }
 return(r_skew)
}
#------------------------------------------------------------------------------

#===================================================
#PLOTS
#===================================================

#Function to make GA plots
#------------------------------------------------------------------------------
plotGA <- function(field_name,
                   indicator_name,
                   unnormalized_granularities,
                   frac_granularities,
                   geometric_granularities,
                   geometric_lim5_granularities,
                   directed_frac_granularities,
                   directed_geo_granularities,
                   unnormalized_values, 
                   frac_values, 
                   geometric_values,
                   geometric_lim5_values,
                   directed_frac_values,
                   directed_geo_values){
 
 #Output file
 #png(paste0(outputFolder, "/",field_name, "_", indicator_name,".png"), width=22, height=12, unit="cm", res=600)
 
 #Set margins
 #par(mar=c(5.1, 4.1, 4.1, 14), xpd=TRUE)
 
 #Get max and min values
 values <- c(unnormalized_values,frac_values,geometric_values,geometric_lim5_values, directed_frac_values, directed_geo_values)
 values[values == Inf] <- NA
 ylim_max = max(values, na.rm=T)
 ylim_min = min(values, na.rm=T)
 xlim_max = max(c(unnormalized_granularities,frac_granularities,geometric_granularities,geometric_lim5_granularities, directed_frac_granularities, directed_geo_granularities))
 xlim_min = min(c(unnormalized_granularities,frac_granularities,geometric_granularities,geometric_lim5_granularities, directed_frac_granularities, directed_geo_granularities))
 
 #Create the plot
 plot(NA, 
      type="l", 
      log="x", 
      xlim=c(xlim_min,xlim_max), 
      ylim = c(ylim_min, ylim_max), 
      col="white", 
      xlab="Granularity", 
      ylab=indicator_name,
      cex.lab=1.5, cex.axis=1 # cex.main=2, cex.sub=1.5
 )
 xspline(unnormalized_granularities[!is.infinite(unnormalized_values)], unnormalized_values[!is.infinite(unnormalized_values)], lty=1, shape = 1, border=safe_colorblind_palette[1], pch=19)
 points(unnormalized_granularities[!is.infinite(unnormalized_values)], unnormalized_values[!is.infinite(unnormalized_values)], col=safe_colorblind_palette[1])
 xspline(frac_granularities[!is.infinite(frac_values)], frac_values[!is.infinite(frac_values)], lty=1, shape = 1, border=safe_colorblind_palette[2])
 points(frac_granularities[!is.infinite(frac_values)], frac_values[!is.infinite(frac_values)], col=safe_colorblind_palette[2])
 xspline(geometric_granularities[!is.infinite(geometric_values)], geometric_values[!is.infinite(geometric_values)], lty=1, shape = 1, border=safe_colorblind_palette[3])
 points(geometric_granularities[!is.infinite(geometric_values)], geometric_values[!is.infinite(geometric_values)], col=safe_colorblind_palette[3])
 xspline(geometric_lim5_granularities[!is.infinite(geometric_lim5_values)], geometric_lim5_values[!is.infinite(geometric_lim5_values)], lty=1, shape = 1, border=safe_colorblind_palette[4])
 points(geometric_lim5_granularities[!is.infinite(geometric_lim5_values)], geometric_lim5_values[!is.infinite(geometric_lim5_values)], col=safe_colorblind_palette[4])
 xspline(directed_frac_granularities[!is.infinite(directed_frac_values)], directed_frac_values[!is.infinite(directed_frac_values)], lty=1, shape = 1, border=safe_colorblind_palette[5])
 points(directed_frac_granularities[!is.infinite(directed_frac_values)], directed_frac_values[!is.infinite(directed_frac_values)], col=safe_colorblind_palette[5])
 xspline(directed_geo_granularities[!is.infinite(directed_geo_values)], directed_geo_values[!is.infinite(directed_geo_values)], lty=1, shape = 1, border=safe_colorblind_palette[6])
 points(directed_geo_granularities[!is.infinite(directed_geo_values)], directed_geo_values[!is.infinite(directed_geo_values)], col=safe_colorblind_palette[6])
 title(main=field_name, adj=0, font.main=1, cex.main=2)
 
 #dev.off()
 
}

gaPlotToFile <- function(fileName, measure, field_names, rs, approaches, rs_names){
 
 #Output file
 png(paste0(outputFolder, "/", fileName), width=22, height=25, unit="cm", res=600)
 
 layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE), heights=c(4, 4, 2))
 
 par(mai=c(0.7, 0.7, 0.5, 0.1), oma=c(1,1,1,1))
 
 for(i in 1:length(field_names)){
  plotGA(field_names[i],
         indicator_name=measure,
         unnormalized_granularities=rs[[i]]$unnormalized_granularities,
         frac_granularities=rs[[i]]$frac_granularities,
         geometric_granularities=rs[[i]]$geometric_granularities,
         geometric_lim5_granularities=rs[[i]]$geometric_lim5_granularities,
         directed_frac_granularities=rs[[i]]$directed_frac_granularities,
         directed_geo_granularities=rs[[i]]$directed_geo_granularities,
         unnormalized_values=rs[[i]][rs_names[1]][[1]], 
         frac_values=rs[[i]][rs_names[2]][[1]],
         geometric_values=rs[[i]][rs_names[3]][[1]],
         geometric_lim5_values=rs[[i]][rs_names[4]][[1]],
         directed_frac_values=rs[[i]][rs_names[5]][[1]],
         directed_geo_values=rs[[i]][rs_names[6]][[1]])
 }
 par(mai=c(0,0,0,0))
 plot.new()
 legend("center", inset=c(0,0), approaches,
        lty = c(1,1,1,1,1,1),
        col = safe_colorblind_palette[1:6],
        cex = 1.5)
 dev.off()
}

createAllPlots <- function(field_names, rs, approaches){
 
 #Plot histograms over the number of relations per publication
 for(i in 1:length(field_names)){
  plot_hist_relations(field_names[i], rs[[i]])
 }
 
 #Plot Silhouette width
 rs_names <- c("unnormalized_silhouttes_means",
               "frac_silhouttes_means",
               "geometric_silhouttes_means",
               "geometric_lim5_silhouttes_means",
               "directed_frac_silhouttes_means",
               "directed_geo_silhouttes_means")
 measure = "Silhouette width"
 gaPlotToFile("silhouette_width.png", measure, field_names, rs, approaches, rs_names)
 
 #Internal citations
 rs_names <- c("unnormalized_internal_citations",
               "frac_internal_citations",
               "geometric_internal_citations",
               "geometric_lim5_internal_citations",
               "directed_frac_internal_citations",
               "directed_geo_internal_citations")
 gaPlotToFile("internal_citations.png", measure="# Internal citations", field_names, rs, approaches, rs_names)
 
 # Probably inaccurate
 for(i in 1:length(field_names)){
  rs[[i]]$unnormalized_PIA <- getInaccurate(rs[[i]]$unnormalized_item_rel_dfs)
  rs[[i]]$frac_PIA <- getInaccurate(rs[[i]]$frac_item_rel_dfs)
  rs[[i]]$geometric_PIA <- getInaccurate(rs[[i]]$geometric_item_rel_dfs)
  rs[[i]]$geometric_lim5_PIA <- getInaccurate(rs[[i]]$geometric_lim5_item_rel_dfs)
  rs[[i]]$directed_frac_PIA <- getInaccurate(rs[[i]]$directed_frac_item_rel_dfs)
  rs[[i]]$directed_geo_PIA <- getInaccurate(rs[[i]]$directed_geo_item_rel_dfs)
 }
 rs_names <- c("unnormalized_PIA",
               "frac_PIA",
               "geometric_PIA",
               "geometric_lim5_PIA",
               "directed_frac_PIA",
               "directed_geo_PIA")
 gaPlotToFile("probably_inaccurate.png", measure="PIA", field_names, rs, approaches, rs_names)
 
 #ARI
 r_ari <- list()
 for(i in 1:length(field_names)){
  #Get ARI for each field
  r_ari[[i]] <- getARI(rs[[i]])
  
  rs[[i]]$unnormalized_ARI <- r_ari[[i]]$r_unnormalized_ARI
  rs[[i]]$frac_ARI <- r_ari[[i]]$r_frac_ARI
  rs[[i]]$geometric_ARI <- r_ari[[i]]$r_geometric_ARI
  rs[[i]]$geometric_lim5_ARI <- r_ari[[i]]$r_geometric_lim5_ARI
  rs[[i]]$directed_frac_ARI <- r_ari[[i]]$r_directed_frac_ARI
  rs[[i]]$directed_geo_ARI <- r_ari[[i]]$r_directed_geo_ARI
 }
 rs_names <- c("unnormalized_ARI",
               "frac_ARI",
               "geometric_ARI",
               "geometric_lim5_ARI",
               "directed_frac_ARI",
               "directed_geo_ARI")
 gaPlotToFile("ari.png", measure="ARI", field_names, rs, approaches, rs_names)
 
 # Plot skewness
 for(i in 1:length(field_names)){
  rs[[i]]$unnormalized_skew <- getSkewness(rs[[i]]$unnormalized_item_rel_dfs)
  rs[[i]]$frac_skew <- getSkewness(rs[[i]]$frac_item_rel_dfs)
  rs[[i]]$geometric_skew <- getSkewness(rs[[i]]$geometric_item_rel_dfs)
  rs[[i]]$geometric_lim5_skew <- getSkewness(rs[[i]]$geometric_lim5_item_rel_dfs)
  rs[[i]]$directed_frac_skew <- getSkewness(rs[[i]]$directed_frac_item_rel_dfs)
  rs[[i]]$directed_geo_skew <- getSkewness(rs[[i]]$directed_geo_item_rel_dfs)
 }
 
 rs_names <- c("unnormalized_skew",
               "frac_skew",
               "geometric_skew",
               "geometric_lim5_skew",
               "directed_frac_skew",
               "directed_geo_skew")
 gaPlotToFile("skewness.png", measure="Skewness", field_names, rs, approaches, rs_names)
 
 #return all values plotted
 #=========================
 df <- data.frame("Field"=character(), "Approach"=character(), "Resolution"=numeric(), "Granularity"=numeric(), "Silhouette width"=numeric(), "Internal citations"=integer(), "PIA"=integer(), "ARI"=numeric(), "Skewness"=numeric())
 
 for(i in 1:length(field_names)){
  df_temp <- data.frame("Field"=field_names[i], 
                        "Approach"="Unnormalized",
                        "Resolution" = resolutions, 
                        "Granularity"=rs[[i]]$unnormalized_granularities, 
                        "Silhouette.width"=rs[[i]]$unnormalized_silhouttes_means,
                        "Internal.citations"=rs[[i]]$unnormalized_internal_citations,
                        "PIA"=rs[[i]]$unnormalized_PIA,
                        "ARI"=rs[[i]]$unnormalized_ARI,
                        "Skewness"=rs[[i]]$unnormalized_skew
  )   
  df <- rbind(df, df_temp)
 }
 for(i in 1:length(field_names)){
  df_temp <- data.frame("Field"=field_names[i], 
                        "Approach"="Fractional",
                        "Resolution" = resolutions, 
                        "Granularity"=rs[[i]]$frac_granularities, 
                        "Silhouette.width"=rs[[i]]$frac_silhouttes_means,
                        "Internal.citations"=rs[[i]]$frac_internal_citations,
                        "PIA"=rs[[i]]$frac_PIA,
                        "ARI"=rs[[i]]$frac_ARI,
                        "Skewness"=rs[[i]]$frac_skew
  )   
  df <- rbind(df, df_temp)
 }
 for(i in 1:length(field_names)){
  df_temp <- data.frame("Field"=field_names[i], 
                        "Approach"="Geometric",
                        "Resolution" = resolutions, 
                        "Granularity"=rs[[i]]$geometric_granularities, 
                        "Silhouette.width"=rs[[i]]$geometric_silhouttes_means,
                        "Internal.citations"=rs[[i]]$geometric_internal_citations,
                        "PIA"=rs[[i]]$geometric_PIA,
                        "ARI"=rs[[i]]$geometric_ARI,
                        "Skewness"=rs[[i]]$geometric_skew
  )   
  df <- rbind(df, df_temp)
 }
 for(i in 1:length(field_names)){
  df_temp <- data.frame("Field"=field_names[i], 
                        "Approach"="Geometric",
                        "Resolution" = resolutions, 
                        "Granularity"=rs[[i]]$geometric_lim5_granularities, 
                        "Silhouette.width"=rs[[i]]$geometric_lim5_silhouttes_means,
                        "Internal.citations"=rs[[i]]$geometric_lim5_internal_citations,
                        "PIA"=rs[[i]]$geometric_lim5_PIA,
                        "ARI"=rs[[i]]$geometric_lim5_ARI,
                        "Skewness"=rs[[i]]$geometric_lim5_skew
  )   
  df <- rbind(df, df_temp)
 }
 for(i in 1:length(field_names)){
  df_temp <- data.frame("Field"=field_names[i], 
                        "Approach"="unnormalized",
                        "Resolution" = resolutions, 
                        "Granularity"=rs[[i]]$directed_frac_granularities, 
                        "Silhouette.width"=rs[[i]]$directed_frac_silhouttes_means,
                        "Internal.citations"=rs[[i]]$directed_frac_internal_citations,
                        "PIA"=rs[[i]]$directed_frac_PIA,
                        "ARI"=rs[[i]]$directed_frac_ARI,
                        "Skewness"=rs[[i]]$directed_frac_skew
  )   
  df <- rbind(df, df_temp)
 }
 for(i in 1:length(field_names)){
  df_temp <- data.frame("Field"=field_names[i], 
                        "Approach"="unnormalized",
                        "Resolution" = resolutions, 
                        "Granularity"=rs[[i]]$directed_geo_granularities, 
                        "Silhouette.width"=rs[[i]]$directed_geo_silhouttes_means,
                        "Internal.citations"=rs[[i]]$directed_geo_internal_citations,
                        "PIA"=rs[[i]]$directed_geo_PIA,
                        "ARI"=rs[[i]]$directed_geo_ARI,
                        "Skewness"=rs[[i]]$directed_geo_skew
  )   
  df <- rbind(df, df_temp)
 }
 
 return(df)
}

#------------------------------------------------------------------------------

#Function to run all plots
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#Function to plot histogram of number of relations per publication
#------------------------------------------------------------------------------
plot_hist_relations <- function(field_name, r){
 
 pubs_cnt1 <- r$pubs[r$pubs$cnt_relations > 0 & r$pubs$cnt_relations <= 100,]
 
 mainname=paste0(field_name, " - Histogram of publication relation counts\n(restricted to 100 relations or less)")
 
 #get bins
 n_bins <- ceiling(max(pubs_cnt1$cnt_relations))
 
 outputfilename=paste0(field_name, "_hist_relations.png")
 #Output to file
 png(paste0(outputFolder, "/",outputfilename), width=22, height=12, unit="cm", res=600)
 
 image <- ggplot(data.frame(pubs_cnt1), aes(cnt_relations)) +               
  geom_histogram(bins = n_bins, color = "#000000", fill="#DCDCDC") +
  #scale_x_log10() + 
  theme_classic() + 
  scale_color_manual(values=c("#FFF000")) +
  ggtitle(mainname) + 
  xlab("# Relations") + 
  ylab("Count")
 
 print(image)
 
 dev.off()
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
if(getNewResults){
 #--------------------------------------- Social psychology
 r_sp <- createOutput(
  pmids_file = paste0(dataFolder, "/social_psycology_pmids_1995_2021_2023-08-31.csv"),
  relations_file = paste0(dataFolder, "/social_psycology_relations_1995_2021_2023-08-31.csv"),
  field_name = "Social psychology",
  consoleFile = paste0(outputFolder, "/console_sp.txt")
 )
 saveRDS(r_sp, paste0(outputFolder,"/r_sp.rdata"))
 #--------------------------------------- Autoimmune diseases
 r_aid <- createOutput(
  pmids_file = paste0(dataFolder, "/autoimmune_diseases_pmids_1995_2021_2023-08-31.csv"),
  relations_file = paste0(dataFolder, "/autoimmune_diseases_relations_1995_2021_2023-08-31.csv"),
  field_name = "Autoimmune diseases",
  consoleFile = paste0(outputFolder, "/console_aid.txt")
 )
 saveRDS(r_aid, paste0(outputFolder,"/r_aid.rdata"))
 #--------------------------------------- Metabolism
 r_mb <- createOutput(
  pmids_file = paste0(dataFolder, "/metabolism_pmids_1995_2021_2023-08-31.csv"),
  relations_file = paste0(dataFolder, "/metabolism_relations_1995_2021_2023-08-31.csv"),
  field_name = "Metabolism",
  consoleFile = paste0(outputFolder, "/console_mb.txt")
 )
 saveRDS(r_mb, paste0(outputFolder,"/r_mb.rdata"))
 #--------------------------------------- Stem cells
 r_sc <- createOutput(
  pmids_file = paste0(dataFolder, "/Stem Cells_pmids_1995_2021_2023-08-31.csv"),
  relations_file = paste0(dataFolder, "/Stem Cells_relations_1995_2021_2023-08-31.csv"),
  field_name = "Stem cells",
  consoleFile = paste0(outputFolder, "/console_sc.txt")
 )
 saveRDS(r_sc, paste0(outputFolder,"/r_sc.rdata"))
 #---------------------------------------
 
} else {
 r_sp <- readRDS(paste0(readFromFolder,"/r_sp.rdata"))
 r_aid <- readRDS(paste0(readFromFolder,"/r_aid.rdata"))
 r_mb <- readRDS(paste0(readFromFolder,"/r_mb.rdata"))
 r_sc <- readRDS(paste0(readFromFolder,"/r_sc.rdata"))
}

#Print console to file
con <- file(paste0(outputFolder, "/console_plots.txt"))
write(' ', file=con, append=F)
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

# Create plots
#=============
field_names <- c("Social psychology", "Autoimmune diseases", "Metabolism", "Stem cells")
rs <- list(r_sp, r_aid, r_mb, r_sc)
approaches <- c("Unnormalized", "Fractional", "Geometric mean", "Geometric mean-limit5", "Directional-fractional", "Directional-geometric")

plot_data <- createAllPlots(field_names, rs, approaches)

write.csv(plot_data, paste0(outputFolder, "/plot_data.csv"))

# Plot violine of relations distribution
#=======================================

#Get data
data <- rbind(data.frame(name="Social psycology",value=r_sp$pubs$cnt_relations),
              data.frame(name="Autoimmune diseases",value=r_aid$pubs$cnt_relations),
              data.frame(name="Metabolism",value=r_mb$pubs$cnt_relations),
              data.frame(name="Stem cells",value=r_sc$pubs$cnt_relations))

max_relations = 100

#Sample size - for printing

sample_size_1to10 = data[data$value >= 1 & data$value <= 10,] %>% group_by(name) %>% summarize("1_10"=n())
sample_size_11to100 = data[data$value >= 11 & data$value <= 100,] %>% group_by(name) %>% summarize("11_100"=n())
sample_size_101to1000 = data[data$value >= 101 & data$value <= 1000,] %>% group_by(name) %>% summarize("101_1000"=n())
sample_size_1001 = data[data$value > 1000,] %>% group_by(name) %>% summarize(">1000"=n())

rel_stats <- data.frame(sample_size_1to10$name, 
                        sample_size_1to10$"1_10", 
                        sample_size_11to100$"11_100", 
                        sample_size_101to1000$"101_1000", 
                        sample_size_1001$`>1000`)

write.csv(rel_stats, paste0(outputFolder, "/rel_stats.csv"), row.names = FALSE)

sample_size = data[data$value <= 100,] %>% group_by(name) %>% summarize(num=n())

png(paste0(outputFolder, "/violine_plot_num_relations.png"), width=22, height=12, unit="cm", res=600)

#Plot (restrict to n relations or less)
data[data$value <= max_relations,] %>%
 left_join(sample_size) %>%
 mutate(myaxis = paste0(name, "\n", "n=", num)) %>%
 ggplot(aes(x=myaxis, y=value, fill=name)) +
 geom_violin(width=1.2) +
 geom_boxplot(width=0.1, color="black", alpha=0) +
 scale_fill_manual(values=c("#EBECF0","#EBECF0","#EBECF0","#EBECF0")) +
 theme_classic() +
 theme(
  legend.position="none",
  plot.title = element_text(size=16),
  axis.text = element_text(size=14),
  axis.title = element_text(size=14)
 ) +
 ggtitle("") +
 xlab("") + ylab("# Relations")

dev.off()

# Restore output to console
sink() 
sink(type="message")
#------------------------------------------------------------------------------

# print elapsed time
endTimeStamp <- Sys.time() - startTimeStamp 
print(endTimeStamp) 
