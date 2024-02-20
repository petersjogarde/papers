getRAfromSeeds <- function(seeds, sysSleepTime=0.4){

 ra_score_list <- vector("list", length = length(seeds))
 temp_ra_score <- data.frame()
 
 ra_score_dfs <- vector("list", length = length(seeds))
 
 for(i in 1:length(seeds)){
  seed <- seeds[i]
  
  #Try to get results 5 times 
  r <- NULL 
  q <- ""
  q <- paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pubmed&id=', seed, '&cmd=neighbor_score&retmode=json')
  
  attempt <- 1
  while(is.null(r) && attempt <= 5) {
   print(paste0("Getting RA. Try: ", attempt, ". Seed: ", seed))
   
   check = FALSE
   tryCatch({
    
    r <- GET(q)
    
    #Check if content is correct - otherwise stop
    rContent <- jsonlite::fromJSON(content(r, type="text"), simplifyVector = FALSE)
    
    if(rContent$linksets[[1]][[3]][[1]]$linkname == "pubmed_pubmed") check = TRUE
    
   }, error = function(e) {
     message('Nothing retrieved from PubMed')
   })
   
   if(!check){
    r <- NULL
   }
   
   attempt <- attempt + 1
   Sys.sleep(3) #test every second
  } 
  
  if(check){
   t <- rContent$linksets[[1]][[3]][[1]]$links
   
   related_pmid <- vector("integer", length(t))
   score <- vector("integer", length(t))
   for(j in 1:length(t)){
    related_pmid[j] <- as.integer(t[[j]]$id)
    score[j] <- t[[j]]$score
   }
   
   ra_score_dfs[[i]] <- data.frame(related_pmid, score)
  }
 }

 Sys.sleep(sysSleepTime)
 
 # Bind all data frames
 ra_score_df <- do.call("rbind", ra_score_dfs)
 
 # Sum scores 
 ra_score <- ra_score_df %>% group_by(related_pmid) %>% summarise(score=sum(score))

 #Delmit to years
 raPmids <- paste0(ra_score_df$related_pmid, collapse=",")
 raPmidsdelimited <- dbGetQuery(db, paste0("
SELECT
 omp.pmid
FROM ", nihOccMetadataTable, " omp
WHERE omp.pmid IN (",raPmids, ")
 AND omp.year BETWEEN ", refFirstYear," AND ", refLastYear,"
"))
 
 ra_score <- ra_score[ra_score$related_pmid %in% raPmidsdelimited$pmid,]
 
 ra_score <- ra_score[!ra_score$related_pmid %in% seeds,]

 print("Done getting RA")
 return(ra_score) 
}
