#Direct citation connections - both cited and citing publications
#==============================
#Input: 
# seedPmids_str = a comma separated string of pmids. 
# removePmids = a vector of pmids that should not be regarded as citing/cited. 
getDCfromSeeds <- function(seedPmids_str, removePmids=NULL){

 sql <- paste0("
  WITH seeds AS (
    SELECT
    unnest(array[",seedPmids_str,"]) AS seed_pmid
  )
  -- pmids citing the seeds
  SELECT 
   seeds.seed_pmid,
   cp.cited AS related_pmid
  FROM seeds
  INNER JOIN ", nihOccCitationsTable, " cp ON (cp.citing = seeds.seed_pmid)
  INNER JOIN ", nihOccMetadataTable, " ac ON (ac.pmid = cp.cited) -- meta data of citing set 
  WHERE 
   ac.year BETWEEN ", refFirstYear," AND ", refLastYear,"
  -- pmids cited by the seeds 
  UNION ALL
  SELECT 
   seeds.seed_pmid,
   cp.citing AS related_pmid
  FROM seeds
  INNER JOIN ", nihOccCitationsTable, " cp ON (cp.cited = seeds.seed_pmid)
  INNER JOIN ", nihOccMetadataTable, " ac ON (ac.pmid = cp.citing)
  WHERE 
   ac.year BETWEEN ", refFirstYear," AND ", refLastYear,"
 ")
 
 #Execulte sql 
 dc_df <- dbGetQuery(db, sql)
 
 if(!is.null(removePmids)){
  dc_df <- dc_df[!dc_df$related_pmid %in% removePmids, ]
 }
 
 #Count DC
 dc_score <- dc_df %>% 
  group_by(related_pmid) %>% 
  summarise(score = n())
 
 #Order DC 
 dc_score <- dc_score[order(-dc_score$score),]
 
 return(dc_score)
}
#==============================