getCCfromSeeds <- function(seedPmids_str, removePmids=NULL, minCCscore){
 
 remove_pmids_str <- paste0(removePmids, collapse=",")

 sql <- paste0("WITH seeds AS (
 SELECT
   unnest(array[",seedPmids_str,"]) AS seed_pmid
 )
 , t AS (
 	SELECT DISTINCT
 	 seeds.seed_pmid,
 	 cp2.cited AS related_pmid,
 	 cp.citing
 	FROM seeds 
 	INNER JOIN ", nihOccCitationsTable, " cp ON (cp.cited = seeds.seed_pmid AND cp.citing NOT IN (", remove_pmids_str, ")) -- seeds are in reference lists of cp.citing
 	INNER JOIN ", nihOccCitationsTable, " cp2 ON (cp2.citing = cp.citing AND cp2.cited != cp.cited) -- cp2.cited are also in reference list of cp.citing
 )
 , y AS (
  SELECT 
   pmid, 
   ac.year 
  FROM ", nihOccMetadataTable, " ac 
  WHERE pmid IN (SELECT DISTINCT related_pmid FROM t) AND year BETWEEN ", refFirstYear," AND ", refLastYear,"
 )
 SELECT 
  t.related_pmid,
  count(*) AS score
 FROM t
 WHERE t.related_pmid IN (SELECT pmid FROM y)
 GROUP BY 1
 HAVING count(*) >= ", minCCscore, "
 ORDER BY score DESC
 ")
 
 df <- dbGetQuery(db, sql)
  
 return(df)
}
