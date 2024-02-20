getBCfromSeeds <- function(seedPmids_str, minBCscore){
 
 sql <- paste0("WITH seeds AS (
 SELECT
   unnest(array[",seedPmids_str,"]) AS seed_pmid
 )
 , t AS (
 	SELECT DISTINCT
 	 seeds.seed_pmid,
 	 cp2.citing AS related_pmid,
 	 cp.cited
 	FROM seeds 
 	INNER JOIN ", nihOccCitationsTable, " cp ON (cp.citing = seeds.seed_pmid) -- cp.cited includes refs from seeds
 	INNER JOIN ", nihOccCitationsTable, " cp2 ON (cp2.cited = cp.cited AND cp2.citing NOT IN (", seedPmids_str, ")) -- cp2.citing includes pmids citing the same refs (seeds should not be connected to seeds by other seeds or it self)
 ) 
 , y AS (
  SELECT 
   pmid 
  FROM ", nihOccMetadataTable, " 
  WHERE pmid IN (SELECT DISTINCT related_pmid FROM t) AND year BETWEEN ", refFirstYear," AND ", refLastYear,"
 )
 SELECT 
  t.related_pmid,
  count(*) AS score
 FROM t 
 GROUP BY 1
 HAVING count(*) >= ", minBCscore, "
 ORDER BY score DESC
 ")
 
 df <- dbGetQuery(db, sql)
 
 return(df)
}
