#============================ Function to get precision and recall from a data frame
#inputDf - should include columns: pmid, related_pmid, score
getPrecisionRecall <- function(df, pmidsToFind, maxRank){
 
 #Mark hits 
 df$hit <- 0
 df$hit[df$related_pmid %in% pmidsToFind] <- 1
 
 #Rank
 df$rand <- sample(1:nrow(df)) #Add a random variable
 df <- df[order(-df$score, df$rand),]
 df$rnk <- 1:nrow(df)

 #Total precision and recall
 total_recall <- sum(df$hit) / length(pmidsToFind) * 100
 total_precision <- mean(df$hit) * 100
 
 #Total number of retrieved publications
 total_count <- nrow(df)
 
 #Restrict to maxRank
 df <- df[df$rnk <= maxRank, ]
 
 #Precision and recall per rank (up to maxRank)
 df$hit_sum <- cumsum(df$hit) #Get cumulative sum
 if(total_recall == 0){
  df$rnk_recall <- 0
 } else{
  df$rnk_recall <- df$hit_sum / length(pmidsToFind) * 100
 }
 df$rnk_precision <- df$hit_sum / 1:nrow(df) * 100

 df <- df[, c("related_pmid", "rnk", "rnk_recall", "rnk_precision", "hit", "hit_sum", "score")]
 return(list("df"=df, "total_recall"=total_recall, "total_precision"=total_precision, "total_count"=total_count))
}
