R_square <- function(preded_Y, true_Y) {
  colnames(preded_Y) <- 'Predited_Y'
  Pre_and_True <- merge(preded_Y, true_Y, by = "row.names")
  colnames(Pre_and_True) <- c('samples', 'preds', 'ture')
  Pre_and_True$samples <- NULL
  Pre_and_True <- as.matrix(Pre_and_True)
  rss <- sum((Pre_and_True[,1] - Pre_and_True[,2]) ^ 2)  ## residual sum of squares
  tss <- sum((Pre_and_True[,2] - mean(Pre_and_True[,2])) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  print('The R-square of predicted value to ture value is')
  print(rsq)
  return(rsq)
}
