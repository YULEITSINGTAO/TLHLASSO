cross_lasso <- function(X,Y, Group) {
  suppressPackageStartupMessages({library(glmnet)})
  print("Doing the cross validation of Lasso Regression")
  #Perform lasso on 10 groups
  Pre_Y <- list()
  for(i in 1: length (Group)){
    testX <- X[Group[[i]],]
    testY <- Y[Group[[i]], ]
    train_names <- setdiff(rownames(X), Group[[i]])
    trainX <- X[train_names, ]
    trainY <- Y[train_names, ]

    #CV Lasso
    ################
    # Fit lasso model on training data
    cat("Doing the ", i, "th", " group training", "\n")
    cv.out = cv.glmnet(trainX, trainY, type.measure = 'mse', nfolds = 10)
    cat(i, "th", " group training finished", "\n")
    response <- predict(cv.out, testX, s= "lambda.min")
    rownames(response) <- Group[[i]]
    Pre_Y[[i]] <- response
  }
  Predicted_Y <- do.call(rbind, Pre_Y)
  return(Predicted_Y)
}
