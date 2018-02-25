########################################################################################
#               Part I: Use KNN classified method on original data                     #
#                                                              by Shiou-Shiou Deng     #
########################################################################################

# Load data
data_numeric <- read.csv("train_numeric.csv", header = TRUE, sep = ",", nrows = 10000)
data_categorical <- read.csv("train_categorical.csv", header = TRUE, sep = ",", nrows = 10000, stringsAsFactors = FALSE)

# Data processing
# Replace all Na value to 0 in numerica data
data_numeric[is.na(data_numeric)] <- 0

# Make a data.frame to represent the categorical data
categorical <- matrix(0, nrow = nrow(data_categorical), ncol = ncol(data_categorical))
for (i in 2:ncol(data_categorical)){
  if (is.logical(data_categorical[,i]) == FALSE){
    indexw <- which(data_categorical[,i] %in% "" == FALSE)
    categorical[,i][indexw] <- 1
  }
}
categorical <- as.data.frame(categorical)

# Divide the data into 2 part (80% training and validation, 20% testing)
set.seed(1234)
index <- sample(2, nrow(data_numeric), replace=TRUE, prob=c(0.8, 0.2))
train <- data_numeric[index==1,]
trainC <- categorical[index==1,]
test <- data_numeric[index==2,] 
testC <- categorical[index==2,]

# Make a function to Normalize the data
normalize <- function(x) {
  denominator <- sqrt(sum(x^2))
  if (denominator!=0){
    norm2 <- x/denominator
  } else {
    norm2 <- x
  }
  return(norm2)
}
# Divide the data to the training data and the response, and then nornomalize the data.
train_response <- train[,970]
train <- train[,2:969]
trainC <- trainC[,2:ncol(trainC)]
trainC <- matrix((sapply(trainC, normalize)), nrow = nrow(trainC))
train <- matrix((sapply(train, normalize)), nrow = nrow(train))
test_response <- test[,970]
test <- test[,2:(ncol(test)-1)]
testC <- testC[,2:ncol(testC)]
testC <- matrix(sapply(testC, normalize), nrow = nrow(testC))
test <- matrix(sapply(test, normalize), nrow = nrow(test)) 

# Make a assement function to evaluate the results
Assement <- function(TP, FP, FN, TN){
  Accuracy = (TP+TN)/ (TP+TN+FP+FN)
  MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  Sensitivity = TP / (TP+FN) 
  Precision = TP/(TP+FP)
  F.measure = 2*TP/(2*TP+FN+FP)
  return(c(Accuracy, Sensitivity, Precision, F.measure, MCC))
}
# Setting the k value
k <- 20
# Make a function to perform knn model on combination data
knn.comb <- function(train, trainC, valid, validC, train_response, k, alpha = 1, beta = 0){
  predict_response <- numeric(nrow(valid))
  # Compute the distances from the each valid data to all train data
  distancesC <- beta*(1 - validC %*% t(trainC))
  distancesN <- alpha*(1 - valid %*% t(train))
  distances <- distancesC + distancesN
  # predict the new variable
  for (m in 1:nrow(valid)) {
    neighbors_response <- numeric(k)
    ind <- order(distances[m,])
    neighbors_response <- train_response[ind][1:k]
    if (sum(neighbors_response) == (k/2)) {
      predict_response[m] <- ifelse(runif(1,0,1)>0.5, 1, 0)
    } else if (sum(neighbors_response) > (k/2)){ 
      predict_response[m] <- 1
    } else {
      predict_response[m] <- 0
    }
  }
  return(predict_response)
}

# Perform KNN on train data using different k, alpha, beta
alpha <- seq(0, 1, by = 0.1)
beta <- 1 - alpha
totalresults <- matrix(0, nrow = 0, ncol = 5)
# Use different combinations of weights for numerical data and categorical data (alpha, beta)
for (j in 1:length(alpha)){
  set.seed(1234)
  # Generate the numbers for 10 folds
  folds_i <- sample(rep(1:10, length.out = nrow(train)))
  # Store the results after applying models
  resultsB <- matrix(0, nrow = k, ncol = 5)
  colnames(resultsB) <- c("Accuracy", "Sensitivity", "Precision", "F-measure","MCC" )
  rownames(resultsB) <- rownames(resultsB, do.NULL = FALSE, prefix = "k = ")
  for (i in 1:k){
    cvresult <- matrix(0, nrow = 10, ncol = 5)
    # Use 10 fold cross-validation
    for (n in 1:10) {
      valid_i <- which(folds_i == n)
      train_dataN <- train[-valid_i, ]
      train_dataC <- trainC[-valid_i, ]
      trainlabel <- train_response[-valid_i]
      valid_dataN <- train[valid_i, ]
      valid_dataC <- trainC[valid_i, ]
      validlabel <- train_response[valid_i]
      cat("alpha= ", alpha[j],"\n")
      cat("k= ", i,"\n")
      cat("fold= ", n,"\n")
      # Perform kNN method
      model <- knn.comb(train_dataN, train_dataC, valid_dataN, valid_dataC, trainlabel, k=i, alpha = alpha[j], beta = beta[j])
      confusion.matrix <- table(factor(validlabel), factor(model, levels = c(0,1))) 
      cvresult[n,] <- Assement(confusion.matrix[4], confusion.matrix[3], confusion.matrix[2], confusion.matrix[1])
      resultsB[i,] <- colMeans(cvresult)
    }
  }
  totalresults <- rbind(totalresults, resultsB)
}

# Apply the optimatized models on test data
# Perform KNN on train data using different k, alpha, beta
totalresultsT <- matrix(0, nrow = 0, ncol = 5)
# Store the results after applying models
k <- 2
# Use different combinations of weights for numerical data and categorical data (alpha, beta)
for (j in 1:length(alpha)){
  # Store the results after applying models
  resultsT <- matrix(0, nrow = k, ncol = 5)
  colnames(resultsT) <- c("Accuracy", "Sensitivity", "Precision", "F-measure","MCC" )
  rownames(resultsT) <- rownames(resultsT, do.NULL = FALSE, prefix = "k = ")
  for (i in 1:k){
      cat("alpha= ", alpha[j],"\n")
      cat("k= ", i,"\n")
      # Perform knn method
      model <- knn.comb(train, trainC, test, testC, train_response, k=i, alpha = alpha[j], beta = beta[j])
      confusion.matrix <- table(factor(test_response), factor(model, levels = c(0,1))) 
      resultsT[i,] <- Assement(confusion.matrix[4], confusion.matrix[3], confusion.matrix[2], confusion.matrix[1])
  }
  totalresultsT <- rbind(totalresultsT, resultsT)
}

# Plot the validation parts to visualize the results
# Set the plots' parameters 
par(mfrow=c(1,1))
opts = c("p","l","o","b","c") 
colors <- rainbow(5) 
linetype <- c(1:5) 
plotchar <- seq(18,18+5,1)

# Plot the results using different combination data
# alpha = 0.1, beta = 0.9
plot(0, 0, xlim = c(0, 20), ylim = c(0,1), type="n", xlab="K", ylab="Percentage" ) 
title(main = "Combination Data (alpha=0.1, beta=0.9)", col.lab="black", font.lab = 4, cex.lab= 1.2)
for (i in 1:5){
  lines(1:20, totalresults[21:40,i], col=colors[i], lwd=2, lty=linetype[i], type="b", pch=plotchar[i])
}
legend("topright", legend = c("A", "S", "P", "F", "MCC"), cex=0.8, col=colors, pch=plotchar, lty=linetype)

# alpha = 0.5, beta = 0.5
plot(0, 0, xlim = c(0, 20), ylim = c(0,1), type="n", xlab="K", ylab="Percentage" ) 
title(main = "Combination Data (alpha=0.5, beta=0.5)", col.lab="black", font.lab = 4, cex.lab= 1.2)
for (i in 1:5){
  lines(1:20, totalresults[101:120,i], col=colors[i], lwd=2, lty=linetype[i], type="b", pch=plotchar[i])
}
legend("topright", legend = c("A", "S", "P", "F", "MCC"), cex=0.8, col=colors, pch=plotchar, lty=linetype)

# alpha = 0.9, beta = 0.1
plot(0, 0, xlim = c(0, 20), ylim = c(0,1), type="n", xlab="K", ylab="Percentage" ) 
title(main = "Combination Data (alpha=0.9, beta=0.1)", col.lab="black", font.lab = 4, cex.lab= 1.2)
for (i in 1:5){
  lines(1:20, totalresults[181:200,i], col=colors[i], lwd=2, lty=linetype[i], type="b", pch=plotchar[i])
}
legend("topright", legend = c("A", "S", "P", "F", "MCC"), cex=0.8, col=colors, pch=plotchar, lty=linetype)

# Plot the results using different combination data given k = 1
plot(0, 0, xlim = c(0, 1), ylim = c(0,1), type="n", xlab="alpha", ylab="Percentage" ) 
title(main = "Combination Data (K = 1)", col.lab="black", font.lab = 4, cex.lab= 1.2)
for (i in 1:5){
  lines(seq(0,1,by=0.1), totalresults[1+20*(0:10),i], col=colors[i], lwd=2, lty=linetype[i], type="b", pch=plotchar[i])
}
legend("topright", legend = c("A", "S", "P", "F", "MCC"), cex=0.8, col=colors, pch=plotchar, lty=linetype)

# Plot the results using different combination data given k = 2
plot(0, 0, xlim = c(0, 1), ylim = c(0,1), type="n", xlab="alpha", ylab="Percentage" ) 
title(main = "Combination Data (K = 2)", col.lab="black", font.lab = 4, cex.lab= 1.2)
for (i in 1:5){
  lines(seq(0,1,by=0.1), totalresults[2+20*(0:10),i], col=colors[i], lwd=2, lty=linetype[i], type="b", pch=plotchar[i])
}
legend("topright", legend = c("A", "S", "P", "F", "MCC"), cex=0.8, col=colors, pch=plotchar, lty=linetype)

# Plot the test parts to visualize the results

# Given k = 1, comparing the results in different combinations
plot(0, 0, xlim = c(0, 1), ylim = c(0,1), type="n", xlab="alpha", ylab="Percentage" ) 
title(main = "Combination Data (K = 1)", col.lab="black", font.lab = 4, cex.lab= 1.2)
for (i in 1:5){
  lines(seq(0,1,by=0.1), totalresultsT[1+2*(0:10),i], col=colors[i], lwd=2, lty=linetype[i], type="b", pch=plotchar[i])
}
legend("topright", legend = c("A", "S", "P", "F", "MCC"), cex=0.8, col=colors, pch=plotchar, lty=linetype)

# Given k = 2, comparing the results in different combinations
plot(0, 0, xlim = c(0, 1), ylim = c(0,1), type="n", xlab="alpha", ylab="Percentage" ) 
title(main = "Combination Data (K = 2)", col.lab="black", font.lab = 4, cex.lab= 1.2)
for (i in 1:5){
  lines(seq(0,1,by=0.1), totalresultsT[2+2*(0:10),i], col=colors[i], lwd=2, lty=linetype[i], type="b", pch=plotchar[i])
}
legend("topright", legend = c("A", "S", "P", "F", "MCC"), cex=0.8, col=colors, pch=plotchar, lty=linetype)

