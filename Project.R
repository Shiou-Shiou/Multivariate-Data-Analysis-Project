rm(list = ls())
graphics.off()

# set working direction
setwd("C:/Users/daughter/Documents/Classes_SJSU/Math 257/Project")

# Install the packages #####

library(andrews) # for plot andrews curves
library(ggplot2)
library(GGally) # for parallel plot
library(rgl) # for 3D plot
library(car) # for MANONA
library(rpart) # for decision tree
library(rpart.plot) # for tree plot
library(randomForest) # for random forest
library(e1071) # for SVM

# Function made for the project##########
# Make a assement function to evaluate the results 
Assement <- function(TP, FP, FN, TN){
  Accuracy = (TP+TN)/ (TP+TN+FP+FN)
  Sensitivity = TP / (TP+FN) 
  Precision = TP/(TP+FP)
  F.measure = 2*TP/(2*TP+FN+FP)
  return(c(Accuracy, Sensitivity, Precision, F.measure))
}


# Data Preprocessing #####
# Load the data 
breast.cancer <- read.csv("breast_cancer.csv", header = TRUE)
# Delete personal imformation, such as "ID"
breast.cancer <- breast.cancer[,-1]


# Generate the numbers for 10 folds
set.seed(12345)
folds_i <- sample(rep(1:10, length.out = nrow(breast.cancer)))
# Store the 10 folder results
results <- matrix(0, nrow = 8, ncol = 4)
colnames(results) <- c("Accuracy", "Sensitivity", "Precision", "F-measure")

for (a in 1:10){
  valid_i <- which(folds_i == a)
  # Divide the data into 2 part (90% training and validation, 10% testing)
  train <- breast.cancer[-valid_i, ]
  test <- breast.cancer[valid_i,]
  # Scale the data exclude the label
  train.s <- data.frame(scale(train[,-1], center = T, scale= T))
  train.s <- cbind(train$diagnosis, train.s)
  colnames(train.s)[1] <- c("diagnosis")
  test.s <- data.frame(scale(test[,-1], center = T, scale = T))
  test.s <- cbind(test$diagnosis, test.s)
  colnames(test.s)[1] <- c("diagnosis")
  if (a == 1){
    summary(train)
    summary(train.s)
    # Display the data ######
    
    ######## Andrew plot #####
    andrews(train.s[,-1], type = 1, ymax = 13, main = "Andrews Plot for the breast cancer Data") # no label
    andrews(train.s, type = 1, clr = 1, ymax = 13, main = "Andrews Plot for the breast cancer Data") # with labels
    ######## Parallel corrdinate plot #####
    p1 <- ggparcoord(data = train, columns = 2:31, groupColumn = 1, showPoints = TRUE, alphaLines = 0.3) + theme(axis.title.x=element_blank(),
                                                                                                               axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle("Parallel Coordinate Plot") + theme(plot.title = element_text(hjust = 0.5))
    p1
    
    p2 <- ggparcoord(data = train.s, columns = 2:31, groupColumn = 1, showPoints = TRUE, alphaLines = 0.3) + theme(axis.title.x=element_blank(),
                                                                                                                   axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle("Parallel Coordinate Plot") + theme(plot.title = element_text(hjust = 0.5))
    p2
    ######## SVD plot #####
    svd <- svd(as.matrix(train[,-1]))
    plot3d(svd$u[, 1], svd$u[, 2], svd$u[, 3], main = "SVD", col = as.integer(train$diagnosis))
    svd.s <- svd(as.matrix(train.s[,-1]))
    plot3d(svd.s$u[, 1], svd.s$u[, 2], svd.s$u[, 3], main = "SVD", col = as.integer(train.s$diagnosis))
    ######## scatter plot ##### 
    # define the function to calculate the correlation between variables
    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- abs(cor(x, y))
      txt <- format(c(r, 0.123456789), digits = digits)[1]
      txt <- paste0(prefix, txt)
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex = cex.cor * r *1.4)
    }
    # define the function to draw the histograms
    panel.hist <- function(x, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y,col = "cyan", ...)
    }
    # unstandarized
    pairs(train[2:11], pch=21, bg = c("red", "blue")[unclass(train$diagnosis)],
          lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist,
          main="Breast Cancer Scatterplot Matrix for 1-10 features")
    pairs(train[12:21], pch=21, bg = c("red", "blue")[unclass(train$diagnosis)],
          lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist,
          main="Breast Cancer Scatterplot Matrix for 11-20 features")
    pairs(train[22:31], pch=21, bg = c("red", "blue")[unclass(train$diagnosis)],
          lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist,
          main="Breast Cancer Scatterplot Matrix for 21-30 features")
    # standarized
    pairs(train.s[2:11], pch=21, bg = c("red", "blue")[unclass(train.s$diagnosis)],
          lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist,
          main="Breast Cancer Scatterplot Matrix for 1-10 features (Standarized)")
    pairs(train.s[12:21], pch=21, bg = c("red", "blue")[unclass(train.s$diagnosis)],
          lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist,
          main="Breast Cancer Scatterplot Matrix for 11-20 features (Standarized)")
    pairs(train.s[22:31], pch=21, bg = c("red", "blue")[unclass(train.s$diagnosis)],
          lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist,
          main="Breast Cancer Scatterplot Matrix for 21-30 features (Standarized)")
    
    
    # Detect the potential outliers ########
    
    # Finding the statistical distances 
    trains.mat <- as.matrix(train.s[,-1])
    n <- nrow(trains.mat)
    p <- ncol(trains.mat)
    xbar <- colMeans(trains.mat)
    Si <- solve(var(trains.mat))
    # Then, calculate all Mahalanobis distances
    distance <- sapply(1:n, function(k) (t(trains.mat[k,] - xbar)%*% Si %*%(trains.mat[k,] - xbar)))
    # Compute quantiles of a chi-square distribution
    q1 <- qchisq((1:n-0.5)/n,p)
    #  Create the chi-squared probability plot
    qqplot(q1,distance,xlab = "Chi-square quantiles", ylab = "Sample statistical distances", main = "Chi-square Probability Plot")
    lines(q1,q1)
    # Delete the 21 outliers
    #q2 <- qchisq((1:(n-21)-0.5)/(n-21),p)
    #qqplot(q2,d[d<100],xlab = "Chi-square quantiles", ylab = "Sample statistical distances", main = "Chi-square Probability Plot")
    #lines(q2,q2)
    
    index <- distance > 100
    observations <- 1:nrow(train.s)
    outliers <- observations[index]
    cat("Outliers are", outliers,"\n")

    
  }
  
  # Compare the mean vectors ########
  # Compute sample mean vector and sample covariance matrix.
  xbar1 <- colMeans(train.s[train.s$diagnosis == "B" ,2:31])
  xbar2 <- colMeans(train.s[train.s$diagnosis== "M" ,2:31])
  S1 <- var(train.s[train.s$diagnosis == "B" ,2:31])
  S2 <- var(train.s[train.s$diagnosis== "M" ,2:31])
  difference <- xbar1-xbar2
  ## Barlett's test for equality of covariance matrices
  k <- 2
  n1 = nrow(train.s[train.s$diagnosis=="B",])
  n2 = nrow(train.s[train.s$diagnosis=="M",])
  Spool <- S1*(n1-1)/(n1+n2-2)+S2*(n2-1)/(n1+n2-2)
  # Calculating our test statistic
  (M <- (n1-1)*log(det(Spool))+(n2-1)*log(det(Spool))-(n1-1)*log(det(S1))-(n2-1)*log(det(S2)))
  (Cinv <- 1-((2*p^2+3*p-1)/(6*(p+1)*(k-1)))*(1/(n1-1)+1/(n2-1)-1/(n1+n2-2)))
  df <- 1/2*(k-1)*(p+1)*p
  # Comparing to a chisquare
  Barlett.result <- M*Cinv > qchisq(0.05,df,lower.tail = FALSE)
  cat("Rejecting the hypothesis of equal variance is", Barlett.result,"\n")
  ## Testing the difference between two means, different variances/large-sample approximation
  mu0 <- numeric(30)
  Sj <- 1/n1*S1 + 1/n2*S2
  # Use pool variance
  T2 <- t(difference - mu0) %*% solve((1/n1+1/n2)*Sj) %*% (difference - mu0)
  # Compare to critical value
  critical.value <- (n-2)*p/(n-p-1)*qf(0.05,30, n-p-1,lower.tail=FALSE)
  cat("Rejecting the hypothesis of equal mean vectors between groups is", T2 > critical.value,"\n")
  
  # Feature Selection ##### 
  ######## Method 1 : MANOVA ####
  
  # Reconstruct the data with the standarized data
  features <- c("radius","texture","perimeter","area","smoothness","compactness","concavity","concave points","symmetry","fractal dimension")
  group <- as.factor(rep(features, each = nrow(train.s)))
  
  group.data <- train.s[,c(2,10+2, 20+2)]
  colnames(group.data) <- c("mean", "se","worst")
  for (i in 3:11){
    sub.data <- train.s[,c(i,10+i, 20+i)]
    colnames(sub.data) <- c("mean", "se","worst")
    group.data <- rbind(group.data, sub.data)
  }
  group.data <- data.frame(cbind(group, group.data))
  # Run MANOVA #
  summary(group.data)
  
  # Specify the linear relationship in MANOVA
  fit.lm <- lm(cbind(mean, se, worst)~group, data = group.data) # we need to do this because Manova takes as its input a model from lm, glm, or multinom
  # Run the Manova
  fit.manova <- Manova(fit.lm)
  # See results for each of the tests
  summary(fit.manova)
  
  ######## Method 2 : Using CIs ######
  ## Check the CIs of each features
  cis <- diag(30)
  CIs.2groups <- matrix(NA, nrow = 30, ncol = 2)
  colnames(CIs.2groups) <- c("Lower","Upper")
  for (i in 1:30){
    CIs.2groups[i,] <- round(c(t(cis[i,])%*%difference - sqrt(qchisq(0.05,2,lower.tail = FALSE))*sqrt(t(cis[i,])%*%Sj%*%cis[i,]),t(cis[i,])%*%difference + sqrt(qchisq(0.05,2,lower.tail = FALSE))*sqrt(t(cis[i,])%*%Sj%*%cis[i,])),2)
    
  }
  # Check the CIs do not contain 0
  ind <- !(CIs.2groups[,1] < 0 & CIs.2groups[,2] > 0)
  
  # Find the new data through CIs
  train.data2 <- data.frame(diagnosis = train.s$diagnosis, train.s[,c(FALSE,ind)])
  #transform test by same dimensions
  test.data2 <- test.s[,c(FALSE,ind)]
  test.data2 <- as.data.frame(test.data2)
  
  
  ######## Method 3 : PCA #####
  
  # Find the pricipal components
  trains.pc <- prcomp(train.s[,-1])
  cumvar <- cumsum(trains.pc$sdev^2)/sum(trains.pc$sdev^2)
  propvar <- trains.pc$sdev^2/sum(trains.pc$sdev^2)
  if(a==1){
    # Plot proportion of variance explained by each component & cumulative proportion of variance
    plot(propvar[1:15], ylim=c(0,1), xaxt = "n", main = "Proportion of variance explained by each PC", xlab = "principal component", ylab = "Proportion of variance explained", pch = 16, bty = "n")
    axis(1, at = c(1:15), labels = c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]), expression(lambda[4]), expression(lambda[5]), 
                                     expression(lambda[6]), expression(lambda[7]), expression(lambda[8]), expression(lambda[9]), expression(lambda[10]), 
                                     expression(lambda[11]), expression(lambda[12]), expression(lambda[13]), expression(lambda[14]), expression(lambda[15])))
    lines(propvar[1:15], lty = 19)
    lines(cumvar[1:15], lty = 20)
    points(cumvar[1:15], pch = 21, bg = "white")
    legend(7, 0.6, legend=c("Proportion of variance", "Cumulative variance"), lty = c(19, 20), pch = c(16, 21), pt.bg = c(NA,"white"), cex = 0.8, bty = "n")
    #  plot component scores
    par(pch=5, fin=c(3,3))
    pairs(trains.pc$x[,c(1,2,3)],labels=c("PC1","PC2","PC3"), col = as.integer(train.s$diagnosis))
    plot3d(trains.pc$x[, 1], trains.pc$x[, 2], trains.pc$x[, 3], main = "PC", col = as.integer(train.s$diagnosis))
    
  }
  
  # Find the new data through PCA
  train.data3 <- data.frame(diagnosis = train.s$diagnosis, trains.pc$x)
  #we are interested in those PCAs which have total variance smaller than 96%
  num <- length(cumvar[cumvar < 0.96])
  train.data3 <- train.data3[,1:(num+1)]
  #transform test into PCA
  test.data3 <- predict(trains.pc, newdata = test.s[,-1])
  test.data3 <- as.data.frame(test.data3)
  #select the first 10 components
  test.data3 <- test.data3[,1:num]
  
  
  # Classification ###########
  
  
  ######## Method 1 : Decision Tree ########
  
  # Use Original
  
  #run a decision tree
  set.seed(12345)
  rpart.model <- rpart(diagnosis ~ .,data = train.data2, method = "class")
  if (a == 1){
    printcp(rpart.model) # display the results 
    plotcp(rpart.model) # visualize cross-validation results 
    summary(rpart.model) # detailed summary of splits
    # plot tree 
    plot(rpart.model, uniform=TRUE, 
         main="Classification Tree for Breast Cancer")
    text(rpart.model, use.n=TRUE, all=TRUE, cex=0.8)
    rpart.plot(rpart.model)
    # prune the tree 
    pfit<- prune(rpart.model, cp=rpart.model$cptable[which.min(rpart.model$cptable[,"xerror"]),"CP"])
    # plot the pruned tree 
    plot(pfit, uniform=TRUE, 
         main="Pruned Classification Tree for Breast Cancer (pruned)")
    text(pfit, use.n=TRUE, all=TRUE, cex=.8)
    rpart.plot(pfit)
  }
  
  #make prediction on test data
  rpart.prediction <- predict(rpart.model, test.data2, type = "class")
  confusion.map <- table(test.s$diagnosis, rpart.prediction)
  # Store in odd rows
  results[1,] <- results[1,] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
  # Use PCA
  
  #run a decision tree
  set.seed(12345)
  rpart.model <- rpart(diagnosis ~ .,data = train.data3, method = "class")
  if(a==1){
    printcp(rpart.model) # display the results 
    plotcp(rpart.model) # visualize cross-validation results 
    summary(rpart.model) # detailed summary of splits
    # plot tree 
    plot(rpart.model, uniform=TRUE, 
         main="Classification Tree for Breast Cancer")
    text(rpart.model, use.n=TRUE, all=TRUE, cex=0.8)
    rpart.plot(rpart.model)
    
    # prune the tree 
    pfit<- prune(rpart.model, cp=rpart.model$cptable[which.min(rpart.model$cptable[,"xerror"]),"CP"])
    # plot the pruned tree 
    plot(pfit, uniform=TRUE, 
         main="Pruned Classification Tree for Breast Cancer (pruned)")
    text(pfit, use.n=TRUE, all=TRUE, cex=.8)
    rpart.plot(pfit)
    
  }
  
  #make prediction on test data
  rpart.prediction <- predict(rpart.model, test.data3, type = "class")
  confusion.map <- table(test.s$diagnosis, rpart.prediction)
  results[2,] <- results[2,] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
  ######## Method 2 : Random Forest #########
  
  # Original 
  
  # Run Random Forest 
  set.seed(12345)
  RF.model <- randomForest(diagnosis ~ ., data=train.data2)
  
  if(a==1){
    print(RF.model) # view results 
    importance(RF.model) # importance of each predictor
    # Plot the error of tree
    layout(matrix(c(1,2),nrow=1),width=c(4,1)) 
    par(mar=c(5,4,4,0)) #No margin on the right side
    plot(RF.model, main ="Random Forest Model")
    par(mar=c(5,0,4,2)) #No margin on the left side
    plot(c(0,1),type="n", axes=F, xlab="", ylab="")
    legend("top", colnames(RF.model$err.rate),col=1:4,cex=0.8,fill=1:4)
    layout(matrix(c(1,1)))
    par(mar=c(5,4,4,4))
    varImpPlot(RF.model)
    
  }
  
  #make prediction on test data
  RF.prediction <- predict(RF.model, newdata = test.data2)
  confusion.map <- table(test.s$diagnosis, RF.prediction)
  results[3,] <- results[3,] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
  
  # PCA
  
  # Run Random Forest 
  set.seed(12345)
  RF.model <- randomForest(diagnosis ~ ., data=train.data3)
  if(a==1){
    print(RF.model) # view results 
    importance(RF.model) # importance of each predictor
    # Plot the error of tree
    layout(matrix(c(1,2),nrow=1),width=c(4,1)) 
    par(mar=c(5,4,4,0)) #No margin on the right side
    plot(RF.model, main ="Random Forest Model(PCA)")
    par(mar=c(5,0,4,2)) #No margin on the left side
    plot(c(0,1),type="n", axes=F, xlab="", ylab="")
    legend("top", colnames(RF.model$err.rate),col=1:4,cex=0.8,fill=1:4)
    layout(matrix(c(1,1)))
    par(mar=c(5,4,4,4))
    varImpPlot(RF.model)
    
  }
  
  #make prediction on test data
  RF.prediction <- predict(RF.model, newdata = test.data3)
  confusion.map <- table(test.s$diagnosis, RF.prediction)
  results[4, ]  <- results[4, ] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
  ######## Method 3 : Logistic Regression (Not finished) #####
  
  # Original 
  
  #logistic regression model
  set.seed(12345)
  LR.model <- glm(diagnosis~ .-diagnosis,family=binomial(link='logit'), data = train.data2)
  if(a==1){
    table(train.data2$diagnosis,predict(LR.model,type='response')>=0.5)
    summary(LR.model)
    anova(LR.model, test="Chisq")
    
    slope <- coef(LR.model)[2]/(-coef(LR.model)[3])
    intercept <- coef(LR.model)[1]/(-coef(LR.model)[3]) 
    
    library(lattice) # for decision boundary
    xyplot( texture_mean  ~ radius_mean  , data = train.data2, groups = diagnosis,
            panel=function(...){
              panel.xyplot(...)
              panel.abline(intercept , slope)
              panel.grid(...)
            })
    
  }
  # apply to test data
  fitted.probs <- predict(LR.model,test.data2,type='response')
  LR.prediction <- ifelse(fitted.probs > 0.5,"M","B")
  confusion.map <- table(test.s$diagnosis, LR.prediction)
  results[5,] <- results[5,] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
  # PCA
  
  #logistic regression model
  set.seed(12345)
  LR.model <- glm(diagnosis~ .-diagnosis,family=binomial(link='logit'), data = train.data3)
  if(a==1){
    table(train.data3$diagnosis,predict(LR.model,type='response')>=0.5)
    summary(LR.model)
    anova(LR.model, test="Chisq")
    
    slope <- coef(LR.model)[2]/(-coef(LR.model)[3])
    intercept <- coef(LR.model)[1]/(-coef(LR.model)[3]) 
    library(lattice) # for decision boundary
    xyplot( PC2 ~ PC1 , data = train.data3, groups = diagnosis,
            panel=function(...){
              panel.xyplot(...)
              panel.abline(intercept , slope)
              panel.grid(...)
            })
  }
  
  # apply to test data
  fitted.probs <- predict(LR.model,test.data3,type='response')
  LR.prediction <- ifelse(fitted.probs > 0.5,"M","B")
  confusion.map <- table(test.s$diagnosis, LR.prediction)
  results[6,] <- results[6,] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
  ######## Method 4 : Surport Vector Machine #######
  
  # Original
  
  # Build the model
  SVM.model <- svm(diagnosis ~ .-diagnosis, data=train.data2)
  summary(SVM.model)
  if(a==1){
    # plot the model
    plot(SVM.model, train.data2, perimeter_worst ~ concave.points_worst)
  }
  
  # Predict the label
  SVM.prediction <- predict(SVM.model, newdata = test.data2)
  confusion.map <- table(test.s$diagnosis, SVM.prediction)
  results[7,] <- results[7,] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
  # PCA
  
  # Build the model
  SVM.model <- svm(diagnosis ~ .-diagnosis, data=train.data3)
  summary(SVM.model)
  if(a==1){
    # plot the model
    plot(SVM.model, train.data3, PC1 ~ PC2)
  }

  # Predict the label
  SVM.prediction <- predict(SVM.model, newdata = test.data3)
  confusion.map <- table(test.s$diagnosis, SVM.prediction)
  results[8,] <- results[8,] + Assement(confusion.map[1], confusion.map[3], confusion.map[2], confusion.map[4])
  
}
final.result <- results/10 


