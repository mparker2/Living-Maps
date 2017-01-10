
# Create subset to split training points for training and testing
random.subset <- function(df, frac.training) {
  frac <- floor((nrow(df) * frac.training))
  subset <- sample(nrow(df), size=frac)
  return (subset)
}

# Produces user / producer accuracy
# depending on whether rowSums or colSums is passed
user.producer.accuracy <- function(i, confusion.table, f=rowSums) {
  tots <- f(confusion.table)
  # replace NaNs with 0 if there is no data for the row/col
  if (tots[i] == 0) {
    return(0.0)
  }
  else {
    return (confusion.table[i,i]/tots[i] * 100)
  }
}

# produce confusion matrix and user / produce accuracies
confusion.matrix <- function(user, producer) {
  user <- as.factor(user)
  producer <- as.factor(producer)
  user <- as.factor(as.character(user))
  levels(producer) <- levels(user)
  # if user and producer are different lengths, something is wrong!
  stopifnot(length(user) == length(producer))
  
  n <- length(levels(user))
  confusion.table <- table(producer, user)
  
  range <- as.matrix(1:n)
  # total accuracy is the sum of the diagonal divided by the sum of the whole
  total.accuracy <- sum(diag(confusion.table)) / sum(confusion.table) * 100
  user.accuracy <- apply(range, MARGIN=1, FUN=user.producer.accuracy,
                         confusion.table=confusion.table, f=colSums)
  producer.accuracy <- apply(range, MARGIN=1, FUN=user.producer.accuracy,
                             confusion.table=confusion.table, f=rowSums)
  names(user.accuracy) <- colnames(confusion.table)
  names(producer.accuracy) <- rownames(confusion.table)
  
  return (list(table=confusion.table,
               user.accuracy=user.accuracy,
               producer.accuracy=producer.accuracy,
               total.accuracy=total.accuracy))
}
# 
# library(randomForest)
# library(e1071)
# library(monmlp)
# 
# setwd('C:\\Users\\Matthew\\Documents\\GitHub\\Living-Maps\\Living_Maps_Segmentation_Dartmoor')
# source('../glmulti2.R')
# 
# training.data <- read.table("training_data.txt", sep="\t", header=T)
# # Use only Tier 1 data and complete cases, as svm / random forest cannot deal with missing values
# training.data.clean <- training.data[complete.cases(training.data) & training.data$Tier == 1,
#                                      c(3, 5:(ncol(training.data) - 3))]
# 
# # Use 80% of data for training and 20% for testing
# subset <- random.subset(training.data.clean, 0.8)
# training.data.input <- training.data.clean[subset,]
# training.data.test <- training.data.clean[-subset,]
# 
# M.glm <- glmulti2("Feature_Ty ~", data=training.data.test,
#                   variables=colnames(training.data.test)[2:ncol(training.data.test)])
# p.glm <- predict(M.glm, newdata=training.data.test)
# confusion.matrix(training.data.test$Feature_Ty, p.glm)
# 
# M.rf <- randomForest(Feature_Ty ~ ., data=training.data.input, ntree=200)
# p.rf <- predict(M.rf, newdata=training.data.test)
# confusion.matrix(training.data.test$Feature_Ty, p.rf)
# 
# M.svm <- svm(Feature_Ty ~ ., data=training.data.input)
# p.svm <- predict(M.svm, newdata=training.data.test)
# confusion.matrix(training.data.test$Feature_Ty, p.svm)
