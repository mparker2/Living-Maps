library(reshape2)
library(ggplot2)

# Create subset to split training points for training and testing
random.subset <- function(df, frac.training, seed=NA) 
{
   if (!is.na(seed))
   {
      set.seed(1001)    
   }
  
   frac <- floor((nrow(df) * frac.training))
   subset <- sample(nrow(df), size=frac)
   return (subset)
}

# produce confusion matrix and user / produce accuracies
confusion.matrix <- function(user, producer) {
  user <- as.factor(user)
  producer <- as.factor(producer)
  # if user and producer are different lengths, something is wrong!
  stopifnot(length(user) == length(producer))
  
  n <- length(levels(user))
  confusion.table <- table(producer, user)
  
  range <- as.matrix(1:n)
  # total accuracy is the sum of the diagonal divided by the sum of the whole
  total.accuracy <- sum(diag(confusion.table)) / sum(confusion.table) * 100
  #user accuracy is the diagonal divided by colSums
  user.accuracy <- diag(confusion.table) / colSums(confusion.table) * 100
  #producer accuracy is the diagonal divided by rowSums
  producer.accuracy <- diag(confusion.table) / rowSums(confusion.table) * 100

  names(user.accuracy) <- colnames(confusion.table)
  names(producer.accuracy) <- rownames(confusion.table)
  
  return (list(table=confusion.table,
               user.accuracy=user.accuracy,
               producer.accuracy=producer.accuracy,
               total.accuracy=total.accuracy))
}

broadclass.confusion.matrix <- function(user, broaduser, producer) {
  user <- as.factor(user)
  broaduser <- as.factor(broaduser)
  producer <- as.factor(producer)
  # if user and producer are different lengths, something is wrong!
  stopifnot(length(user) == length(producer))
  stopifnot(length(user) == length(broaduser))
  
  # make a vector mapping subclasses to broad classes
  subclasses <- levels(user)
  n.subclass <- length(subclasses)
  broadclasses <- character(n.subclass)
  for (i in 1:n.subclass) {
    broadclasses[i] <- as.character(
      broaduser[user == subclasses[i]][1]
    )
  }
  names(broadclasses) <- subclasses
  
  # get predictions at broader class level
  broadproducer <- factor(
    unlist(lapply(producer,
                  FUN=function(x) {broadclasses[[x]]})),
    levels=levels(broaduser)
  )
  return (confusion.matrix(broaduser, broadproducer))
}

barplot.confusion.matrix <- function(confusion, plot_user=TRUE) {
  cm <- melt(confusion$table)
  colnames(cm) <- c('User', 'Producer', 'Count')
  if (plot_user) {
    xticklabs <- paste(names(confusion$user.accuracy),
                       sprintf('%.1f%%', confusion$user.accuracy))
    xvar <- 'Producer'
    fillvar <- 'User'
  }
  else {
    xticklabs <- paste(names(confusion$producer.accuracy),
                       sprintf('%.1f%%', confusion$producer.accuracy))
    xvar <- 'User'
    fillvar <- 'Producer'
  }
  
  g <- ggplot(data=cm, aes_string(x=xvar, y='Count', fill=fillvar)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels=xticklabs) +
    ggtitle(sprintf('Total Accuracy - %.2f%%', confusion$total.accuracy))
  return(g)
}