#### Predicting the class based on a 60% training set #######

library(mlbench)
library(caret)
testClust <- data.frame(cbind(RA.pca$x,RAexp$Treatment))
colnames(testClust)[60] <- "Treatment"

pcaColNumber <- 15
pcacols <- testClust[,1:pcaColNumber]

train1 <- cbind(pcacols, treatment = as.factor(testClust$Treatment))
set.seed(34)
rows <- sample(nrow(train1))
train1 <- train1[rows, ]
split <- round(nrow(train1) * .60)
train <- train1[1:split, ]
test <- train1[(split + 1):nrow(train1), ]

metric = "Accuracy"
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = T)
set.seed(33)
fit <- train(treatment ~., data = train, method = "knn", metric = metric, 
             preProc = c("center", "scale"), trControl = control)
fit$results
preds <- predict(fit$finalModel, train1[,1:pcaColNumber]) 
# do the pca on the things you are training  the model with, and then apply the PCA on the ones you want to test

### This is the part that only keeps the predicted samples that have been predicted with > 0.8 (or whatever
# value is chosen condifence; this can also be removed but has been added for more "safety" of the 
#classes that are predicted)
keepOrNot <- rep(0,1,nrow(preds))
for(i in 1:nrow(preds)){
  for(j in 1:ncol(preds)){
    if(preds[i,j] >=0.8){
      keepOrNot[i] <- j
    }
  }
}
cbind(preds,keepOrNot)

predsClass <- apply(preds, 1, which.max) # for each row find which is the max, and gives you the column number
#colnames(preds)[predsClass]

toPlot <- data.frame(PC1 = train1$PC1, PC2 = train1$PC2, 
                     prediction = predsClass, treatment = train1$treatment)

gpred <- ggplot(toPlot, aes(x = PC1, y = PC2, group = as.factor(keepOrNot))) + 
  geom_point(aes(shape = as.factor(treatment), colour = as.factor(keepOrNot))) + 
  ggtitle("Predicted Groups vs. Treatment Condition, RA-treated hGlds")
#stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = as.factor(keepOrNot)), level = 0.95) 

gpred +theme_minimal()+ theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) + 
  labs(shape = "Treatment", fill = "Prediction") + 
  scale_color_discrete(name = "Prediction", labels = c("NA", "0.5 Chi RA","0.5 Chi control"))+
  scale_shape_discrete(labels=c( "0.5 Chi RA","0.5 Chi control"))


accuracy <- data.frame(fit$resample) 

##### Additional plots and exploration ########
#### Look at which PCs may be important 
featurePlot(x = train[, 76:85],
y = train$treatment,
plot = "box",
strip=strip.custom(par.strip.text=list(cex=.6)),
scales = list(x = list(relation="free"),
y = list(relation="free")))

### 
featurePlot(x = train[, 76:85], 
            y = train$treatment, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

############## Plot the accuracies of different k values ##############
ggplot(fit) + xlab("# Neighbours")+
  ggtitle("Averaged Model Accuracies with Knn Classifier") + theme_minimal() + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0)))
  
######### Get the confusion matrix and different model metrics #################
xtab <- confusionMatrix(as.factor(predsClass), train1$treatment)
xtab
xtabmat <- as.matrix(xtab, "classes")
xtabmat
predsClass
test$treatment
