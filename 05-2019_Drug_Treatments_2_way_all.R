#all plots for the individual drug treatments 48h 
library(gridExtra)
# PCAs

VPA.pca <- prcomp(VPAexp[, -c(5,6,23,37,86,96,97)], center = TRUE, scale. = TRUE)
RA.pca <- prcomp(RAexp[, -c(5,6,23,37,86,96,97)], center = TRUE, scale. = TRUE)
TD.pca <- prcomp(TDexp[, -c(5,6,23,37,86,96,97)], center = TRUE, scale. = TRUE)

TDpcaplot <- fviz_pca_ind(TD.pca, axes = c(1,2), label="none", 
                       geom.ind = "point",
                       habillage = RAexp$Treatment,
                       pointshape = 1, pointsize =1,
                       palette = "Paired",
                       addEllipses = T,
                       title = "PCA: 0.5 uM Chi + 25 mM Thalidomide",ellipse.level = 0.95) + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) 

VPApcaplot <- fviz_pca_ind(VPA.pca, axes = c(1,2), label="none", 
                                              geom.ind = "point",
                                              habillage = VPAexp$Treatment,
                                              pointshape = 1, pointsize =1,
                                              palette = "Paired",
                                              addEllipses = T,
                                              title = "PCA: 0.5 uM Chi + 60 mM Valproic Acid",ellipse.level = 0.95) + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) 

RApcaplot <- fviz_pca_ind(RA.pca, axes = c(1,2), label="none", 
                       geom.ind = "point",
                       habillage = RAforpca$Treatment,
                       pointshape = 1, pointsize =1,
                       palette = "Paired",
                       addEllipses = T,
                       title = "PCA: 0.5 uM Chi + 50 mM Retinoic Acid",ellipse.level = 0.95)+ 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) 
RApcaplot

grid.arrange(TDpcaplot,VPApcaplot,RApcaplot, nrow=1)

######## Clustering ##########

######## Additional Clustering #############
dfTD <- TDexp[-c(5,6,23,37,86,96,97)]
dfVPA <- VPAexp[-c(5,6,23,37,86,96,97)]
dfRA <- RAexp[-c(5,6,23,37,86,96,97)]
dfchi <- ChiPT[ , -c(5,6,23,37,86,96,97,101)]
# As we don’t want the clustering algorithm to depend to an arbitrary variable 
#unit, we start by scaling/standardizing the data using the R function 
dfTD <- scale(dfTD)
dfVPA <- scale(dfVPA)
dfRA <- scale(dfRA)
dfTD <- scale(dfchi)
############## K-means clustering ##############
set.seed(1)

#### Distance Matrix #####
distance1 <- get_dist(dfTD, method="euclidean") # try out different methods here
dp1<- fviz_dist(distance1, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) + 
  ggtitle("Distance Matrix: TD Treatment") + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) 
dp1
  
  
distance2 <- get_dist(dfVPA, method="euclidean") # try out different methods here
dp2 <- fviz_dist(distance2, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
distance3 <- get_dist(dfRA, method="euclidean") # try out different methods here
dp3 <- fviz_dist(distance3, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


####### Kmeans clustering ############
#### Determining optimal number of clusters 
# https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
set.seed(123)
# Elbow method
elbowtd <- fviz_nbclust(dfTD, kmeans, method = "wss")
elbowvpa <- fviz_nbclust(dfVPA, kmeans, method = "wss")
elbowra <- fviz_nbclust(dfRA, kmeans, method = "wss")

grid.arrange(elbowtd,elbowvpa,elbowra, nrow=1)

#Avg silhouette method
siltd <- fviz_nbclust(dfTD, kmeans, method = "silhouette")
silvpa <- fviz_nbclust(dfVPA, kmeans, method = "silhouette")
silra <- fviz_nbclust(dfRA, kmeans, method = "silhouette")

siltd + ggtitle("Optimal number of clusters", subtitle = "0.5 uM Chi + 50 uM Retinoic Acid") +
  theme(plot.title = element_text(size=11, face="bold",margin = margin(0, 0, 7, 0))) +theme(plot.subtitle = element_text(size=9, face="bold",margin = margin(0, 0, 0, 0))) 

grid.arrange(siltd,silvpa,silra, nrow=1)

# Gap statistic method 
gaptd <- fviz_nbclust(dfTD, kmeans, nstart = 25,  method = "gap_stat") + ggtitle("Optimal number of clusters", subtitle = "0.5 uM Chi + 25 mM Thalidomide") +
  theme(plot.title = element_text(size=11, face="bold",margin = margin(0, 0, 5, 0))) +theme(plot.subtitle = element_text(size=9, face="bold",margin = margin(0, 0, 0, 0))) 
gapvpa <- fviz_nbclust(dfVPA, kmeans, nstart = 25,  method = "gap_stat") + ggtitle("Optimal number of clusters",subtitle = "0.5 uM Chi + 60 mM Valproic Acid")+
  theme(plot.title = element_text(size=11, face="bold",margin = margin(0, 0, 5, 0))) +theme(plot.subtitle = element_text(size=9, face="bold",margin = margin(0, 0, 0, 0))) 
gapra <- fviz_nbclust(dfRA, kmeans, nstart = 25,  method = "gap_stat") + ggtitle("Optimal number of clusters",subtitle = "0.5 uM Chi + 50uM Retinoic Acid") +
  theme(plot.title = element_text(size=11, face="bold",margin = margin(0, 0, 5, 0))) +theme(plot.subtitle = element_text(size=9, face="bold",margin = margin(0, 0, 0, 0))) 

grid.arrange(gaptd,gapvpa,gapra, nrow=1)

# NbClust Package to determine the best number of clusters, uses 30 indices to compute and provides the best
# solution 
nbcl <- NbClust::NbClust(dfTD, diss=NULL, distance="euclidean", method="mcquitty", min.nc=2, max.nc=7, index="kl")
nbcl

k2 <- kmeans(dfTD, centers = 2, nstart = 25)
k3 <- kmeans(dfVPA, centers = 2, nstart = 25)
k4 <- kmeans(dfRA, centers = 2, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)
ddk <- cbind(df, cluster = k2$cluster)

# plots to compare
p1 <- fviz_cluster(k2, data = dfTD, geom = "",  palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("K-Means Clusters TD: k = 2") + geom_point(aes(shape = TDexp$Treatment)) +
  labs(shape = "Treatment")
p2 <- fviz_cluster(k3, geom = "",data = dfVPA, palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("K-Means Clusters VPA: k = 2") + geom_point(aes(shape = VPAexp$Treatment)) +
  labs(shape = "Treatment")
p3 <- fviz_cluster(k4, geom = "",data = dfRA, palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("K-Means Clusters RA: k = 2") + geom_point(aes(shape = RAexp$Treatment)) +
  labs(shape = "Treatment")
p4 <- fviz_cluster(k5, geom = "",data = ddk, palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("k = 5") + geom_point(aes(shape = VPAexp$Treatment))+
  labs(shape = "Treatment")

grid.arrange(p1, p2, p3, nrow = 1)

# Add the clusters to the initial data 
drugexpts2 %>%
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

######## PAM clustering - less sensitive to outliers ############

pam.res1 <- pam(dfTD, 2)
pam.res2 <- pam(dfVPA, 2)
pam.res3 <- pam(dfRA, 2)
pam.res4 <- pam(df, 7)
dd <- cbind(df, cluster = pam.res$cluster)

#head(dd, n = 3)
pam.res1$medoids[,38]

pamp1<- fviz_cluster(pam.res1, data = dfTD, geom = "", palette="Paired", ggtheme =theme_minimal()) + ggtitle("PAM TD: k=2") + geom_point(aes(shape = VPAexp$Treatment)) +
  labs(shape = 'Treatment')
pamp2<- fviz_cluster(pam.res2, data = dfVPA, geom = "", palette="Paired", ggtheme =theme_minimal()) + ggtitle("PAM VPA: k=2") + geom_point(aes(shape = VPAexp$Treatment)) +
  labs(shape = 'Treatment')
pamp3<- fviz_cluster(pam.res3, data = dfRA, geom = "", palette="Paired", ggtheme =theme_minimal()) + ggtitle("PAM RA: k=2") + geom_point(aes(shape = RAexp$Treatment)) +
  labs(shape = 'Treatment')
pamp1

grid.arrange(pamp1,pamp2, pamp3, nrow = 1)

######### Compare them ###########

######## Predict ################
library(caret)
TDexp$Treatment <- as.numeric(TDexp$Treatment)
TDexp$Treatment <- as.factor(TDexp$Treatment)
testClust <- data.frame(cbind(TD.pca$x,TDexp$Treatment))

colnames(testClust)[60] <- "Treatment"

pcaColNumber <- 16
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

fit <- train(treatment ~., data = train, method = "knn", metric = metric, 
             preProc = c("center", "scale"), trControl = control)

fit$results
preds <- predict(fit$finalModel, test[,1:pcaColNumber]) # do the pca on the things you are training 
#the model with, and then apply the PCA on the ones you want to test 

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

toPlot <- data.frame(PC1 = test$PC1, PC2 = test$PC2, 
                     prediction = predsClass, treatment = test$treatment)

gpred <- ggplot(toPlot, aes(x = PC1, y = PC2, group = as.factor(keepOrNot))) + 
  geom_point(aes(shape = as.factor(treatment), colour = as.factor(keepOrNot))) + 
  ggtitle("Predicted Groups vs. Treatment Condition, TD-treated hGlds")
#stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = as.factor(keepOrNot)), level = 0.95) 

gpred + theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) + 
  labs(shape = "Treatment", fill = "Prediction") + 
  scale_color_discrete(name = "Prediction", labels = c("NA", "0.5 Chi TD", "0.5 Chi control"))+
  scale_shape_discrete(labels=c("0.5 Chi TD", "0.5 Chi control"))


xtab <- confusionMatrix(as.factor(predsClass), test$treatment)
xtab
xtabmat <- as.matrix(xtab, "classes")
ddd <- as.matrix(xtab, "overall")
rapredfeatures <- rbind(xtabmat, ddd)

write.csv(rapredfeatures, file = "/Users/alexandrabaranowski/Desktop/Images for report /Drug Expt/TD Prediction Features for test v2.csv", row.names = TRUE)

predsClass
test$treatment

####### Prediction using Optimal Weighted Nearest Neighbour Classifier ##########
library(FNN)
TDexp$Treatment <- as.numeric(TDexp$Treatment)
TDexp$Treatment <- as.factor(TDexp$Treatment)
testClust <- data.frame(cbind(TD.pca$x,TDexp$Treatment))

colnames(testClust)[60] <- "Treatment"

pcaColNumber <- 16
pcacols <- testClust[,1:pcaColNumber]

train1 <- cbind(pcacols, treatment = as.factor(testClust$Treatment))

rows <- sample(nrow(train1))
train1 <- train1[rows, ]
split <- round(nrow(train1) * .60)
train <- train1[1:split, ]
test <- train1[(split + 1):nrow(train1), ]

cl <- factor(c(rep("TD",), rep("Control",25), rep("v",25)))
testcl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
out <- ownn(train, test, cl, testcl)
out
