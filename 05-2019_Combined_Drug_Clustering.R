library(gridExtra)
library(cluster)
library(tidyverse)
###### Clustering Analysis #######
testClust <- data.frame(cbind(drug2.pca$x,drugexpts2$Treatment))
colnames(testClust)[90] <- "Treatment"

unique(testClust$Treatment)
ttt <- testClust[,1:16]

#### Determining optimal number of clusters 
set.seed(123)
# Elbow method
fviz_nbclust(ttt, kmeans, method = "wss")

#Avg silhouette method
fviz_nbclust(ttt, kmeans, method = "silhouette")

# Gap statistic method 
fviz_nbclust(ttt, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)

nbcl <- NbClust::NbClust(ttt, diss=NULL, distance="canberra", method="ward.D2", min.nc=2, max.nc=7, index="all")
## Plot a histogram showing what best # clusters all the different methods gave 
bestnclust <- data.frame(nbcl$Best.nc)
bestnclust <- t(bestnclust)
ggplot(bestnclust, aes(x=Number_clusters)) + geom_histogram(binwidth=1, color="cadetblue", fill = "lightcyan3") +
  theme_minimal() + xlab("Calculated Optimal Number of Clusters") + ylab("Count") +
  ggtitle("Optimal Number of Clusters as determined by 23 Separate Methods") + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(0, 0, 10, 0)))

distance <- get_dist(ttt, method="spearman") # try out different methods here
fviz_dist(distance, gradient = list(low = "#FC4E07", mid = "white", high = "#00AFBB"))
# low = low dissimilarity, high = high dissimilarity
#################################

distm <- dist(ttt, method = "canberra")
hclustttt <- hclust(distm, method = "ward.D2")

#hclustttt <- hclust(distm, method = "complete")
plot(hclustttt, xlab="Distance \nDistance calculation: Canberra, Clustering Method: Complete")
rect.hclust(hclustttt, k =4, border = 2:6)

cuttt <- cutree(hclustttt, k = 4)

afterclust <- mutate(ttt, cluster = cuttt)
#count(seeds_df_cl,cluster)

afterclust <- cbind(afterclust, treatment = testClust$Treatment)
#View(afterclust[,11:12])

g1 <- ggplot(afterclust, aes(x = PC1, y = PC2, group = as.factor(cluster))) + 
  geom_point(aes(color = as.factor(treatment))) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = as.factor(cluster)), level = 0.95) +
  ggtitle("Clusters produced from first 16 Principal Components, \nCanberra and Ward.D2")

g1+ theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) + 
  labs(color = "Treatment", fill = "Cluster") + scale_color_discrete(labels = c("TD", "RA", 
                                                                                "VPA", "Control"))

plot3d(afterclust,col=as.factor(cuttt), shape=drugexpts2[,"Treatment"], type ='p')

ggplot(afterclust, aes(x = PC1, y = PC2, group = as.factor(cluster))) + 
  geom_point(aes(shape = as.factor(Rallh$Pre.treatment), color = as.factor(Rallh$Pre.treatment))) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = as.factor(cluster)), level = 0.95)

hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")
mycol <- colorpanel(40, "darkblue", "yellow", "white")
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(mycl))

######## Additional Clustering #############
df <- drugexpts2[-c(5,6,23,37,86,96,97)]
# As we don’t want the clustering algorithm to depend to an arbitrary variable 
#unit, we start by scaling/standardizing the data using the R function 
df <- scale(df)
############## K-means clustering ##############
set.seed(1)



k2 <- kmeans(df, centers = 2, nstart = 25)
k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 7, nstart = 25)
ddk <- cbind(df, cluster = k2$cluster)

# plots to compare
p1 <- fviz_cluster(k2, data = ddk, geom = "",  palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("k = 2") + geom_point(aes(shape = drugexpts2$Treatment))
p2 <- fviz_cluster(k3, geom = "",data = ddk, palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("k = 3") + geom_point(aes(shape = drugexpts2$Treatment)) 
p3 <- fviz_cluster(k4, geom = "",data = ddk, palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("k = 4") + geom_point(aes(shape = drugexpts2$Treatment))
p4 <- fviz_cluster(k5, geom = "",data = ddk, palette = "Paired", 
                   ggtheme =theme_minimal()) + ggtitle("k = 5") + geom_point(aes(shape = drugexpts2$Treatment))+
  labs(shape = "Treatment")
p4<- fviz_cluster(k5, geom = "",data = ddk, palette = "Paired", 
                   ggtheme =theme_minimal()) 
grid.arrange(p1, p2, p3, p4, nrow = 2)

# Add the clusters to the initial data 
drugexpts2 %>%
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

######## PAM clustering - less sensitive to outliers ############

pam.res1 <- pam(df, 2)
pam.res2 <- pam(df, 3)
pam.res3 <- pam(df, 4)
pam.res4 <- pam(df, 7)
dd <- cbind(df, cluster = pam.res$cluster)

#head(dd, n = 3)
pam.res1$medoids[,38]

pamp1<- fviz_cluster(pam.res1, data = dd, geom = "", palette="Paired", ggtheme =theme_minimal()) + ggtitle("PAM: k=2") + geom_point(aes(shape = drugexpts2$Treatment))
pamp2<- fviz_cluster(pam.res2, data = dd, geom = "", palette="Paired", ggtheme =theme_minimal()) + ggtitle("PAM: k=3") + geom_point(aes(shape = drugexpts2$Treatment))
pamp3<- fviz_cluster(pam.res3, data = dd, geom = "", palette="Paired", ggtheme =theme_minimal()) + ggtitle("PAM: k=4") + geom_point(aes(shape = drugexpts2$Treatment))
pamp4<- fviz_cluster(pam.res4, data = dd, geom = "", palette="Paired", ggtheme =theme_minimal()) + ggtitle("PAM: k=9") + geom_point(aes(shape = drugexpts2$Treatment))
pamp3

grid.arrange(pamp1,pamp2, pamp3, pamp4, nrow = 2)

######### Compare them ###########

grid.arrange(p1, p3, p4, pamp1, pamp3, pamp4, nrow=2)

#################################### Heatmaps ##############################################
## Row- and column-wise clustering 
druge2forclust <- drugexpts2[-c(5,6,23,37,86,96,97)]
y<- scale(druge2forclust)
heatmap.2(y)

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
## Plot heatmap 
mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 

pheatmap(y, color=brewer.pal(9,"Blues"))
# new heatmap
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")
mycol <- colorpanel(40, "darkblue", "yellow", "white")
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(mycl))

### Correlation matrix
c <- cor(t(y), method="pearson") 
as.matrix(c)[1:89,1:89]
# Correlation-based distance matrix 
d <- as.dist(1-c)
as.matrix(d)[1:89,1:89]
hr <- hclust(d, method = "complete", members=NULL)
names(hr)
par(mfrow = c(1, 2)); plot(hr, hang = 0.1); plot(hr, hang = -1) 

#k-means PAM

pamy <- pam(d, 4)
(kmcol <- pamy$clustering)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(kmcol))
### k-means Fuzzy clustering 

fannyy <- fanny(d, k=4, memb.exp = 1.5)
round(fannyy$membership, 2)[1:89,]

fannyyMA <- round(fannyy$membership, 2) > 0.10 
apply(fannyyMA, 1, function(x) paste(which(x), collapse="_"))

