###### Clustering Analysis #######
testClust <- data.frame(cbind(SB43chi.pca$x,ChiPT$Treatment))
colnames(testClust)[88] <- "Treatment"


unique(testClust$Treatment)
ttt <- testClust[,1:14] # however many PCs you want to include;  those explaining 
# >90 % cumulative variation / based on the scree plot

distm <- dist(ttt, method = "canberra")
hclustttt <- hclust(distm, method = "complete")

#hclustttt <- hclust(distm, method = "complete")
plot(hclustttt, xlab="Distance \nDistance calculation: Canberra, Clustering Method: Complete")
rect.hclust(hclustttt, k = 2, border = 2:6)

cuttt <- cutree(hclustttt, k = 2)

afterclust <- mutate(ttt, cluster = cuttt)
#count(seeds_df_cl,cluster)

afterclust <- cbind(afterclust, treatment = testClust$Treatment)
#View(afterclust[,11:12])

g1 <- ggplot(afterclust, aes(x = PC1, y = PC2, group = as.factor(cluster))) + 
  geom_point(aes(color = as.factor(treatment))) +
stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = as.factor(cluster)), level = 0.95) +
  ggtitle("Clusters produced using first 16 Principal Components of Chi Pre-treated hGlds")

g1+ theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) + 
  labs(color = "Treatment", fill = "Cluster") + scale_color_discrete(labels = c("0.5 Chi", "3 Chi"))
  

ggplot(afterclust, aes(x = PC1, y = PC2, group = as.factor(cluster))) + 
  geom_point(aes(shape = as.factor(Rallh$Pre.treatment), color = as.factor(Rallh$Pre.treatment))) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = as.factor(cluster)), level = 0.95)
