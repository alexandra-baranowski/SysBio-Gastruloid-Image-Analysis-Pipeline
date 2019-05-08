library(plyr)


prop1<- Rallh %>%
  select(Stdev.fluo.Ch1,TreatmentConditions, Date, Hrs.Post.A)
prop1$TreatmentConditions <- factor(prop1$TreatmentConditions, levels = c("Normal", "TD Control", "0.25mM TD", 
                                                                          "RA control", "50uM RA", "VA control", 
                                                                          "60mM VA"))
prop2<- Rallh %>%
  select(Stdev.fluo.Ch2,Mean.fluo.Ch2 ,TreatmentConditions, Date, Hrs.Post.A)
prop3<- Rallh %>%
  select(Stdev.fluo.Ch3,Mean.fluo.Ch3,TreatmentConditions, Date, Hrs.Post.A)

df <- data.frame(
  gp = factor(rep(letters[1:3], each = 10)),
  y = rnorm(30)
)

propplot1_2 <- ggplot(prop1, aes(x=Rallh$Hrs.Post.A, y=Rallh$Stdev.fluo.Ch1)) + geom_boxplot(varwidth = TRUE, aes(colour=Hrs.Post.A)) + xlab("Hrs AA") + ylab("Proportion of gastruloid expressing Sox2")
propplot2 <- ggplot(prop2, aes(x=prop2$Hrs.Post.A, y=prop2$Stdev.fluo.Ch2)) + geom_boxplot()
propplot3 <- ggplot(prop3, aes(x=prop3$Hrs.Post.A, y=prop3$Stdev.fluo.Ch3)) + geom_boxplot()

png("propplot1_2.png")
propplot1_2

ggsave("propplot1_2.png")
dev.off()

propplot2 + (aes(colour=Hrs.Post.A)) + xlab("Hours After Aggregation") + ylab("Channel 2 Standard Deviation in Fluorescence") +
  ggtitle("Standard Deviation in Fluorescence Values of SOX2::mCitrine") + theme(plot.title = element_text(size=11, face="bold",
                                                                                                           margin = margin(10, 0, 10, 0)))
propplot3

propplotcomb <- ggplot(proportions, aes(x=prop1$TreatmentConditions, y=prop1$Proportion.Ch1)) + geom_boxplot(varwidth = TRUE, aes(colour=Controls)) + xlab("Treatment") + ylab("Proportion of gastruloid expressing T/Bra")
propplotcomb

Rallh$Treatment<- factor(Rallh$Treatment, levels = c("N2B27", "DMSO", "Chi"))
Rallh$Pre.treatment<- factor(Rallh$Pre.treatment, levels = c("Control", "Chordin", "FGF2", "Noggin", "DKK", "Wnt3a"))
## Plot comparing different conditions by subgroup
pl <- ggplot(subset(Rallh, Treatment %in% c("N2B27", "DMSO", "Chi") & Hrs.Post.A %in% c("72", "96", "120")),
       aes(x = Pre.treatment, y = Rallh$Proportion.Ch1,  colour = Area.Ch0)) + facet_wrap( ~ Hrs.Post.A) +
  geom_point(alpha = 0.3) + geom_boxplot(alpha = 0, colour = "azure4",  scale = "width") +
  ylab("Gastruloid Area") + xlab("Treatment") +
  ggtitle("Chi Treatment Leads to Increased Area Over Time") + theme(plot.title = element_text(size=11, face="bold",
                                                                                                     margin = margin(10, 0, 10, 0)))


pl + labs(color = "T/Bra Expression \n(Proportion of \nGastruloid") +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
  stat_compare_means( label = "p.signif", method = "t.test",
                     ref.group = "N2B27", size = 2)

###### DKK 
## Plot comparing different conditions by subgroup
pl <- ggplot(subset(Rallh, Pre.treatment %in% c("Chordin", "Control", "DKK", "FGF2", "Noggin", "Wnt3a") & Hrs.Post.A %in% c("72", "96", "120")),
             aes(x = Pre.treatment, y = Mean.Fluorescence,  colour = Area.Ch0))  +
  geom_point(alpha = 0.3) + geom_boxplot(alpha = 0, colour = "black",  scale = "width",outlier.shape = NA) +
  ylab("Normalised Mean T/Bra Fluorescence ") + xlab("Treatment") +
  ggtitle("T/Bra Expression over time for Chi vs. DMSO and N2B27") + theme(plot.title = element_text(size=11, face="bold",
                                                                                                     margin = margin(10, 0, 10, 0)))


pl + labs(color = "Gastruloid Area") +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
  stat_compare_means( label = "p.signif", method = "t.test",
                      ref.group = "Control")

#######
## Plot comparing different conditions by subgroup
pl <- ggplot(subset(Rallh, Pre.treatment %in% c("Chordin", "Control", "DKK", "FGF2", "Noggin") & Hrs.Post.A %in% c("72", "96", "120")),
             aes(x = Pre.treatment, y = Proportion.Ch1,  colour = Area.Ch0)) + facet_wrap(~ Hrs.Post.A) +
  geom_point(alpha = 0.6) + geom_boxplot(alpha = 0, color = "azure4", outlier.shape = NA) +
  ylab("Normalised Mean T/Bra Fluorescence") + xlab("Pre-treatment Condition") +
  ggtitle("DKK Pre-treatment Delays the onset of T/Bra expression") 


pl + labs(color = "Gastruloid Area") +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
  stat_compare_means( label = "p.signif", method = "t.test",
                      ref.group = "Control") + theme(plot.title = element_text(size=11, face="bold",
                                                                               margin = margin(10, 0, 10, 0))) +
  theme(axis.text.x=element_text(angle=40, size=9, vjust=0.5))
 

#### Do anova 

anova <- compare_means(Proportion.Ch1 ~ Hrs.Post.A,  data = Rallh, method = "anova")

anova

pl2 <- ggplot(subset(Rallh, Treatment %in% c("DMSO", "N2B27", "Chi") & Hrs.Post.A %in% c("72")),
              aes(x = Treatment, y = Proportion.Ch1,  colour = Proportion.Ch1)) +
   geom_boxplot(alpha = 0, colour = "black") + geom_point(alpha = 0.5)

my_comparisons <- list( c("Chi", "DMSO"), c("DMSO", "N2B27"), c("Chi", "N2B27") )
pl2 + stat_compare_means(method = "anova") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "N2B27")


my_comparisons <- list( c("72", "96"), c("96", "120"), c("72", "120") ) 

ggplot(subset(Rallh, Treatment %in% c("Chi", "DMSO", "N2B27") & Hrs.Post.A %in% c("72", "96", "120")),
       aes(x = Hrs.Post.A, y = Circularity.Ch0,  colour = Area.Ch0)) + facet_wrap( ~ Treatment) +
  geom_point(alpha = 0.5, position = "jitter") + geom_violin(alpha = 0, colour = "grey") 
################### Group by ###############

###### plotting fluorescences #########
# Making a heatmap of avg max, avg avg, and avg stdev 

forheatmap <- Rallh %>%
  dplyr::select("Max.fluo.Ch1",  "Max.fluo.Ch2", "Max.fluo.Ch3", "Mean.fluo.Ch1", 
         "Mean.fluo.Ch2", "Mean.fluo.Ch3", "Stdev.fluo.Ch1", "Stdev.fluo.Ch2", "Stdev.fluo.Ch3")
subgroup_means <- aggregate(. ~ Rallh$Hrs.Post.A, Rallh[, c(71,72,73,74,75,76,85,86,87)], mean)
head(subgroup_means)
head(forheatmap)
subgroup_means <- as.data.frame(subgroup_means)
# need to reshape the dataframe 
library(tidyverse)
sbm2 <- subgroup_means %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
  
head(sbm2)

ggplot(means2, aes(x = colname, y = rowname, fill = value)) +
  geom_tile()

###### doing the three separately so they have separate fill values 

#1. st dev of the fluorescences 
stdevhm <- Rallh %>%
  dplyr::select("Stdev.fluo.Ch1", "Stdev.fluo.Ch2", "Stdev.fluo.Ch3")
stdevhm2 <- stdevhm %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(stdevhm2)

ggplot(stdevhm2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()

#2. Mean fluorescence 
meanfluohm <- Rallh %>%
  select("Mean.fluo.Ch1", "Mean.fluo.Ch2", "Mean.fluo.Ch3")
meanfluohm2 <- meanfluohm %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(meanfluohm2)

ggplot(meanfluohm2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()

#3. Max fluorescence 
maxfluohm <- Rallh %>%
  select("Max.fluo.Ch1", "Max.fluo.Ch2", "Max.fluo.Ch3")
head(maxfluohm)
maxfluohm2 <- maxfluohm %>%
  rownames_to_column() %>%
  gather(colname, value, -(rowname))
head(maxfluohm2)

ggplot(maxfluohm2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()

# Correlation matrix
maxfluohm <- Rallh %>%
  select("Max.fluo.Ch1", "Max.fluo.Ch2", "Max.fluo.Ch3",
         "Mean.fluo.Ch1", "Mean.fluo.Ch2", "Mean.fluo.Ch3")
cormat <- round(cor(maxfluohm),2)
head(cormat)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri((cormat))
upper_tri

# Melt the correlation matrix

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#### reorder the matrix 
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

# to add the correlation coefficients 
ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
