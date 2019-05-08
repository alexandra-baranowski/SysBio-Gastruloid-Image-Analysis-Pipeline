library(devtools)
library(ggbiplot)
library(ggfortify)
library(dplyr)
library(ggplot2)
library(tidyr)
library(factoextra)
library(ggpubr)


######### Looking at individual parameters ########################
library(matrixStats)
library(ggplot2)

################ Do some summary statistics first ##################

# Finding the means of all the columns 
colMeans1 <- colMeans(ChiPT[, -c(5,6,24,38,98,108,109) ]) # only use the numeric ones
colMeans1
# since they differ by subgroups, calculate the mean for each:

# Find the means by subgroup
subgroup_means <- aggregate(. ~ ChiPT$Treatment, ChiPT[, c()], mean)
subgroup_means
ggplot(subgroup_means[,1:12])
# Calculate the variance - Covariance matrix of the dataframe 
varcovar <- var(ChiPT[, -c(5,6,24,38,98,108,109)])
varcovar
# Calculate a correlation matrix 
cormatrix <- cor(ChiPT[, -c(5,6,24,38,98,108,109)])

# Visualise this using corrplot 
library(corrplot)

res1 <- cor.mtest(ChiPT[, -c(5,6,24,38,98,108,109)], conf.level = 0.95)

corrplot(cor(ChiPT[, -c(5,6,24,38,98,108,109)]), method = "square", title = "Correlation Matrix of All Variables",
         order = "FPC", tl.cex = 0.5, tl.col = 'black', p.mat = res1$p, insig = "pch", pch.cex = 0.05)
# the tilt of the ellipse specifies the sign of the correlation 
# how elliptical it is specifies the magnitude 

# Plot pairwise scatter plot of all the variables in the dataset; not great though bc too many variables 
#can also add log = 'x' or y or both there 
scatter <- pairs(RAPT[, c(1,2,11,12,13,19,20)], cex.labels = 0.8, cex = 0.3,
                 col = ChiPT$Treatment, main = "Pairwise Correlations for Selected Parameters")


library(GGally)
# this shows scatter plots, smoothed densities for each variable by group, and overall correlations
ggpairs(data = ChiPT, columns = c(1,2,3,4,7,8), mapping = aes(color = Treatment))

# Can also make 3D plots 
library(scatterplot3d)
scatterplot3d(Rallh[, c(83, 4, 75)], color = as.numeric(Rallh$Pre.treatment)) #area ch0, propch2


########### Do MDS ##############################
library(MASS)

swiss.x <- as.matrix(Rallh[, -c(5,6,23,37,86,96,97,101)])
swiss.x <- as.matrix(Rallh[, -c(4,5,11,25,63,66,67)])
swiss.dist <- dist(swiss.x)
swiss.mds <- isoMDS(swiss.dist)
ggplot(swiss.mds, aes(x=swiss.mds$points[,1], y = swiss.mds$points[,2])) + geom_point(aes(colour=Rallh$Pre.treatment)) +
  scale_colour_brewer(palette = "Paired") + scale_x_continuous(name="MDS Points Dim1") +
  scale_y_continuous(name="MDS Points Dim2") +theme_minimal()+ theme(legend.position="right") 

plot(swiss.mds$points, type = "n")
text(swiss.mds$points, labels = as.factor(Rallh$Pre.treatment))
swiss.sh <- Shepard(swiss.dist, swiss.mds$points)
ggplot(swiss.sh, aes(x=swiss.sh$x, y=swiss.sh$y)) + geom_point() +
  scale_colour_brewer(palette = "Paired")

plot(swiss.sh, pch = ".")
lines(swiss.sh$x, swiss.sh$yf, type = "S")

################ Do PCA ############################################

#Rallh05 <- bind_rows(ChiActAPT05, ChiPT, ChiSB4305, SB4305, Wnt3a05, RAPT05)
#Without LogHu of all the channels
#Rallh.pca <- prcomp(ChiPT[ , -c(5,6,23,37,77,87,88)], center = TRUE, scale. = TRUE)
SB43chi.pca <- prcomp(ChiPT[ , -c(5,6,23,37,86,96,97,101)], center = TRUE, scale. = TRUE) # remove the
# columns that have discrete / non-numerical values and only keep the active variables for ALL PCAs

Rallh.pca <- prcomp(Rallh[ , -c(5,6,23,37,86,96,97,101)], center = TRUE, scale. = TRUE)

str(Rallh.pca)

###
library(rgl)
library(pca3d)

plot3d(Rallh.pca$scores[,1:3], col=Rallh$Treatment)

pca3d(Rallh.pca, group = Rallh[, 101], show.group.labels = FALSE,
      bg = 'white', show.ellipses = F)

snapshotPCA3d(file="first_3d.png")

individuals_plot <- fviz_pca_ind(Rallh.pca, col.ind = "cos2") + labs(title ="PCA", x = "PC1", y = "PC2") + scale_color_gradient2(low="white", mid="blue",
                                                                                                                                 high="red", midpoint=0.6)
individuals_plot
###### this is the one where normal and treated are separated and coloured by group ! ############
Rallh2 <- fviz_pca_ind(Rallh.pca, axes = c(1,2), label="none", 
                       habillage = R2$Pre.treatment, geom="point",
                       palette = "Paired",
                       addEllipses = T,
                       title = "PCA on Chi Pre-treated Gastruloids",ellipse.level = 0.95
) # scale_color_brewer(palette="Paired")#,
#addEllipses=TRUE, ellipse.level=0.95)

#png("Rallh2")

Rallh2 + theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) 


### selecting and visualising individuals with cos > 0.96 
individuals_plot_3 <- fviz_pca_ind(SB43chi.pca, select.ind = list(cos2 = 0.7), habillage = ChiPT$Treatment,
                                   ggtheme = theme_minimal(), title = "Signalling Experiment 72h PCA, cos2 > 0.5")
individuals_plot_3 + theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0)))


hi <- fviz_pca_biplot(Rallh.pca, axes = c(1,2), geom = "point", select.var = list(contrib = 10), size.var = 4, geom.var = c("arrow"), 
                      repel = TRUE, habillage = R2$Pre.treatment, col.var = "darkslategrey",addEllipses = F,
                      title = "PCA Biplot Showing Chi Pre-Treated hGlds and Top Contributing features") 

hi + scale_color_brewer(palette="Paired") + theme(plot.title = element_text(size=11, face="bold",margin = margin(, 0, 5, 0))) +
  theme(legend.position="right")

## Apply the PCA onto unknown gastruloids 

treatment <- ChiPT[, 86]
pretreatment <- ChiPT[, 76]
hi <- ggbiplot::ggbiplot(SB43chi.pca, geom = "point", alpha=0.5,
                         groups = ChiPT$Treatment, ellipse = TRUE,
                         ellipse.prob = 0.95, scale = 1, obs.scale = 1, var.scale = 1,
                         var.axes = F,varname.abbrev = F
                         ) + scale_color_discrete(name = '') +theme_minimal() #+ ggtitle("PCA Applied on Chi Pre-treated hGlds and Projected to other PTs") +
hi # this is the one used to then project the PCA onto it after 

### Remove the same columns as have been removed from the "Normal" Condition for the application of the PCA 
## only do this if required, in some cases all the columns will be the same already 
cols_to_keep <- intersect(colnames(ChiPT),colnames(ChiActAPT))

ChiActAPT <- ChiActAPT[,cols_to_keep, drop=FALSE]
RAPT <- ChiPT[,cols_to_keep, drop=FALSE]
ChiSB43PT <- ChiSB43PT[, cols_to_keep, drop = FALSE]
SB43PT <- SB43PT[, cols_to_keep, drop = FALSE]
Wnt3a <- Wnt3a[, cols_to_keep, drop = FALSE]
ChiPT <- ChiPT[, cols_to_keep, drop = FALSE]

## Keep only the 0.5 Treatment 
ChiActAPT05 <- subset(ChiActAPT, Treatment == "0.5 Chi")
RAPT05 <- subset(RAPT, Treatment == "0.5 Chi")
ChiSB4305 <- subset(ChiSB43PT, Treatment == "0.5 Chi")
SB4305 <- subset(SB43PT, Treatment == "0.5 Chi")
Wnt3a05 <- subset(Wnt3a, Treatment == "0.5 Chi")

## Alternatively, eep only the 3 Treatment 
ChiActAPT3 <- subset(ChiActAPT, Treatment == "3 Chi")
RAPT3 <- subset(RAPT, Treatment == "3 Chi")
ChiSB433 <- subset(ChiSB43PT, Treatment == "3 Chi")
SB433 <- subset(SB43PT, Treatment == "3 Chi")
Wnt3a3 <- subset(Wnt3a, Treatment == "3 Chi")

# Keep only the active variables 
ChiActAPT2 <- ChiActAPT05 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
SB43PT2 <- SB4305 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
Wnt3a2 <- Wnt3a05 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
ChiSB43PT2 <- ChiSB4305 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
RAPT2 <- RAPT05 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))

ChiActAPT33 <- ChiActAPT3 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
SB43PT33 <- SB433 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
Wnt3a33 <- Wnt3a3 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
ChiSB43PT3 <- ChiSB433 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))
RAPT33 <- RAPT3 %>%
  dplyr::select(-c(Cell.line,Cell.number,Date,Hrs.Post.A,Pre.treatment,Treatment,TreatmentConditions, PT.Treatment))

set.seed(123)
project.ChiActA <- predict(SB43chi.pca, ChiActAPT2)
project.SB43 <- predict(SB43chi.pca, SB43PT2)
project.Wnt3a <- predict(SB43chi.pca, Wnt3a2)
project.ChiSB43 <- predict(SB43chi.pca, ChiSB43PT2)
project.RA <- predict(SB43chi.pca, RAPT2)

set.seed(123)
project.ChiActA3 <- predict(SB43chi.pca, ChiActAPT33)
project.SB433 <- predict(SB43chi.pca, SB43PT33)
project.Wnt3a3 <- predict(SB43chi.pca, Wnt3a33)
project.ChiSB433 <- predict(SB43chi.pca, ChiSB43PT3)
project.RA3 <- predict(SB43chi.pca, RAPT33)

# Map it back onto the original PCA

hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.ChiActA[, "PC1"],
                   yvar = project.ChiActA[, "PC2"],
                   groups = "Chi ActA 0.5 Chi" 
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.ChiActA3[, "PC1"],
                   yvar = project.ChiActA3[, "PC2"],
                   groups = "Chi ActA 3 Chi" 
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.SB43[, "PC1"],
                   yvar = project.SB43[, "PC2"],
                   groups = "SB43 0.5 Chi"
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.SB433[, "PC1"],
                   yvar = project.SB433[, "PC2"],
                   groups = "SB43 3 Chi"
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.Wnt3a[, "PC1"],
                   yvar = project.Wnt3a[, "PC2"],
                   groups = "Wnt3a 0.5 Chi"
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.Wnt3a3[, "PC1"],
                   yvar = project.Wnt3a3[, "PC2"],
                   groups = "Wnt3a 3 Chi"
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.ChiSB43[, "PC1"],
                   yvar = project.ChiSB43[, "PC2"],
                   groups = "ChiSB43 0.5 Chi"
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.ChiSB433[, "PC1"],
                   yvar = project.ChiSB433[, "PC2"],
                   groups = "ChiSB43 3 Chi"
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.RA[, "PC1"],
                   yvar = project.RA[, "PC2"],
                   groups = "RA 0.5 Chi"
                 )
)
hi$data <- rbind(hi$data, 
                 data.frame(
                   xvar = project.RA3[, "PC1"],
                   yvar = project.RA3[, "PC2"],
                   groups = "RA 3 Chi"
                 )
)
print(hi) 

#### do QQ plot to check for any outliers that need to be removed 
library(car)
plot(Rallh.pca$x[,1], pch = 20, col = ChiPT$Treatment)  
qqPlot(Rallh.pca$x[,1],pch = 20, col = Rallh$Treatment)

plot(Rallh, pch = 20, col = Rallh$Pre.treatment)  
qqPlot(Rallh$Max.Fluorescence,pch = 20, col = Rallh$Pre.treatment)

### Remove the outliers and Re-run the prediction
hi$data<- hi$data[-c(423), ]


print(hi) + scale_color_brewer(palette = "Paired") + ggtitle("PCA Performed on Chi Pre-treated hGlds \nand Projected to other PTs") +
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0)))

##### Looking at the PCA results and the properties #########

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
# get eigenvalues and variance retained by each component 
eig.val <- get_eigenvalue(Rallh.pca)
(eig.val)


# do scree plot of the pca, rule of thumb don't use any PCs beyond the 'knee' if doing elbow curve
# or can also use eigenvalues

fviz_screeplot(SB43chi.pca, ncp=10, addlabels=TRUE) # % of variance explained
fviz_screeplot(SB43chi.pca, ncp=10, choice="eigenvalue") # eigenvalue of the PC; eigenvalue of >1 is considered 
#to be the cut-off as it means that the PC accounts for more variance than the one of the original variables  

# or, can plot the cumulative proportion of variation explained 
# 1. variance explained 
pc.var <- Rallh.pca$sdev^2
# 2. proportion of variation 
pc.pvar <- pc.var / sum(pc.var)
# 3. cumulative proportion 
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9) # where 90% of the variation is explained 


####
# look at how much each variable contributes to each PC 
attributes(Rallh.pca)
Rallh.pca$rotation

# Looking at correlations and variances with factoextra; extracting results for the active variables 
var <- get_pca_var(Rallh.pca)


# graphing a correlation circle, showing the variables that contribute the most 
head(var$coord, 6)
fviz_pca_var(SB43chi.pca, col.var = "contrib", alpha = 0.6,
             repel = TRUE, select.var = list(contrib = 16)) + theme_minimal()
#Positively correlated variables are grouped together.
#Negatively correlated variables are positioned on opposite sides of the plot origin (opposed quadrants).
#The distance between variables and the origine measures the quality of the variables on the factor 
#map. Variables that are away from the origin are well represented on the factor map.

# Quality of representation 
head(var$cos2, 4)
library(corrplot)

#is.corr=false means that the input matrix is not a correlation matrix
corrplot(var$cos2, is.corr=FALSE, method = "circle", tl.cex = 0.4)
# if the variable is perfectly represented by the first 2 PCs, then that value will be 1 and it will lie 
# on the circumference of the above circe 

#we can also make a bar plot of total cos2 of the variables on dim1 and dim 2
fviz_cos2(Rallh.pca, choice = "var", axes = 1:2)

### we can then return to the correlation circle and colour by cos2 values for each variable 
fviz_pca_var(Rallh.pca, col.var = "cos2", select.var = list(contrib = 7),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

####### looking at the contributions of different variables to the PCs
head(var$contrib, 4)

# look at which are the most contributing 
corrplot(var$contrib, is.corr=FALSE, tl.cex = 0.4)

# we can also show this in a bar plot 
# Contributions of variables to PC1
fviz_contrib(Rallh.pca, choice = "var", axes = 1, top = 20)

# Contributions of variables to PC2
fviz_contrib(Rallh.pca, choice = "var", axes = 2, top = 20, fill = "lightcyan3", color = "cadetblue") +
  ggtitle("Contribution of Features to PC1") +
  theme(plot.title = element_text(size=12, face="bold",margin = margin(0, 10, 5, 10))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1))

# to PC1 and 2 combined 
fviz_contrib(Rallh.pca, choice = "var", axes = 1:2, top = 15)
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
# The red dashed line on the graph above indicates the expected average contribution. If the contribution of the 
#variables were uniform, the expected value would be 1/length(variables) = 1/10 = 10%. For a given component, a 
#variable with a contribution larger than this cutoff could be considered as important in contributing to the component

#looking at the most contributing factors 
fviz_pca_var(Rallh.pca, col.var = "contrib", select.var = list(contrib = 35),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE,
)

############################# Signalling plots ###########################################
### There are many different types of plot below, each used as required, with the variables able to 
# be changed ###

library(plyr)
library(tidyr)

prop1 <- ChiActAPT %>%
  dplyr::select(Proportion.Ch1, Proportion.Ch2, Proportion.Ch3, St.Dev.Fluorescence.Ch1, St.Dev.Fluorescence.Ch2,
                St.Dev.Fluorescence.Ch3, Ch.1.3.Overlap...Gld,Overlap.of.1.3..1.3area,Treatment, 
                Mean.Fluorescence.Ch3,Date, Hrs.Post.A, Pre.treatment, Area.Ch0, Ch1.3.Overlap...Ch3.Fluo,
                Contour.Area.Ch1)

prop2<- ChiSB43PT %>%
  dplyr::select(Proportion.Ch1, Proportion.Ch2, Proportion.Ch3, St.Dev.Fluorescence.Ch1, St.Dev.Fluorescence.Ch2,
                St.Dev.Fluorescence.Ch3, Ch.1.3.Overlap...Gld,Overlap.of.1.3..1.3area,Treatment, 
                Mean.Fluorescence.Ch3,Date, Hrs.Post.A, Pre.treatment, Area.Ch0, Ch1.3.Overlap...Ch3.Fluo,
                Contour.Area.Ch1)
prop3<- ChiPT %>%
  dplyr::select(Proportion.Ch1, Proportion.Ch2, Proportion.Ch3, St.Dev.Fluorescence.Ch1, St.Dev.Fluorescence.Ch2,
                St.Dev.Fluorescence.Ch3, Ch.1.3.Overlap...Gld,Overlap.of.1.3..1.3area,Treatment, 
                Mean.Fluorescence.Ch3,Date, Hrs.Post.A, Pre.treatment, Area.Ch0, Ch1.3.Overlap...Ch3.Fluo,
                Contour.Area.Ch1)
prop4 <- SB43PT %>%
  dplyr::select(Proportion.Ch1, Proportion.Ch2, Proportion.Ch3, St.Dev.Fluorescence.Ch1, St.Dev.Fluorescence.Ch2,
                St.Dev.Fluorescence.Ch3, Ch.1.3.Overlap...Gld,Overlap.of.1.3..1.3area,Treatment, 
                Mean.Fluorescence.Ch3,Date, Hrs.Post.A, Pre.treatment, Area.Ch0, Ch1.3.Overlap...Ch3.Fluo,
                Contour.Area.Ch1)
prop5<- Wnt3a %>%
  dplyr::select(Proportion.Ch1, Proportion.Ch2, Proportion.Ch3, St.Dev.Fluorescence.Ch1, St.Dev.Fluorescence.Ch2,
                St.Dev.Fluorescence.Ch3, Ch.1.3.Overlap...Gld,Overlap.of.1.3..1.3area,Treatment, 
                Mean.Fluorescence.Ch3,Date, Hrs.Post.A, Pre.treatment, Area.Ch0, Ch1.3.Overlap...Ch3.Fluo,
                Contour.Area.Ch1)
prop6<- RAPT %>%
  dplyr::select(Proportion.Ch1, Proportion.Ch2, Proportion.Ch3, St.Dev.Fluorescence.Ch1, St.Dev.Fluorescence.Ch2,
                St.Dev.Fluorescence.Ch3, Ch.1.3.Overlap...Gld,Overlap.of.1.3..1.3area,Treatment, 
                Mean.Fluorescence.Ch3,Date, Hrs.Post.A, Pre.treatment, Area.Ch0, Ch1.3.Overlap...Ch3.Fluo,
                Contour.Area.Ch1)

together <- rbind(prop1, prop2, prop3, prop4, prop5, prop6)


library(ggridges)
h <- ggplot(Rallh, aes(x = Rallh$Max.Fluorescence.Ch3, y = Rallh$Pre.treatment)) +
  geom_density_ridges(aes(fill = Rallh$Treatment, alpha =0.9),scale =0.9) +
  scale_fill_manual(values = c("cadetblue", "cadetblue2")) + scale_y_discrete(expand = c(0.14, 0)) +
  theme_minimal()

ggplot(Rallh, aes(x = Rallh$Overlap.of.1.3..1.3area, y = Rallh$Overlap.of.2.3..2.3area, size=Rallh$Proportion.Ch3)) +
  geom_point(aes(colour = Rallh$PT.Treatment)) + theme_minimal() +scale_color_brewer(palette="Paired")

ggplot(Rallh, aes(x = Rallh$ch1, y = Rallh$Proportion.Ch2, size=Rallh$Proportion.Ch1, 
                  colour=Rallh$PT.Treatment)) +
  geom_point()  +theme_minimal() +scale_color_brewer(palette="Paired") + geom_smooth(method="lm", se=F)

#scale_fill_manual(values = c("cadetblue", "cadetblue1", "azure"))
h +ggtitle("Distribution of Gastruloid Proportion Expressing Sox2 for each PT") + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) + 
  xlab("Proportion of Gastruloid") + ylab("Pre-Treatment") + labs(fill = "Aggregation \nConditions (uM)") +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9)) 

##### Make it into a faceted plot to show T/Bra Expression Levels as well 

pl <- ggplot(subset(together, Pre.treatment %in% c("Chi", "Chi ActA", "Chi SB43", "SB43", "RA", "Wnt3a") & Treatment %in% c("0.5 Chi", "3 Chi")),
             aes(x = Treatment, y = Proportion.Ch1,  colour = Proportion.Ch3, alpha =0.3, cmap='viridis')) + facet_wrap(~ Pre.treatment) +
  geom_point(alpha = 0.6) + geom_boxplot(alpha = 0, color = "azure4") +
  ylab("Proportion of Gastruloid Expressing SOX2") + xlab("Aggregation Condition (uM)") +
  ggtitle("Expression of SOX2 and T/BRA for Gastruloids \nof Different Pre-Treatments")

pl + theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) + 
  scale_fill_continuous(name = "T/BRA Proportion \n of Gastruloid")

######
ggplot(together, aes(x = together$Ch1.3.Overlap...Ch3.Fluo, y = prop6$Treatment)) +
  geom_density_ridges_gradient(aes(fill = prop6$Treatment)) +
  scale_fill_viridis_c(name = "Overlap", option = "C") + 
  scale_fill_manual(values = c("cadetblue", "cadetblue1", "azure")) +
  labs(title = 'title')

quartiles <- ggplot(prop6, aes(x=Ch1.3.Overlap...Ch3.Fluo, y=Treatment, fill=factor(..quantile..))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) +
  scale_fill_viridis_d( name = "Quartiles")
quartiles


pl + labs(color = "Gastruloid Area") +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
  stat_compare_means( label = "p.signif", method = "t.test",
                      ref.group = "Control")

############################## Clustering Analysis ###############################################

testClust <- data.frame(cbind(SB43chi.pca$x,ChiPT$Treatment))
colnames(testClust)[88] <- "Treatment"


unique(testClust$Treatment)
ttt <- testClust[,1:14]

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
