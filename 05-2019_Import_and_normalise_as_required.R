# import the csv file 
ChiActAPT <- read.csv("~/Desktop/Project Code/20192904/2019-03-04 Chi ActA 72h 2904.csv", row.names = 1)
ChiPT <- read.csv("~/Desktop/Project Code/20192904/2019-03-04 Chi 72h 2904.csv", row.names = 1)
ChiSB43PT <- read.csv("~/Desktop/Project Code/20192904/2019-03-04 Chi SB43 72h 2904.csv", row.names = 1)
Wnt3a <- read.csv("~/Desktop/Project Code/20192904/2019-03-04 Wnt3a 72h 2904.csv", row.names = 1)
#Wnt3a <- Wnt3a[-c(51, 74,87),]
SB43PT <- read.csv("~/Desktop/Project Code/20192904/2019-03-04 SB43 72h 2904.csv", row.names = 1)
#SB43PT <- SB43PT[-c(23,87),]
RAPT <- read.csv("~/Desktop/Project Code/20192904/2019-03-04 RA 72h 2904.csv", row.names = 1)

#######
# Check to see whether log Hu or Hu is better, and remove the worse ones 
Hupairs1 <- pairs(ChiPT[, c(39,42,45,48,4,16,102,103,104)], cex.labels = 0.8, cex = 0.3,
                  col = ChiPT$Treatment, main = "Pairwise Correlations for Selected Parameters")
LogHupairs1 <- pairs(tester[, c(63,66,69,72,4,16,102,103,104)], cex.labels = 0.8, cex = 0.3,
                     col = Rallh$Treatment, main = "Pairwise Correlations for Selected Parameters")

Hupairs2 <- pairs(ChiPT[, c(51,54,57,4,16,102,103,104)], cex.labels = 0.8, cex = 0.3,
                  col = ChiPT$Treatment, main = "Pairwise Correlations for Selected Parameters")

LogHupairs <- pairs(tester[, c(74,77,80,4,16,102,103,104)], cex.labels = 0.8, cex = 0.3,
                    col = Rallh$Treatment, main = "Pairwise Correlations for Selected Parameters")

### Keep Hu Here 
RAPT <- RAPT %>%
  dplyr::select(-Log.Hu.I1.Ch0, -Log.Hu.I2.Ch0, -Log.Hu.I3.Ch0, -Log.Hu.I4.Ch0, -Log.Hu.I5.Ch0, 
                -Log.Hu.I6.Ch0,-Log.Hu.I7.Ch0, -Log.Hu.I1.Ch1,
                -Log.Hu.I2.Ch1, -Log.Hu.I3.Ch1, -Log.Hu.I4.Ch1, -Log.Hu.I5.Ch1, -Log.Hu.I6.Ch1,-Log.Hu.I7.Ch1,
                -Log.Hu.I1.Ch3,
                -Log.Hu.I2.Ch3, -Log.Hu.I3.Ch3, -Log.Hu.I5.Ch3, -Log.Hu.I6.Ch3,-Log.Hu.I7.Ch3)

# Extract only the active variables, and remove anything with arrays in one column
RAPT <- RAPT %>%
  dplyr::select(-Min.rot.rect.Ch0, -Min.rot.rect.Ch1, -Min.rot.rect.Ch3, -Ch2.3.Overlap...Ch3.Fluo)

# Make sure all have the same number of columns 
cols_to_keep <- intersect(colnames(RAPT),colnames(ChiPT))

ChiActAPT <- ChiActAPT[,cols_to_keep, drop=FALSE]
ChiActAPT <- ChiActAPT[-c(22),]
ChiPT <- ChiPT[,cols_to_keep, drop=FALSE]
ChiSB43PT <- ChiSB43PT[, cols_to_keep, drop = FALSE]
SB43PT <- SB43PT[, cols_to_keep, drop = FALSE]
Wnt3a <- Wnt3a[, cols_to_keep, drop = FALSE]
RAPT <- RAPT[, cols_to_keep, drop = FALSE]


# Replace NaN and inf values with 0 
is.na(ChiPT)<-sapply(ChiPT, is.infinite)
ChiPT[is.na(ChiPT)]<-0

is.na(ChiSB43PT)<-sapply(ChiSB43PT, is.infinite)
ChiSB43PT[is.na(ChiSB43PT)]<-0

is.na(RAPT)<-sapply(RAPT, is.infinite)
RAPT[is.na(RAPT)]<-0

is.na(ChiActAPT)<-sapply(ChiActAPT, is.infinite)
ChiActAPT[is.na(ChiActAPT)]<-0

is.na(Wnt3a)<-sapply(Wnt3a, is.infinite)
Wnt3a[is.na(Wnt3a)]<-0

is.na(SB43PT)<-sapply(SB43PT, is.infinite)
SB43PT[is.na(SB43PT)]<-0

# Convert all the discrete variables to factors 
ChiActAPT[,'Hrs.Post.A']<-factor(ChiActAPT[,'Hrs.Post.A'])
ChiActAPT[, 'Cell.number'] <-factor(ChiActAPT[,'Cell.number'])
ChiActAPT[, 'Treatment'] <-factor(ChiActAPT[,'Treatment'])
ChiActAPT[, 'Pre.treatment'] <-factor(ChiActAPT[,'Pre.treatment'])

ChiSB43PT[,'Hrs.Post.A']<-factor(ChiSB43PT[,'Hrs.Post.A'])
ChiSB43PT[, 'Cell.number'] <-factor(ChiSB43PT[,'Cell.number'])
ChiSB43PT[, 'Treatment'] <-factor(ChiSB43PT[,'Treatment'])
ChiSB43PT[, 'Pre.treatment'] <-factor(ChiSB43PT[,'Pre.treatment'])

RAPT[,'Hrs.Post.A']<-factor(RAPT[,'Hrs.Post.A'])
RAPT[, 'Cell.number'] <-factor(RAPT[,'Cell.number'])
RAPT[, 'Treatment'] <-factor(RAPT[,'Treatment'])
RAPT[, 'Pre.treatment'] <-factor(RAPT[,'Pre.treatment'])

SB43PT[,'Hrs.Post.A']<-factor(SB43PT[,'Hrs.Post.A'])
SB43PT[, 'Cell.number'] <-factor(SB43PT[,'Cell.number'])
SB43PT[, 'Treatment'] <-factor(SB43PT[,'Treatment'])
SB43PT[, 'Pre.treatment'] <-factor(SB43PT[,'Pre.treatment'])

ChiPT[,'Hrs.Post.A']<-factor(ChiPT[,'Hrs.Post.A'])
ChiPT[, 'Cell.number'] <-factor(ChiPT[,'Cell.number'])
ChiPT[, 'Treatment'] <-factor(ChiPT[,'Treatment'])
ChiPT[, 'Pre.treatment'] <-factor(ChiPT[,'Pre.treatment'])

Wnt3a[,'Hrs.Post.A']<-factor(Wnt3a[,'Hrs.Post.A'])
Wnt3a[, 'Cell.number'] <-factor(Wnt3a[,'Cell.number'])
Wnt3a[, 'Treatment'] <-factor(Wnt3a[,'Treatment'])
Wnt3a[, 'Pre.treatment'] <-factor(Wnt3a[,'Pre.treatment'])



##### Normalise all the variables that need normalising 

############### Normalise the channel shape features ###################
## ChiPT
ChiPT$PT.Treatment <- paste(ChiPT$Pre.treatment, ChiPT$Treatment)

ChiPT$Total.Area.Ch1 <- ChiPT$Total.Area.Ch1 / max(ChiPT$Total.Area.Ch1)
ChiPT$Total.Area.Ch3 <- ChiPT$Total.Area.Ch3 / max(ChiPT$Total.Area.Ch3)
ChiPT$Total.Area.Ch2 <- ChiPT$Total.Area.Ch2 / max(ChiPT$Total.Area.Ch2)

ChiPT$Contour.Area.Ch1 <- ChiPT$Contour.Area.Ch1 / max(ChiPT$Contour.Area.Ch1)
ChiPT$Contour.Area.Ch3 <- ChiPT$Contour.Area.Ch1 / max(ChiPT$Contour.Area.Ch3)

ChiPT$Perimeter.Ch1 <- ChiPT$Perimeter.Ch1 / max(ChiPT$Perimeter.Ch1)
ChiPT$Perimeter.Ch3 <- ChiPT$Perimeter.Ch3 / max(ChiPT$Perimeter.Ch3)

ChiPT$Union.Ch1.Ch2 <- ChiPT$Union.Ch1.Ch2 / max(ChiPT$Union.Ch1.Ch2)
ChiPT$Union.Ch1.Ch3 <- ChiPT$Union.Ch1.Ch3 / max(ChiPT$Union.Ch1.Ch3)
ChiPT$Union.Ch2.Ch3 <- ChiPT$Union.Ch2.Ch3 / max(ChiPT$Union.Ch2.Ch3)

ChiPT$MECradius.Ch1 <- ChiPT$MECradius.Ch1 / max(ChiPT$MECradius.Ch1)
ChiPT$MECradius.Ch3 <- ChiPT$MECradius.Ch3 / max(ChiPT$MECradius.Ch3)

ChiPT$Hull.Area.Ch1 <- ChiPT$Hull.Area.Ch1 / max(ChiPT$Hull.Area.Ch1)
ChiPT$Hull.Area.Ch3 <- ChiPT$Hull.Area.Ch3 / max(ChiPT$Hull.Area.Ch3)

## ChiActAPT
ChiActAPT$PT.Treatment <- paste(ChiActAPT$Pre.treatment, ChiActAPT$Treatment)

ChiActAPT$Total.Area.Ch1 <- ChiActAPT$Total.Area.Ch1 / max(ChiActAPT$Total.Area.Ch1)
ChiActAPT$Total.Area.Ch3 <- ChiActAPT$Total.Area.Ch3 / max(ChiActAPT$Total.Area.Ch3)
ChiActAPT$Total.Area.Ch2 <- ChiActAPT$Total.Area.Ch2 / max(ChiActAPT$Total.Area.Ch2)

ChiActAPT$Contour.Area.Ch1 <- ChiActAPT$Contour.Area.Ch1 / max(ChiActAPT$Contour.Area.Ch1)
ChiActAPT$Contour.Area.Ch3 <- ChiActAPT$Contour.Area.Ch1 / max(ChiActAPT$Contour.Area.Ch3)

ChiActAPT$Perimeter.Ch1 <- ChiActAPT$Perimeter.Ch1 / max(ChiActAPT$Perimeter.Ch1)
ChiActAPT$Perimeter.Ch3 <- ChiActAPT$Perimeter.Ch3 / max(ChiActAPT$Perimeter.Ch3)

ChiActAPT$Union.Ch1.Ch2 <- ChiActAPT$Union.Ch1.Ch2 / max(ChiActAPT$Union.Ch1.Ch2)
ChiActAPT$Union.Ch1.Ch3 <- ChiActAPT$Union.Ch1.Ch3 / max(ChiActAPT$Union.Ch1.Ch3)
ChiActAPT$Union.Ch2.Ch3 <- ChiActAPT$Union.Ch2.Ch3 / max(ChiActAPT$Union.Ch2.Ch3)

ChiActAPT$MECradius.Ch1 <- ChiActAPT$MECradius.Ch1 / max(ChiActAPT$MECradius.Ch1)
ChiActAPT$MECradius.Ch3 <- ChiActAPT$MECradius.Ch3 / max(ChiActAPT$MECradius.Ch3)

ChiActAPT$Hull.Area.Ch1 <- ChiActAPT$Hull.Area.Ch1 / max(ChiActAPT$Hull.Area.Ch1)
ChiActAPT$Hull.Area.Ch3 <- ChiActAPT$Hull.Area.Ch3 / max(ChiActAPT$Hull.Area.Ch3)

### Chi SB43
ChiSB43PT$PT.Treatment <- paste(ChiSB43PT$Pre.treatment, ChiSB43PT$Treatment)

ChiSB43PT$Total.Area.Ch1 <- ChiSB43PT$Total.Area.Ch1 / max(ChiSB43PT$Total.Area.Ch1)
ChiSB43PT$Total.Area.Ch3 <- ChiSB43PT$Total.Area.Ch3 / max(ChiSB43PT$Total.Area.Ch3)
ChiSB43PT$Total.Area.Ch2 <- ChiSB43PT$Total.Area.Ch2 / max(ChiSB43PT$Total.Area.Ch2)

ChiSB43PT$Contour.Area.Ch1 <- ChiSB43PT$Contour.Area.Ch1 / max(ChiSB43PT$Contour.Area.Ch1)
ChiSB43PT$Contour.Area.Ch3 <- ChiSB43PT$Contour.Area.Ch1 / max(ChiSB43PT$Contour.Area.Ch3)

ChiSB43PT$Perimeter.Ch1 <- ChiSB43PT$Perimeter.Ch1 / max(ChiSB43PT$Perimeter.Ch1)
ChiSB43PT$Perimeter.Ch3 <- ChiSB43PT$Perimeter.Ch3 / max(ChiSB43PT$Perimeter.Ch3)

ChiSB43PT$Union.Ch1.Ch2 <- ChiSB43PT$Union.Ch1.Ch2 / max(ChiSB43PT$Union.Ch1.Ch2)
ChiSB43PT$Union.Ch1.Ch3 <- ChiSB43PT$Union.Ch1.Ch3 / max(ChiSB43PT$Union.Ch1.Ch3)
ChiSB43PT$Union.Ch2.Ch3 <- ChiSB43PT$Union.Ch2.Ch3 / max(ChiSB43PT$Union.Ch2.Ch3)

ChiSB43PT$MECradius.Ch1 <- ChiSB43PT$MECradius.Ch1 / max(ChiSB43PT$MECradius.Ch1)
ChiSB43PT$MECradius.Ch3 <- ChiSB43PT$MECradius.Ch3 / max(ChiSB43PT$MECradius.Ch3)

ChiSB43PT$Hull.Area.Ch1 <- ChiSB43PT$Hull.Area.Ch1 / max(ChiSB43PT$Hull.Area.Ch1)
ChiSB43PT$Hull.Area.Ch3 <- ChiSB43PT$Hull.Area.Ch3 / max(ChiSB43PT$Hull.Area.Ch3)

## SB43
SB43PT$PT.Treatment <- paste(SB43PT$Pre.treatment, SB43PT$Treatment)

SB43PT$Total.Area.Ch1 <- SB43PT$Total.Area.Ch1 / max(SB43PT$Total.Area.Ch1)
SB43PT$Total.Area.Ch3 <- SB43PT$Total.Area.Ch3 / max(SB43PT$Total.Area.Ch3)
SB43PT$Total.Area.Ch2 <- SB43PT$Total.Area.Ch2 / max(SB43PT$Total.Area.Ch2)

SB43PT$Contour.Area.Ch1 <- SB43PT$Contour.Area.Ch1 / max(SB43PT$Contour.Area.Ch1)
SB43PT$Contour.Area.Ch3 <- SB43PT$Contour.Area.Ch1 / max(SB43PT$Contour.Area.Ch3)

SB43PT$Perimeter.Ch1 <- SB43PT$Perimeter.Ch1 / max(SB43PT$Perimeter.Ch1)
SB43PT$Perimeter.Ch3 <- SB43PT$Perimeter.Ch3 / max(SB43PT$Perimeter.Ch3)

SB43PT$Union.Ch1.Ch2 <- SB43PT$Union.Ch1.Ch2 / max(SB43PT$Union.Ch1.Ch2)
SB43PT$Union.Ch1.Ch3 <- SB43PT$Union.Ch1.Ch3 / max(SB43PT$Union.Ch1.Ch3)
SB43PT$Union.Ch2.Ch3 <- SB43PT$Union.Ch2.Ch3 / max(SB43PT$Union.Ch2.Ch3)

SB43PT$MECradius.Ch1 <- SB43PT$MECradius.Ch1 / max(SB43PT$MECradius.Ch1)
SB43PT$MECradius.Ch3 <- SB43PT$MECradius.Ch3 / max(SB43PT$MECradius.Ch3)

SB43PT$Hull.Area.Ch1 <- SB43PT$Hull.Area.Ch1 / max(SB43PT$Hull.Area.Ch1)
SB43PT$Hull.Area.Ch3 <- SB43PT$Hull.Area.Ch3 / max(SB43PT$Hull.Area.Ch3)

### RAPT 
RAPT$PT.Treatment <- paste(RAPT$Pre.treatment, RAPT$Treatment)

RAPT$Total.Area.Ch1 <- RAPT$Total.Area.Ch1 / max(RAPT$Total.Area.Ch1)
RAPT$Total.Area.Ch3 <- RAPT$Total.Area.Ch3 / max(RAPT$Total.Area.Ch3)
RAPT$Total.Area.Ch2 <- RAPT$Total.Area.Ch2 / max(RAPT$Total.Area.Ch2)

RAPT$Contour.Area.Ch1 <- RAPT$Contour.Area.Ch1 / max(RAPT$Contour.Area.Ch1)
RAPT$Contour.Area.Ch3 <- RAPT$Contour.Area.Ch1 / max(RAPT$Contour.Area.Ch3)

RAPT$Perimeter.Ch1 <- RAPT$Perimeter.Ch1 / max(RAPT$Perimeter.Ch1)
RAPT$Perimeter.Ch3 <- RAPT$Perimeter.Ch3 / max(RAPT$Perimeter.Ch3)

RAPT$Union.Ch1.Ch2 <- RAPT$Union.Ch1.Ch2 / max(RAPT$Union.Ch1.Ch2)
RAPT$Union.Ch1.Ch3 <- RAPT$Union.Ch1.Ch3 / max(RAPT$Union.Ch1.Ch3)
RAPT$Union.Ch2.Ch3 <- RAPT$Union.Ch2.Ch3 / max(RAPT$Union.Ch2.Ch3)

RAPT$MECradius.Ch1 <- RAPT$MECradius.Ch1 / max(RAPT$MECradius.Ch1)
RAPT$MECradius.Ch3 <- RAPT$MECradius.Ch3 / max(RAPT$MECradius.Ch3)

RAPT$Hull.Area.Ch1 <- RAPT$Hull.Area.Ch1 / max(RAPT$Hull.Area.Ch1)
RAPT$Hull.Area.Ch3 <- RAPT$Hull.Area.Ch3 / max(RAPT$Hull.Area.Ch3)

#Wnt3a 
Wnt3a$PT.Treatment <- paste(Wnt3a$Pre.treatment, Wnt3a$Treatment)

Wnt3a$Total.Area.Ch1 <- Wnt3a$Total.Area.Ch1 / max(Wnt3a$Total.Area.Ch1)
Wnt3a$Total.Area.Ch3 <- Wnt3a$Total.Area.Ch3 / max(Wnt3a$Total.Area.Ch3)
Wnt3a$Total.Area.Ch2 <- Wnt3a$Total.Area.Ch2 / max(Wnt3a$Total.Area.Ch2)

Wnt3a$Contour.Area.Ch1 <- Wnt3a$Contour.Area.Ch1 / max(Wnt3a$Contour.Area.Ch1)
Wnt3a$Contour.Area.Ch3 <- Wnt3a$Contour.Area.Ch1 / max(Wnt3a$Contour.Area.Ch3)

Wnt3a$Perimeter.Ch1 <- Wnt3a$Perimeter.Ch1 / max(Wnt3a$Perimeter.Ch1)
Wnt3a$Perimeter.Ch3 <- Wnt3a$Perimeter.Ch3 / max(Wnt3a$Perimeter.Ch3)

Wnt3a$Union.Ch1.Ch2 <- Wnt3a$Union.Ch1.Ch2 / max(Wnt3a$Union.Ch1.Ch2)
Wnt3a$Union.Ch1.Ch3 <- Wnt3a$Union.Ch1.Ch3 / max(Wnt3a$Union.Ch1.Ch3)
Wnt3a$Union.Ch2.Ch3 <- Wnt3a$Union.Ch2.Ch3 / max(Wnt3a$Union.Ch2.Ch3)

Wnt3a$MECradius.Ch1 <- Wnt3a$MECradius.Ch1 / max(Wnt3a$MECradius.Ch1)
Wnt3a$MECradius.Ch3 <- Wnt3a$MECradius.Ch3 / max(Wnt3a$MECradius.Ch3)

Wnt3a$Hull.Area.Ch1 <- Wnt3a$Hull.Area.Ch1 / max(Wnt3a$Hull.Area.Ch1)
Wnt3a$Hull.Area.Ch3 <- Wnt3a$Hull.Area.Ch3 / max(Wnt3a$Hull.Area.Ch3)

###########

# Put the dfs together - not always required
Rallh <- bind_rows(ChiActAPT, ChiPT, ChiSB43PT, RAPT, Wnt3a, SB43PT)
Rallh$PT.Treatment <- paste(Rallh$Pre.treatment, Rallh$Treatment)
#############

is.na(Rallh)<-sapply(Rallh, is.infinite)
Rallh[is.na(Rallh)]<-0

Rallh[,'Hrs.Post.A']<-factor(Rallh[,'Hrs.Post.A'])
Rallh[, 'Cell.number'] <-factor(Rallh[,'Cell.number'])
Rallh[, 'Treatment'] <-factor(Rallh[,'Treatment'])
Rallh[, 'Pre.treatment'] <-factor(Rallh[,'Pre.treatment'])

