### F score Plots

signFscores <- read.csv("~/Desktop/Signalling Expts F scores.csv")

precision<- ggplot(signFscores, aes(x=Condition, y = Precision, group = signFscores$Segmentation.Method)) + geom_point() +
  theme_minimal() + scale_color_brewer(palette="Paired") + expand_limits(y=c(.5, 1)) + 
  geom_line(aes(colour=Segmentation.Method)) + ylab("Precision Score") +
  ggtitle("Contour Precision Score") + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) +
  theme(axis.text.x = element_text(angle = 0,
                                   size = 8))

recall <- ggplot(signFscores, aes(x=Condition, y = Recall, group = signFscores$Segmentation.Method)) + geom_point() +
  theme_minimal() + scale_color_brewer(palette="Paired") + expand_limits(y=c(.5, 1)) + 
  geom_line(aes(colour=Segmentation.Method)) + ylab("Recall Score") +
  ggtitle("Contour Recall Score") + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) +
  theme(axis.text.x = element_text(angle = 0,
                                   size = 8))
f1<- ggplot(signFscores, aes(x=Condition, y = F1.Score, group = signFscores$Segmentation.Method)) + geom_point() +
  theme_minimal() + scale_color_brewer(palette="Paired") + expand_limits(y=c(.5, 1)) + 
  geom_line(aes(colour=Segmentation.Method)) + ylab("F1 Score") +
  ggtitle("Contour F1 Score for Three Different Segmentation Methods") + 
  theme(plot.title = element_text(size=11, face="bold",margin = margin(10, 0, 10, 0))) +
  theme(axis.text.x = element_text(angle = 0,
                                   size = 8))

grid.arrange(precision, recall, f1, nrow = 1)

