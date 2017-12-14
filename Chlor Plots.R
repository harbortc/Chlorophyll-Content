#Chlorophyll content plots master version
# 14.12.17
# CJ Harbort



#Set the working directory 

# Used internet example for stats on my data without understanding it
# https://stackoverflow.com/questions/18771516/is-there-a-function-to-add-aov-post-hoc-testing-results-to-ggplot2-boxplot
# This script gives me what I want, though I don't quite understand how yet. Works for now...

rm(list=ls())
library(plyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(multcompView)

# First, set directory in new folder for experiment and upload dataframe
setwd()
data.df <- read.delim()

# Rename Col-0 removing dash just in case
data.df$Genotype <- revalue(data.df$Genotype, c("Col-0" = "Col0"))


##########################
# Plot Fresh Weight input and stats
FW.aov <- aov(FW ~ Genotype, data = data.df)

FW.tHSD <- TukeyHSD(FW.aov, ordered = FALSE, conf.level = 0.95)


# this is a little voodoo to get stat group labels

generate_label_FW <- function(HSD, group){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[group]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(data.df, group, function (x) max(fivenum(x$FW))*1.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = group, sort = FALSE)
  
  return(labels.df)
}

FW.plot <- ggplot(data.df, aes(x= Genotype, y=FW)) + geom_boxplot() +
  geom_text(data = generate_label_FW(FW.tHSD, "Genotype"), aes(x = plot.labels, y = max(V1), label = labels)) +
  geom_jitter() +
  scale_y_continuous(name = "Fresh Weight (g)") + #limits = c(0, 3)) +
  ggtitle("Shoot Input Weight") +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

FW.plot


#######################################################
# Chlorophyll calculations and plots

# formula for chlor.total calculation
# Chlor.total (mg/L) = A652*1000/34.5
# Chlor.A (mg/L) = 20.2*A663
# Chlor.B (mg/L) = 8.02*A645

data.df$Chlor.total_perL <- (data.df$A652 * 1000/34.5)
data.df$Chlor.A_perL <- (data.df$A663 * 20.2)
data.df$Chlor.B_perL <- (data.df$A645 * 8.02)

# define a function that calculates mgChlor/gFW 
# extr.vol is volume of DMSO used for extraction; default = 1

Chlor.mggFW <- function(x, gFW, extr.vol = 1) {(x/1000)*(extr.vol/gFW)}

# Calculate mgChlor/gFW and add to DF
data.df$Chlor.total.mgpergFW <- Chlor.mggFW(data.df$Chlor.total_perL, data.df$FW) # add extr.vol if not = 1
data.df$Chlor.A.mgpergFW <- Chlor.mggFW(data.df$Chlor.A_perL, data.df$FW) # add extr.vol if not = 1
data.df$Chlor.B.mgpergFW <- Chlor.mggFW(data.df$Chlor.B_perL, data.df$FW) # add extr.vol if not = 1


# Stats and plot Chlor.total

Chlor.total.aov <- aov(Chlor.total.mgpergFW ~ Genotype, data = data.df)

Chlor.total.tHSD <- TukeyHSD(Chlor.total.aov, ordered = FALSE, conf.level = 0.95)



generate_label_ChlT <- function(HSD, group){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[group]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(data.df, group, function (x) max(fivenum(x$Chlor.total.mgpergFW))*1.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = group, sort = FALSE)
  
  return(labels.df)
}

Chlor.total.plot <- ggplot(data.df, aes(x= Genotype, y=Chlor.total.mgpergFW)) + 
  geom_boxplot() +
  geom_text(data = generate_label_ChlT(Chlor.total.tHSD, "Genotype"), aes(x = plot.labels, y = max(V1), label = labels)) +
  geom_jitter() +
  scale_y_continuous(name = "Total Chlorophyll\n(mg/gFW)") + #, limits = c(0, 3)) +
  ggtitle("Total Chlorophyll Content") +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

Chlor.total.plot



##############################
# Chlor A content
# Stats and plot Chlor.total

Chlor.A.aov <- aov(Chlor.A.mgpergFW ~ Genotype, data = data.df)

Chlor.A.tHSD <- TukeyHSD(Chlor.A.aov, ordered = FALSE, conf.level = 0.95)



generate_label_ChlA <- function(HSD, group){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[group]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(data.df, group, function (x) max(fivenum(x$Chlor.A.mgpergFW))*1.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = group, sort = FALSE)
  
  return(labels.df)
}

Chlor.A.plot <- ggplot(data.df, aes(x= Genotype, y=Chlor.A.mgpergFW)) + 
  geom_boxplot() +
  geom_text(data = generate_label_ChlA(Chlor.A.tHSD, "Genotype"), aes(x = plot.labels, y = max(V1), label = labels)) +
  geom_jitter() +
  scale_y_continuous(name = "Chlorophyll A\n(mg/gFW)") + #, limits = c(0, 3)) +
  ggtitle("Chlorophyll A Content") +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

Chlor.A.plot

##########################################
# Chlor B content

Chlor.B.aov <- aov(Chlor.B.mgpergFW ~ Genotype, data = data.df)

Chlor.B.tHSD <- TukeyHSD(Chlor.B.aov, ordered = FALSE, conf.level = 0.95)



generate_label_ChlB <- function(HSD, group){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[group]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(data.df, group, function (x) max(fivenum(x$Chlor.B.mgpergFW))*1.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = group, sort = FALSE)
  
  return(labels.df)
}

Chlor.B.plot <- ggplot(data.df, aes(x= Genotype, y=Chlor.B.mgpergFW)) + 
  geom_boxplot() +
  geom_text(data = generate_label_ChlB(Chlor.B.tHSD, "Genotype"), aes(x = plot.labels, y = max(V1), label = labels)) +
  geom_jitter() +
  scale_y_continuous(name = "Chlorophyll B\n(mg/gFW)") + #, limits = c(0, 3)) +
  ggtitle("Chlorophyll B Content") +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

Chlor.B.plot

######################################
# Ratio Chlor A to B 

data.df$RatioAB <- data.df$Chlor.A.mgpergFW/data.df$Chlor.B.mgpergFW


Chlor.rat.aov <- aov(RatioAB ~ Genotype, data = data.df)

Chlor.rat.tHSD <- TukeyHSD(Chlor.rat.aov, ordered = FALSE, conf.level = 0.95)


generate_label_Chl.rat <- function(HSD, group){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[group]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(data.df, group, function (x) max(fivenum(x$RatioAB))*1.05)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = group, sort = FALSE)
  
  return(labels.df)
}

Chlor.rat.plot <- ggplot(data.df, aes(x= Genotype, y=RatioAB)) + 
  geom_boxplot() +
  geom_text(data = generate_label_Chl.rat(Chlor.rat.tHSD, "Genotype"), aes(x = plot.labels, y = max(V1), label = labels)) +
  geom_jitter() +
  scale_y_continuous(name = "Chlorophyll A/B Ratio", limits = c(0, 11.5)) +
  ggtitle("Chlorophyll A/B Ratio") +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))

Chlor.rat.plot



# Save files as pdf
ggsave("Shoot Fresh Weight.pdf", plot = SFW.plot)
ggsave("Total Chlorophyll.pdf", plot = Chlor.total.plot)
ggsave("Chlorophyll A.pdf", plot = Chlor.A.plot)
ggsave("Chlorophyll B.pdf", plot = Chlor.B.plot)
ggsave("Chlorophyl AB ratio.pdf", plot = Chlor.rat.plot)

