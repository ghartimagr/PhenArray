#setting the working directory
getwd()
setwd("/home/mamata/Desktop/project/code/")
#loading the required libraries
library("xlsx")
library("reshape2")
library("lattice")
library("dplyr")
library("tibble")
library("limma")
source("Scaling.R") #contains all my functions for the project
library(contrast)

#STEP 1: PLOTTING THE RAW DATA: BOX PLOTS

#CONTROL PLATE
#loading the table from the excel sheet
pmcontrol = read.xlsx("../Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheetName = "control plate")
head(pmcontrol)
#melt works like stack: stacking LL and HL columns into one for the box plot
newpmcontrol = melt(pmcontrol, id = "position")
head(newpmcontrol)
#using bwplot from lattice library: values in y axis and positions in x axis
#splitting the axes into LL and HL, coloring by the positions
bwplot(newpmcontrol$value~newpmcontrol$position|newpmcontrol$variable, 
       col = as.numeric(as.factor(newpmcontrol$position)), ylab= "Values", xlab= "Position")

#we see the position effect in the HL samples, the values on the left hand side
# of the  plate seem to be higher, it might be because of the position of the light source
#they might have positioned the light in the left hand corner that is why the samples there
#received more light

#PM1 PLATE
#reading the data from the excel sheet
pm1 = read.xlsx("../Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheetName = "PM1")
head(pm1)
pm1reduced <- pm1[,c(4:9)] #extracting the numeric part of the table LLs and HLs
head(pm1reduced)
Tpm1 <- t(pm1reduced) # transposing pm2reduced so that wen can have metabolites as columns
Tpm1<- as.data.frame(Tpm1)
head(Tpm1)
#setting the column names as the metabolite names and rownames as the light conditions
colnames(Tpm1) <- pm1$metabolite
rownames(Tpm1) <- c("LL1", "HL1", "LL2", "HL2", "LL3", "HL3")
head(Tpm1)

#to change rownames into legit column, we use the tibble library command
Tpm1 <- tibble::rownames_to_column(Tpm1, "light_condition")
head(Tpm1)
#stacking all the metabolites into one column and their respective values into another column
#we keep the first column of rownames as it is 
Tpm1stacked <- data.frame(Tpm1[,1], stack(Tpm1[2:97])) 
head(Tpm1stacked)
#setting the column names of our newly stacked column
colnames(Tpm1stacked) <- c("light_condition", "values", "metabolite")
head(Tpm1stacked)
#using my Stacked function from Scaling.R file to set HL1, HL2, HL3 into HL
#and LL1, LL2 and LL3 into LL so that we can group the variables by HL and LL in boxplots
pm1stacked<- Stacked(Tpm1stacked)
head(pm1stacked)
#boxplot(Tpm1stacked$values~Tpm1stacked$sample , col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")
bwplot(pm1stacked$values~pm1stacked$metabolite|pm1stacked$light_condition, col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")

# PM2 PLATE
#reading the PM2A sheet using xlsx library
pm2 = read.xlsx("../Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheetName = "PM2A")
head(pm2)
#again extracting the columns for LLs and HLs
pm2reduced <- pm2[,c(4:9)]
head(pm2reduced)
#transposing the dataframe pm2reduced so that we have metabolites as the columns 
Tpm2 <- t(pm2reduced)
Tpm2<- as.data.frame(Tpm2)
head(Tpm2)
#setting the column names as the metabolite names and row names as the light conditions
colnames(Tpm2) <- pm2$metabolite
rownames(Tpm2) <- c("LL1", "HL1", "LL2", "HL2", "LL3", "HL3")
head(Tpm2)
#to change rownames into legit column
Tpm2 <- tibble::rownames_to_column(Tpm2, "light_condition")
head(Tpm2)
##stacking the categorical and numeric values into 2 separate columns
#the categorical variables are the metabolites and the numeric values are their respctive expression values
Tpm2stacked <- data.frame(Tpm2[,1], stack(Tpm2[2:97])) #2:97 because we 97 columns
head(Tpm2stacked)
#set the column names again
colnames(Tpm2stacked) <- c("light_condition", "values", "metabolite")
#using my Stacked function from Scaling.R file to set HL1, HL2, HL3 into HL
#and LL1, LL2 and LL3 into LL so that we can group the variables by HL and LL in boxplots
pm2stacked <- Stacked(Tpm2stacked)
head(pm2stacked)
#using bwplot from the lattice library
#values of the metabolites are in the y axis and metabolites are in the x-axis
#the values are grouped together by HL and LL
bwplot(pm2stacked$values~pm2stacked$metabolite|pm2stacked$light_condition, col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")

#STEP2 : SCALING AND LOG TRANSFORMATION
#SCALING THE DATA TO REMOVE THE POSITION EFFECT SO THAT ONLY THE METABOLITE AND LIGHT EFFECTS REMAIN
#TO REMOVE THE RIGHT SKEWNESS OF OUR DATA, WE USE LOG TRANSFORMATION


#first we calculate the scaling factor
#Scaling factor of a data point in the HL column is the value divided by the mean of whole HL col 
#and the same holds for the LL col.

#we initialize sftab to store the scaling factors then we iterate through rows and cols of the 
#control table, and divide each value by mean of the whole column
mat = as.matrix(pmcontrol)
sftab = pmcontrol
for (i in 2:ncol(pmcontrol))
{
  for ( j in 1:nrow(pmcontrol))
  {
    sftab[j,i] = pmcontrol[j,i]/mean(pmcontrol[,i])
    sftab = as.data.frame(sftab)
  }
}

head(sftab) #this returns the scaling factor

#SCALING THE CONTROL TABLE

#CtableUpdate is a function to scale the whole table
#each datapoint is divided by the scaling factor of that position
spmcontrol= CTableUpdate(pmcontrol)
head(spmcontrol)
#taking the LL and HL columns (2,3)  and log transforming
logcontrol <- log(spmcontrol[, c(2,3)])
#binding the position column from the original table and the new log transformed table
snewpmcontrol <- cbind(spmcontrol[,1], logcontrol)
colnames(snewpmcontrol) <- c("position", "LL", "HL")
snewpmcontrol = melt(snewpmcontrol, id = "position") #sort of rotating the data
head(snewpmcontrol)
#boxplot for the scaled values gives a straight line for HLs and LLs
bwplot(snewpmcontrol$value~snewpmcontrol$position|snewpmcontrol$variable, 
       col = as.numeric(as.factor(snewpmcontrol$position)), ylab= "Values", xlab= "Position")

#SCALING THE PM1 TABLE
#extracting the scaling factors for LL and HL
LLsf <-  as.vector(sftab[,'LL']) #taking out the low light scaling factor
head(LLsf)
HLsf <- as.vector(sftab[,"HL"]) #taking out the high light scaling factor
head(HLsf)
head(sftab)
#dividing each position with respective scaling factor and binding each column together
# And merging everything together as in the orignal table
scaledpm1 <- cbind(pm1[,c(1,2,3)],pm1[,4]/LLsf, pm1[,5]/HLsf, pm1[,6]/LLsf, pm1[,7]/HLsf, pm1[,8]/LLsf, pm1[,9]/HLsf)
head(scaledpm1) #this is the scaled table
#setting the colnames
colnames(scaledpm1) <- c("strain", "position", "metabolite", "LL1", "HL1", "LL2", "HL2", "LL3", "HL3" )
head(scaledpm1)

#Like before, taking out the LL and HL columns and setting the col and row names
#and transposing the dataframe so that we have metabolites as the columns
#then we log transform in the same step
pm1new <- log(t(scaledpm1[, c(4:9)])) 
#to make sure the class of pm1new remains as dataframe 
pm1new <- as.data.frame(pm1new)
head(pm1new)
colnames(pm1new) <- pm1$metabolite
rownames(pm1new) <- c("LL1", "HL1", "LL2", "HL2", "LL3", "HL3")
head(pm1new)
pm1new <- tibble::rownames_to_column(pm1new, "light_condition")
head(pm1new)
#again Stacking the metabolites together for the boxplots
spm1stacked <- data.frame(pm1new[,1], stack(pm1new[2:97])) 
head(spm1stacked)
colnames(spm1stacked) <- c("light_condition", "values", "metabolite")
head(spm1stacked)
spm1 <- Stacked(spm1stacked)
head(spm1)
#plotting values vs. metabolites grouped by HL and LL
bwplot(spm1$values~spm1$metabolite | spm1$light_condition, col= rainbow(ncol(pm1new)), ylab = "values", xlab = "metabolites")

# SCALING PM2 TABLE

scaledpm2 <- cbind(pm2[,c(1,2,3)],pm2[,4]/LLsf, pm2[,5]/HLsf, pm2[,6]/LLsf, pm2[,7]/HLsf, pm2[,8]/LLsf, pm2[,9]/HLsf)
head(scaledpm2)
colnames(scaledpm2) <- c("strain", "position", "metabolite", "LL1", "HL1", "LL2", "HL2", "LL3", "HL3" )
head(scaledpm2)

#again extracting the HL and LL cols and transposing the data to have metabolites as columns
pm2new <-log(t(scaledpm2[, c(4:9)]))
pm2new <- as.data.frame(pm2new)
head(pm2new)
colnames(pm2new) <- pm2$metabolite
pm2new <- tibble::rownames_to_column(pm2new, "light_condition")
head(pm2new)

spm2stacked <- data.frame(pm2new[,1], stack(pm2new[2:97])) #stacking the categorical and numeric values into 2 cols
head(spm2stacked)
colnames(spm2stacked) <- c("light_condition", "values", "metabolite")
head(spm2stacked)
spm2 <- Stacked(spm2stacked)
head(spm2)
#Values vs. Metabolite grouped by light condition 
bwplot(spm2$values~spm2$metabolite | spm2$light_condition, col= rainbow(ncol(pm1new)), ylab = "values", xlab = "metabolites")

#STEP 3: FITTING THE LINEAR MODELS AFTER SCALING AND TRANSFORMING THE DATA

#LINEAR MODEL FOR CONTROL PLATE
# TRY model with intercept

X <- model.matrix(~variable+position, data = snewpmcontrol)
fitpmc <- lm(value~variable+position, data =snewpmcontrol)
summary(fitpmc)
plot(fitpmc)

#for model without intercept
#forming a group row in the dataframe with position and light conditions (variable) columns combined
snewpmcontrol$group <- factor(paste0(snewpmcontrol$position, snewpmcontrol$variable))
head(snewpmcontrol)
colnames(X)
head(X)
X <- model.matrix(~0+group, data = snewpmcontrol)
fitpmc <- lm(value~0+group, data = snewpmcontrol)
summary(fitpmc)
plot(fitpmc)


# LINEAR MODEL FOR PM1

#this our scaled pm1 table with values , light condition and 
#metabolites stacked in their repective columns and we use it to fit the model now
head(spm1stacked)
head(spm1) 

#grouping the metabolite and light condition together
spm1$group <- factor(paste0(spm1$metabolite, spm1$light_condition))
head(spm1)
Xpm1<-model.matrix(~0+group, data =spm1)
head(Xpm1)

#using default lm function
fitpm1 <- lm(values~0+group, data = spm1)
summarypm1 <- (summary(fitpm1))
coefs <- as.data.frame(summarypm1$coefficients)
coefs <- tibble::rownames_to_column(coefs)
head(coefs)
plot(fitpm1)


#using the limma package
# fitpm1 <- lmFit(spm1$values, Xpm1)
# limma.pm1 <- eBayes(fitpm1)
# pm1.limma.model <-topTable(limma.pm1)
# head(pm1.limma.model)
# results <- decideTests(limma.pm1)
# vennDiagram(results)

#LINEAR MODELS FOR PM2 PLATE
head(spm2) 
#grouping the metabolite and light condition together
spm2$group <- factor(paste0(spm2$metabolite, spm2$light_condition))
head(spm2)
Xpm2<-model.matrix(~0+group, data =spm2)
head(Xpm1)
fitpm2 <- lm(values~0+group, data = spm2)
summarypm2 <- summary(fitpm2)
plot(fitpm2)
coefs2 <- as.data.frame(summarypm2$coefficients)
coefs2 <- tibble::rownames_to_column(coefs2)
head(coefs2)


#STEP 4 : CONTRAST MATRICES TO GET THE FOLLOWING:
# LLcontrol - LLmetabolite
# HLcontrol - HLmetabolite
# (LLcontrol - LLmetabolite) - (HLcontrol - HLmetabolites)
 
#plate pm1
#All of our HL metabolite combos are in odd rows
#All of our LL metabolites are in the even rows

coefs[c(165,166), ] # the columns for te negative controls
contrast1 <- Contrast(coefs, LLc=166, HLc=165)
head(contrast1)
#putting back the metabolites names in the contrast- matrix
contrastbind = cbind(coefs[,1], contrast)
head(contrastbind)
#diffcontrast gives us (LLcontrol - LLmetabolite) - (HLcontrol - HLmetabolites) value
difcontrast1 <- Difference(contrast1)
head(difcontrast1)


#plate pm2
coefs2[c(161,162),] #te columns for te negative controls
contrast2 <- Contrast(coefs2, LLc=162, HLc=161)
head(contrast2)
#putting back the metabolites names in the contrast matrix
contrastbind2 = cbind(coefs2[,1], contrast2)
head(contrastbind2)
difcontrast2 <- Difference(contrast2)
head(difcontrast2)












