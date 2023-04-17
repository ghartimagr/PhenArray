setwd("/home/mamata/Desktop/project/code/")
#setwd("/media/mamata/3306c0a9-648d-4863-ab8f-6d482a3f647a/home/mamata/Desktop/project/code/")
library("readxl")
library("reshape2")
library("lattice")
library("dplyr")
library("tibble")
library("limma")
source("Scaling.R") #contains all my functions for the project
library("multcomp")
library("ggplot2")
library("ggrepel")
library("stringr")
library ("caret")
library("PRROC")
library("pROC")
library("yardstick")

#STEP 1: PLOTTING THE RAW DATA: BOX PLOTS
#CONTROL PLATE
pmcontrol = read_excel("/home/mamata/Desktop/project/code/Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheet = "control plate")
head(pmcontrol)
#melt works like stack
newpmcontrol = melt(pmcontrol, id = "position")
head(newpmcontrol)
save(newpmcontrol, file = "newpmcontrol.RData")
bwplot(newpmcontrol$value~newpmcontrol$position|newpmcontrol$variable, 
       col = as.numeric(as.factor(newpmcontrol$position)), ylab= "luminescence values", 
       xlab= "Position", scales=list(x=list(draw=FALSE)))

#we see the position effect in the samples, the values on the left hand side
# of the  plate seem to be higher, it might be because of the position of the light source
#they might have positioned the light in the left hand corner that is why the samples there received more light
#this is the position effect, because the samples show different expression values although there is
#the same condition of light , and there is no metabolite present.

#PM1 PLATE
#reading the data from the excel sheet
pm1 = read_excel("Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheet = "PM1")
head(pm1)
#convert metabolite names to conform R variable name syntax
pm1$metabolite=conv_metnames(pm1$metabolite)
pm1reduced <- pm1[,c(4:9)] #extracting the numeric part of the table LLs and HLs
head(pm1reduced)
Tpm1 <- t(pm1reduced) # transposing pm2reduced so that wen can have metabolites as columns
Tpm1<- as.data.frame(Tpm1)
head(Tpm1)
#setting the column names as the metabolite names and rownames as the light conditions
colnames(Tpm1) <- pm1$metabolite
rownames(Tpm1) <- c("LL1", "HL1", "LL2", "HL2", "LL3", "HL3")
head(Tpm1)
#to change rownames into legit column
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
#pm1stacked<- Stacked(Tpm1stacked)
#head(pm1stacked)
#save(pm1stacked, Tpm1, file ="pm1stacked.RData")
Tpm1stacked$light_condition <- gsub("[0-9]", "", Tpm1stacked$light_condition)
#boxplot(Tpm1stacked$values~Tpm1stacked$sample , col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")
bwplot(Tpm1stacked$values~Tpm1stacked$metabolite|Tpm1stacked$light_condition, 
       col= rainbow(ncol(Tpm1)), ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

# PM2 PLATE
#reading the PM2A sheet using xlsx library
pm2 = read_excel("Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheet = "PM2A")
head(pm2)
pm2$metabolite=conv_metnames(pm2$metabolite)
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
#pm2stacked <- Stacked(Tpm2stacked)
#head(pm2stacked)
#save(pm2stacked, Tpm2, file ="pm2stacked.RData" )
#using bwplot from the lattice library
#values of the metabolites are in the y axis and metabolites are in the x-axis
#the values are grouped together by HL and LL
#pm2stacked$values <- format(pm2stacked$values, scientific = TRUE)
Tpm2stacked$light_condition <- gsub("[0-9]", "", Tpm2stacked$light_condition)

bwplot(Tpm2stacked$values~Tpm2stacked$metabolite|Tpm2stacked$light_condition, 
       col= rainbow(ncol(Tpm2)), ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE), y = list(tick.number =4)))

#STEP2 : SCALING AND LOG TRANSFORMATION
#SCALING THE DATA TO REMOVE THE POSITION EFFECT SO THAT ONLY THE METABOLITE AND LIGHT EFFECTS REMAIN
#TO REMOVE THE RIGHT SKEWNESS OF OUR DATA, WE USE LOG TRANSFORMATION
#first we calculate the scaling factor
#Scaling factor for HL :  data point in the HL column is  divided by the mean of whole HL column and the same holds for the LL col.
#we initialize sftab to store the scaling factors then we iterate through rows and cols of the 
#control table, and divide each value by mean of the whole column
# mat = as.matrix(pmcontrol)
# sftab = as.matrix(pmcontrol)
# 
# for (i in 2:ncol(pmcontrol))
# {
#   for ( j in 1:nrow(pmcontrol))
#   {
#     sftab[j,i] = pmcontrol[j,i]/mean(pmcontrol[,i])
#     sftab = as.data.frame(sftab)
#    
#   }
# }
sftab <-as.data.frame(apply(pmcontrol[,2:ncol(pmcontrol)],2,function(x){x/sum(x)}))
sftab <- cbind(pmcontrol[,1], sftab)
#SCALING THE CONTROL TABLE
#CtableUpdate is a function to scale the whole table
#each datapoint is divided by the scaling factor of that position
spmcontrol= CTableUpdate(pmcontrol)
head(spmcontrol)
#taking the LL and HL columns (2,3)  and log transforming
logcontrol <- log2(spmcontrol[, c(2,3)])
#binding the position column from the original table and the new log transformed table
snewpmcontrol <- cbind(spmcontrol[,1], logcontrol)
colnames(snewpmcontrol) <- c("position", "LL", "HL")
snewpmcontrol = melt(snewpmcontrol, id = "position") #sort of rotating the data
head(snewpmcontrol)

#SAVING THE DATA FOR PLOTTING CONTROL DATA BEFORE AND AFTER POSITION CORRECTION
save(snewpmcontrol, file = "snewpmcontrol.RData")

control.ranked <- data.frame(cbind(rank(newpmcontrol$value, ties.method = 'first'),
                                rank(snewpmcontrol$value, ties.method = 'first')))
corcontrol <- cor(control.ranked$X1, control.ranked$X2, method = "spearman")
corcontrol
save(corcontrol, file = "cor.RData")

#boxplot for the scaled values gives a straight line for HLs and LLs
bwplot(snewpmcontrol$value~snewpmcontrol$position|snewpmcontrol$variable, 
       col = as.numeric(as.factor(snewpmcontrol$position)), ylab= "luminescence values", xlab= "Position", scales=list(x=list(draw=FALSE)))#
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

#position correction
pm1new <- t(scaledpm1[, c(4:9)])
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
spm1_notlog <- spm1stacked
#plotting values vs. metabolites grouped by HL and LL
spm1_notlog$light_condition <- gsub("[0-9]", "", spm1_notlog$light_condition)

bwplot(spm1_notlog$values~spm1_notlog$metabolite | spm1_notlog$light_condition, col= rainbow(ncol(pm1new)), ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))
#for spearman'S corelation
pm1ranked <- transform(Tpm1stacked,rank=ave(1:nrow(Tpm1stacked),metabolite,
                                           FUN=function(x) order(values[x],decreasing=TRUE)))

spm1ranked <- transform(spm1_notlog,rank=ave(1:nrow(spm1_notlog),metabolite,
                                             FUN=function(x) order(values[x],decreasing=TRUE)))
corpm1 <- cor(pm1ranked$rank, spm1ranked$rank, method = "spearman")
corpm1

#SAVING THE DATA FOR PLOTTING PM1 BEFORE AND AFTER POSITION CORRECTION
#save(pm1stacked, Tpm1, spm1_notlog, pm1new, corpm1, file= "spm1.RData")

# SCALING PM2 TABLE
#Like before, taking out the LL and HL columns and setting the col and row names
#and transposing the dataframe so that we have metabolites as the columns
#then we log transform in the same step
pm1new <- log2(t(scaledpm1[, c(4:9)]))


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
spm1 <- spm1stacked
#head(spm1)
#save(spm1, pm1new, file = "spm1log.RData")
spm1$light_condition <- gsub("[0-9]", "", spm1$light_condition)

#save(pm1new, file = "pm1new.RData")
#plotting values vs. metabolites grouped by HL and LL
bwplot(spm1$values~spm1$metabolite | spm1$light_condition, col= rainbow(ncol(pm1new)), ylab = "log2 luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

#to check the heteroscedacity
#Heteroscedasticity means unequal scatter. 
#In regression analysis, we talk about heteroscedasticity in the context of the residuals or error term
# To check for heteroscedasticity, you need to assess the residuals by fitted value plots specifically.
#Typically, the telltale pattern for heteroscedasticity is that as the fitted values increases, 
#the variance of the residuals also increases.
#As I mentioned earlier, linear regression assumes that the spread of the residuals is constant across the plot. 
#Anytime that you violate an assumption, there is a chance that you can’t trust the statistical results.

#plot(lm(values~light_condition*metabolite, data=spm1_notlog))
#plot(lm(values~light_condition*metabolite, data=spm1))
# SCALING PM2 TABLE

#POSITION CORRECTION FOR PM2
scaledpm2 <- cbind(pm2[,c(1,2,3)],pm2[,4]/LLsf, pm2[,5]/HLsf, pm2[,6]/LLsf, pm2[,7]/HLsf, pm2[,8]/LLsf, pm2[,9]/HLsf)
head(scaledpm2)
colnames(scaledpm2) <- c("strain", "position", "metabolite", "LL1", "HL1", "LL2", "HL2", "LL3", "HL3" )
head(scaledpm2)
#again extracting the HL and LL cols and transposing the data to have metabolites as columns
pm2new <-t(scaledpm2[, c(4:9)])
pm2new <- as.data.frame(pm2new)
head(pm2new)
colnames(pm2new) <- pm2$metabolite
pm2new <- tibble::rownames_to_column(pm2new, "light_condition")
head(pm2new)
spm2stacked <- data.frame(pm2new[,1], stack(pm2new[2:97])) #stacking the categorical and numeric values into 2 cols
head(spm2stacked)
colnames(spm2stacked) <- c("light_condition", "values", "metabolite")
head(spm2stacked)
spm2_notlog <- spm2stacked
head(spm2_notlog)

spm2_notlog$light_condition <- gsub("[0-9]", "", spm2_notlog$light_condition)

#Values vs. Metabolite grouped by light condition 
bwplot(spm2_notlog$values~spm2_notlog$metabolite | spm2_notlog$light_condition, col= rainbow(ncol(pm2new)), 
       ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))


scaledpm2 <- cbind(pm2[,c(1,2,3)],pm2[,4]/LLsf, pm2[,5]/HLsf, pm2[,6]/LLsf, pm2[,7]/HLsf, pm2[,8]/LLsf, pm2[,9]/HLsf)
head(scaledpm2)
colnames(scaledpm2) <- c("strain", "position", "metabolite", "LL1", "HL1", "LL2", "HL2", "LL3", "HL3" )
head(scaledpm2)
#again extracting the HL and LL cols and transposing the data to have metabolites as columns
pm2new <-log2(t(scaledpm2[, c(4:9)]))
pm2new <- as.data.frame(pm2new)
head(pm2new)
colnames(pm2new) <- pm2$metabolite
pm2new <- tibble::rownames_to_column(pm2new, "light_condition")
head(pm2new)
spm2stacked <- data.frame(pm2new[,1], stack(pm2new[2:97])) #stacking the categorical and numeric values into 2 cols
head(spm2stacked)
colnames(spm2stacked) <- c("light_condition", "values", "metabolite")
head(spm2stacked)
spm2 <- spm2stacked
head(spm2)
spm2$light_condition <- gsub("[0-9]", "", spm2$light_condition)

#Values vs. Metabolite grouped by light condition 
bwplot(spm2$values~spm2$metabolite | spm2$light_condition, col= rainbow(ncol(pm2new)), 
       ylab = "log2 luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

#for spearman'S corelation
pm2ranked <- transform(Tpm2stacked,rank=ave(1:nrow(Tpm2stacked),metabolite,
                                           FUN=function(x) order(values[x],decreasing=TRUE)))
spm2ranked <- transform(spm2_notlog,rank=ave(1:nrow(spm2_notlog),metabolite,
                                             FUN=function(x) order(values[x],decreasing=TRUE)))
corpm2 <- cor(pm2ranked$rank, spm2ranked$rank, method = "spearman")
corpm2

#plot(lm(values~light_condition*metabolite, data=spm2_notlog))
#plot(lm(values~light_condition*metabolite, data=spm2))
#save(pm2stacked, Tpm2, spm2_notlog, pm2new, spm2, corpm2,  file = "spm2.RData")

#STEP 3: FITTING THE LINEAR MODELS AFTER SCALING AND TRANSFORMING THE DATA
spm1$group <- factor(paste0(spm1$metabolite, spm1$light_condition))
spm2$group <- factor(paste0(spm2$metabolite, spm2$light_condition))
#combining spm1 and spm2
#data <- rbind(spm1, spm2)
#head(data)
#building design matrix for the model
Design1<-model.matrix(~0+group, data =spm1)
colnames(Design1)=gsub("^group", "", colnames(Design1))
col <- gsub("group", "", colnames(Design1))
#rownames(Design1) <- spm1$group
################################################################################
# CODE for presentation for dummy design and contrast matrices
which(Design1[,1]==1) #194 196 198 
spm1[c(194, 196, 198),]
which(grepl("MaltoseHL", rownames(Design1))==TRUE)
which(grepl("Negative.ControlHL", rownames(Design1))==TRUE)
#[1] 200 202 204
spm1[c(200, 202, 204),]
#So desin matrix gives a score of 1 whenever it finds the sample in the rows
#ie. matrix gets 1 score if rows and columns match
dummyspm1 <-as.data.frame(spm1[c(5:6,203:204),])
dummyDesign <- model.matrix(~0+group, data = dummyspm1)
colnames(dummyDesign)=gsub("^group", "", colnames(dummyDesign))
rownames(dummyDesign) <- dummyspm1$group
dummyDesign
# colnames(dummyDesign)[colSums(dummyDesign)!=0]
# dummyDesign<- dummyDesign[dummyDesign != "0"]
#example
contpm=makeContrasts(contrasts =c("MaltoseHL-Negative.ControlHL", 
                                  "MaltoseLL-Negative.ControlLL", 
                                  "(MaltoseHL-Negative.ControlHL)-(MaltoseLL-Negative.ControlLL)"),
                                  levels = Design1)
#this is the end of code for dummy design and contrast matrices for presentation
################################################################################

vecLL <- col[which(grepl("LL$", col))]

vecHL <- col[which(grepl("HL$", col))]
vecLL <- paste(vecLL, "Negative.ControlLL", sep = "-")
vecHL <- paste(vecHL, "Negative.ControlHL", sep = "-" )
vecLL = setdiff(vecLL, c("Negative.ControlLL-Negative.ControlLL"))
vecHL = setdiff(vecHL, c("Negative.ControlHL-Negative.ControlHL"))
metabolitevecLL <- paste( "(", vecLL, ")")
metabolitevecHL <- paste( "(", vecHL, ")")
metabolitevec <-  paste(metabolitevecHL, metabolitevecLL, sep ="-")
contrast1=makeContrasts(contrasts =c(vecLL, vecHL, metabolitevec), levels=Design1)
model1<- lm(values~0+group, data = spm1)

mod1Summary <- summary(model1)
mod1coef <- as.data.frame(mod1Summary$coefficients)

#hypothesis testing
confit1=glht(model1, t(contrast1))
summary1= summary(confit1, test = adjusted(type = "fdr"))
hist(summary1$test$pvalues)
pval =as.data.frame(summary1$test$pvalues) #pvalues
lfc= as.data.frame(summary1$test$coefficients)

lm1 <- as.data.frame(cbind(lfc, pval))
lm1<- tibble::rownames_to_column(lm1, "contrasts") 
colnames(lm1) = c( "contrasts", "lfc", "p.value")
head(lm1)
p <- ggplot(data=lm1, aes(x=lfc, y=-log10(p.value))) + geom_point() + theme_minimal()
p
lm1$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
lm1$diffexpressed[lm1$lfc > 0.6 & lm1$p.value< 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
lm1$diffexpressed[lm1$lfc < -0.6 & lm1$p.value< 0.05] <- "DOWN"
p <- ggplot(data=lm1, aes(x=lfc, y=-log10(p.value), col=diffexpressed)) + geom_point() + theme_minimal()
p

#checking pvalues with limma
spm1limma=data.frame(t(spm1$values))
colnames(spm1limma)=spm1$group
#generate canonical linear model
limfit = lmFit(spm1limma, Design1)
#extract contrasts
limcont = contrasts.fit(limfit, contrast1)
#calculate empirical bayes moderated t statistics
eb=eBayes(limcont)
limmaLFC=eb$coefficients
#pval = eb$p.value
#Benjamini Hochberg correction 
limmapval <- p.adjust(eb$p.value,method="fdr")
#check for NAS
p= as.data.frame(limmapval)
lf = as.data.frame(limmaLFC)
any(is.na(p))
which(is.na(p)==TRUE) #taking out the positions of NAs
any(is.na(lf))
#find out which components are NAs
lf[,c(83, 179, 275)]
#this can be removed if i remove these comparisons in lines 278-279 by using set diff
#making a dataframe of limmaLFC ie. estimate values and limmapval ie. p values for volcano plot
limmaLFC = t(limmaLFC)
dfpm1 <- as.data.frame(cbind(limmaLFC, limmapval))
dfpm1 <- tibble::rownames_to_column(dfpm1, "contrasts") 
#dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(dfpm1) = c( "contrasts", "limmaLFC", "limmapval")
head(dfpm1)
save(dfpm1, file = "dfpm1.RData")
# Volcano plots indicate the fold change (either positive or negative) in the x axis 
#and a significance value (such as the p-value or the adjusted p-value, i.e. limmapval) in the y axis
#The ‘limmapval’ columns contains the corrected pvalues; 
#these must be converted to the negative of their logarithm base 10 before plotting, 
#i.e -log10(p-value) or -log10(limmapval).
#Since volacno plots are scatter plots, we can use geom_point() to generate one with ggplot2
#the higher the position of a point, the more significant its value is (y axis).
#Points with positive fold change values (to the right) are up-regulated and 
#points with negative fold change values (to the left) are down-regulated (x axis).
p <- ggplot(data=dfpm1, aes(x=limmaLFC, y=-log10(limmapval))) + geom_point() + theme_minimal()
p
dfpm1$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
dfpm1$diffexpressed[dfpm1$limmaLFC > 0.6 & dfpm1$limmapval < 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
dfpm1$diffexpressed[dfpm1$limmaLFC < -0.6 & dfpm1$limmapval < 0.05] <- "DOWN"
p <- ggplot(data=dfpm1, aes(x=limmaLFC, y=-log10(limmapval), col=diffexpressed)) + geom_point() + theme_minimal()
p

#make heatmaps of the glht function
lfc1 =as.data.frame(cbind(summary1$test$coefficients, summary1$test$pvalues))
lfc1 <- tibble::rownames_to_column(lfc1, "contrasts")
colnames(lfc1) = c("contrasts", "lfc", "pval" )
head(lfc1)
for ( i in 1:nrow(lfc1))
{
  if (lfc1$pval[i] >= 0.05)
  {
    lfc1$pval[i]= NA
    lfc1$lfc[i] = NA
  }
}
#before kicking out the groups we have 285 entries in the groups, and wach contrast column
# ie. LL, HL and HL-LL contrasts had 95 entries but ow we have 260 entries in the group
# and the no. of groups in each contrast is not the same 
# so we cannot use this: heat_pm2 = as.data.frame(cbind(lfcpm2[1:95,], lfcpm2[96:190,], lfcpm2[191:285,]))
#loop to find out how many LL, HL and HL-LL contrast groups are there after removal of certain groups

idxLLgroups <- which(grepl("LL-Negative.ControlLL$", lfc1$contrasts)==TRUE)
LLgroups <- lfc1[idxLLgroups,]
idxHLgroups <- which(grepl("HL-Negative.ControlHL$", lfc1$contrasts)==TRUE)
HLgroups <- lfc1[idxHLgroups,]
HL_LL_diff <-lfc1[which(grepl("^\\(.+\\)$", lfc1$contrasts)==TRUE),]
heat1 = as.data.frame(cbind(LLgroups, HLgroups, HL_LL_diff))
metNames<-gsub("LL-.+", "", heat1[,1])
heat1 <- heat1[,c(2,5,8)]

rownames(heat1) <- metNames
heat1<- heat1[rowSums(is.na(heat1)) != ncol(heat1), ]
colnames(heat1) <- c("LL contrast", "HL contrast", "HL-LL contrast")
heat1Zeros <- heat1
heat1Zeros#s[is.na(heat1Zeros)] = 0

library("superheat")
superheat(heat1,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "white", "red"),
          left.label.text.size = 5, heat.lim = c(-7, 6), heat.pal.values = c(0,0.5,1))


#data2 <- rbind(spm1, spm2_modified)
Design2<-model.matrix(~0+group, data =spm2)
colnames(Design2)=gsub("^group", "", colnames(Design2))
col <- gsub("group", "", colnames(Design2))
vecLL <- col[which(grepl("LL$", col))]
vecHL <- col[which(grepl("HL$", col))]
vecLL <- paste(vecLL, "Negative.ControlLL", sep = "-")
vecHL <- paste(vecHL, "Negative.ControlHL", sep = "-" )
vecLL = setdiff(vecLL, c("Negative.ControlLL-Negative.ControlLL"))
vecHL = setdiff(vecHL, c("Negative.ControlHL-Negative.ControlHL"))
metabolitevecLL <- paste( "(", vecLL, ")")
metabolitevecHL <- paste( "(", vecHL, ")")
metabolitevec <-  paste(metabolitevecHL, metabolitevecLL, sep ="-")
contrast2=makeContrasts(contrasts =c(vecLL, vecHL, metabolitevec), levels=Design2)
model2<- lm(values~0+group, data = spm2)

#hypothesis testing
confit2=glht(model2, t(contrast2))
summary2= summary(confit2, test = adjusted(type = "fdr"))
pval =as.data.frame(summary2$test$pvalues) #pvalues
lfc= as.data.frame(summary2$test$coefficients)

lm2 <- as.data.frame(cbind(lfc, pval))
lm2<- tibble::rownames_to_column(lm2, "contrasts") 
#dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(lm2) = c( "contrasts", "lfc", "p.value")
head(lm2)
p <- ggplot(data=lm2, aes(x=lfc, y=-log10(p.value))) + geom_point() + theme_minimal()
p
lm2$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
lm2$diffexpressed[lm2$lfc > 0.6 & lm2$p.value< 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
lm2$diffexpressed[lm2$lfc < -0.6 & lm2$p.value < 0.05] <- "DOWN"
p <- ggplot(data=lm2, aes(x=lfc, y=-log2(p.value), col=diffexpressed)) + geom_point() + theme_minimal()
p

#limma implementation
spm2limma=data.frame(t(spm2$values))
colnames(spm2limma)=spm2$group
#generate canonical linear model
limfit = lmFit(spm2limma,Design2)
#extract contrasts
limcont = contrasts.fit(limfit, contrast2)
#calculate empirical bayes moderated t statistics
eb=eBayes(limcont)
limmaLFC=eb$coefficients
#pval = eb$p.value
#Benjamini Hochberg correction 
limmapval <- p.adjust(eb$p.value,method="fdr")
#making a dataframe of limmaLFC ie. estimate values and limmapval ie. p values for volcano plot
limmaLFC = t(limmaLFC)
dfpm2 <- as.data.frame(cbind(limmaLFC, limmapval))
dfpm2 <- tibble::rownames_to_column(dfpm2, "contrasts")     
#dfpm2 <- tibble::rownames_to_column(dfpm2, "Numbers") 
colnames(dfpm2) = c("contrasts", "limmaLFC", "limmapval")
head(dfpm2)
p <- ggplot(data=dfpm2, aes(x=limmaLFC, y=-log10(limmapval))) + geom_point() + theme_minimal()
p
dfpm2$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
dfpm2$diffexpressed[dfpm1$limmaLFC > 0.6 & dfpm2$limmapval < 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
dfpm2$diffexpressed[dfpm2$limmaLFC < -0.6 & dfpm2$limmapval < 0.05] <- "DOWN"
save(dfpm2, file = "dfpm2.RData")
p <- ggplot(data=dfpm2, aes(x=limmaLFC, y=-log10(limmapval), col=diffexpressed)) + geom_point() + theme_minimal()
p


lfc2 =as.data.frame(cbind(summary2$test$coefficients, summary2$test$pvalues))
lfc2 <- tibble::rownames_to_column(lfc2, "contrasts")
colnames(lfc2) = c("contrasts", "lfc", "pval" )
head(lfc2)
for ( i in 1:nrow(lfc2))
{
  if (lfc2$pval[i] >= 0.05)
  {
    lfc2$pval[i]= NA
    lfc2$lfc[i] = NA
  }
}

idxLLgroups <- which(grepl("LL-Negative.ControlLL$", lfc2$contrasts)==TRUE)
LLgroups <- lfc2[idxLLgroups,]
idxHLgroups <- which(grepl("HL-Negative.ControlHL$", lfc2$contrasts)==TRUE)
HLgroups <- lfc2[idxHLgroups,]
HL_LL_diff <-lfc2[which(grepl("^\\(.+\\)$", lfc2$contrasts)==TRUE),]
heat2 = as.data.frame(cbind(LLgroups, HLgroups, HL_LL_diff))
metNames<-gsub("LL-.+", "", heat2[,1])
heat2 <- heat2[,c(2,5,8)]
rownames(heat2) <- metNames
heat2<- heat2[rowSums(is.na(heat2)) != ncol(heat2), ]
colnames(heat2) <- c("LL contrast", "HL contrast", "HL-LL contrast")
library("superheat")
superheat(heat2,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "white", "red"), 
          left.label.text.size = 5, heat.lim = c(-9,4), heat.pal.values = c(0,0.7,1))

#checking variances for anova
pm12var=data.frame(pm1var=aggregate(values~group, data=spm1, var), pm2var=aggregate(values~group, data=spm2, var))
#barlett.test null hypothesis is that the variances are equal for all
#test if sample variances are equal for pm1
pm1vartest=bartlett.test(values~group, data=spm1)
# p-value = 0.9767, means that the p value is insignificant, we accept the null hypothesis
#that the variances for the groups in pm1 dataset are equal
#test if sample variances are equal for pm2
pm2vartest=bartlett.test(values~group, data=spm2)
#p-value = 2.308e-12, we have to reject the null hypothesis, the variances are not equal
#plot boxplot of sample variances with test results
#plotting varinces of pm1 sample and pm2 samples as a boxplot
#each point in the boxplot is a group
#par(mfrow=c(1,0))
save(pm12var, pm1vartest, pm2vartest, file = "pm12var.RData")
limits=boxplot(pm12var[, c("pm1var.values", "pm2var.values")], col ="bisque", 
               names=c("PM1", "PM2"), ylab="Var", ylim=c(0,max(c(pm12var$pm1var.values, pm12var$pm2var.values))*1.2))
ypos= y=apply(pm12var[, c("pm1var.values", "pm2var.values")], 2, max)
text( 
  x=c(1:2), 
  y=ypos + 0.1*max(ypos), 
  paste("p = ",signif(c(pm1vartest$p.value, pm2vartest$p.value),3),sep="")
  #shows p values in the boxplot
)
#how to deal with the variance?
# Remove the outliers ie. the groups with the high variances one by one until we have non significant variances
#in our case we will remove all the variances above 4 in the figure
# then we rerun everything, means we have to now adjust the contrast matrix

outliers <- boxplot(pm12var$pm2var.values, plot=FALSE)$out #taking out the outliers
outliers
idx_outliers <- which(pm12var$pm2var.values %in% outliers) #which ones are the outliers
pm2_outliers <- pm12var[idx_outliers,]
pm2_outliers <- as.data.frame(pm2_outliers[, c(3,4)])
pm2_outliers
# ordering the outliers from largest values to lowest, so that the largest ones get removed first
pm2_outliers <- pm2_outliers[order(pm2_outliers$pm2var.values, decreasing = TRUE),] #ordering 
pm2_outliers

#new trial for removing the high variance groups
# look in the old code from Modified.R for older version of the removal of groups
var_mets <- pm2_outliers$pm2var.group
var_mets <-gsub ("LL$", "", var_mets)
var_mets <-gsub ("HL$", "", var_mets)
var_mets
spm2after <- spm2

for ( i in 1:length(var_mets))
{
  barlet <- bartlett.test(values~group, data=spm2after)
  pvalue <- barlet$p.value
  if (pvalue <= 0.05)
  {
    #idx= grep(paste("^", var_mets[i], "$", sep = ""), spm2after$metabolite)
    idx <- grep(var_mets[i], spm2after$metabolite)
    expel <-  spm2after[idx,]
    spm2after <- setdiff(spm2after, expel)
    barlet<- bartlett.test(values~group, data=spm2after)
    pvalue <- barlet$p.value
    print(pvalue)
    print(i)
    i=i+1
  }
  else
  {
    break
  }
}

pm2var=aggregate(values~group, data=spm2after, var)
pm2vartest=bartlett.test(values~group, data=spm2after)
pm2vartest$p.value
setdiff(spm2$metabolite, spm2after$metabolite) #expelled
dim(spm2after)
dim(spm2)

boxplot(pm2var[, c("values")], xlab =paste("p = ",signif(c(pm2vartest$p.value),3)), 
        main = "PM2", labels = paste("p = ",signif(c(pm2vartest$p.value),3)), col = 
          "bisque", ylim=c(0,max(c(pm2var$values))*1.2))

write.csv(spm2after, file = "spm2after.csv")
spm2after <-read.csv("/home/mamata/Desktop/project/code/spm2after.csv")

Design2<-model.matrix(~0+group, data =spm2after)
colnames(Design2)=gsub("^group", "", colnames(Design2))
col <- gsub("group", "", colnames(Design2))

vecLL <- col[which(grepl("LL$", col))]
vecHL <- col[which(grepl("HL$", col))]
vecLL <- paste(vecLL, "Negative.ControlLL", sep = "-")
vecHL <- paste(vecHL, "Negative.ControlHL", sep = "-" )
vecLL = setdiff(vecLL, c("Negative.ControlLL-Negative.ControlLL"))
vecHL = setdiff(vecHL, c("Negative.ControlHL-Negative.ControlHL"))
metabolitevecLL <- paste( "(", vecLL, ")")
metabolitevecHL <- paste( "(", vecHL, ")")
metabolitevec <-  paste(metabolitevecHL, metabolitevecLL, sep ="-")
contrast2=makeContrasts(contrasts =c(vecLL, vecHL, metabolitevec), levels=Design2)
model2<- lm(values~0+group, data = spm2after)
#hypothesis testing
confit2=glht(model2, t(contrast2))
# summary2= summary(confit2, test = adjusted(type = "fdr"))
summarypm2= summary(confit2, test = adjusted(type = "fdr"))
pval =as.data.frame(summarypm2$test$pvalues) #pvalues
lfc= as.data.frame(summarypm2$test$coefficients)

#limma implementation
spm2limma=data.frame(t(spm2after$values))
colnames(spm2limma)=spm2after$group
#generate canonical linear model
limfit = lmFit(spm2limma,Design2)
#extract contrasts
limcont = contrasts.fit(limfit, contrast2)
#calculate empirical bayes moderated t statistics
eb=eBayes(limcont)
limmaLFC=eb$coefficients
#pval = eb$p.value
#Benjamini Hochberg correction 
limmapval <- p.adjust(eb$p.value,method="fdr")
#making a dataframe of limmaLFC ie. estimate values and limmapval ie. p values for volcano plot
limmaLFC = t(limmaLFC)
dfpm2 <- as.data.frame(cbind(limmaLFC, limmapval))
dfpm2 <- tibble::rownames_to_column(dfpm2, "contrasts")     
#dfpm2 <- tibble::rownames_to_column(dfpm2, "Numbers") 
colnames(dfpm2) = c("contrasts", "limmaLFC", "limmapval")
head(dfpm2)
p <- ggplot(data=dfpm2, aes(x=limmaLFC, y=-log10(limmapval))) + geom_point() + theme_minimal()
p
dfpm2$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
dfpm2$diffexpressed[dfpm1$limmaLFC > 0.6 & dfpm2$limmapval < 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
dfpm2$diffexpressed[dfpm2$limmaLFC < -0.6 & dfpm2$limmapval < 0.05] <- "DOWN"
save(dfpm2, file = "dfpm2.RData")
p <- ggplot(data=dfpm2, aes(x=limmaLFC, y=-log10(limmapval), col=diffexpressed)) + geom_point() + theme_minimal()
p

#glht implementation 
pval =as.data.frame(summarypm2$test$pvalues) #pvalues
lfc= as.data.frame(summarypm2$test$coefficients) #coefficients
lmpm2 <- as.data.frame(cbind(lfc, pval))
lmpm2<- tibble::rownames_to_column(lmpm2, "contrasts") 
#dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(lmpm2) = c( "contrasts", "lmlfc", "lmpval")
head(lmpm2)
p <- ggplot(data=lmpm2, aes(x=lmlfc, y=-log2(lmpval))) + geom_point() + theme_minimal()
p
lmpm2$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
lmpm2$diffexpressed[lmpm2$lmlfc > 0.6 & lmpm2$lmpval< 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
lmpm2$diffexpressed[lmpm2$lmlfc < -0.6 & lmpm2$lmpval < 0.05] <- "DOWN"

save(lmpm2, file =  "lmpm2after.RData")
p <- ggplot(data=lmpm2, aes(x=lmlfc, y=-log2(lmpval), col=diffexpressed)) + geom_point() + theme_minimal()
p
lfcpm2 =as.data.frame(cbind(summarypm2$test$coefficients, summarypm2$test$pvalues))
lfcpm2 <- tibble::rownames_to_column(lfcpm2, "contrasts")
colnames(lfcpm2) = c("contrasts", "lfc", "pval" )
head(lfcpm2)
for ( i in 1:nrow(lfcpm2))
{
  if (lfcpm2$pval[i] >= 0.05)
  {
    lfcpm2$pval[i]= NA
    lfcpm2$lfc[i] = NA
  }
}
dim(lfcpm2) # 246instead of 285
#before kicking out the groups we have 285 entries in the groups, and wach contrast column
# ie. LL, HL and HL-LL contrasts had 95 entries but ow we have 260 entries in the group
# and the no. of groups in each contrast is not the same 
# so we cannot use this: heat_pm2 = as.data.frame(cbind(lfcpm2[1:95,], lfcpm2[96:190,], lfcpm2[191:285,]))
#loop to find out how many LL, HL and HL-LL contrast groups are there after removal of certain groups

idxLLgroups <- which(grepl("LL-Negative.ControlLL$", lfcpm2$contrasts)==TRUE)
LLgroups <- lfcpm2[idxLLgroups,]
idxHLgroups <- which(grepl("HL-Negative.ControlHL$", lfcpm2$contrasts)==TRUE)
HLgroups <- lfcpm2[idxHLgroups,]
HL_LL_diff <-lfcpm2[which(grepl("^\\(.+\\)$", lfcpm2$contrasts)==TRUE),]

heat_pm2 = as.data.frame(cbind(LLgroups, HLgroups, HL_LL_diff))
expelledgroups <- setdiff( spm2$group, spm2after$group)
metNamespm2 <-gsub("LL-.+", "", heat_pm2[,1])
heat_pm2 <- heat_pm2[,c(2,5,8)]
rownames(heat_pm2) <- metNamespm2
colnames(heat_pm2) <- c("LL contrast", "HL contrast", "HL-LL contrast")

heat_pm2<- heat_pm2[rowSums(is.na(heat_pm2)) != ncol(heat_pm2), ]
library("superheat")
superheat(heat_pm2,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "white", "red"), left.label.text.size = 2.5, heat.lim = c(-9,5))

dim(heat2) # 36
dim(heat_pm2) # 58
#so after the removal of the high variance, more metabolites turned significant 
# in pm2 dataset. ie., 58-36 = 22 metabolites

#calculating jaccards coefficient
head(pm1new)
head(spm1)
pm1_anova_0 <- aov(values~light_condition*metabolite, data = spm1) 
summary(pm1_anova_0)
thsd0=TukeyHSD(pm1_anova_0)
group_cont <-as.data.frame(thsd0$`light_condition:metabolite`)
group_cont <-  tibble::rownames_to_column(group_cont, "MetNames")
#group_cont <- group_cont[group_cont$`p adj`<0.05,]
#write.csv(group_cont, "/home/mamata/Desktop/project/group_cont.csv", row.names=FALSE)
sig_metLL<-group_cont[which(grepl("LL:.+-LL:Negative.Control", group_cont$MetNames)== TRUE),]
length(which(heat1$`LL contrast`<=0)) #12
dim(sig_metLL) #13
#sig_metLL
#View(heat_pm1)
sig_metHL <- group_cont[which(grepl("HL:.+-HL:Negative.Control", group_cont$MetNames, ignore.case = TRUE)== TRUE),]
dim(sig_metHL)
length(which(heat1$`HL contrasts`<=0)) #40 metHL-Negative.ControlHL values are significant in lm
dim(heat1) #compare to the no. of significant groups 
head(summary1$test$coefficients) # the coefficients from the glht function


for ( i in 1:nrow(sig_metLL))
  
{
  if (sig_metLL$`p adj`[i] >= 0.05)
  {
    sig_metLL$`p adj`[i]= NA
    sig_metLL$diff[i] =NA
    i =i+1
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm1LL = as.data.frame(sig_metLL[, c(2,5)])

#taking out the metabolite names 
metNames_LL= gsub("LL:Negative.Control", "", sig_metLL[, "MetNames"])
metNames_LL = gsub("LL:", "", metNames_LL)
metNames_LL = gsub("-$", "", metNames_LL)
#taking out only lfcs


rownames(heatpm1LL) = metNames_LL
heatpm1LL <- tibble::rownames_to_column(heatpm1LL, "MetNames")

heatpm1LL

for ( i in 1:nrow(sig_metHL))
{
  if (sig_metHL$`p adj`[i] >= 0.05)
  {
    sig_metHL$`p adj`[i]= NA
    sig_metHL$diff[i] =NA
    i=i+1
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm1HL = as.data.frame(sig_metHL[c(2,5)])

#taking out the metabolite names 
metNames_HL= gsub("HL:Negative.Control", "", sig_metHL[, "MetNames"])
metNames_HL= gsub("HL:", "", metNames_HL)
metNames_HL= gsub("-$", "", metNames_HL)
#taking out only lfcs

rownames(heatpm1HL) = metNames_HL
heatpm1HL <- tibble::rownames_to_column(heatpm1HL, "MetNames")


heatpm1anova <- cbind (heatpm1LL[,2], heatpm1HL[,2])
rownames(heatpm1anova)= heatpm1HL[,"MetNames"]
colnames(heatpm1anova) = c("LL contrasts" ,"HL contrasts" )

#Delete rows with complete NAs
heatpm1anova<- as.data.frame(heatpm1anova[rowSums(is.na(heatpm1anova)) != ncol(heatpm1anova), ])

superheat(heatpm1anova,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "white"), left.label.text.size = 2.6)



pm2_anova_0 <- aov(values~light_condition*metabolite, data = spm2) 
summary(pm2_anova_0)
thsd0=TukeyHSD(pm2_anova_0)
group_cont <-as.data.frame(thsd0$`light_condition:metabolite`)
group_cont <-  tibble::rownames_to_column(group_cont, "MetNames")
#group_cont <- group_cont[group_cont$`p adj`<0.05,]
#write.csv(group_cont, "/home/mamata/Desktop/project/group_cont.csv", row.names=FALSE)
sig_metLL<-group_cont[which(grepl("LL:.+-LL:Negative.Control", group_cont$MetNames)== TRUE),]
length(which(heat_pm2$`LL contrast`<=0)) #52
dim(sig_metLL) #52
#sig_metLL
#View(heat_pm1)
sig_metHL <- group_cont[which(grepl("HL:.+-HL:Negative.Control", group_cont$MetNames, ignore.case = TRUE)== TRUE),]
dim(sig_metHL)
length(which(heat_pm2$`HL contrasts`<=0)) #40 metHL-Negative.ControlHL values are significant in lm
dim(heat_pm2) #compare to the no. of significant groups 
head(summarypm2$test$coefficients) # the coefficients from the glht function


for ( i in 1:nrow(sig_metLL))
  
{
  if (sig_metLL$`p adj`[i] >= 0.05)
  {
    sig_metLL$`p adj`[i]= NA
    sig_metLL$diff[i] =NA
    i*i+1
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm2LL = as.data.frame(sig_metLL[, c(2,5)])

#taking out the metabolite names 
metNames_LL= gsub("LL:Negative.Control", "", sig_metLL[, "MetNames"])
metNames_LL = gsub("LL:", "", metNames_LL)
metNames_LL = gsub("-$", "", metNames_LL)
#taking out only lfcs


rownames(heatpm2LL) = metNames_LL
heatpm2LL <- tibble::rownames_to_column(heatpm2LL, "MetNames")


for ( i in 1:nrow(sig_metHL))
{
  if (sig_metHL$`p adj`[i] >= 0.05)
  {
    sig_metHL$`p adj`[i]= NA
    sig_metHL$diff[i] =NA
    i=i+1
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm2HL = as.data.frame(sig_metHL[c(2,5)])
#taking out the metabolite names 
metNames_HL= gsub("HL:Negative.Control", "", sig_metHL[, "MetNames"])
metNames_HL= gsub("HL:", "", metNames_HL)
metNames_HL= gsub("-$", "", metNames_HL)
#taking out only lfcs
rownames(heatpm2HL) = metNames_HL
heatpm2HL <- tibble::rownames_to_column(heatpm2HL, "MetNames")
heatpm2HL

heatpm2anova <- cbind (heatpm2LL[,2], heatpm2HL[,2])
rownames(heatpm2anova)= heatpm2HL[,"MetNames"]
colnames(heatpm2anova) = c("LL contrasts" ,"HL contrasts" )
#Delete rows with complete NAs
heatpm2anova<- as.data.frame(heatpm2anova[rowSums(is.na(heatpm2anova)) != ncol(heatpm2anova), ])
superheat(heatpm2anova,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "white"), left.label.text.size = 2.6)


#jaccard index to see if anova and lm model have similar results
head(heat1)
head(heatpm1anova)
heat1 <- tibble::rownames_to_column(heat1, "Metabolites")
heatpm1anova <- tibble::rownames_to_column(heatpm1anova, "Metabolites")
jaccard1 <- jaccard(heat1$Metabolites, heatpm1anova$Metabolites)
jaccard1
#[1] 0.3636364
head(heat_pm2)
head(heatpm2anova)
heat_pm2<- tibble::rownames_to_column(heat_pm2, "Metabolites")
heatpm2anova <- tibble::rownames_to_column(heatpm2anova, "Metabolites")
heatpm2anova$Metabolites<- gsub("-$", "", heatpm2anova$Metabolites)
jaccard2 <- jaccard(heat_pm2$Metabolites, heatpm2anova$Metabolites)
jaccard2
#[1] 0.5555556


pm2_anova_0 <- aov(values~light_condition*metabolite, data = spm2after) 
summary(pm2_anova_0)
thsd0=TukeyHSD(pm2_anova_0)
group_cont <-as.data.frame(thsd0$`light_condition:metabolite`)
group_cont <-  tibble::rownames_to_column(group_cont, "MetNames")
#group_cont <- group_cont[group_cont$`p adj`<0.05,]
#write.csv(group_cont, "/home/mamata/Desktop/project/group_cont.csv", row.names=FALSE)
sig_metLL<-group_cont[which(grepl("LL:.+-LL:Negative.Control", group_cont$MetNames)== TRUE),]
length(which(heat_pm2$`LL contrast`<=0)) #14
dim(sig_metLL) 
#sig_metLL
#View(heat_pm1)
sig_metHL <- group_cont[which(grepl("HL:.+-HL:Negative.Control", group_cont$MetNames, ignore.case = TRUE)== TRUE),]
dim(sig_metHL)
length(which(heat_pm2$`HL contrasts`<=0)) #40 metHL-Negative.ControlHL values are significant in lm
dim(heat_pm2) #compare to the no. of significant groups 
head(summarypm2$test$coefficients) # the coefficients from the glht function


for ( i in 1:nrow(sig_metLL))
  
{
  if (sig_metLL$`p adj`[i] > 0.05)
  {
    sig_metLL$`p adj`[i]= NA
    sig_metLL$diff[i] =NA
    i =i+1
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm2LL = as.data.frame(sig_metLL[, c(2,5)])

#taking out the metabolite names 
metNames_LL= gsub("LL:Negative.Control", "", sig_metLL[, "MetNames"])
metNames_LL = gsub("LL:", "", metNames_LL)
metNames_LL = gsub("-$", "", metNames_LL)
#taking out only lfcs
rownames(heatpm2LL) = metNames_LL
heatpm2LL <- tibble::rownames_to_column(heatpm2LL, "MetNames")
heatpm2LL

for ( i in 1:nrow(sig_metHL))
{
  if (sig_metHL$`p adj`[i] > 0.05)
  {
    sig_metHL$`p adj`[i]= NA
    sig_metHL$diff[i] =NA
    i =i+1
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm2HL = as.data.frame(sig_metHL[c(2,5)])
#taking out the metabolite names 
metNames_HL= gsub("HL:Negative.Control", "", sig_metHL[, "MetNames"])
metNames_HL= gsub("HL:", "", metNames_HL)
metNames_HL= gsub("-$", "", metNames_HL)
metNames_HL
#taking out only lfcs
rownames(heatpm2HL) = metNames_HL
heatpm2HL <- tibble::rownames_to_column(heatpm2HL, "MetNames")
heatpm2anova <- cbind (heatpm2LL[,2], heatpm2HL[,2])
rownames(heatpm2anova)= heatpm2HL[,"MetNames"]
colnames(heatpm2anova) = c("LL contrasts" ,"HL contrasts" )
#Delete rows with complete NAs
heatpm2anova<- as.data.frame(heatpm2anova[rowSums(is.na(heatpm2anova)) != ncol(heatpm2anova), ])
superheat(heatpm2anova,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "white"), left.label.text.size = 2.6)


#for the pm2 dataset after removal of high variance groups


#ENRICHMENT SIMPLE
pm1compounds = read_excel("Data/PM01_PM2A_compound list.xlsx", sheet = "PM01", col_names = FALSE) # header=FALSE)
colnames(pm1compounds) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7")
pm1compounds$X3 = conv_metnames(pm1compounds$X3)
pm1compounds$X4 = conv_metnames(pm1compounds$X4)
head (lm1)

#separating metabolites that have effect and that do not
effect1 = as.vector(c(which((lm1$diffexpressed=="DOWN" & lm1$p.value <0.05)==TRUE), which((lm1$diffexpressed=="UP" & lm1$p.value <0.05)==TRUE)))
length(effect1) #68 
eff= data.frame(matrix(NA, nrow =68, ncol = 4))
j = 1
for ( i in effect1 )
{
  eff[j,] = lm1[i,]
  j =j+1
}
eff <- as.data.frame(eff)
head(eff)
#taking out the metabolite Names´
effectmets1 <- gsub("L-Negative.Control.+", "", eff$X1)
effectmets1 <- gsub("\\( ", "", effectmets1)
effectmets1 <- gsub("[A-Z]$", "", effectmets1)
#take unique cuz some metabolites are repeated
effectmets1 <-  unique(effectmets1)
Csource1 <- CsourceExtract(effectmets1)

noeffect1 = as.vector(which((lm1$diffexpressed=="NO")==TRUE))
length(noeffect1) #68 
NOeff= data.frame(matrix(NA, nrow =217, ncol = 4))
j = 1
for ( i in noeffect1 )
{
  NOeff[j,] = lm1[i,]
  j =j+1
}
NOeff <- as.data.frame(NOeff)
head(NOeff)
#taking out the metabolite Names´
NOeffectmets1 <- gsub("L-Negative.Control.+", "", NOeff$X1)
NOeffectmets1 <- gsub("\\( ", "", NOeffectmets1)
NOeffectmets1 <- gsub("[A-Z]$", "", NOeffectmets1)
#take unique cuz some metabolites are repeated
NOeffectmets1 <-  unique(NOeffectmets1)
Csource1no <- CsourceExtract(NOeffectmets1)

#contigency matrix
table_effect <- as.data.frame(table(Csource1$Csource))
table_noeffect <- as.data.frame(table(Csource1no$Csource))
table_effect
table_noeffect

contingency_mat1 <- matrix(NA, nrow=8, ncol= 2)
colnames(contingency_mat1) = c("Effect", "No_Effect")
rownames (contingency_mat1) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "Fatty acid")
contingency_mat1[,1] <- c(0, table_effect$Freq)
contingency_mat1[,2] <- table_noeffect$Freq
contingency_mat1
c1 <-as.data.frame(contingency_mat1)

probabilities1 = data.frame()
probs = data.frame()
for (class in row.names(contingency_mat1)){
  #goes from 0 to total no. of upregulation in that class
  q = contingency_mat1[class, "Effect"]
  p = phyper(
    q = q-1,
    m = sum(contingency_mat1[class, c("Effect", "No_Effect")]),
    n = sum(contingency_mat1[which(row.names(contingency_mat1) != class),
                            c("No_Effect", "No_Effect")]),
    k = sum(contingency_mat1[, "Effect"]), lower.tail = FALSE,
  )
  probs = data.frame("No. of significant Metabolites" = q,
                     "p-value" = p)
  #print (probs)
  #probabilities[[class]] = probs
  probs = data.frame(q,p)
  probabilities1 = rbind(probabilities1, probs)
  
}
colnames(probabilities1) <- c("No. of Metabolites", "Effect pvalue")
rownames (probabilities1) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "Fatty acid")
probabilities1


probabilitiesno = data.frame()
probs = data.frame()
for (class in row.names(contingency_mat1)){
  #goes from 0 to total no. of upregulation in that class
  q = contingency_mat1[class, "No_Effect"]
  p = phyper(
    q = q-1,
    m = sum(contingency_mat1[class, c("Effect", "No_Effect")]),
    n = sum(contingency_mat1[which(row.names(contingency_mat1) != class),
                             c("No_Effect", "No_Effect")]),
    k = sum(contingency_mat1[, "No_Effect"]), lower.tail = FALSE,
  )
  # probs = data.frame("No. of significant Metabolites" = q,
  #                    "p-value" = p)
  #print (probs)
  #probabilities[[class]] = probs
  probs = data.frame(q,p)
  probabilitiesno= rbind(probabilitiesno, probs)
  
}
probabilitiesno

colnames(probabilitiesno) <- c("No. of Metabolites", "No_effect pvalue")
rownames (probabilitiesno) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "Fatty acid")
probabilitiesno


#for presentation
enrich1 <- cbind(contingency_mat1[,1], probabilities1[,2],contingency_mat1[,2], probabilitiesno[,2])
colnames(enrich1) <- c("Effect", "p-value", "No effect", "p-value")
enrich1effect <- cbind(contingency_mat1[,1], probabilities1[,2])
colnames(enrich1effect) <- c("Effect", "p-value")
enrich1noeffect <- cbind(contingency_mat1[,2], probabilitiesno[,2])
colnames(enrich1noeffect) <- c("No Effect", "p-value")

pm2compounds = read_excel("/home/mamata/Desktop/project/code/Data/PM01_PM2A_compound list.xlsx", sheet = "PM2A", col_names = F)
colnames(pm2compounds) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7")
pm2compounds$X3 = conv_metnames(pm2compounds$X3)
pm2compounds$X4 = conv_metnames(pm2compounds$X4)
#pm2compounds <- pm2compounds[, -c(5,6,7,8, 11)]
dim (lm2) # before removal of groups #285
dim(lmpm2) #after removal of groups #246
# pm2compounds$Negative.Control= conv_metnames(pm2compounds$Negative.Control)
# pm2compounds$C.Source..negative.control= conv_metnames(pm2compounds$C.Source..negative.control)
# head (lm2)
#separating metabolites that have effect and that do not
effect2 = as.vector(c(which((lm2$diffexpressed=="DOWN" & lm2$p.value <0.05)==TRUE), which((lm2$diffexpressed=="UP" & lm2$p.value <0.05)==TRUE)))
length(effect2) #52
eff= data.frame(matrix(NA, nrow =52, ncol = 4))
j = 1
for ( i in effect2 )
{
  eff[j,] = lm2[i,]
  j =j+1
}
eff <- as.data.frame(eff)
head(eff)
#taking out the metabolite Names´
effectmets2 <- gsub("L-Negative.Control.+", "", eff$X1)
effectmets2 <- gsub("\\( ", "", effectmets2)
effectmets2 <- gsub("[A-Z]$", "", effectmets2)
#take unique cuz some metabolites are repeated
effectmets2 <-  unique(effectmets2)
Csource2<- CsourceExtractpm2(effectmets2)

noeffect2= as.vector(which((lm2$diffexpressed=="NO")==TRUE))
length(noeffect2) #233
NOeff= data.frame(matrix(NA, nrow =233, ncol=4))
j = 1
for ( i in noeffect2 )
{
  NOeff[j,] = lm2[i,]
  j =j+1
}
NOeff <- as.data.frame(NOeff)
head(NOeff)

#taking out the metabolite Names´
NOeffectmets2<- gsub("L-Negative.Control.+", "", NOeff$X1)
NOeffectmets2 <- gsub("\\( ", "", NOeffectmets2)
NOeffectmets2 <- gsub("[A-Z]$", "", NOeffectmets2)
#take unique cuz some metabolites are repeated
NOeffectmets2 <-  unique(NOeffectmets2)
Csource2no<- CsourceExtractpm2(NOeffectmets2)

#contigency matrix
table_effect <- as.data.frame(table(Csource2$Csource))
table_noeffect <- as.data.frame(table(Csource2no$Csource))
table_effect
table_noeffect

contingency_mat2<- matrix(NA, nrow=8, ncol= 2)
colnames(contingency_mat2) = c("Effect", "No_Effect")
rownames (contingency_mat2) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer")
contingency_mat2[,1] <- c(4, 0, 2,6,10,9,0,5)
contingency_mat2[,2] <- table_noeffect$Freq
contingency_mat2

probabilities2 = data.frame()
for (class in row.names(contingency_mat2)){
  #goes from 0 to total no. of upregulation in that class
  q = contingency_mat2[class, "Effect"]
  p = phyper(
    q = q-1,
    m = sum(contingency_mat2[class, c("Effect", "No_Effect")]),
    n = sum(contingency_mat2[which(row.names(contingency_mat2) != class),
                             c("No_Effect", "No_Effect")]),
    k = sum(contingency_mat2[, "Effect"]), lower.tail = FALSE,
  )
  probs = data.frame(q,p)
  #print (probs)
  probabilities2 = rbind(probabilities2, probs)
}
probabilities2
colnames(probabilities2) <- c("No. of Metabolites", "Effect pvalue")
rownames (probabilities2) =c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer")
probabilities2


probabilitiesno = data.frame()
probs = data.frame()
for (class in row.names(contingency_mat2)){
  #goes from 0 to total no. of upregulation in that class
  q = contingency_mat2[class, "No_Effect"]
  p = phyper(
    q = q-1,
    m = sum(contingency_mat2[class, c("Effect", "No_Effect")]),
    n = sum(contingency_mat2[which(row.names(contingency_mat2) != class),
                             c("No_Effect", "No_Effect")]),
    k = sum(contingency_mat2[, "No_Effect"]), lower.tail = FALSE,
  )
  # probs = data.frame("No. of significant Metabolites" = q,
  #                    "p-value" = p)
  #print (probs)
  #probabilities[[class]] = probs
  probs = data.frame(q,p)
  probabilitiesno= rbind(probabilitiesno, probs)
  
}
probabilitiesno

colnames(probabilitiesno) <- c("No. of Metabolites", "No_effect pvalue")
rownames (probabilitiesno) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer")
probabilitiesno

#for presentation
enrich2 <- cbind(contingency_mat2[,1], probabilities2[,2],contingency_mat2[,2], probabilitiesno[,2])
colnames(enrich2) <- c("Effect", "p-value", "No effect", "p-value")
enrich2effect <- cbind(contingency_mat2[,1], probabilities2[,2])
colnames(enrich2effect) <- c("Effect", "p-value")
enrich2noeffect <- cbind(contingency_mat2[,2], probabilitiesno[,2])
colnames(enrich2noeffect) <- c("No Effect", "p-value")

#Enrichment for spm2 after the removalof the high variance groups
head (lmpm2)
expelledMets <- setdiff(spm2$metabolite, spm2after$metabolite) # these metabolites are excluded from the pm2 dataset after variance removal
# now we expell these groups from the compounds list as well for the enrichment analysis
#taking out the indices of the compounds
indices = vector(length = length(expelledMets))
for ( i in 1: length(expelledMets))
{
  indices[i] <- grep (paste("^", expelledMets[i], "$", sep = ""), pm2compounds$X3)
  i = i+1
}
indices
pm2compounds[indices,]
dim(pm2compounds)
pm2compounds <- pm2compounds[-indices,]
dim(pm2compounds)

#separating metabolites that have effect and that do not
effect2 = as.vector(c(which((lmpm2$diffexpressed=="DOWN" & lmpm2$lmpval <0.05)==TRUE), which((lmpm2$diffexpressed=="UP" & lmpm2$lmpval <0.05)==TRUE)))
length(effect2) #78
eff= data.frame(matrix(NA, nrow =78, ncol = 4))
j = 1
for ( i in effect2 )
{
  eff[j,] = lmpm2[i,]
  j =j+1
}
eff <- as.data.frame(eff)
head(eff)
#taking out the metabolite Names´
effectmets2 <- gsub("L-Negative.Control.+", "", eff$X1)
effectmets2 <- gsub("\\( ", "", effectmets2)
effectmets2 <- gsub("[A-Z]$", "", effectmets2)
#take unique cuz some metabolites are repeated
effectmets2 <-  unique(effectmets2)
Csource2<- CsourceExtractpm2(effectmets2)

noeffect2= as.vector(which((lmpm2$diffexpressed=="NO")==TRUE))
length(noeffect2) #168
NOeff= data.frame(matrix(NA, nrow =168, ncol=4))
j = 1
for ( i in noeffect2 )
{
  NOeff[j,] = lmpm2[i,]
  j =j+1
}
NOeff <- as.data.frame(NOeff)
head(NOeff)

#taking out the metabolite Names´
NOeffectmets2<- gsub("L-Negative.Control.+", "", NOeff$X1)
NOeffectmets2 <- gsub("\\( ", "", NOeffectmets2)
NOeffectmets2 <- gsub("[A-Z]$", "", NOeffectmets2)
#take unique cuz some metabolites are repeated
NOeffectmets2 <-  unique(NOeffectmets2)
Csource2no<- CsourceExtractpm2(NOeffectmets2)

#contigency matrix
table_effect <- as.data.frame(table(Csource2$Csource))
table_noeffect <- as.data.frame(table(Csource2no$Csource))
table_effect
table_noeffect

contingency_mat2<- matrix(NA, nrow=8, ncol= 2)
colnames(contingency_mat2) = c("Effect", "No_Effect")
rownames (contingency_mat2) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer")
contingency_mat2[,1] <- c(4, 0, 2,6,10,9,0,5)
contingency_mat2[,2] <- table_noeffect$Freq
contingency_mat2

probabilities2 = data.frame()
for (class in row.names(contingency_mat2)){
  #goes from 0 to total no. of upregulation in that class
  q = contingency_mat2[class, "Effect"]
  p = phyper(
    q = q-1,
    m = sum(contingency_mat2[class, c("Effect", "No_Effect")]),
    n = sum(contingency_mat2[which(row.names(contingency_mat2) != class),
                             c("No_Effect", "No_Effect")]),
    k = sum(contingency_mat2[, "Effect"]), lower.tail = FALSE,
  )
  probs = data.frame(q,p)
  #print (probs)
  probabilities2 = rbind(probabilities2, probs)
}
probabilities2

colnames(probabilities2) <- c("No. of Metabolites", "Effect pvalue")
rownames (probabilities2) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer")
probabilities2

probabilitiesno = data.frame()
probs = data.frame()
for (class in row.names(contingency_mat2)){
  #goes from 0 to total no. of upregulation in that class
  q = contingency_mat2[class, "No_Effect"]
  p = phyper(
    q = q-1,
    m = sum(contingency_mat2[class, c("Effect", "No_Effect")]),
    n = sum(contingency_mat2[which(row.names(contingency_mat2) != class),
                             c("No_Effect", "No_Effect")]),
    k = sum(contingency_mat2[, "No_Effect"]), lower.tail = FALSE,
  )
  # probs = data.frame("No. of significant Metabolites" = q,
  #                    "p-value" = p)
  #print (probs)
  #probabilities[[class]] = probs
  probs = data.frame(q,p)
  probabilitiesno= rbind(probabilitiesno, probs)
  
}
probabilitiesno

colnames(probabilitiesno) <- c("No. of Metabolites", "No_effect pvalue")
rownames (probabilitiesno) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer")
probabilitiesno

#for presentation
enrich2 <- cbind(contingency_mat2[,1], probabilities2[,2],contingency_mat2[,2], probabilitiesno[,2])
colnames(enrich2) <- c("Effect", "p-value", "No effect", "p-value")
enrich2effect <- cbind(contingency_mat2[,1], probabilities2[,2])
colnames(enrich2effect) <- c("Effect", "p-value")
enrich2noeffect <- cbind(contingency_mat2[,2], probabilitiesno[,2])
colnames(enrich2noeffect) <- c("No Effect", "p-value")

# STEP 5 : BUILD SUPPORT VECTOR MACHINES AND RANDOM FOREST CLASSIFIERS TO SEPARATE
#THE METABOLITES WHICH HAVE EFFECT FROM THE ONES WHICH DONOT
#fingerprints
fingerprint1 <- read.csv("/home/mamata/Desktop/project/code/fingerprintPM1.csv")
fingerprint2 <- read.csv("/home/mamata/Desktop/project/code/fingerprintPM2.csv")
#removing columns with all zero entries cuz they give me all NAS in correlation matrix
fingerprint1<- subset(fingerprint1, select = -c(nF, tbonds))
namesfingerprint1 <- setdiff(pm1$metabolite, c("Gly.Asp", "Negative.Control" ))

fingerprint2<- subset(fingerprint2, select = -c(nF, tbonds))
namesfingerprint2 <- read.csv("/home/mamata/Desktop/project/code/Data/metspm2")

#giving the metabolite names to the fingerprints
fingerprint1 <- cbind(namesfingerprint1, fingerprint1)
fingerprint2 <-cbind(namesfingerprint2, fingerprint2)
colnames(fingerprint1) <- c( "names", "cid", "abonds" ,
                             "atoms", "bonds", "dbonds","HBA1","HBA2" ,"HBD","logP", "MR",  "MW","sbonds","TPSA")
#combine the finger prints
colnames(fingerprint2) <- colnames(fingerprint1)
fingerprint <- rbind(fingerprint1, fingerprint2)
fingerprint$names <- conv_metnames(fingerprint$names)
#correlation matrix of the features for the presentation 
corfeatures <-cor(fingerprint[,3 : ncol(fingerprint)])
library(corrplot)
corrplot(as.matrix(corfeatures), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
#modify the correlation matrix # setting evering in the upper triangle to zero
cormodified <- corfeatures                  
cormodified[upper.tri(cormodified)] <- 0
diag(cormodified) <- 0
cormodified
#Now, we can use this updated correlation matrix to remove all variables
#with a correlation above a certain threshold (i.e. > 0.9)
# Remove highly correlated variables
#chcks in columns  in cor1modified if there are any entries greater than 0.9
correlated <- apply(cormodified, 2, function(x) any(x > 0.6))
uncor_features<- fingerprint[ ,!correlated]
#these are the features for RF, SVM, Ridge
# cluster features and use one represnetative method
# random forest or svm  ( are good for classification problem)
# ridge regression 
#classifier for just classifying metabolites that effect and not effect regardless of up/down regulation
#using uncorrelated features for svm

str(uncor_features)
uncor_features <- uncor_features[, 2:5]
#take out the significant metabolites from heat1 and heat_pm2
# heat_pm2 is the one after removal of significant groups

met1 <- heat1$Metabolites
#after the removal of variance
met2before <- rownames(heat2)
length(met2before) #36

met2 <- heat_pm2$Metabolites
length(met2) #58
sigmets <- c(met1, met2)
sigmets
length(sigmets)

#initialize zeros
expression <- rep(0, nrow(fingerprint))
for ( i in 1: nrow(fingerprint))
{
  if ( fingerprint$names[i] %in% sigmets)
  {
    expression[i] <- 1
    i = i+1
  }
}

uncor_features <- cbind(expression, uncor_features)
head(uncor_features)
table(uncor_features$expression)

#gamma parameter defines how far the influence of a single training example reaches, 
#with low values meaning ‘far’ and high values meaning ‘close’. 
#The gamma parameters can be seen as the inverse of the radius of influence of 
#samples selected by the model as support vectors
# Soft margin classification : we modify the optimization problem to optimize both 
#the fit of the line to data and penalizing the amount of samples inside the margin at the same time, 
#where C defines the weight of how much samples inside the margin contribute to the overall error. 
#Consequently, with C you can adjust how hard or soft your large margin classification should be.
#With a low C, samples inside the margins are penalized less than with a higher C. 
#With a C of 0, samples inside the margins are not penalized anymore 
#- which is the one possible extreme of disabling the large margin classification. 
#With an infinite C you have the other possible extreme of hard margins.

#SPLITTING THE DATA ONLY ONCE FOR ALL THREE CLASSIFIERS
#SVM with e1071 library
library(caTools)
#scaling the features: z-score
features <- fingerprint[, 3: ncol(fingerprint)]
head(features)
features <- z_score(features)
features <- cbind(expression, features)
features[["expression"]] = factor(features[["expression"]])
set.seed(1000)
intrainrf <- createDataPartition(y =features$expression, p= 0.8, list = FALSE)
trainingrf <- features[intrainrf,]
testingrf <- features[-intrainrf,]
training_set <- uncorrelated(trainingrf)
test_set<- uncorrelated(testingrf)
library(e1071)
set.seed(100)
tune_svm <- tune(svm, expression ~ ., data =training_set,ranges = 
                   list(gamma=exp(seq(from=log(1/nrow(uncor_features)), to=2,len=20)), cost = seq(0.1, 10, len=20)), 
                 tunecontrol =tune.control(sampling = "cross", nrepeat=5, cross = 5, best.model = T, performances = T), probability=T)
bestmod <- tune_svm$best.model
min(tune_svm$performances$error) # the best model choses the parameter corresponding to it

plot(tune_svm, type = "contour")
xlabel <-paste("gamma=", as.character(round(tune_svm$best.parameters$gamma,3)), ", C=", as.character(round(tune_svm$best.parameters$cost, 3)))
#taking out the performances
performanceSVM <- tune_svm$performances
library(plotly)
p <- ggplot(performanceSVM, aes(gamma, cost, z= error)) +
  stat_contour(geom="polygon",aes(fill=stat(level))) +
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  geom_vline(xintercept = tune_svm$best.parameters$gamma, linetype = "dashed", color = "red", size= 0.25)+
  annotate("text", tune_svm$best.parameters$gamma-1.5, tune_svm$best.parameters$cost+0.1, label=xlabel, size=3)+
  geom_hline(yintercept = tune_svm$best.parameters$cost,  linetype = "dashed", color = "red", size = 0.25) 
ggplotly(p)
#took out the best model and did prediction on the training and the test data sets
library("raster")
#to calculate the probabilities
training_pred <-raster:: predict (bestmod, newdata = training_set, probability = TRUE)
training_pred
head(attributes(training_pred)$probabilities)
prob_svmTrain <- attr(training_pred, "probabilities")
prob_svmTrain
library("labelled")
#removing the attributes
training_pred <- remove_attributes(training_pred, "probabilities")
training_pred
#training_pred <- ifelse(training_pred>0.5, 1, 0)
confusion_training <- confusionMatrix(as.factor(training_pred), as.factor(training_set$expression), mode = "prec_recall")
confusion_training
test_pred <- raster::predict(bestmod, newdata = test_set, probability = T)
test_pred
prob_svmTest <- attr(test_pred, "probabilities")
prob_svmTest
test_pred <- remove_attributes(test_pred, "probabilities")
test_pred
#test_pred <- ifelse(test_pred>0.5, 1, 0)
#test_pred
confusion_test <- confusionMatrix(table(as.factor(test_pred), as.factor(test_set$expression)), mode = "everything")

confusion_test
prob_svmTest
roc_svmTrain <- pROC::roc(as.factor(training_set[,1]), as.numeric(prob_svmTrain[,2]), plot = TRUE )#
pr_svmTrain <- pROC::coords(roc_svmTrain, "all", ret = c("threshold", "recall", "precision"), transpose = FALSE)

roc_svmTest<- pROC::roc(as.factor(test_set[,1]), as.numeric(prob_svmTest[,2]), plot = TRUE )
pr_svmTest <- pROC::coords(roc_svmTest, "all", ret = c("threshold", "recall", "precision"), transpose = FALSE)

ggplot(pr_svmTrain, aes(recall, precision)) + 
  geom_path(aes(recall, precision), colour="salmon")+ geom_path(data = pr_svmTest, aes(recall, precision), colour="blue")

rocSVM<- pROC::roc(as.factor(training_set[,1]), as.numeric(prob_lasoTest), plot = TRUE )
pr <- pROC::coords(roclassoTest, "all", ret = c("threshold", "recall", "precision"), transpose = FALSE)

Train <- as.data.frame(confusion_training$byClass)
F1svmTrain<- Train[7,]
Test<- as.data.frame(confusion_test$byClass)
F1svmTest<- Test[7,]

library(randomForest)
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes
# mtry <- sqrt(ncol(trainingrf))
# tunegrid <- expand.grid(.mtry=mtry)
# Define Random Forest model and hyperparameters to tune
control <- trainControl(method = "cv", number = 5)
tuneGrid <- expand.grid(.mtry = seq(ceiling(sqrt(ncol(trainingrf))), ncol(trainingrf), len = 10), .ntree=c(1000, 1500, 2000, 2500))
#tuneGrid <- expand.grid(mtry = seq(0.1, 12, by = 0.1), ntree = seq(0.1, 500, by = 10))
#tuneGrid <- expand.grid(.mtry = seq(ceiling(sqrt(ncol(trainingrf))), ncol(trainingrf), 1))#, length.out = 50))
set.seed(1111)
rf_model <-caret::train(expression ~ ., data = trainingrf, method = customRF,
                        trControl = control, tuneGrid = tuneGrid)

rfResults <- rf_model$results
rfResults
summary(rf_model)
max(rf_model$results$Accuracy) # best tune corresponds to this value

ggplot(data=rfResults, aes(x=mtry, y=Accuracy, group=ntree))+
  geom_line(aes(color=ntree))+
  geom_point(aes(color=ntree))

xlabel <-paste("mtry=", as.character(rf_model$bestTune$mtry), ", ntree=", as.character( rf_model$bestTune$ntre))
p <- ggplot(rfResults, aes(mtry, ntree, z= Accuracy)) +
  stat_contour(geom="polygon",aes(fill=stat(level))) +
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  geom_vline(xintercept = rf_model$bestTune$mtry, linetype = "dashed", color = "red", size= 0.25)+
  annotate("text", rf_model$bestTune$mtry+1.5, rf_model$bestTune$ntree+5, label=xlabel, size=3)+
  geom_hline(yintercept = rf_model$bestTune$ntree,  linetype = "dashed", color = "red", size = 0.25) 
ggplotly(p)


# ggplotly(p)
RF_model <- randomForest(formula = expression~., data = trainingrf, ntree = rf_model$bestTune$ntree,
                         mtry = rf_model$bestTune$mtry, replace = TRUE, importance=TRUE)

randomForest::varImpPlot(RF_model)

important <- importance(RF_model, type=1 )  #Author DataFlair
Important_Features <- data.frame(Feature = row.names(important), Importance = important[, 1])
plot_ <- ggplot(Important_Features, aes(x= reorder(Feature, Importance) , y = Importance), size = 2 ) +
  geom_bar(stat = "identity", fill = "#800080") + coord_flip()+ theme_light(base_size = 10) + xlab("") +  ylab("Importance")+
  ggtitle("") + theme(plot.title = element_text(size=8))

ggsave("important_features.png",  plot_)
plot_
#Conditional=True, adjusts for correlations between predictors.
i_scores <- varImp(RF_model, conditional=TRUE)

pred_RF_train <- predict (object = RF_model, newdata = trainingrf, type = "class")
pred_RF_train
cm_rf <-confusionMatrix(table(trainingrf$expression, pred_RF_train), mode = "everything")
cm_rf


pred_RF_test <- predict (object = RF_model, newdata = testingrf, type = "class")
pred_RF_test
cm_rftest <-confusionMatrix(table(testingrf$expression, pred_RF_test), mode = "everything")
cm_rftest

Train <- as.data.frame(cm_rf$byClass)
F1rfTrain<- Train[7,]
F1rfTrain
Test<- as.data.frame(cm_rftest$byClass)
F1rfTest<- Test[7,]
F1rfTest

#predicting the probabilites for PR curve
rf_probTrain <- predict (object = RF_model, newdata = trainingrf, type = "prob")
rf_probTrain


rf_probTest <- predict (object = RF_model, newdata = testingrf, type = "prob")
rf_probTest


roc_rftrain <- pROC::roc(as.factor(trainingrf[,1]), as.numeric(rf_probTrain[,1]), plot = TRUE )
pr_rfTrain <- pROC::coords(roc_rftrain, "all", ret = c("threshold", "recall", "precision"), transpose = FALSE)
# ggplot(prcoords, aes(recall, precision)) + 
#   geom_path(aes(recall, precision), colour="salmon") + # needed to connect points in order of appearance +
#   theme(aspect.ratio = 1) +
#   theme(panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"))
roc_rfTest<- pROC::roc(as.factor(testingrf[,1]), as.numeric(rf_probTest[,1]), plot = TRUE )
pr_rfTest <- pROC::coords(roc_rfTest, "all", ret = c("threshold", "recall", "precision"), transpose = FALSE)

library(ROCR)
predobj <- prediction(as.vector(pred_RF_test), testingrf[,1])
perf <- performance(predobj,"prec", "rec")
plot(perf)

ggplot(pr_rfTrain, aes(recall, precision)) + 
  geom_path(aes(recall, precision), colour="salmon")+ geom_path(data = pr_rfTest, aes(recall, precision), colour="blue") 

#lasso regression
library(glmnet)   
set.seed(100)
traininglasso <- trainingrf
testinglasso <- testingrf
str(traininglasso)
train_x= traininglasso[, 2: ncol(traininglasso)]
train_y=traininglasso$expression
test_x= testinglasso[, 2: ncol(testinglasso)]
test_y=testinglasso$expression
set.seed(143)
train_matrix <- model.matrix(train_y ~ ., data = train_x)
test_matrix <- model.matrix(test_y ~ ., data = test_x)

#scaling and considering all variables as numeric(!!)
# CROSS VALIDATION WITH DEFAULT METERS
# ALPHA=1
# LAMBDAS = 10^seq(2, -3, by = -.1)
# set.seed(456)
# CV=cv.glmnet(train_matrix, train_y,alpha=ALPHA,lambda = LAMBDAS, type.measure = "mse", nfolds = 5, family = "binomial")
# #class give misclassification error for the binomial 
# plot(CV)
# bestLAMBDA=CV$lambda.1se
# bestLAMBDA
# lasso.mod=glmnet(x,y,alpha = 1,lambda = bestLAMBDA, family = "binomial")
# coef(lasso.mod)
# 
# #CROSS VALIDATION
str(traininglasso)
LAMBDA <-10^seq(2, -3, by = -.1)
#LAMBDA <- expand.grid(lambda=LAMBDA)
# 5-fold cross validation and we wanna save all thw prediction
ctrl <-trainControl("cv", number=5, savePredictions = "all" )
# specify lasso regression model to be estimated using training data
set.seed(1234)
lasso <- caret::train(expression~., data = traininglasso,
                      method="glmnet", tuneGrid = expand.grid(alpha=1,lambda = LAMBDA),
                      trControl = ctrl)#, family = "binomial")
lassoResults <- lasso$results
head(lassoResults)
lasso$bestTune
bestLambda <- lasso$bestTune$lambda
plot(lasso)
lasso$bestTune$lambda
varImp(lasso)
#lasso model parameter estimates
coef(lasso$finalModel, lasso$bestTune$lambda)

#other way to find the best lambda with highest accurcy
idx  <- which.max(lassoResults$Accuracy)
idx
maxParam <- lassoResults[idx,]
maxParam
bestLambda <- maxParam$lambda
bestLambda

xlabel = paste("lambda=", round(maxParam$lambda,4))
ggplot(lassoResults, aes(Accuracy, lambda)) + geom_point()+
  geom_point(aes(x = maxParam$Accuracy, y = maxParam$lambda),color = "red")+
  geom_hline(yintercept = maxParam$lambda, linetype = "dashed", color = "red", size= 0.25)+
  geom_vline(xintercept = maxParam$Accuracy,  linetype = "dashed", color = "red", size= 0.25)+
  annotate("text",maxParam$Accuracy-0.005, maxParam$lambda+3, label=xlabel, size=3)

#col = ifelse(1:nrow(lassoResults) ==lasso$bestTune$lambda, "red", "blue"), size = 2)

bestlasso =  glmnet(train_matrix, train_y, family = "binomial", alpha=1, lambda = maxParam$lambda)
coefLasso <- coef(bestlasso)
coefLasso
#coefLasso <- tibble::rownames_to_column(coefLasso, "Features")

#model prediction on the training data set
pred_train_lasso <- predict (bestlasso, s = maxParam$lambda,  
                             newx=train_matrix, type="response")
pred_train_lasso
prob_lassoTrain <- pred_train_lasso

pred_train_lasso <- ifelse(pred_train_lasso>0.5, 1, 0)
table <- table(Predicted = pred_train_lasso, Actual = traininglasso[,1])
table
Cmlassotrain<- confusionMatrix(table, mode = "everything")
Cmlassotrain

#model prediction on the test data
# testinglasso <- as.matrix(testinglasso)
# newy <- as.numeric(testinglasso[,1])
# newx =as.matrix(testinglasso[, 2: ncol(testinglasso)])
predLasso <- predict(bestlasso, newx=test_matrix, s=maxParam$lambda, type = "response")
predLasso
prob_lasoTest<- predLasso
predLasso <- ifelse(predLasso>0.5, 1, 0)
table <- table(Predicted = predLasso, Actual = testinglasso$expression)
table
Cmlassotest <- confusionMatrix(table, mode = "everything")
Cmlassotest

Train <- as.data.frame(Cmlassotrain$byClass)
F1lassoTrain<- Train[7,]
F1lassoTrain
Test<- as.data.frame(Cmlassotest$byClass)
F1lassoTest<- Test[7,]
F1lassoTest

roclassoTrain <- pROC::roc(as.factor(traininglasso[,1]), as.numeric(prob_lassoTrain), plot = TRUE )
pr_lassoTrain <- pROC::coords(roclassoTrain, "all", ret = c("threshold", "recall", "precision"), transpose = FALSE)
# ggplot(prcoords, aes(recall, precision)) + 
#   geom_path(aes(recall, precision), colour="salmon") + # needed to connect points in order of appearance +
#   theme(aspect.ratio = 1) +
#   theme(panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"))
roclassoTest<- pROC::roc(as.factor(testinglasso[,1]), as.numeric(prob_lasoTest), plot = TRUE )
pr_lassoTest <- pROC::coords(roclassoTest, "all", ret = c("threshold", "recall", "precision"), transpose = FALSE)

ggplot(pr_lassoTrain, aes(recall, precision)) + 
  geom_path(aes(recall, precision), colour="salmon")+ geom_path(data = pr_lassoTest, aes(recall, precision), colour="blue") # + # need+ed to connect points in order of appearance +
  # theme(aspect.ratio = 1) +
  # theme(panel.background = element_blank(), 
  #       axis.line = element_line(colour = "black"))

svmtrain <- prediction(prob_svmTrain[,1], training_set[,1])
pred1train <- performance(svmtrain, "prec", "rec")
plot(pred1train, colorize=T, lwd=2, main="Precision-Recall Curve")

rfCM <- cm_rf$table
rfCM
rftrain <- prediction(rf_probTrain[,2], trainingrf$expression)
pred2train <- performance(rftrain, "prec", "rec")
plot(pred2train, colorize=T, lwd=2, main="Precision-Recall Curve")

Cmlassotrain$table
lassoTrain <- prediction (prob_lasoTrain[,1], traininglasso[,1])
pred3train <- performance(lassoTrain, "prec", "rec")
plot(pred3train, colorize=T, lwd=2, main="Precision-Recall Curve")

#now do confusionMatrix(pred, actual, mode = "everything", positive="1") to get the f1 score
#plot the F1 score for the test and training data as a barplot

#extract f1 scores for all the methods

# Store precision and recall scores at different cutoffs



vec_F1 <- c(F1svmTrain, F1svmTest, F1rfTrain, F1rfTest, F1lassoTrain, F1lassoTest )
vec_names <- c("svm Train", "svm Test", "rf Train", "rf Test", "lasso Train", "lasso Test" )

F1 <- data.frame(Classifier = c("svm Train", "svm Test", "rf Train", "rf Test", "lasso Train", "lasso Test" ), 
                 F1_score = c(F1svmTrain, F1svmTest, F1rfTrain, F1rfTest, F1lassoTrain, F1lassoTest ) )

match =c ("svm", "rf", "lasso")

# Assigning default and different colors to bar plot
perf <-ggplot(data = F1, aes(x= Classifier, y=F1_score, fill=Classifier))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("yellow",
                             "yellow",
                             "red",
                             "red",
                             "blue",
                             "blue"))
perf

#(You can include the crossvalidation mean and standard error as error bars for the training data set)

#Since you have unbalanced data sets I would rather use Precission recall curves than receiver-operator curves

#Precision recall curve for test set for all methods
ggplot()+ geom_path(data=pr_svmTest,  aes(recall, precision, color = "SVM"))+
  geom_path(data = pr_rfTest, aes(recall, precision, colour="RF"))+
  geom_path(data = pr_lassoTest, aes(recall, precision, colour="Lasso"))+
  scale_color_manual(name = "Classifiers", values = c("SVM" = "darkblue", "RF" = "red", "Lasso"= "yellow"))

ggplot()+ geom_path(data=pr_svmTrain,  aes(recall, precision, color = "SVM"))+
  geom_path(data = pr_rfTrain, aes(recall, precision, colour="RF"))+
  geom_path(data = pr_lassoTrain, aes(recall, precision, colour="Lasso"))+
  scale_color_manual(name = "Classifiers", values = c("SVM" = "darkblue", "RF" = "red", "Lasso"= "yellow"))

F1rfTest
F1rfTrain

F1svmTest
F1svmTrain

F1lassoTest
F1svmTrain

#The next step would be to extract feature importance from the best perfoming classifier.    



