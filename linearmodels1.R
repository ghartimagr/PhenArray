#setting the working directory
setwd("/home/mamata/Desktop/project/code/")
#loading the required libraries
library("xlsx")
library("reshape2")
library("lattice")
library("dplyr")
library("tibble")
library("limma")
source("Scaling.R") #contains all my functions for the project
library("contrast")
library("multcomp")
library("ggplot2")
library("ggrepel")
conv_metnames=function(metn) {
  #Function to convert an character vector or metabolite names into valid
  #R variable names
  
  #substitute all invalid characters with . 
  metn=gsub(" |-|,", ".",metn)
  #quotes are removed
  metn=gsub("`", "", metn)
  #Add a 'met' as prefix to all metabolite names starting with a number
  metn[grep("^[[:digit:]]", metn)]=paste("met", metn[grep("^[[:digit:]]", metn)], sep="")
  return(metn)
}

#STEP 1: PLOTTING THE RAW DATA: BOX PLOTS
#CONTROL PLATE
#loading the table from the excel sheet
pmcontrol = read.xlsx("Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheetName = "control plate")
head(pmcontrol)
#melt works like stack: stacking LL and HL columns into one for the box plot
newpmcontrol = melt(pmcontrol, id = "position")
head(newpmcontrol)
save(newpmcontrol, file = "newpmcontrol.RData")
bwplot(newpmcontrol$value~newpmcontrol$position|newpmcontrol$variable, 
       col = as.numeric(as.factor(newpmcontrol$position)), ylab= "luminescence values", xlab= "Position",scales=list(x=list(draw=FALSE))) #scales =list(newpmcontrol$position, cex = 0.5, rot = 90))
#we see the position effect in the samples, the values on the left hand side
# of the  plate seem to be higher, it might be because of the position of the light source
#they might have positioned the light in the left hand corner that is why the samples there received more light

#PM1 PLATE
#reading the data from the excel sheet
pm1 = read.xlsx("Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheetName = "PM1")
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
pm1stacked<- Stacked(Tpm1stacked)
head(pm1stacked)
save(pm1stacked, Tpm1, file ="pm1stacked.RData")
#boxplot(Tpm1stacked$values~Tpm1stacked$sample , col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")
bwplot(pm1stacked$values~pm1stacked$metabolite|pm1stacked$light_condition, 
       col= rainbow(ncol(Tpm1)), ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))


# PM2 PLATE
#reading the PM2A sheet using xlsx library
pm2 = read.xlsx("Data/three experiments_biolog PM1-PM2A_LL-HL.xlsx", sheetName = "PM2A")
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
pm2stacked <- Stacked(Tpm2stacked)
head(pm2stacked)
save(pm2stacked, Tpm2, file ="pm2stacked.RData" )
#using bwplot from the lattice library
#values of the metabolites are in the y axis and metabolites are in the x-axis
#the values are grouped together by HL and LL
bwplot(pm2stacked$values~pm2stacked$metabolite|pm2stacked$light_condition, 
       col= rainbow(ncol(Tpm2)), ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))
#STEP2 : SCALING AND LOG TRANSFORMATION
#SCALING THE DATA TO REMOVE THE POSITION EFFECT SO THAT ONLY THE METABOLITE AND LIGHT EFFECTS REMAIN
#TO REMOVE THE RIGHT SKEWNESS OF OUR DATA, WE USE LOG TRANSFORMATION
#first we calculate the scaling factor
#Scaling factor for HL :  data point in the HL column is  divided by the mean of whole HL column and the same holds for the LL col.

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
logcontrol <- log2(spmcontrol[, c(2,3)])
#binding the position column from the original table and the new log transformed table
snewpmcontrol <- cbind(spmcontrol[,1], logcontrol)
colnames(snewpmcontrol) <- c("position", "LL", "HL")
snewpmcontrol = melt(snewpmcontrol, id = "position") #sort of rotating the data
head(snewpmcontrol)

#SAVING THE DATA FOR PLOTTING CONTROL DATA BEFORE AND AFTER POSITION CORRECTION
save(snewpmcontrol, file = "snewpmcontrol.RData")

corcontrol <- cor(newpmcontrol$value, snewpmcontrol$value, method = "spearman")
save(corcontrol, file = "cor.RData")

#boxplot for the scaled values gives a straight line for HLs and LLs
bwplot(snewpmcontrol$value~snewpmcontrol$position|snewpmcontrol$variable, 
       col = as.numeric(as.factor(snewpmcontrol$position)), ylab= "luminescence values", xlab= "Position", scales=list(x=list(draw=FALSE)))
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
spm1_notlog <- Stacked(spm1stacked)
#plotting values vs. metabolites grouped by HL and LL
bwplot(spm1_notlog$values~spm1_notlog$metabolite | spm1_notlog$light_condition, col= rainbow(ncol(pm1new)), ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

# pm1ranked <- pm1stacked %>%
#   mutate(good_ranks = order(order(values, decreasing=TRUE)))
# 
# spm1ranked <- spm1_notlog %>%
#   mutate(good_ranks = order(order(values, decreasing=TRUE)))
# corpm1 <- cor(pm1ranked$good_ranks, spm1ranked$good_ranks, method = "spearman")
# corpm1

#for spearman'S corelation
pm1ranked <- transform(pm1stacked,rank=ave(1:nrow(pm1stacked),metabolite,
                              FUN=function(x) order(values[x],decreasing=TRUE)))

spm1ranked <- transform(spm1_notlog,rank=ave(1:nrow(spm1_notlog),metabolite,
                               FUN=function(x) order(values[x],decreasing=TRUE)))
corpm1 <- cor(pm1ranked$rank, spm1ranked$rank, method = "spearman")
corpm1

#SAVING THE DATA FOR PLOTTING PM1 BEFORE AND AFTER POSITION CORRECTION
save(pm1stacked, Tpm1, spm1_notlog, pm1new, corpm1, file= "spm1.RData")

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
spm1 <- Stacked(spm1stacked)
head(spm1)
save(spm1, pm1new, file = "spm1log.RData")
#save(pm1new, file = "pm1new.RData")
#plotting values vs. metabolites grouped by HL and LL
bwplot(spm1$values~spm1$metabolite | spm1$light_condition, col= rainbow(ncol(pm1new)), ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))


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
spm2_notlog <- Stacked(spm2stacked)
head(spm2_notlog)


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
spm2 <- Stacked(spm2stacked)
head(spm2)
#Values vs. Metabolite grouped by light condition 
bwplot(spm2$values~spm2$metabolite | spm2$light_condition, col= rainbow(ncol(pm2new)), 
       ylab = "luminescence values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

#for spearman'S corelation
pm2ranked <- transform(pm2stacked,rank=ave(1:nrow(pm2stacked),metabolite,
                                           FUN=function(x) order(values[x],decreasing=TRUE)))

spm2ranked <- transform(spm2_notlog,rank=ave(1:nrow(spm2_notlog),metabolite,
                                             FUN=function(x) order(values[x],decreasing=TRUE)))
corpm2 <- cor(pm2ranked$rank, spm2ranked$rank, method = "spearman")
corpm2

save(pm2stacked, Tpm2, spm2_notlog, pm2new, spm2, corpm2,  file = "spm2.RData")
#STEP 3: FITTING THE LINEAR MODELS AFTER SCALING AND TRANSFORMING THE DATA
#LINEAR MODEL FOR CONTROL PLATE
#for model without intercept
#forming a group row in the dataframe with position and light conditions (variable) columns combined
snewpmcontrol$group <- factor(paste0(snewpmcontrol$position, snewpmcontrol$variable))
head(snewpmcontrol)
colnames(X)
head(X)
X <- model.matrix(~0+group, data = snewpmcontrol)
fitpmc <- lm(value~0+group, data = snewpmcontrol)
summary(fitpmc)
#plot(fitpmc)


# LINEAR MODEL FOR PM1 - 
#this our scaled pm1 table with values , light condition and 
#metabolites stacked in their repective columns and we use it to fit the model now
#head(spm1stacked)
head(spm1) 
#grouping the metabolite and light condition together
spm1$group <- factor(paste0(spm1$metabolite, spm1$light_condition))
head(spm1)
Xpm1<-model.matrix(~0+group, data =spm1)
### R base LM IMPLEMENTATION
#remove group prefix - here is one example on how to get contrast which you can then use to extract inference statistics please write a loop that extracts all other contrasts
colnames(Xpm1)=gsub("^group", "", colnames(Xpm1))
#taking out the column names of the Xpm1 matrix for the metabolite names
cHL= Xpm1[,165] #negative control HL 
cLL= Xpm1[,166] # negative controlLL
#contpm=makeContrasts(contrasts =c("MaltoseHL-Negative.ControlHL", "MaltoseLL-Negative.ControlLL", "(MaltoseHL-Negative.ControlHL)-(MaltoseLL-Negative.ControlLL)"), levels=Xpm1)
#head(contpm)
#the following code is to create a character vector of the column names to calculate the
#differences and differences of the differences using the makeContrasts 
vec <- vector(mode = "character", length = ncol(Xpm1)/2)
vec2 <- vector(mode = "character", length = ncol(Xpm1)/2)
for ( i in 1:ncol(Xpm1))
{
  if (i%%2==0)
  {
    vec[i/2] = paste(colnames(Xpm1)[i], "Negative.ControlLL", sep = "-")

  }
  else
  {
    vec2[i%/%2+1] = paste(colnames(Xpm1)[i], "Negative.ControlHL", sep = "-")
  }
}


#removing "Negative.ControlHL-Negative.ControlHL", "Negative.ControlLL-Negative.ControlLL", 
#"Negative.ControlHL-Negative.ControlHL-Negative.ControlLL-Negative.ControlLL"
#so that they do not introduce NA's in pvalue calculation 
#summary function literally fails if i have NAs
vec = setdiff(vec, c("Negative.ControlLL-Negative.ControlLL"))
vec2 = setdiff(vec2, c("Negative.ControlHL-Negative.ControlHL"))
metabolitevec1 = vector()
for ( i in 1: length(vec))
{
  metabolitevec1[i] = paste( c("("), vec[i], c(")"))
}
metabolitevec2 = vector()
for ( i in 1: length(vec2))
{
  metabolitevec2[i] = paste( c("("), vec2[i], c(")"))
}
metabolitevec <- paste(metabolitevec2, metabolitevec1, sep ="-")
#metabolitevec <- setdiff(metabolitevec, c("( Negative.ControlHL-Negative.ControlHL )-( Negative.ControlLL-Negative.ControlLL )"))
contpm1=makeContrasts(contrasts =c(vec2, vec, metabolitevec), levels=Xpm1)
#using default lm function
fitpm1 <- lm(values~0+group, data = spm1)
# glht already adjust pvalues
confit=glht(fitpm1, t(contpm1))
summarypm1= summary(confit, test = adjusted(type = "fdr"))
#plot pvalues
hist(summarypm1$test$pvalues)

#checking if the p values are actually zero as shown in the values
#no olny 7 are zero and the rest are not
table(summarypm1$test$pvalues ==0) 

pval =as.data.frame(summarypm1$test$pvalues) #pvalues
lfc= as.data.frame(summarypm1$test$coefficients) #coefficients
lmpm1 <- as.data.frame(cbind(lfc, pval))
lmpm1 <- tibble::rownames_to_column(lmpm1, "contrasts") 
#dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(lmpm1) = c( "contrasts", "lmlfc", "lmpval")
head(lmpm1)
save(lmpm1, file = "lmpm1.RData")
# Volcano plot
p <- ggplot(data=lmpm1, aes(x=lmlfc, y=-log2(lmpval))) + geom_point() + theme_minimal()
p
lmpm1$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
lmpm1$diffexpressed[lmpm1$lmlfc > 0.6 & lmpm1$lmpval < 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
lmpm1$diffexpressed[lmpm1$lmlfc < -0.6 & lmpm1$lmpval< 0.05] <- "DOWN"
p <- ggplot(data=lmpm1, aes(x=lmlfc, y=-log2(lmpval), col=diffexpressed)) + geom_point() + theme_minimal()
p

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially 
#expressed (NA in case they are not)
#commenting the labels. we do not need it
# lmpm1$label <- NA
# lmpm1$label[lmpm1$diffexpressed != "NO"] <- lmpm1$contrasts[lmpm1$diffexpressed != "NO"]
# 
# ggplot(data=lmpm1, aes(x=lmlfc, y=-log10(lmpval), col=diffexpressed)) + 
#   geom_point() + 
#   theme_minimal() +
#   geom_text(data = lmpm1, aes(label = label), check_overlap = TRUE, size = 3, hjust = 0.5)

### LIMMA IMPLEMENTATION
#We trasnfrom the data here into limma fromat for 1 gene
spm1limma=data.frame(t(spm1$values))
colnames(spm1limma)=spm1$group
#generate canonical linear model
limfit = lmFit(spm1limma,Xpm1)

#extract contrasts
limcont = contrasts.fit(limfit, contpm1)
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

#which genes are down regulated, upregulated or none ?
downregulated = as.vector(which((dfpm1$diffexpressed=="DOWN")==TRUE))
#length(downregulated) : 52 dowwnregulated 

# down = data.frame(matrix(NA, nrow =52, ncol = 4))
# j = 1
# for ( i in downregulated)
# {
#   down[j,] = dfpm1[i,]
#   j =j+1
# }
# down
# upregulated = as.vector(which((dfpm1$diffexpressed=="UP")==TRUE))
# length(upregulated)
# up = data.frame(matrix(NA, nrow =16, ncol = 4))
# j = 1
# for ( i in upregulated)
# {
#   up[j,] = dfpm1[i,]
#   j =j+1
# }
# up
#binding coefficients and p values 
lfcpm1 =as.data.frame(cbind(summarypm1$test$coefficients, summarypm1$test$pvalues))
lfcpm1 <- tibble::rownames_to_column(lfcpm1, "contrasts")


colnames(lfcpm1) = c("contrasts", "lfc", "pval" )
head(lfcpm1)
#setting all the insignificant p values to NA
for ( i in 1:nrow(lfcpm1))
{
  if (lfcpm1$pval[i] >= 0.05)
  {
    lfcpm1$pval[i]= NA
    lfcpm1$lfc[i] = NA
  }
}

#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heat_pm1 = as.data.frame(cbind(lfcpm1[1:95,], lfcpm1[96:190,], lfcpm1[191:285,]))
colnames(heat_pm1) = c("HLgroups", "HL_lfc", "HLpvalues", "LLgroups", "LL_lfc","LLpvalues" , "HL_LL_dif_groups", "HL_LL_dif_lfc", "HL_LL_dif_pvalues")
#taking out the metabolite names 
metNames_pm1= gsub("HL.+", "", heat_pm1[, "HLgroups"])
#taking out only lfcs
heat_pm1 = as.data.frame(cbind(heat_pm1$LL_lfc, heat_pm1$HL_lfc, heat_pm1$HL_LL_dif_lfc)) # only lfcs
rownames(heat_pm1) = metNames_pm1

heat_pm1 <- heat_pm1[rowSums(is.na(heat_pm1)) != ncol(heat_pm1), ]
colnames(heat_pm1) = c("LL contrasts", "HL contrasts", "HL-LL contrasts")

save(heat_pm1, file= "heat_pm1.RData")

# llcontrast = as.data.frame(lmpm1$lmlfc[which(as.numeric(lmpm1$indices)<=95)])
# #take out the signifcant columns and 
# hlcontrast = which(as.numeric(lmpm1$indices)>95 & as.numeric(lmpm1$indices)<=190)
# hl_ll_diffcontrast = which(as.numeric(lmpm1$indices)>190)
# significant_pm1 = as.data.frame(c(llcontrast, hlcontrast, hl_ll_diffcontrast))
library("superheat")
superheat(heat_pm1,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "red"), left.label.text.size = 3) 

# superheat(heat_pm1,
#           # scale the matrix
#           # change color of missing values
#           heat.na.col = "white", heat.col.scheme =  "red")

#ENRICHMENT ANALYSIS: to find the categories of compounds inhibiting/accelerating the gene expression
#partition the metabolite by 2 criterias: 
#class 1 : what class of compounds they belong to 
#class 2 : if they are over or under-expressed
# we then get a matrix class 1 in rows and class 2 in columns
#we fill every element in the matrix by  no. of metabolite in the class x, thus we obtain a contigency matrix
#then we can test whether the particular entry of contingency matrix is statistically significant

#importing the compound list
pm1compounds = read.xlsx("Data/PM01_PM2A_compound list.xlsx", sheetName = "PM01", header=FALSE)
pm1compounds$X3 = conv_metnames(pm1compounds$X3)
pm1compounds$X4 = conv_metnames(pm1compounds$X4)
head (lmpm1)
#which genes are down regulated, upregulated or none ?
downregulated = as.vector(which((lmpm1$diffexpressed=="DOWN" & lmpm1$lmpval <0.05)==TRUE))
length(downregulated) #: 52 dowwnregulated 
down = data.frame(matrix(NA, nrow =52, ncol = 4))
j = 1
for ( i in downregulated)
{
  down[j,] = lmpm1[i,]
  j =j+1
}
head(down)
dim(down) #51

downMets <-downdiff(down)
# CsourceDownLL<- CsourceExtract(downMets$LLmets)
# CsourceDownHL <- CsourceExtract(downMets$HLmets)
df_downLL= CsourceExtract(downLLmets)
df_downHL= CsourceExtract(downHLmets)

#Upregulation
upregulated = as.vector(which((lmpm1$diffexpressed=="UP" & lmpm1$lmpval <0.05)==TRUE))
length(upregulated)
up = data.frame(matrix(NA, nrow =16, ncol = 4))
j = 1
for ( i in upregulated)
{
  up[j,] = lmpm1[i,]
  j =j+1
}
up
# up consists of ohly those metabolites in the third colum of the heat map ie. the interaction component
upMets <- gsub("HL(-.+).+", "", up$X1)
upMets <- gsub(".+ ", "", upMets)
#upMets <- gsub('.+([a-z]+)HL.+', '(\\1)', up$X1)
CsourceUp = vector(length = length(upMets), mode = "character")
df_up <- CsourceExtract(upMets)


Noreg = as.vector(which((lmpm1$diffexpressed=="NO")==TRUE))
length(Noreg)
No= data.frame(matrix(NA, nrow =217, ncol = 4))
j = 1
for ( i in Noreg)
{
  No[j,] = lmpm1[i,]
  j =j+1
}
No = as.data.frame(No)
head(No)
dim(No) #217
NoMets <- Nodiff(No)
df_noLL= CsourceExtract(NoMets$LLmets)
df_noHL= CsourceExtract(NoMets$HLmets)
df_noHL_LL = CsourceExtract(NoMets$HL_LLmets)


#no regulation for seaprated ones
tableHL_LL <- as.data.frame(table(df_noHL_LL$Csource))
tablenoLL <-  as.data.frame(table(df_noLL$Csource))
tablenoHL <- as.data.frame(table(df_noHL$Csource))
#ureg for separated ones
tableup <- as.data.frame(table(df_up$Csource))
#downreg for the separated ones
tabledownLL <- as.data.frame(table(df_downLL$Csource))
tabledownHL <- as.data.frame(table(df_downHL$Csource))
#creating contingecy matrix
contingency_mat <- matrix(NA, nrow=8, ncol= 6)
colnames(contingency_mat) = c( "DownregulationLL", "DownregulationHL", "UpregulationHL_LL",
                               "NoRegulationLL", "NoRegulationHL", "NoregulationHL_LL")
rownames (contingency_mat) = c("Alcohol", "Amide", "Amine", "Amino acid", 
                                "Carbohydrate", "Carboxylic acid", "Ester", "Fatty acid")
contingency_mat[,1] = c(0, tabledownLL$Freq)
contingency_mat[,2] =  c(0, 1,1,0, 4,3,0,3)
contingency_mat[,3] =  c(0,0, 0, 1, 7, 8, 0, 0)
contingency_mat[,4] = c(2, 0,0,12,25,16,0,0)
contingency_mat[,5] = c(2,0,1,16,34,29,1,0)
contingency_mat[,6] =tableHL_LL$Freq

contingency_mat
hyperup <- phyperfunction(contingency_mat, "UpregulationHL_LL")
hyperdownLL <- phyperfunction(contingency_mat, "DownregulationLL")
hyperdownHL <- phyperfunction(contingency_mat, "DownregulationHL")
hypernoLL <- phyperfunction(contingency_mat, "NoRegulationLL")
hypernoHL <-  phyperfunction(contingency_mat, "NoRegulationHL")
hyperno <- phyperfunction(contingency_mat, "NoregulationHL_LL")

library(ggpubr)

plotfunction(hyperup)
plotfunction(hyperdownHL)
plotfunction(hyperdownLL)

#check contigencyMat.R for the plots without no_regulation
#q : vector of quantiles representing the number of white balls drawn without 
#replacement from a set which contains both black and white balls.
#m: the number of white balls in the set
#n: the number of black balls in the set
#k: the number of balls drawn from the set, hence must be in 0,1,\dots, m+n0,1,…,m+n.

#from the contingency matrix
#m = rowsum of the entry we are interested in
# n = grandtotal - the rowsum of the entry
# k = colsum of the column the entry belongs to

# EXAMPLES ##############
# first contingency table
####            Up  Down
#### Amide       0    2
#### not amide  16   50
#phyper(q=c(0:0), m=2, n=66, k=16)

# carbohydrates contingency table
####            Up  Down
#### Carb       7    17
#### not carb   9    35
# phyper(q=c(0:7), m=24, n=44, k=16)

# LINEAR MODEL FOR PM2 - 
#this our scaled pm1 table with values , light condition and 
#metabolites stacked in their repective columns and we use it to fit the model now
#before removal of variannce

# LINEAR MODEL FOR PM2 - 
#this our scaled pm1 table with values , light condition and 
#metabolites stacked in their repective columns and we use it to fit the model now

#before removal of var
head(spm2) 
#grouping the metabolite and light condition together
spm2$group <- factor(paste0(spm2$metabolite, spm2$light_condition))
head(spm2)
Xpm2<-model.matrix(~0+group, data =spm2)

### R base LM IMPLEMENTATION
#remove group prefix - here is one example on how to get contrast which you can then use to extract inference statistics please write a loop that extracts all other contrasts
colnames(Xpm2)=gsub("^group", "", colnames(Xpm2))

#taking out the column names of the Xpm1 matrix for the metabolite names
cHL= Xpm2[,161] #negative control HL 
cLL= Xpm2[,162] # negative controlLL
#the following code is to create a character vector of the column names to calculate the
#differences and differences of the differences using the makeContrasts 
vec <- vector(mode = "character", length = ncol(Xpm2)/2)
vec2 <- vector(mode = "character", length = ncol(Xpm2)/2)
for ( i in 1:ncol(Xpm2))
{
  if (i%%2==0)
  {
    vec[i/2] = paste(colnames(Xpm2)[i], "Negative.ControlLL", sep = "-")
  }
  else
  {
    vec2[i%/%2+1] = paste(colnames(Xpm2)[i], "Negative.ControlHL", sep = "-")
  }
}
#removing "Negative.ControlHL-Negative.ControlHL", "Negative.ControlLL-Negative.ControlLL", 
#"Negative.ControlHL-Negative.ControlHL-Negative.ControlLL-Negative.ControlLL"
#so that they do not introduce NA's in pvalue calculation 
#summary function literally fails if i have NAs
vec = setdiff(vec, c("Negative.ControlLL-Negative.ControlLL"))
vec2 = setdiff(vec2, c("Negative.ControlHL-Negative.ControlHL"))

metabolitevec1 = vector()
for ( i in 1: length(vec))
{
  metabolitevec1[i] = paste( c("("), vec[i], c(")"))
}
metabolitevec2 = vector()

for ( i in 1: length(vec2))
{
  metabolitevec2[i] = paste( c("("), vec2[i], c(")"))
}
metabolitevec <- paste(metabolitevec2, metabolitevec1, sep ="-")

#using default lm function
contpm2=makeContrasts(contrasts =c(vec2, vec, metabolitevec), levels=Xpm2)
head(contpm2)
#using default lm function
fitpm2 <- lm(values~0+group, data = spm2)
summary(fitpm2)
# glht already adjust pvalues
confit=glht(fitpm2, t(contpm2))

summarypm2= summary(confit, test = adjusted(type = "fdr"))
#to see the p value distribution
hist(summarypm2$test$pvalues)

pval =as.data.frame(summarypm2$test$pvalues) #pvalues
lfc= as.data.frame(summarypm2$test$coefficients) #coefficients

lmpm2 <- as.data.frame(cbind(lfc, pval))
lmpm2 <- tibble::rownames_to_column(lmpm2, "contrasts") 
#dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(lmpm2) = c( "contrasts", "lmlfc", "lmpval")
head(lmpm2)
# Volcano plot
save(lmpm2, file = "lmpm2.RData")
p <- ggplot(data=lmpm2, aes(x=lmlfc, y=-log10(lmpval))) + geom_point() + theme_minimal()
p
lmpm2$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
# refers to 1.5 nominal change
lmpm2$diffexpressed[lmpm2$lmlfc > 0.6 & lmpm2$lmpval < 0.05] <- "UP" 
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
lmpm2$diffexpressed[lmpm2$lmlfc < -0.6 & lmpm2$lmpval< 0.05] <- "DOWN"
p <- ggplot(data=lmpm2, aes(x=lmlfc, y=-log10(lmpval), col=diffexpressed)) + geom_point() + theme_minimal()
p

### LIMMA IMPLEMENTATION
#We trasnfrom the data here into limma fromat for 1 gene
spm2limma=data.frame(t(spm2$values))
colnames(spm2limma)=spm2$group
#generate canonical linear model
limfit = lmFit(spm2limma,Xpm2)
#extract contrasts
limcont = contrasts.fit(limfit, contpm2)
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
#which metabolites are down regulated, upregulated or none ?
#which genes are down regulated, upregulated or none ?
downregulated = as.vector(which((dfpm2$diffexpressed=="DOWN")==TRUE))
length(downregulated) #50

down = data.frame(matrix(NA, nrow =50, ncol = 4))
j = 1
for ( i in downregulated)
{
  down[j,] = dfpm2[i,]
  j =j+1
}
down
upregulated = as.vector(which((dfpm2$diffexpressed=="UP")==TRUE))
length(upregulated) #1 only one metabolite is upregulated i.e. 266


dfpm2[226,]

#STEP 4 : HEATMAPS
#For the results obtained from glht, prepare a heat map of the log fold changes 
#in the three contrasts (columns are the contrasts and rows are metabolites (n_met*3) cells,
#do NOT scale the values,
#only cluster the rows, 
#only show log fold change that is linked to a significant p-value 
#set all other LFC values to
lfcpm2 =as.data.frame(cbind(summarypm2$test$coefficients, summarypm2$test$pvalues))
lfcpm2 <- tibble::rownames_to_column(lfcpm2, "contrasts")
colnames(lfcpm2) = c("contrasts", "lfc", "pval" )

for ( i in 1:nrow(lfcpm2))
{
  if (lfcpm2$pval[i] >= 0.05)
  {
    lfcpm2$pval[i]= NA
    lfcpm2$lfc[i] = NA
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heat_pm2 = as.data.frame(cbind(lfcpm2[1:95,], lfcpm2[96:190,], lfcpm2[191:285,]))
colnames(heat_pm2) = c("HLgroups", "HL_lfc", "HLpvalues", "LLgroups", "LL_lfc","LLpvalues" , "HL_LL_dif_groups", "HL_LL_dif_lfc", "HL_LL_dif_pvalues")
#taking out the metabolite names 
metNames_pm2= gsub("HL.+", "", heat_pm2[, "HLgroups"])
#taking out only lfcs
heat_pm2 = as.data.frame(cbind(heat_pm2$LL_lfc, heat_pm2$HL_lfc, heat_pm2$HL_LL_dif_lfc))
rownames(heat_pm2) = metNames_pm2
#Delete rows with complete NAs
heat_pm2 <- heat_pm2[rowSums(is.na(heat_pm2)) != ncol(heat_pm2), ]
colnames(heat_pm2) = c("LL contrasts","HL contrasts", "HL-LL contrasts")


# saving data for volcano plots and heatmaps for pm2
save(lmpm2, dfpm2, heat_pm2, file ="heat_volcano_data.RData")

library("superheat")
superheat(heat_pm2,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "red"), left.label.text.size = 3)
# superheat(heat_pm2,
#           # scale the matrix
#           # change color of missing values
#           heat.na.col = "white", heat.col.scheme =  "red")

#lm implementation after outlier removal

# REPEATING THE ANALYSIS AFTER REMOVAL OF HIGH VARIANCE DATA
# we will only use glht function and ANOVA for this part
#since we kicked yout several groups from the pm2 dataset, we have to also kick them our from the contrast matrix
head(spm2) 
#grouping the metabolite and light condition together
spm2$group <- factor(paste0(spm2$metabolite, spm2$light_condition))
head(spm2)
Xpm2<-model.matrix(~0+group, data =spm2)

### R base LM IMPLEMENTATION
#remove group prefix - here is one example on how to get contrast which you can then use to extract inference statistics please write a loop that extracts all other contrasts
colnames(Xpm2)=gsub("^group", "", colnames(Xpm2))
###Compare the variances
#Create table of sample variances
#calculates the variances of each group in pm1 and pm2 samples are joins them together to form a dataframe
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
save(pm12var, pm1vartest, pm2vartest, file = "pm12var.RData")
limits=boxplot(pm12var[, c("pm1var.values", "pm2var.values")], col ="bisque", names=c("PM1", "PM2"), ylab="Var", ylim=c(0,max(c(pm12var$pm1var.values, pm12var$pm2var.values))*1.2))
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
library(stringr)
#outliers in the pm2 dataset
#smallvariances <-which(pm12var$pm2var.values<2.202078) #170 from her and 22 from outliers. total = 192
outliers <- boxplot(pm12var$pm2var.values, plot=FALSE)$out #taking out the outliers
idx_outliers <- which(pm12var$pm2var.values %in% outliers) #which ones are the outliers
pm2_outliers <- pm12var[idx_outliers,]
pm2_outliers <- as.data.frame(pm2_outliers[, c(3,4)])
#pm2_outliers<- tibble::rownames_to_column(pm2_outliers, "met_indices") 
#sorting the data from highest to lowest variance
pm2_outliers <- pm2_outliers[order(pm2_outliers$pm2var.values, decreasing = TRUE),] #ordering 
#pm2_outliers<- tibble::rownames_to_column(pm2_outliers, "met_indices") 
var_mets <- as.character(pm2_outliers$pm2var.group)
spm2_modified <- spm2
met_spm2 <- as.character(spm2_modified$group)
for ( i in 1:length(var_mets))
{
  idxes <- which(grepl(var_mets[i], spm2_modified$group)==TRUE)
  print(idxes)
  # 1, 8 and 11th entries give 2 groups 
}
var_mets [1]
idx1 <- which(grepl(var_mets[1], spm2_modified$group)==TRUE)
spm2_modified[idx1,] # need to excluse the first 3
var_mets [8]
idx8 <- which(grepl(var_mets[8], spm2_modified$group)==TRUE)
spm2_modified[idx8,] # need to excluse the last 3
var_mets [11]
idx11 <- which(grepl(var_mets[11], spm2_modified$group)==TRUE)
spm2_modified[idx11,]
#need to exclude the first 3

for ( i in 1:length(var_mets))
{
  barlet <- bartlett.test(values~group, data=spm2_modified)
  pvalue <- barlet$p.value
  if (pvalue <= 0.05)
  {
    if ( i %in% c(1, 11))
    {
      #Grepl() return TRUE if the given pattern is present in the vector. Otherwise, it return FALSE
      #grepl('pattern', vector)
      #searching var_met(i) ie each group in the ooutlier in the spm2_modified groups
      # returns the indices where  we find it
      idx= which(grepl(var_mets[1], spm2_modified$group)==TRUE)
      #Butyric.acidLL is the  var_mets[1]
      # since it returns 3 entries for the other group g.Amino.N.Butyric.acidLL, that has Butyric.acidLL in them
      # we do the following to exlude first 3 indices
      idx <- idx[(length(idx)-2):length(idx)] # exclusding 271 273 275 indces
      
      # now removing the Butyric.acidLL group 
      spm2_modified = spm2_modified[-c(idx),]
      # doing the barlett test and taking out the p values
      barlet<- bartlett.test(values~group, data=spm2_modified)
      pvalue <- barlet$p.value
    }
    else if ( i ==8)
    {
      idx= which(grepl(var_mets[8], spm2_modified$group)==TRUE)
      idx <- idx[1:3] # exclusding last 3 indices
      # now removing the Butyric.acidLL group 
      spm2_modified = spm2_modified[-c(idx),]
      # doing the barlett test and taking out the p values
      barlet<- bartlett.test(values~group, data=spm2_modified)
      pvalue <- barlet$p.value
    }
    else
    {
      idx= which(grepl(var_mets[i], spm2_modified$group)==TRUE)
      spm2_modified = spm2_modified[-c(idx),]
      barlet<- bartlett.test(values~group, data=spm2_modified)
      pvalue <- barlet$p.value
    }
    
  }
  else
  {
    break
  }
}
pm1var=aggregate(values~group, data=spm1, var)
pm2var=aggregate(values~group, data=spm2_modified, var)
pm1vartest=bartlett.test(values~group, data=spm1)
pm1vartest$p.value
pm2vartest=bartlett.test(values~group, data=spm2_modified)
pm2vartest$p.value
par(mfrow=c(1,2))
save(pm1var, pm1vartest, pm2var, pm2vartest, file = "pm1pm2var.RData")

boxplot(pm1var[, c("values")], xlab =paste("p = ",signif(c(pm1vartest$p.value),3)), 
        main = "PM1", col = "bisque", ylim=c(0,max(c(pm1var$values))*1.2))
# ggplot(pm1var, aes(x = values)) +  # ggplot function
#   geom_boxplot()+
#   labs(title="PM1",x=paste("p = ",signif(c(pm1vartest$p.value),3)), y = "Var")+
#   theme_classic()
#ypos= y=apply(pm1var[, c("values")], 2, max)
boxplot(pm2var[, c("values")], xlab =paste("p = ",signif(c(pm2vartest$p.value),3)), 
        main = "PM2", labels = paste("p = ",signif(c(pm2vartest$p.value),3)), col = 
          "bisque", ylim=c(0,max(c(pm1var$values))*1.2))


limits=boxplot(pm2var[, c("values")], col ="bisque", main = paste("p = ",signif(c(pm2vartest$p.value),3),sep=""), names=c("PM2"), ylab="Var", xlab = "Figure 15: pm2 variances after outlier removal", ylim=c(0,max(c(pm2var$values))*1.2))
#ypos= y=apply(pm2var[, c("values")], 2, max)
#text(x=c(1:2), y=ypos + 0.1*max(ypos), paste("p = ",signif(c(pm2vartest$p.value),3),sep=""))

c = colnames(Xpm2)
c[161]
c[162]
#taking out the column names of the Xpm2 matrix for the metabolite names
cHL= Xpm2[,161] #negative control HL 
cLL= Xpm2[,162] # negative controlLL
#to extract the groups that are excluded from spm2 dataset
setdif<- setdiff(spm2, spm2_modified)
expelled_mets <- setdif$group
length(expelled_mets)
#54
#means 54 entries were deleted from spm2_modified which means 54/3 = 18 groups
'%ni%'<- Negate('%in%')
vec <- vector(mode = "character", length = ncol(Xpm2)/2)
vec2 <- vector(mode = "character", length = ncol(Xpm2)/2)

for ( i in 1:ncol(Xpm2))
{
  if (i%%2==0)
  {
    if (colnames(Xpm2)[i] %ni% expelled_mets)
    {
      vec[i/2] = paste(colnames(Xpm2)[i], "Negative.ControlLL", sep = "-")
    }
  }
  else
  {
    if (colnames(Xpm2)[i] %ni% expelled_mets)
    {
      vec2[i%/%2+1] = paste(colnames(Xpm2)[i], "Negative.ControlHL", sep = "-")
    }
  }
}
length(vec)
# the empty strings are the expelled groups and we remove them now from our  contrast
# vec = vec[vec != ""]
# length(vec)
# length(vec2)
# vec2 = vec2[vec2 != ""]
# length(vec2)
# removing so Negative.Control entries so that we dont get NAs which actually stops glht function from running
vec = setdiff(vec, c("Negative.ControlLL-Negative.ControlLL"))
vec2 = setdiff(vec2, c("Negative.ControlHL-Negative.ControlHL"))
#adding the brackets 
metabolitevec1 = vector()
for ( i in 1: length(vec))
{
  metabolitevec1[i] = paste( c("("), vec[i], c(")"))
}

metabolitevec2 = vector()

for ( i in 1: length(vec2))
{
  metabolitevec2[i] = paste( c("("), vec2[i], c(")"))
}
metabolitevec <- paste(metabolitevec2, metabolitevec1, sep ="-")

idx <-which(grepl("\\(  \\)", metabolitevec))  # returns the indices of the deleted entries

deleted <- metabolitevec[idx]
metabolitevec<- setdiff(metabolitevec, deleted)
metabolitevec <- setdiff(metabolitevec, c("Negative.ControlHL-Negative.ControlHL-Negative.ControlLL-Negative.ControlLL"))
vec = vec[vec != ""]
length(vec)
length(vec2)
vec2 = vec2[vec2 != ""]
length(vec2)
contpm2=makeContrasts(contrasts =c(vec2, vec, metabolitevec), levels=Xpm2)
#using default lm function
fitpm2 <- lm(values~0+group, data = spm2)
# glht already adjust pvalues
confit2=glht(fitpm2, t(contpm2))
summarypm2= summary(confit2, test = adjusted(type = "hochberg"))
pval =as.data.frame(summarypm2$test$pvalues) #pvalues
lfc= as.data.frame(summarypm2$test$coefficients) #coefficients

lmpm2 <- as.data.frame(cbind(lfc, pval))

lmpm2<- tibble::rownames_to_column(lmpm2, "contrasts") 
#dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(lmpm2) = c( "contrasts", "lmlfc", "lmpval")
head(lmpm2)
p <- ggplot(data=lmpm2, aes(x=lmlfc, y=-log10(lmpval))) + geom_point() + theme_minimal()
p
lmpm2$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
lmpm2$diffexpressed[lmpm2$lmlfc > 0.6 & lmpm2$lmpval< 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
lmpm2$diffexpressed[lmpm2$lmlfc < -0.6 & lmpm2$lmpval < 0.05] <- "DOWN"

save(lmpm2, file =  "lmpm2after.RData")
p <- ggplot(data=lmpm2, aes(x=lmlfc, y=-log10(lmpval), col=diffexpressed)) + geom_point() + theme_minimal()
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
dim(lfcpm2) # 260 instead of 285
#before kicking out the groups we have 285 entries in the groups, and wach contrast column
# ie. LL, HL and HL-LL contrasts had 95 entries but ow we have 260 entries in the group
# and the no. of groups in each contrast is not the same 
# so we cannot use this: heat_pm2 = as.data.frame(cbind(lfcpm2[1:95,], lfcpm2[96:190,], lfcpm2[191:285,]))
#loop to find out how many LL, HL and HL-LL contrast groups are there after removal of certain groups
LL = 0
HL = 0
HL_LL_diff = 0
for ( i in 1:nrow(lfcpm2))
{
  if (grepl("\\LL$", lfcpm2$contrasts[i])==TRUE)
  {
    LL = LL+1
  }
  else if  (grepl("\\HL$", lfcpm2$contrasts[i])==TRUE)
  {
    HL = HL+1
  }
  else
  {
    HL_LL_diff = HL_LL_diff+1
  }
}
LL
HL
HL_LL_diff
HLgroups <- lfcpm2[1:HL,]
LLgroups<-lfcpm2[85:172,]
HL_LL_diff_groups <- lfcpm2[173:260,]
# lets add four rows at the bottom so that we can bind them together
x = matrix(data=NA,nrow=4,ncol=3)
colnames(x) = c("contrasts", "lfc", "pval" )
HLgroups <- rbind(HLgroups, x)
#lets check if the dimensions of the contrast columns are equal
dim(HLgroups) == dim(LLgroups) 
dim(HLgroups)== dim(HL_LL_diff_groups)
heat_pm2 = as.data.frame(cbind(HLgroups, LLgroups, HL_LL_diff_groups))
#for rownames 
metabolites = pm2$metabolite #consists of all the metabolites
# but we have expelled some of them 
length(metabolites)
dim(heat_pm2) 
expelled_mets
#looking at the expelled mets we see that all 6 entries ( ie. 3 samples fom HL and 3 from LL) 
#of b.D.Allose, L.Glucose, Caproic.acid are included in the expelled_mets, so we can expell them altogether 
#from our metabolite vector´
#this gives us the 8 meatabolite names that were expelled including NegativeControlLL and NegativeControlHL
metabolites = setdiff(metabolites, c("Negative.Control", "b.D.Allose", "L.Glucose", "Caproic.acid" ))
#rownames(heat_pm2) = metabolites
colnames(heat_pm2) = c("HLgroups", "HL_lfc", "HLpvalues", "LLgroups", "LL_lfc","LLpvalues" , "HL_LL_dif_groups", "HL_LL_dif_lfc", "HL_LL_dif_pvalues")

head(heat_pm2)
# we can also get the metabolite ames after exclusion from the spm2_modified dataframe # i forgot about it 
met = unique(spm2_modified$metabolite)
met = setdiff(met, c("Negative.Control"))
met == metabolites #true for all
#met = as.data.frame(met)
length(met) # 92, thus there should be 92 entries in the contrast matrix but we have 88
#now need to compare the metabolite in each row, if they dont match then i set NAs to the positions where they dont

# okay i work with the LLcontrast column because that is the column that differs from the other 2
# i will extrast the metabolite name from the each row and check if the corresponding rows in 
#HLcontrast and HL-LL contrast columns have this metabolite, if not i add NAs to those rows

#s you can see the metabolite groups for LLgroups differ from HLgroups and HL_LL_Diff_groups. 
#i need to come up with function or whatever that extracts the metabolite name excluding "LL" 
#from "LLgroups" column and see if "HLgroups" and "HL_LL_Diff_groups"  have this metabolite 
#in the corresponding rows if not i need to add NAs to the those rows so that i have like 
#for example for second row :  NA   NA   NA   a.Keto.Valeric.acidLL-Negative.ControlLL  -4.299579   1.391480e-04    NA NA NA AND SO ON
# correct metabolite mismatches in rows
d = data.frame() # required for rbind
HLs = gsub("HL.+", "", heat_pm2[, "HLgroups"]) # extract met
LLs = gsub("LL.+", "", heat_pm2[, "LLgroups"])

for (i in seq_len(length(HLs))) { # for every met in HLs
  idx = which(LLs==HLs[i]) # search for same met in LLs
  
  if(!identical(idx,integer(0))){ # false if metabolite in HL contrast  is not in LL contrast
    d=rbind(d, c(heat_pm2[i,1:3], heat_pm2[idx, 4:6], heat_pm2[i,7:9]))
    HLs[i]=NA # mask met if found
    LLs[idx]=NA
  }
}

HLs = as.data.frame(HLs)
LLs = as.data.frame(LLs)
# create new rows for leftover mets only found in either HLs or LL
lonelyHLs=cbind(heat_pm2[which(!is.na(HLs)),1:3], matrix(NA, nrow=length(which(!is.na(HLs))), ncol=3), heat_pm2[which(!is.na(HLs)),7:9])
lonelyLLs = cbind(matrix(NA, nrow=length(which(!is.na(LLs))), ncol=3), heat_pm2[which(!is.na(LLs)),4:6], matrix(NA, nrow=length(which(!is.na(LLs))), ncol=3))
lonelyHLs = as.data.frame(lonelyHLs)
lonelyLLs = as.data.frame(lonelyLLs)
d=as.data.frame(rbind(as.matrix(d), as.matrix(lonelyHLs), as.matrix(lonelyLLs)))
# add mets as rownames
contrast_values=d[,c(5,2,8)]
rownames=gsub("HL.+", "", d[,1]) # taking out the metabolites names again
rownames[which(is.na(rownames))]=gsub("LL.+","",d[(nrow(d)-length(which(is.na(rownames)))+1):nrow(d),4])
rownames(contrast_values)=rownames
# transform into superheat compatible format
#colnames(contrast_values)=colnames(heat_pm2) 
heat_pm2_precursor <- contrast_values[rowSums(is.na(contrast_values)) != ncol(contrast_values), ] # remove rows with only NAs
heat_pm2 = apply(heat_pm2_precursor,2,as.numeric) # convert <NA> into NA
rownames(heat_pm2)=rownames(heat_pm2_precursor) 
colnames(heat_pm2) = c("LL contrasts", "HL contrasts", "HL-LL contrasts")# add lost rownames
save(heat_pm2, file ="heatpm2after.RData")
library("superheat")
superheat(heat_pm2,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "red"), left.label.text.size = 2.5)

 #enrichment test
#importing the compound list
pm2compounds = read.xlsx("Data/PM01_PM2A_compound list.xlsx", sheetName = "PM2A", header=FALSE)
pm2compounds$X3 = conv_metnames(pm2compounds$X3)
pm2compounds$X4 = conv_metnames(pm2compounds$X4)
head (lmpm2)

downregulated = as.vector(which((lmpm2$diffexpressed=="DOWN")==TRUE))
length(downregulated) #: 27

down = data.frame(matrix(NA, nrow =27, ncol = 4))
j = 1
for ( i in downregulated)
{
  down[j,] = lmpm2[i,]
  j =j+1
}
head(down)
dim(down) #51
downMets <-downdiff(down)


df_downLL= CsourceExtract(downLLmets)
df_downHL= CsourceExtract(downHLmets)

#Upregulation
upregulated = as.vector(which((lmpm1$diffexpressed=="UP" & lmpm1$lmpval <0.05)==TRUE))
length(upregulated)
up = data.frame(matrix(NA, nrow =16, ncol = 4))
j = 1
for ( i in upregulated)
{
  up[j,] = lmpm1[i,]
  j =j+1
}
up
# up consists of ohly those metabolites in the third colum of the heat map ie. the interaction component
upMets <- gsub("HL(-.+).+", "", up$X1)
upMets <- gsub(".+ ", "", upMets)
#upMets <- gsub('.+([a-z]+)HL.+', '(\\1)', up$X1)
CsourceUp = vector(length = length(upMets), mode = "character")
df_up <- CsourceExtract(upMets)


Noreg = as.vector(which((lmpm1$diffexpressed=="NO")==TRUE))
length(Noreg)
No= data.frame(matrix(NA, nrow =217, ncol = 4))
j = 1
for ( i in Noreg)
{
  No[j,] = lmpm1[i,]
  j =j+1
}
No = as.data.frame(No)
head(No)
dim(No) #217
NoMets <- Nodiff(No)
df_noLL= CsourceExtract(NoMets$LLmets)
df_noHL= CsourceExtract(NoMets$HLmets)
df_noHL_LL = CsourceExtract(NoMets$HL_LLmets)


#no regulation for seaprated ones
tableHL_LL <- as.data.frame(table(df_noHL_LL$Csource))
tablenoLL <-  as.data.frame(table(df_noLL$Csource))
tablenoHL <- as.data.frame(table(df_noHL$Csource))
#ureg for separated ones
tableup <- as.data.frame(table(df_up$Csource))
#downreg for the separated ones
tabledownLL <- as.data.frame(table(df_downLL$Csource))
tabledownHL <- as.data.frame(table(df_downHL$Csource))
#creating contingecy matrix
contingency_mat <- matrix(NA, nrow=8, ncol= 6)
colnames(contingency_mat) = c( "DownregulationLL", "DownregulationHL", "UpregulationHL_LL",
                               "NoRegulationLL", "NoRegulationHL", "NoregulationHL_LL")
rownames (contingency_mat) = c("Alcohol", "Amide", "Amine", "Amino acid", 
                               "Carbohydrate", "Carboxylic acid", "Ester", "Fatty acid")
contingency_mat[,1] = c(0, tabledownLL$Freq)
contingency_mat[,2] =  c(0, 1,1,0, 4,3,0,3)
contingency_mat[,3] =  c(0,0, 0, 1, 7, 8, 0, 0)
contingency_mat[,4] = c(2, 0,0,12,25,16,0,0)
contingency_mat[,5] = c(2,0,1,16,34,29,1,0)
contingency_mat[,6] =tableHL_LL$Freq

contingency_mat
hyperup <- phyperfunction(contingency_mat, "UpregulationHL_LL")
hyperdownLL <- phyperfunction(contingency_mat, "DownregulationLL")
hyperdownHL <- phyperfunction(contingency_mat, "DownregulationHL")
hypernoLL <- phyperfunction(contingency_mat, "NoRegulationLL")
hypernoHL <-  phyperfunction(contingency_mat, "NoRegulationHL")
hyperno <- phyperfunction(contingency_mat, "NoregulationHL_LL")

library(ggpubr)

plotfunction(hyperup)
plotfunction(hyperdownHL)
plotfunction(hyperdownLL)



#which genes are down regulated, upregulated or none ?
downregulated = as.vector(which((lmpm2$diffexpressed=="DOWN")==TRUE))
length(downregulated) #: 77

down = data.frame(matrix(NA, nrow =27, ncol = 4))
j = 1
for ( i in downregulated)
{
  down[j,] = lmpm2[i,]
  j =j+1
}
head(down)
downMets = gsub("LL-Negative.Control.+", "", down$X1)
downMets = gsub("HL-.+", "", downMets)
downMets = gsub (".+ ", "", downMets)

upregulated = as.vector(which((lmpm2$diffexpressed=="UP")==TRUE))
length(upregulated)#7
up = data.frame(matrix(NA, nrow =7, ncol = 4))
j = 1
for ( i in upregulated)
{
  up[j,] = lmpm2[i,]
  j =j+1
}
head(up)
upMets <- gsub("HL(-.+).+", "", up$X1)
upMets <- gsub(".+ ", "", upMets)
#upMets <- gsub('.+([a-z]+)HL.+', '(\\1)', up$X1)


CsourceUp = vector(length = length(upMets), mode = "character")

for ( i in 1:length(upMets))
{
  if (upMets[i] %in% pm2compounds$X3)
  {
    idx = which(pm2compounds$X3==upMets[i]) # taking out the index of the metabolite pm1compounds
    CsourceUp[i] = pm2compounds$X4[idx]
    i =i+1
  }
}
#CsourceUp
CsourceUp <- gsub("C.Source+\\.\\.", "", CsourceUp)
CsourceUp

df_up <- as.data.frame(cbind(upMets, CsourceUp))

CsourceDown = vector(length = length(downMets), mode = "character")

for ( i in 1:length(downMets))
{
  if (downMets[i] %in% pm2compounds$X3)
  {
    idx = which(pm2compounds$X3==downMets[i]) # taking out the index of the metabolite pm1compounds
    CsourceDown[i] = pm2compounds$X4[idx]
    i =i+1
  }
}
CsourceDown
CsourceDown <- gsub("C.Source+\\.\\.", "", CsourceDown)
CsourceDown
df_down <- as.data.frame(cbind(downMets, CsourceDown))

Noreg = as.vector(which((lmpm2$diffexpressed=="NO")==TRUE))
length(Noreg) #226

No= data.frame(matrix(NA, nrow =226, ncol = 4))
j = 1
for ( i in Noreg)
{
  No[j,] = lmpm2[i,]
  j =j+1
}
No = as.data.frame(No)
head(No)

tableup <- as.data.frame(table(df_up$CsourceUp))
tableup
tabledown <- as.data.frame(table(df_down$CsourceDown))
tabledown
tableno <- as.data.frame(table(df_No$CsourceNo))
tableno
#creating contingecy matrix
unique(df_down$CsourceDown) #7
unique(df_up$CsourceUp) #3
#length(unique(df_No$CsourceNo)) #8

# contingency_mat <- matrix(NA, nrow=4, ncol= 2)
# 
# colnames(contingency_mat) <- c("Upregulation", "Downregulation" )
# rownames(contingency_mat) = c("Alcohol", "Carbohydrate", "Carboxylic acid", "Polymer")
# contingency_mat [,2]<- tabledown$Freq
# contingency_mat[,1] = c(0, tableup$Freq)
# contingency_mat

# probabilitiesUp = list()
# q = list()
# for (class in row.names(contingency_mat)){
#   #goes from 0 to total no. of upregulation in that class
#   q = c(0:contingency_mat[class, "Upregulation"]) 
#   p = phyper(
#     q = q-1,
#     m = sum(contingency_mat[class, c("Upregulation", "Downregulation")]), 
#     n = sum(contingency_mat[which(row.names(contingency_mat) != class),
#                             c("Upregulation", "Downregulation")]),
#     k = sum(contingency_mat[, "Upregulation"]), lower.tail = FALSE,
#   )
#   probs = data.frame("number of upregulations" = q,
#                      "probability" = p)
#   #print (probs)
#   probabilitiesUp[[class]] = probs
# }
# probabilitiesUp
# 
# 
# plotlist = list()
# for ( i in 1:4)
# {
#   data <- as.data.frame(probabilitiesUp[[i]][c(1,2)])
#   data <- data.frame( x = data[,1], y = data[,2])
#   g = ggplot(data, aes(x=factor(x), y=y)) +
#     theme(axis.text=element_text(size=14),
#           axis.title=element_text(size=18,face="bold"),
#           axis.title.x=element_text(margin=margin(20,0,0,0)),
#           axis.title.y=element_text(margin=margin(0,20,0,0))
#     ) +
#     geom_bar(stat="identity", fill=ifelse(data$x < 12,
#                                           rgb(52, 73, 94, maxColorValue=255),
#                                           rgb(231, 76, 60, maxColorValue=255)),
#              colour="black") + labs(x = "No. of upregulation", y = "Probability")
#   plotlist[[i]] = g
#   i = i+1
#   
# }
# library(ggpubr)
# figure <- ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]]
#                     , labels = c("Alcohol", "Carbohydrate", "Carboxylic acid", "Polymer"), ncol = 2, nrow = 2)
# figure
# 
# probabilitiesDown = list()
# for (class in row.names(contingency_mat)){
#   #goes from 0 to total no. of upregulation in that class
#   q = c(0:contingency_mat[class, "Downregulation"]) 
#   p = phyper(
#     q = q-1,
#     
#     m = sum(contingency_mat[class, c("Upregulation", "Downregulation")]), 
#     n = sum(contingency_mat[which(row.names(contingency_mat) != class),
#                             c("Upregulation", "Downregulation")]),
#     
#     k = sum(contingency_mat[, "Downregulation"]), lower.tail = FALSE
#     #sum of uppregulation column
#     
#   )
#   probs = data.frame("number of Downregulation" = q,
#                      "probability" = p)
#   #print (probs)
#   probabilitiesDown[[class]] = probs
# }
# probabilitiesDown
# 
# #example
# #data <- data.frame( x = probabilitiesDown$Carbohydrate[,1], y = probabilitiesDown$Carbohydrate[,2] )
# 
# plotlist = list()
# for ( i in 1:4)
# {
#   data <- as.data.frame(probabilitiesDown[[i]][c(1,2)])
#   data <- data.frame( x = data[,1], y = data[,2])
#   g = ggplot(data, aes(x=factor(x), y=y)) +
#     theme(axis.text=element_text(size=14),
#           axis.title=element_text(size=18,face="bold"),
#           axis.title.x=element_text(margin=margin(20,0,0,0)),
#           axis.title.y=element_text(margin=margin(0,20,0,0))
#     ) +
#     geom_bar(stat="identity", fill=ifelse(data$x < 12,
#                                           rgb(52, 73, 94, maxColorValue=255),
#                                           rgb(231, 76, 60, maxColorValue=255)),
#              colour="black") + labs(x = "No. of Downrregulation", y = "P-value")
#   plotlist[[i]] = g
#   i = i+1
#   
# }
# figure <- ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]]
#                     , labels = c("Alcohol", "Carbohydrate", "Carboxylic acid", "Polymer"), ncol = 2, nrow = 2)
# figure
# 
# 
# colnames(contingency_mat) = c("Upregulation", "Downregulation", "NoRegulation")
# rownames (contingency_mat) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "Fatty acid")
# contingency_mat[,3] = tableno$Freq
# contingency_mat[,2] =  c(0, tabledown$Freq)
# contingency_mat[,1] = c(0,0, 0, 1, 7, 8, 0, 0) # the entry from table down


contingency_mat1 <- matrix(NA, nrow=8, ncol= 3)
colnames(contingency_mat1) = c("Upregulation", "Downregulation", "NoRegulation")
rownames (contingency_mat1) = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer")
contingency_mat1[,3] = tableno$Freq
contingency_mat1[,2] =  c(7, 0, 0,0, 8, 8, 0, 4)
contingency_mat1[,1] = c(0,0, 0,0, 3, 3,  0, 1) # the entry from table down
contingency_mat1

fisher <-fisher.test(contingency_mat, simulate.p.value = TRUE)

probabilitiesUp1 = list()
q = list()
for (class in row.names(contingency_mat1)){
  #goes from 0 to total no. of upregulation in that class
  q = c(0:contingency_mat1[class, "Upregulation"]) 
  p = phyper(
    q = q,
    m = sum(contingency_mat1[class, c("Upregulation", "Downregulation", "NoRegulation")]), 
    n = sum(contingency_mat1[which(row.names(contingency_mat1) != class),
                             c("Upregulation", "Downregulation", "NoRegulation")]),
    k = sum(contingency_mat1[, "Upregulation"])
  )
  probs = data.frame("number of upregulations" = q,
                     "probability" = p)
  #print (probs)
  probabilitiesUp1[[class]] = probs
}
probabilitiesUp1


plotlist = list()
for ( i in 1:8)
{
  data <- as.data.frame(probabilitiesUp1[[i]][c(1,2)])
  data <- data.frame( x = data[,1], y = data[,2])
  g = ggplot(data, aes(x=factor(x), y=y)) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"),
          axis.title.x=element_text(margin=margin(20,0,0,0)),
          axis.title.y=element_text(margin=margin(0,20,0,0))
    ) +
    geom_bar(stat="identity", fill=ifelse(data$x < 12,
                                          rgb(52, 73, 94, maxColorValue=255),
                                          rgb(231, 76, 60, maxColorValue=255)),
             colour="black") + labs(x = "No. of upregulation", y = "p-values")
  plotlist[[i]] = g
  i = i+1
  
}
library(ggpubr)
figure <- ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]], plotlist[[5]], plotlist[[6]], plotlist[[7]], plotlist[[8]]
                    , labels = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer"), ncol = 2, nrow = 4)
figure

probabilitiesDown1 = list()
for (class in row.names(contingency_mat1)){
  #goes from 0 to total no. of upregulation in that class
  q = c(0:contingency_mat1[class, "Downregulation"]) 
  p = phyper(
    q = q,
    
    m = sum(contingency_mat1[class, c("Upregulation", "Downregulation", "NoRegulation")]), 
    n = sum(contingency_mat1[which(row.names(contingency_mat1) != class),
                             c("Upregulation", "Downregulation", "NoRegulation")]),
    
    k = sum(contingency_mat1[, "Downregulation"])
    #sum of uppregulation column
    
  )
  probs = data.frame("number of Downregulation" = q,
                     "p-value" = p)
  #print (probs)
  probabilitiesDown1[[class]] = probs
}
probabilitiesDown1

#example
#data <- data.frame( x = probabilitiesDown$Carbohydrate[,1], y = probabilitiesDown$Carbohydrate[,2] )

plotlist = list()
for ( i in 1:8)
{
  data <- as.data.frame(probabilitiesDown1[[i]][c(1,2)])
  data <- data.frame( x = data[,1], y = data[,2])
  g = ggplot(data, aes(x=factor(x), y=y)) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"),
          axis.title.x=element_text(margin=margin(20,0,0,0)),
          axis.title.y=element_text(margin=margin(0,20,0,0))
    ) +
    geom_bar(stat="identity", fill=ifelse(data$x < 12,
                                          rgb(52, 73, 94, maxColorValue=255),
                                          rgb(231, 76, 60, maxColorValue=255)),
             colour="black") + labs(x = "No. of Downrregulation", y = "P-value")
  plotlist[[i]] = g
  i = i+1
  
}
library(ggpubr)
figure <- ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]], plotlist[[5]], plotlist[[6]], plotlist[[7]], plotlist[[8]],
                    labels = c("Alcohol","Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "polymer"), ncol = 2, nrow = 4)
figure

fisher <-fisher.test(contingency_mat, simulate.p.value = TRUE)


#do a simple anova with gene-expression as the response and "light_condition" and "metabolite" as the predictors
#Null hypothesis: Means of the sample groups are equal
#Alternative hypothesis: At least the mean of one group is not equal
#If the null hypothesis is rejected, we need to determine which of the 
#group means differ from each other.
head(pm1new)
head(spm1)

#anova with intercept
#pm1_anova <- aov(values~ metabolite*light_condition, data= spm1) 
# thsd=TukeyHSD(pm1_anova)
# met_cont <- thsd$metabolite
# mets <- row.names(met_cont)
# # the following extracts the metabolite-Negativecontrol contrasts
# met_contrasts <- mets[which(grepl("-Negative.Control", mets)== TRUE)]
# met_light_met <- thsd$`metabolite:light_condition

# anova with zero intercept as we did with lm

pm1_anova_0 <- aov(values~light_condition*metabolite, data = spm1) 
summary(pm1_anova_0)
thsd0=TukeyHSD(pm1_anova_0)
group_cont <-as.data.frame(thsd0$`light_condition:metabolite`)
group_cont <-  tibble::rownames_to_column(group_cont, "MetNames")
#group_cont <- group_cont[group_cont$`p adj`<0.05,]
#write.csv(group_cont, "/home/mamata/Desktop/project/group_cont.csv", row.names=FALSE)
sig_metLL<-group_cont[which(grepl("LL:.+-LL:Negative.Control", group_cont$MetNames)== TRUE),]
length(which(heat_pm1$`LL contrasts`<=0)) #12
dim(sig_metLL) #13
#sig_metLL
#View(heat_pm1)
sig_metHL <- group_cont[which(grepl("HL:.+-HL:Negative.Control", group_cont$MetNames, ignore.case = TRUE)== TRUE),]
dim(sig_metHL)
length(which(heat_pm1$`HL contrasts`<=0)) #40 metHL-Negative.ControlHL values are significant in lm
dim(heat_pm1) #compare to the no. of significant groups 
head(summarypm1$test$coefficients) # the coefficients from the glht function


for ( i in 1:nrow(sig_metLL))

  {
  if (sig_metLL$`p adj`[i] >= 0.05)
  {
    sig_metLL$`p adj`[i]= NA
    sig_metLL$diff[i] =NA
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm1LL = as.data.frame(sig_metLL[, c(2,5)])

#taking out the metabolite names 
metNames_LL= gsub("LL:Negative.Control", "", sig_metLL[, "MetNames"])
metNames_LL = gsub("LL:", "", metNames_LL)
#taking out only lfcs


rownames(heatpm1LL) = metNames_LL
heatpm1LL <- tibble::rownames_to_column(heatpm1LL, "MetNames")


for ( i in 1:nrow(sig_metHL))
{
  if (sig_metHL$`p adj`[i] >= 0.05)
  {
    sig_metHL$`p adj`[i]= NA
    sig_metHL$diff[i] =NA
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm1HL = as.data.frame(sig_metHL[c(2,5)])

#taking out the metabolite names 
metNames_HL= gsub("HL:Negative.Control", "", sig_metHL[, "MetNames"])
metNames_HL= gsub("HL:", "", metNames_HL)
#taking out only lfcs

rownames(heatpm1HL) = metNames_HL
heatpm1HL <- tibble::rownames_to_column(heatpm1HL, "MetNames")


heatpm1anova <- cbind (heatpm1LL[,2], heatpm1HL[,2])
rownames(heatpm1anova)= heatpm1HL[,"MetNames"]
colnames(heatpm1anova) = c("LL contrasts" ,"HL contrasts" )

#Delete rows with complete NAs
heatpm1anova<- as.data.frame(heatpm1anova[rowSums(is.na(heatpm1anova)) != ncol(heatpm1anova), ])



# saving data for volcano plots and heatmaps for pm2
save(lmpm2, dfpm2, heat_pm2, file ="heat_volcano_data.RData")

library("superheat")
superheat(heatpm1anova,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "red"), left.label.text.size = 3)



pm2_anova_0 <- aov(values~0+light_condition*metabolite, data = spm2_modified) 
summary(pm2_anova_0)
thsd0=TukeyHSD(pm2_anova_0)
group_cont <-as.data.frame(thsd0$`light_condition:metabolite`)
group_cont <-  tibble::rownames_to_column(group_cont, "MetNames")
  
#group_cont <- group_cont[group_cont$`p adj`<0.05,]
#write.csv(group_cont, "/home/mamata/Desktop/project/group_cont.csv", row.names=FALSE)
sig_metLL<-group_cont[which(grepl("LL:.+-LL:Negative.Control", group_cont$MetNames)== TRUE),]
length(which(heat_pm1$`LL contrasts`<=0)) #12
dim(sig_metLL) #13
#sig_metLL
#View(heat_pm1)
sig_metHL <- group_cont[which(grepl("HL:.+-HL:Negative.Control", group_cont$MetNames, ignore.case = TRUE)== TRUE),]
dim(sig_metHL)

sig_metLL[is.na(sig_metLL)] = 100

for ( i in 1:nrow(sig_metLL))
{
  if (sig_metLL$`p adj`[i] >= 0.05)
  {
    sig_metLL$`p adj`[i]= NA
    sig_metLL$diff[i] =NA
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm2LL = as.data.frame(sig_metLL[, c(2,5)])

#taking out the metabolite names 
metNames_LL= gsub("LL:Negative.Control", "", sig_metLL[, "MetNames"])
metNames_LL = gsub("LL:", "", metNames_LL)
#taking out only lfcs


rownames(heatpm2LL) = metNames_LL
heatpm2LL <- tibble::rownames_to_column(heatpm2LL, "MetNames")
sig_metHL[is.na(sig_metHL)] = 100

for ( i in 1:nrow(sig_metHL))
{
  if (sig_metHL$`p adj`[i] >= 0.05)
  {
    sig_metHL$`p adj`[i]= NA
    sig_metHL$diff[i] =NA
  }
}
#taking out LFCs from the second column and arranging them in 3 columns for LL, HL and HL-LL contrast groups
heatpm2HL = as.data.frame(sig_metHL[c(2,5)])

#taking out the metabolite names 
metNames_HL= gsub("HL:Negative.Control", "", sig_metHL[, "MetNames"])
metNames_HL= gsub("HL:", "", metNames_HL)
#taking out only lfcs

rownames(heatpm2HL) = metNames_HL
heatpm2HL <- tibble::rownames_to_column(heatpm2HL, "MetNames")

heatpm2HL[,1] == heatpm2LL[,1]

heatpm2anova <- cbind (heatpm2LL[,2], heatpm2HL[,2])
rownames(heatpm2anova)= heatpm2HL[,"MetNames"]
colnames(heatpm2anova) = c("LL contrasts" ,"HL contrasts" )

#Delete rows with complete NAs
heatpm2anova<- as.data.frame(heatpm2anova[rowSums(is.na(heatpm2anova)) != ncol(heatpm2anova), ])



# saving data for volcano plots and heatmaps for pm2
save(lmpm2, dfpm2, heat_pm2, file ="heat_volcano_data.RData")

library("superheat")
superheat(heatpm2anova,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white", heat.pal = c("blue", "red"), left.label.text.size = 3)

# STEP 5 : BUILD SUPPORT VE`CTOR MACHINES AND RANDOM FOREST CLASSIFIERS TO SEPARATE
#THE METABOLITES WHICH HAVE EFFECT FROM THE ONES WHICH DONOT


