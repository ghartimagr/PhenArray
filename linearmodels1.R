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
#using bwplot from lattice library: values in y axis and positions in x axis
#splitting the axes into LL and HL, coloring by the positions

#boxplot with labels
#bwplot(newpmcontrol$value~newpmcontrol$position|newpmcontrol$variable, 
 #      col = as.numeric(as.factor(newpmcontrol$position)), ylab= "Values", xlab= "Position",scales =list(newpmcontrol$position, cex = 0.5, rot = 90))

#boxplot without labels
bwplot(newpmcontrol$value~newpmcontrol$position|newpmcontrol$variable, 
       col = as.numeric(as.factor(newpmcontrol$position)), ylab= "Values", xlab= "Position",scales=list(x=list(draw=FALSE))) #scales =list(newpmcontrol$position, cex = 0.5, rot = 90))
#we see the position effect in the HL samples, the values on the left hand side
# of the  plate seem to be higher, it might be because of the position of the light source
#they might have positioned the light in the left hand corner that is why the samples there
#received more light

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

#bpxplot with metabolite names
#bwplot(pm1stacked$values~pm1stacked$metabolite|pm1stacked$light_condition, col= 
#rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites", scales=list(x=list(pm1stacked$metabolite), cex = 0.7, rot = 90))

#boxplot without metabolite names
#boxplot(Tpm1stacked$values~Tpm1stacked$sample , col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")
bwplot(pm1stacked$values~pm1stacked$metabolite|pm1stacked$light_condition, 
       col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

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
#using bwplot from the lattice library
#values of the metabolites are in the y axis and metabolites are in the x-axis
#the values are grouped together by HL and LL
bwplot(pm2stacked$values~pm2stacked$metabolite|pm2stacked$light_condition, 
       col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

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
       col = as.numeric(as.factor(snewpmcontrol$position)), ylab= "Values", xlab= "Position", scales=list(x=list(draw=FALSE)))

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
bwplot(spm1$values~spm1$metabolite | spm1$light_condition, col= rainbow(ncol(pm1new)), ylab = "values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

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
bwplot(spm2$values~spm2$metabolite | spm2$light_condition, col= rainbow(ncol(pm1new)), 
       ylab = "values", xlab = "metabolites", scales=list(x=list(draw=FALSE)))

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
pval =as.data.frame(summarypm1$test$pvalues) #pvalues
lfc= as.data.frame(summarypm1$test$coefficients) #coefficients

lmpm1 <- as.data.frame(cbind(lfc, pval))
lmpm1 <- tibble::rownames_to_column(lmpm1, "contrasts") 
#dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(lmpm1) = c( "contrasts", "lmlfc", "lmpval")
head(lmpm1)
# Volcano pl

p <- ggplot(data=lmpm1, aes(x=lmlfc, y=-log10(lmpval))) + geom_point() + theme_minimal()
p


lmpm1$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
lmpm1$diffexpressed[lmpm1$lmlfc > 0.6 & lmpm1$lmpval < 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
lmpm1$diffexpressed[lmpm1$lmlfc < -0.6 & lmpm1$lmpval< 0.05] <- "DOWN"
p <- ggplot(data=lmpm1, aes(x=lmlfc, y=-log10(lmpval), col=diffexpressed)) + geom_point() + theme_minimal()
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

down = data.frame(matrix(NA, nrow =52, ncol = 4))
j = 1
for ( i in downregulated)
{
  down[j,] = dfpm1[i,]
  j =j+1
}
down

upregulated = as.vector(which((dfpm1$diffexpressed=="UP")==TRUE))
length(upregulated)
up = data.frame(matrix(NA, nrow =16, ncol = 4))
j = 1
for ( i in upregulated)
{
  up[j,] = dfpm1[i,]
  j =j+1
}
up

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially 
#expressed (NA in case they are not)
# 
# dfpm1$label <- NA
# dfpm1$label[dfpm1$diffexpressed != "NO"] <- dfpm1$contrasts[dfpm1$diffexpressed != "NO"]
# 
# ggplot(data=dfpm1, aes(x=limmaLFC, y=-log10(limmapval), col=diffexpressed)) + 
#   geom_point() + 
#   theme_minimal() +
#   geom_text(data = dfpm1, aes(label = label), check_overlap = TRUE, size = 3, hjust = 0.5)

# LINEAR MODEL FOR PM2 - 

#this our scaled pm1 table with values , light condition and 
#metabolites stacked in their repective columns and we use it to fit the model now

#lm implementation
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
pm12var=data.frame(pm1var=aggregate(values~group, data=spm1, var), pm2var=aggregate(values~group, data=spm2, var))
#test if sample variances are equal for pm1
pm1vartest=bartlett.test(values~group, data=spm1)
#test if sample variances are equal for pm2
pm2vartest=bartlett.test(values~group, data=spm2)
#plot boxplot of sample variances with test results
limits=boxplot(pm12var[, c("pm1var.values", "pm2var.values")], names=c("PM1", "PM2"), ylab="Var", ylim=c(0,max(c(pm12var$pm1var.values, pm12var$pm2var.values))*1.2))
ypos= y=apply(pm12var[, c("pm1var.values", "pm2var.values")], 2, max)
text( 
  x=c(1:2), 
  y=ypos + 0.1*max(ypos), 
  paste("p = ",signif(c(pm1vartest$p.value, pm2vartest$p.value),3),sep="")  
)
#taking out the column names of the Xpm2 matrix for the metabolite names
cHL= Xpm2[,161] #negative control HL 
cLL= Xpm2[,162] # negative controlLL

#contpm=makeContrasts(contrasts =c("MaltoseHL-Negative.ControlHL", "MaltoseLL-Negative.ControlLL", "(MaltoseHL-Negative.ControlHL)-(MaltoseLL-Negative.ControlLL)"), levels=Xpm1)
#head(contpm)

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

#metabolitevec <- setdiff(metabolitevec, c("Negative.ControlHL-Negative.ControlHL-Negative.ControlLL-Negative.ControlLL"))

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
p <- ggplot(data=lmpm2, aes(x=lmlfc, y=-log10(lmpval), col=diffexpressed)) + geom_point() + theme_minimal()
p

#which metabolites are down regulated, upregulated or none ?

# lmpm2$label <- NA
# lmpm2$label[lmpm2$diffexpressed != "NO"] <- lmpm2$contrasts[lmpm2$diffexpressed != "NO"]
# 
# ggplot(data=lmpm2, aes(x=lmlfc, y=-log10(lmpval), col=diffexpressed)) + 
#   geom_point() + 
#   theme_minimal() +
#   geom_text(data = lmpm2, aes(label = label), check_overlap = TRUE, size = 3 , vjust =0.5 ,hjust = 0.5 )
# 
# #limma implementation

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

#contpm=makeContrasts(contrasts =c("MaltoseHL-Negative.ControlHL", "MaltoseLL-Negative.ControlLL", "(MaltoseHL-Negative.ControlHL)-(MaltoseLL-Negative.ControlLL)"), levels=Xpm1)
#head(contpm)

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

#metabolitevec <- setdiff(metabolitevec, c("Negative.ControlHL-Negative.ControlHL-Negative.ControlLL-Negative.ControlLL"))

contpm2=makeContrasts(contrasts =c(vec2, vec, metabolitevec), levels=Xpm2)
head(contpm2)
#using default lm function
fitpm2 <- lm(values~0+group, data = spm2)
summary(fitpm2)

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
# Volcano plots indicate the fold change (either positive or negative) in the x axis 
#and a significance value (such as the p-value or the adjusted p-value, i.e. limmapval) in the y axis

#The ‘limmapval’ columns contains the corrected pvalues; 
#these must be converted to the negative of their logarithm base 10 before plotting, 
#i.e -log10(p-value) or -log10(limmapval).
#Since volacno plots are scatter plots, we can use geom_point() to generate one with ggplot2
#the higher the position of a point, the more significant its value is (y axis).
#Points with positive fold change values (to the right) are up-regulated and 
#points with negative fold change values (to the left) are down-regulated (x axis).


p <- ggplot(data=dfpm2, aes(x=limmaLFC, y=-log10(limmapval))) + geom_point() + theme_minimal()
p


dfpm2$diffexpressed <- "NO"
# if limmaLFC > 0.6 and pvalue < 0.05, set as "UP"  for up regualtion
dfpm2$diffexpressed[dfpm1$limmaLFC > 0.6 & dfpm2$limmapval < 0.05] <- "UP"
# if limmaLFC < -0.6 and pvalue < 0.05, set as "DOWN" for down regulation
dfpm2$diffexpressed[dfpm2$limmaLFC < -0.6 & dfpm2$limmapval < 0.05] <- "DOWN"
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


# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially 
#expressed (NA in case they are not)
# 
# dfpm2$label <- NA
# dfpm2$label[dfpm2$diffexpressed != "NO"] <- dfpm2$contrasts[dfpm2$diffexpressed != "NO"]
# 
# ggplot(data=dfpm2, aes(x=limmaLFC, y=-log10(limmapval), col=diffexpressed)) + 
#   geom_point() + 
#   theme_minimal() +
#   geom_text(data = dfpm2, aes(label = label), check_overlap = TRUE, size = 3 , vjust =0.5 ,hjust = 0.5 )


#STEP 4: HEATMAPS
#For the results obtained from glht, prepare a heat map of the log fold changes 
#in the three contrasts (columns are the contrasts and rows are metabolites (n_met*3) cells,
#do NOT scale the values,
#only cluster the rows, 
#only show log fold change that is linked to a significant p-value 
#set all other LFC values to NA and remove rows(metabolites) that only have NA values). 

# PM1
#extracting the log fold chnages and pvalues for pm1

lfcpm1 =as.data.frame(summarypm1$test$coefficients)
heat_pm1 = as.data.frame(cbind(lfcpm1[1:95,], lfcpm1[96:190,], lfcpm1[191:285,]))
metabolites = pm1$metabolite
metabolites = setdiff(metabolites, c("Negative.Control"))
rownames(heat_pm1) = metabolites
colnames(heat_pm1) = c("LLcontrasts", "HLcontrasts", "HL-LL contrasts")
pheatmap(heat_pm1, cluster_rows = TRUE, cluster_cols = FALSE)

lfcpm1 =as.data.frame(cbind(summarypm1$test$coefficients, summarypm1$test$pvalues))
lfcpm1 <- tibble::rownames_to_column(lfcpm1, "contrasts")
colnames(lfcpm1) = c("contrasts", "lfc", "pval" )

for ( i in 1:nrow(lfcpm1))
{
  if (lfcpm1$pval[i] >= 0.05)
  {
    lfcpm1$pval[i]= NA
    lfcpm1$lfc[i] = NA
  }
}

heat_pm1 = as.data.frame(cbind(lfcpm1[1:95,2], lfcpm1[96:190,2], lfcpm1[191:285,2]))
metabolites = pm1$metabolite
metabolites = setdiff(metabolites, c("Negative.Control"))
rownames(heat_pm1) = metabolites
colnames(heat_pm1) = c("LLcontrasts", "HLcontrasts", "HL-LL contrasts")
#heat_pm1 <- na.omit(heat_pm1)
#Delete rows with complete NAs
heat_pm1 <- heat_pm1[rowSums(is.na(heat_pm1)) != ncol(heat_pm1), ]
#heat_pm1<- filter(heat_pm1, rowSums(is.na(heat_pm1)) != ncol(heat_pm1))
#heat_pm1[is.na(heat_pm1)] <- as.double("NA")
#heat_pm1[is.na(heat_pm1)] <- ""
pheatmap(heat_pm1, na_col="white", scale = "none", cluster_cols = FALSE, cluster_rows = TRUE) #  not working

# llcontrast = as.data.frame(lmpm1$lmlfc[which(as.numeric(lmpm1$indices)<=95)])
# #take out the signifcant columns and 
# hlcontrast = which(as.numeric(lmpm1$indices)>95 & as.numeric(lmpm1$indices)<=190)
# hl_ll_diffcontrast = which(as.numeric(lmpm1$indices)>190)
# 
# significant_pm1 = as.data.frame(c(llcontrast, hlcontrast, hl_ll_diffcontrast))
library("superheat")
superheat(heat_pm1,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white")  

# heatmap pm2
lfcpm2 =as.data.frame(summarypm2$test$coefficients)
heat_pm2 = as.data.frame(cbind(lfcpm2[1:95,], lfcpm2[96:190,], lfcpm2[191:285,]))
metabolites = pm2$metabolite
metabolites = setdiff(metabolites, c("Negative.Control"))
rownames(heat_pm2) = metabolites
colnames(heat_pm2) = c("LLcontrasts", "HLcontrasts", "HL-LL contrasts")
pheatmap(heat_pm2, cluster_rows = TRUE, cluster_cols = FALSE)

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

heat_pm2 = as.data.frame(cbind(lfcpm2[1:95,2], lfcpm2[96:190,2], lfcpm2[191:285,2]))
metabolites = pm2$metabolite
metabolites = setdiff(metabolites, c("Negative.Control"))
rownames(heat_pm2) = metabolites
colnames(heat_pm2) = c("LLcontrasts", "HLcontrasts", "HL-LL contrasts")
#heat_pm1 <- na.omit(heat_pm1)
#Delete rows with complete NAs
heat_pm2 <- heat_pm2[rowSums(is.na(heat_pm2)) != ncol(heat_pm2), ]
#heat_pm1<- filter(heat_pm1, rowSums(is.na(heat_pm1)) != ncol(heat_pm1))
#heat_pm1[is.na(heat_pm1)] <- as.double("NA")
pheatmap(heat_pm2, na_col="white", cluster_cols = FALSE, cluster_rows = TRUE)

# llcontrast = as.data.frame(lmpm1$lmlfc[which(as.numeric(lmpm1$indices)<=95)])
# #take out the signifcant columns and 
# hlcontrast = which(as.numeric(lmpm1$indices)>95 & as.numeric(lmpm1$indices)<=190)
# hl_ll_diffcontrast = which(as.numeric(lmpm1$indices)>190)
# 
# significant_pm1 = as.data.frame(c(llcontrast, hlcontrast, hl_ll_diffcontrast))

library("superheat")
superheat(heat_pm2,
          # scale the matrix
          # change color of missing values
          heat.na.col = "white")

