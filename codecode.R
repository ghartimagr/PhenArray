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
<<<<<<< HEAD
library("ggplot2")
library("ggrepel")
=======
>>>>>>> 353ec73864ff4a1752f87a7a205db484b7a91ee8
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
bwplot(newpmcontrol$value~newpmcontrol$position|newpmcontrol$variable, 
       col = as.numeric(as.factor(newpmcontrol$position)), ylab= "Values", xlab= "Position")

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
#boxplot(Tpm1stacked$values~Tpm1stacked$sample , col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")
bwplot(pm1stacked$values~pm1stacked$metabolite|pm1stacked$light_condition, col= rainbow(ncol(Tpm1)), ylab = "values", xlab = "metabolites")

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
#plot(fitpmc)

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
head(spm1stacked)
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
#vec = setdiff(vec, c("Negative.ControlLL-Negative.ControlLL"))
#vec2 = setdiff(vec2, c("Negative.ControlHL-Negative.ControlHL"))

metabolitevec <- paste(vec2,vec, sep ="-")
#metabolitevec <- setdiff(metabolitevec, c("Negative.ControlHL-Negative.ControlHL-Negative.ControlLL-Negative.ControlLL"))

contpm1=makeContrasts(contrasts =c(vec2, vec, metabolitevec), levels=Xpm1)
head(contpm1)
#using default lm function
fitpm1 <- lm(values~0+group, data = spm1)
summary(fitpm1)
#plot(fitpm1)
# glht already adjust pvalues
confit=glht(fitpm1, t(contpm1))
<<<<<<< HEAD
lmcoef = confit$coef
summary(confit) #doesnt work, tried lines 283-285,, still doesnt work
=======
summary(confit)
>>>>>>> 353ec73864ff4a1752f87a7a205db484b7a91ee8

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


#making a dataframe of limmaLFC ie. estimate values and limmapval ie. p values for volcano plot
limmaLFC = t(limmaLFC)
dfpm1 <- as.data.frame(cbind(limmaLFC, limmapval))
dfpm1 <- tibble::rownames_to_column(dfpm1, "contrasts") 
dfpm1 <- tibble::rownames_to_column(dfpm1, "Numbers") 
colnames(dfpm1) = c("Numbers", "contrasts", "limmaLFC", "limmapval")
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

volcanoplot(eb, coef = "contrasts")

p1 <- ggplot(dfpm1, aes(x =limmaLFC, y=-log(limmapval,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"limmaLFC")) + 
  ylab(expression("-log"[10]*"limmapval"))
p1

#Adding color to differentially expressed genes (DEGs)
#Differentially expressed genes (DEGs) are usally considered as those with 
#an absolute fold change greater or equal to 2 and a limmapval value of 0.05 or less. 
#So, we can make our volcano plot a bit more informative if we add some color to the DEGs in the plot. 
#To do so, we’ll add an additional column, named ‘Expression’, 
#indicating whether the expression of a gene is up-regulated, down-regulated, or unchanged

dfpm1 <- dfpm1 %>% 
  mutate(
    Expression = case_when(limmaLFC >= log(2) & limmapval <= 0.05 ~ "Up-regulated",
                           limmaLFC <= -log(2) & limmapval <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(dfpm1) %>%
  knitr_table()

p2 <- ggplot(dfpm1, aes(limmaLFC, -log(limmapval,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"limmaLFC")) + 
  ylab(expression("-log"[10]*"limmapval")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

#If we want to know how many genes are up- or down-regulated, or unchanged,
#we can use dplyr’s count() function.
dfpm1 %>% 
  count(Expression) %>% 
  knitr_table()

#further methods for volcano plots
dfpm1 <- dfpm1 %>%
  mutate(
    Significance = case_when(
      abs(limmaLFC) >= log(2) & limmapval <= 0.05 & limmapval > 0.01 ~ "limmapval 0.05",
      abs(limmaLFC) >= log(2) & limmapval <= 0.01 & limmapval > 0.001 ~ "limmapval 0.01",
      abs(limmaLFC) >= log(2) & limmapval <= 0.001 ~ "limmapval 0.001",
      TRUE ~ "Unchanged")
  )
head(dfpm1) %>%
  knitr_table()

p3 <- ggplot(dfpm1, aes(limmaLFC, -log(limmapval,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"limmapval")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5)))

p3

dfpm1 %>%
  count(Expression, Significance) %>%
  knitr_table()

topup <- 1
topdown <- 163
top_genes <- bind_rows(
  dfpm1 %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(limmapval, desc(abs(limmaLFC))) %>%
    head(topup),
  dfpm1 %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(limmapval, desc(abs(limmaLFC)))%>%
    head(topdown)
)
top_genes %>% 
  knitr_table()
p3 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(limmaLFC, -log(limmapval,10), label = Numbers),
                   size = 2)
p3

# LINEAR MODEL FOR PM2 - 

#this our scaled pm1 table with values , light condition and 
#metabolites stacked in their repective columns and we use it to fit the model now

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
#vec = setdiff(vec, c("Negative.ControlLL-Negative.ControlLL"))
#vec2 = setdiff(vec2, c("Negative.ControlHL-Negative.ControlHL"))

metabolitevec <- paste(vec2,vec, sep ="-")
#metabolitevec <- setdiff(metabolitevec, c("Negative.ControlHL-Negative.ControlHL-Negative.ControlLL-Negative.ControlLL"))

contpm2=makeContrasts(contrasts =c(vec2, vec, metabolitevec), levels=Xpm2)
head(contpm2)
#using default lm function
fitpm2 <- lm(values~0+group, data = spm2)
summary(fitpm2)
# glht already adjust pvalues
confit=glht(fitpm2, t(contpm2))
lmcoef = confit$coef
summary(confit) #doesnt work, tried lines 283-285,, still doesnt work

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
dfpm2 <- tibble::rownames_to_column(dfpm2, "Numbers") 
colnames(dfpm2) = c("Numbers", "contrasts", "limmaLFC", "limmapval")
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

p1 <- ggplot(dfpm2, aes(x =limmaLFC, y=-log(limmapval,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"limmaLFC")) + 
  ylab(expression("-log"[10]*"limmapval"))
p1

#Adding color to differentially expressed genes (DEGs)
#Differentially expressed genes (DEGs) are usally considered as those with 
#an absolute fold change greater or equal to 2 and a limmapval value of 0.05 or less. 
#So, we can make our volcano plot a bit more informative if we add some color to the DEGs in the plot. 
#To do so, we’ll add an additional column, named ‘Expression’, 
#indicating whether the expression of a gene is up-regulated, down-regulated, or unchanged

dfpm2 <- dfpm2 %>% 
  mutate(
    Expression = case_when(limmaLFC >= log(2) & limmapval <= 0.05 ~ "Up-regulated",
                           limmaLFC <= -log(2) & limmapval <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(dfpm2) %>%
  knitr_table()

p2 <- ggplot(dfpm2, aes(limmaLFC, -log(limmapval,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"limmapval")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

#If we want to know how many genes are up- or down-regulated, or unchanged,
#we can use dplyr’s count() function.

dfpm2 %>% 
  count(Expression) %>% 
  knitr_table()

p2 <- ggplot(dfpm1, aes(limmaLFC, -log(limmapval,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"limmaLFC")) + 
  ylab(expression("-log"[10]*"limmapval")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

#If we want to know how many genes are up- or down-regulated, or unchanged,
#we can use dplyr’s count() function.
dfpm2 %>% 
  count(Expression) %>% 
  knitr_table()

#further methods for volcano plots
dfpm2 <- dfpm2 %>%
  mutate(
    Significance = case_when(
      abs(limmaLFC) >= log(2) & limmapval <= 0.05 & limmapval > 0.01 ~ "limmapval 0.05",
      abs(limmaLFC) >= log(2) & limmapval <= 0.01 & limmapval > 0.001 ~ "limmapval 0.01",
      abs(limmaLFC) >= log(2) & limmapval <= 0.001 ~ "limmapval 0.001",
      TRUE ~ "Unchanged")
  )
head(dfpm2) %>%
  knitr_table()

p3 <- ggplot(dfpm2, aes(limmaLFC, -log(limmapval,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"limmapval")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5)))

p3

dfpm2 %>%
  count(Expression, Significance) %>%
  knitr_table()

topup <- 1
topdown <- 163
top_genes <- bind_rows(
  dfpm2 %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(limmapval, desc(abs(limmaLFC))) %>%
    head(topup),
  dfpm2 %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(limmapval, desc(abs(limmaLFC)))%>%
    head(topdown)
)
top_genes %>% 
  knitr_table()
p3 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(limmaLFC, -log(limmapval,10), label = Numbers),
                   size = 2)
p3



