library("limma")
Scaling <- function(table)
{
  #converting dataframe into matrix
  mat = as.matrix(table)
  #initialize aan empty matrix to store the scaling factors
  sfmat = matrix(nrow = nrow(mat), ncol = ncol(mat))
  table2 = table
  for (i in 2:ncol(table))
  {
    for ( j in 1:nrow(table))
    {
      #we divide te each value in the matrix by mean of the column
      sfmat[j,i] = table[j,i]/mean(table[,i])
      sfmat = as.data.frame(sfmat)
      #table 2 actually stores the scaling factors
      table2[j,i] = table[j,i]/sfmat[j,i]
    }
  }
  return(table2)
}##


CTableUpdate <- function(table)
{
  table2 = table
  for (i in 2:ncol(table))
  {
    for ( j in 1:nrow(table))
    {
      table2[j,i] = table[j,i]/sftab[j,i]
    }
  }
  return(table2)
}





Stacked <- function(table)
{
  for ( i in 1:nrow(table))
  {
    #changing LL1, LL2, and LL3 into LL
    if (table$light_condition[i] %in% c("LL1", "LL2", "LL3"))
    {
      table$light_condition[i] = "LL"
    }
    #turning HL1, HL2 and HL3 into HL
    else if (table$light_condition[i] %in% c("HL1", "HL2", "HL3"))
    {
      table$light_condition[i] = "HL"
    }
  }
  return(table)
}




downdiff <- function (down)
{
  listMets = list()
  idxLL<-which(grepl("LL-Negative.ControlLL", down$X1 )==TRUE)
  downLL <- down$X1[idxLL]
  downLLmets <- gsub("LL-Negative.Control.+", "", downLL)
  
  ##metabolites causeing downregulation in high light
  idxHL <- which(grepl("HL-Negative.ControlHL", down$X1 )==TRUE)
  downHL <- down$X1[idxHL]
  downHLmets <- gsub("HL-Negative.Control.+", "", downHL)
  
  listMets$LLmets = downLLmets
  listMets$HLmets = downHLmets

  return(listMets)
  
}


CsourceExtract <- function(Mets)
{
  Csource = vector(length = length(Mets), mode = "character")
  for ( i in 1:length(Mets))
  {
    if (Mets[i] %in% pm1compounds$X3)
    {
      idx = which(pm1compounds$X3==Mets[i])
      Csource[i] = pm1compounds$X4[idx]
      i =i+1
    }
  }
  
  Csource <- gsub("C.Source+\\.\\.", "", Csource)
  df <- as.data.frame(cbind(Mets, Csource))
  return(df)
  
}


CsourceExtractpm2 <- function(Mets)
{
  Csource = vector(length = length(Mets), mode = "character")
  for ( i in 1:length(Mets))
  {
    if (Mets[i] %in% pm2compounds$X3)
    {
      idx = which(pm2compounds$X3==Mets[i])
      Csource[i] = pm2compounds$X4[idx]
      i =i+1
    }
  }
  
  Csource <- gsub("C.Source+\\.\\.", "", Csource)
  df <- as.data.frame(cbind(Mets, Csource))
  return(df)
  
}


Nodiff <- function (No)
{
  list_nodiff <- list()
  #metablotites with LL, HL and HL-LL separation
  # metabolites which casue downregulation in lowlight
  idxNoLL<-which(grepl("LL-Negative.ControlLL$", No$X1 )==TRUE)
  #length(idxNoLL) #55
  NoLL <- No$X1[idxNoLL]
  NoLLmets <- gsub("LL-Negative.Control.+", "", NoLL)
  
  ##metabolites causeing downregulation in high light
  idxNoHL <- which(grepl("HL-Negative.ControlHL$", No$X1 )==TRUE)
  #length(idxNoHL) # 83
  NoHL <- No$X1[idxNoHL]
  NoHLmets <- gsub("HL-Negative.Control.+", "", NoHL)
  
  
  idxno <-  which(grepl("^\\(", No$X1 )==TRUE)
  #length(idxno) #79 # 79+83+55 = 217
  Nomets <<- No$X1[idxno]
  Nomets <- gsub("HL-.+", "", Nomets)
  Nomets <- gsub(".+ ", "", Nomets)
  
  list_nodiff$LLmets <- NoLLmets
  list_nodiff$HLmets <- NoHLmets
  list_nodiff$HL_LLmets <- Nomets
  
  return(list_nodiff)
  
}


phyperfunction <- function(contingency_mat, str = "Upregulation")
{
  probabilitiesUp= list()
  columns <- colnames(contingency_mat)
  q = list()
  
  for (class in row.names(contingency_mat))
  {
    #goes from 0 to total no. of upregulation in that class
    q = c(0:contingency_mat[class, str])
    p = stats::phyper(
      q = q-1,
      m = sum(contingency_mat[class, columns]),
      n = sum(contingency_mat[which(row.names(contingency_mat) != class), columns]),
      k = sum(contingency_mat[, str]), lower.tail = FALSE
    )
    probs = data.frame(str= q,
                       "pvalues" = p)
    #print (probs)
    probabilitiesUp[[class]] = probs
  }
  return(probabilitiesUp)
}




plotfunction <- function(probabilities)
{
  plotlist = list()
  for ( i in 1:8)
  {
    data <- as.data.frame(probabilities[[i]][c(1,2)])
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
               colour="black") + labs(x = "No. of differential Expression", y = "p-values")
    plotlist[[i]] = g
    i = i+1
    
  }
  figure <- ggarrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]], plotlist[[5]], plotlist[[6]], plotlist[[7]], plotlist[[8]]
                      , labels = c("Alcohol", "Amide", "Amine", "Amino acid", "Carbohydrate", "Carboxylic acid", "Ester", "Fatty acid"), ncol = 2, nrow = 4)
  figure
  
  return(figure)
  
}

importKEGG <- function(ids) {                                                                                                                                                                                                                                                                                                 
  sdfset <- SDFset()                                                                                                                                                                                                                                                                                                        
  for(i in ids) {                                                                                                                                                                                                                                                                                                           
    url <- paste0("http://www.genome.jp/dbget-bin/www_bget?-f+m+drug+", i)                                                                                                                                                                                                                                                
    tmp <- as(read.SDFset(url), "SDFset")                                                                                                                                                                                                                                                                                 
    cid(tmp) <- makeUnique(i)                                                                                                                                                                                                                                                                                                         
    sdfset <- c(sdfset, tmp)                                                                                                                                                                                                                                                                                              
  }                                                                                                                                                                                                                                                                                                                         
  sdfset                                                                                                                                                                                                                                                                                                                    
}  


smilesfunc <- function(CasNo)
{
  #NOT function but just a forloop
  # smilespm1 <- list()
  # for ( i in 1: length(pm1CASno))
  # {
  #   ismile <- as.data.frame(cir_query(pm1CASno[i], "smiles"))
  #   smilespm1[[i]] <- ismile
  # }
  
  smiles <- list()
  for ( i in 1:length(CasNo))
  {
    ismile <- as.data.frame(cir_query(CasNo[i], "smiles"))
    smiles[[i]] <- ismile
  }
  return(smiles)
  
}
