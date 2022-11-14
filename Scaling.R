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
}
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



Contrast <- function(coef, LLc, HLc ) #takes the coefficients of the model as input
{
  cf <- coef[, 2:ncol(coef)] #take out te numeric rows
  contrast = matrix(0,nrow = nrow(cf), ncol= ncol(cf))
  for ( i in 1: nrow(cf))
    {
    for (j in 1: ncol(cf))
      {
      if (i%%2==0) #even rows contain LL data
        {
        #substracting the LL values from LL negative control
        #LLcontrol - LLmetabolite
        contrast[i,j] = cf[LLc, j]- cf[i,j]
        }
      else
        {
          #HLcontrol - HLmetabolite
          contrast[i,j] = cf[HLc, j]-cf[i,j]
        }
    }
  }
  #contrast <- cbind(coef[, 1], contrast)
  return(contrast)
}


Difference <- function(contrast)
{
  diffcontrast <- matrix(0, ncol= ncol(contrast), nrow = nrow(contrast)/2)
  for (i in 1: nrow(contrast))
    {
    for (j in 1 :ncol(contrast))
      {
      if (i%%2 ==0) # if the row is even ie- LL difference rows
        #then we enter the loop and subtract the value from HL values from the orevious ro
        #then we update the diffcontrast matrix
        {
        diffcontrast[i/2, j] <- contrast[i,j]-contrast[i-1,j]
        #we divide by2 , cuz for every 2 rows in contrast , we have one row in diffcontrast
      }
    }
  }
  return(diffcontrast)
}

