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




