# simulate_diff_tables --   Differential patterns
# Method to simulate first, second and full order differential tables
# Input : Number of tables and number of conditions
# Output: Sampled tables of K conditions
# Created by: Ruby Sharma and Dr. Joe Song
# Date Created: 7 March 2019
# Modified: 
#   20 December 2019.
#   March 19, 2020. MS. Made changes in function name, input, and return.

# Parameters details
# K              : Number of tables(experimental conditions)[DEFAULT = 2]
# r              : Expected number of rows in the table, should be greater than 1[DEFAULT = 3].
# s              : Expected number of columns in the table , should be greater than 1[DEFAULT = 3].
# n              : Expected sample size [DEFAULT = 100]
# B              : Number of iterations [DEFAULT = 100]
# types include  : first-order: marginals are different across the K tables keeping the interaction same,
#                  second-order: marginals are same across K tables but the interaction is different[DEFAULT] ,
#                  full-order: marginals and interaction are different across K tables. [DEFAULT = second-order]

#' @importFrom stats rexp rmultinom runif


#' @export
simulate_diff_tables <- function(
  K = 2, nrow = 3, ncol = 3, n = 100, B = 100,
  type = c("second-order", "first-order", "full-order")
)
{
  type <- match.arg(type)
  
  if(K<2){
    stop("Number of tables should be greater or equal to 2!")
  }
  if (nrow <  2 || ncol < 2){
    stop("Number of rows and columns should be >=2!")
  }
  if(n < 4){
    stop("Sample size should be greater than 4!")
  }
  
  if(type == "first-order"){
    prob.tables = FirstOrder(K, nrow, ncol, n, B)
    Tabletype = "First-order differential tables"
    
  }else if(type == "second-order"){
    prob.tables = SecondOrder(K, nrow, ncol, n, B)
    Tabletype = "Second-order differential tables"
    
  }else{
    prob.tables = FullOrder(K, nrow, ncol, n, B)
    Tabletype = "Full-order differential tables"
  }
  
  if(0) { # We do not need this generally. Remove from future versions
    # Sampling original table using multinomial distribution
    SamP = lapply(seq(K), function(k){
      return(matrix(rmultinom(1, n, prob.tables$P[[k]]),
                    nrow = nrow, ncol = ncol, byrow= FALSE))
    })
  }
  
  #Sampling modified table using multinomial distribution
  SamDeltaP = lapply(seq(K), function(k){
    return(matrix(rmultinom(1, n, prob.tables$DeltaP[[k]]),
                  nrow = nrow, ncol = ncol, byrow= FALSE))
  })
  
  # P = original probability table, DeltaP = modified probability table, 
  # SampP = sampled original table, SamDeltaP = sampled modified table (differential)
  return(structure(list(
    contingency.tables = SamDeltaP,
    probability.tables = prob.tables$DeltaP, # DeltaP = tables$DeltaP, 
    # P= prob.tables$$P, 
    # SamP = SamP, 
    method = Tabletype
  )))
}

# Method to generate first order differential tables
FirstOrder = function(K, nrow, ncol, n, B, seed.prob.table)
{
  # generating different row and column marinals for K tables
  Marg = GenMarginals(K, nrow, ncol, "D")
  
  # Creating Pk tables using joint distribution 
  P = lapply(seq(K),  function(k){
    return(unlist(Marg$Px[[k]]) %*% t(unlist(Marg$Py[[k]])))
  })
  
  DeltaP = P
  
  iter = 1
  while(iter <= B) {
    # Generate two row and two column indices 
    CellRan = genindex(nrow, ncol)                        
    
    # obtain random 4 cells from all K tables 
    Pij = lapply(seq(K), function(k){
      return(DeltaP[[k]][CellRan$Rin , CellRan$Cin])
    }) 
    
    # Second order adjustment depending on the differential table type 
    # (modifying K DeltaP tables together)
    DeltaP = SecondOrdAdj(K, CellRan, Pij, DeltaP, "first-order")
    
    iter = iter + 1
  }
  
  list( DeltaP = DeltaP, P = P)
  
}

# Method to generate seocond order differential tables
SecondOrder = function(K, nrow, ncol, n, B, seed.prob.table)
{
  # generating same row and column marinals for K tables
  Marg = GenMarginals(K, nrow, ncol, "S")
  
  # Creating Pk tables using joint distribution 
  P = lapply(seq(K),  function(k){
    return(unlist(Marg$Px[[k]]) %*% t(unlist(Marg$Py[[k]])))
  })
  
  DeltaP = P
  
  iter = 1
  # Applying second order adjustment B times
  while(iter <= B)
  { 
    # mofifying K DeltaP tables individually
    lapply(seq(K),  function(k){
      CellRan = genindex(nrow, ncol)
      Pij = DeltaP[[k]][CellRan$Rin , CellRan$Cin]
      DeltaP[[k]] <<- SecondOrdAdj(K, CellRan, Pij, DeltaP[[k]], "second-order")
    })
    iter = iter+1
  }
  
  list(DeltaP = DeltaP, P = P)
  
}

# Method to generate fullorder differential tables 
FullOrder = function(K, nrow, ncol, n, B, seed.prob.table)
{
  # generating different row and column marinals for K tables
  Marg = GenMarginals(K, nrow, ncol, "D")
  
  # Creating Pk tables using joint distribution 
  P = lapply(seq(K),  function(k){
    return(unlist(Marg$Px[[k]]) %*% t(unlist(Marg$Py[[k]])))
  })
  
  DeltaP = P
  
  iter = 1
  
  # Applying second order adjustment B times
  while(iter <= B)
  {
    # mofifying K DeltaP tables individually
    lapply(seq(K),  function(k){
      CellRan = genindex(nrow, ncol)
      Pij = DeltaP[[k]][CellRan$Rin , CellRan$Cin]
      DeltaP[[k]] <<- SecondOrdAdj(K, CellRan, Pij, DeltaP[[k]], "full-order")
    })
    iter = iter+1
  }
  
  list(DeltaP = DeltaP, P = P)
}

# Generate marginals depending on tables type using exponential distrisbution. 
#Firstorder and fullorder tables have different marginals. 
# Secondorder tables have same marginals
GenMarginals = function(K, nrow, ncol, marType)
{
  rowMar <- list ()
  colMar <- list ()
  
  if(marType == "D") {
    lapply(seq(K), function(X)
    {
      Mar <- rexp(n = nrow) 
      Mar <- Mar/sum(Mar)
      rowMar[[X]] <<- Mar
      
      Mar <- rexp(n = ncol) 
      Mar <- Mar/sum(Mar)
      colMar[[X]] <<- Mar
    })
  } else {
    Mar <- rexp(n = nrow) 
    Mar <- Mar/sum(Mar)
    rowMar[[1]] <- Mar
    rowMar <- list(rowMar)[rep(1,K)]
    
    Mar <-  rexp(n = ncol) 
    Mar <- Mar/sum(Mar)
    colMar[[1]] <- Mar
    colMar <- list(colMar)[rep(1,K)]
  }
  list(Px = rowMar, Py = colMar)
  
}

# Generates four random cells 
genindex = function(nrow, ncol)
{
  col.in = sort(sample(seq(ncol), 2, replace = FALSE))
  row.in = sort(sample(seq(nrow), 2, replace = FALSE))
  
  list(Rin = row.in, Cin = col.in)
}


# Obtain indices of non-zero elements 
Zpos <- function(Pij)
{
  inds = which(Pij == 0, arr.ind = T) 
  flag = ifelse(length(which(as.vector(inds)==1))>=3 || length(which(as.vector(inds)==2))>=3, 'NA', 'A')
  if(all(is.na(inds))){inds =  which(Pij != 0, arr.ind = T)}
  
  list(Zin = inds, Fg = flag)
}  

# Generate index, from where addition should start
DeInd <- function(Pij, IndTyp)
{
  addInd = matrix(
    IndTyp$Zin[
      ifelse(IndTyp$Fg == 'A', 
             sample(c(rep(1:nrow(IndTyp$Zin))),1), 'NA' ), ]
  )
  
  return(addInd)
}  


# Add epsilon to modify the four cells 
addEpsi <- function(Pij, ep, addInd)
{
  feasible = TRUE
  if(!is.na(addInd[1])){
    if(paste(as.integer(addInd[,1]), collapse ="") == "11" || paste(as.integer(addInd[,1]), collapse ="") == "22")
    {
      Pij[1,1] = Pij[1,1]+ ep
      Pij[2,2] = Pij[2,2]+ ep
      Pij[2,1] = Pij[2,1]- ep
      Pij[1,2] = Pij[1,2]- ep
    }else{
      Pij[1,1] = Pij[1,1]- ep
      Pij[2,2] = Pij[2,2]- ep
      Pij[2,1] = Pij[2,1]+ ep
      Pij[1,2] = Pij[1,2]+ ep
    }
  }
  
  if(any(unlist(Pij)<0))
  {
    feasible = FALSE
  }
  
  list(Pij = Pij, feasible = feasible)
}  

NZF =  function(SCells, K)
{
  inds <- sapply(seq(K), function(k){
    ind = which(SCells[[k]] == 0, arr.ind = T)
    return(paste(ind[,1], ind[,2], sep=" "))
  })
  
  inds = unique(unlist(inds))
  inds = as.numeric(unlist(strsplit(inds, " ")))
  
  flag = ifelse(length(which(as.vector(inds)==1))>=3 || length(which(as.vector(inds)==2))>=3, 'NA', 'A')
  if(all(is.na(inds))){
    inds =  which(SCells[[1]] != 0, arr.ind = T) 
  }
  
  inds = matrix(as.vector(inds), ncol = 2, byrow =TRUE)
  list(Zin = inds, Fg = flag)
}

# Second order adjustment depending on differential table type: Deciding and adding epsilon to the selected cells
SecondOrdAdj <- function(K, CellRan, Pij, DeltaP, type)
{
  
  nonZero = unlist(Pij)
  
  # getting minimum non zero value from four cells
  min =  min(nonZero[nonZero!=0])
  
  # generate epsilon between 0 and min
  ep = runif(min = 0, max = min, n =1)  
  
  # deciding the starting index for adding epsilon
  if(type == "first-order")
  {
    IndTyp = NZF(Pij, K)
    addInd = DeInd(Pij, IndTyp)
    lapply(seq(K), function(k)
    {
      DeltaPij <- addEpsi(Pij[[k]], ep, addInd)
      if(DeltaPij$feasible == TRUE)
      {
        DeltaP[[k]][CellRan$Rin, CellRan$Cin] <<- DeltaPij$Pij
      }
    })
  } else {
    IndTyp = Zpos(Pij)
    addInd = DeInd(Pij, IndTyp)
    DeltaPij <- addEpsi(Pij, ep, addInd)
    if(DeltaPij$feasible == TRUE)
    {
      DeltaP[CellRan$Rin , CellRan$Cin] = DeltaPij$Pij
    }
  }
  
  return(DeltaP)
}
