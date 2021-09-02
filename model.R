# Importing required libraries....
library(stringr) # for string processing..
library(pROC)
library(ggplot2)

# This program to process data for the recovery group, which has 55 unique V genes and 22 unique J gense. The amino acid sequences for thier CDR3 regions have maximum lenght of 34

# List of standard amino acids
AAlist <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# List of J segment genes
J <- c("TRBJ1-1*01", "TRBJ1-2*01" , "TRBJ1-3*01" , "TRBJ1-4*01",  "TRBJ1-5*01" , "TRBJ1-6*02",  "TRBJ2-1*01" , "TRBJ2-2*01" , "TRBJ2-2P*01" ,"TRBJ2-3*01" , "TRBJ2-4*01" , "TRBJ2-5*01" , "TRBJ2-6*01" , "TRBJ2-7*01" )

# List of V Segment genes.
V <- c("TRBV1*01",  "TRBV10-1*01", "TRBV10-2*01" ,"TRBV10-3*01" ,"TRBV11-1*01", "TRBV11-2*01", "TRBV11-3*01", "TRBV12-1*01" ,"TRBV12-2*01" ,"TRBV12-3*01" ,"TRBV12-4*01", "TRBV12-5*01" ,"TRBV13*01" ,  "TRBV14*01" ,  "TRBV15*01" ,  "TRBV16*02",   "TRBV17*01",   "TRBV18*01",  "TRBV19*01" ,  "TRBV2*01",    "TRBV20-1*01", "TRBV23-1*01", "TRBV24-1*01" ,"TRBV25-1*01", "TRBV27*01",   "TRBV28*01",   "TRBV29-1*01" ,"TRBV3-1*01" , "TRBV3-2*01",  "TRBV30*01" ,  "TRBV4-1*01" , "TRBV4-2*01",  "TRBV4-3*01" , "TRBV5-1*01" , "TRBV5-3*01" , "TRBV5-4*01" ,"TRBV5-5*01",  "TRBV5-6*01" , "TRBV5-7*01",  "TRBV5-8*01",  "TRBV6-1*01",  "TRBV6-2*01",  "TRBV6-4*01" , "TRBV6-5*01",  "TRBV6-6*01",  "TRBV6-7*01",  "TRBV6-8*01",  "TRBV6-9*01",  "TRBV7-1*01",  "TRBV7-2*01",  "TRBV7-3*01" , "TRBV7-4*01",  "TRBV7-6*01" , "TRBV7-7*01" , "TRBV7-8*01" ,"TRBV7-9*01" ,"TRBV9*01")
#******************************************************************************
# A function to construct a parameter vector
parm <- function()
{
  v <- lapply(1:20, function(x) lapply(AAlist,function(y) paste("q(",y ,x,")", sep = "")))
  v <- unlist(v)
  v <- c("Frq", paste("q(", V, ")", sep = ""),    paste("q(", J, ")", sep = ""), v)

  p <-  rep(0, length(v))#rnorm(length(v))
  names(p) <- v
  return(as.data.frame(p))
}
#*****************************************************************************
 # A function to process an input vector for give tripelt (CDR3, v,j)
getInputVector <- function(cdr, v, j, frq)
{
  Q <- parm()
  p <- c(paste("q(", v, ")", sep = ""), paste("q(", j, ")", sep = ""))
  cdr <- unlist(strsplit(toupper(cdr), split = ""))
  for (k in 1:length(cdr))
  {
    p <- c(p, paste("q(", cdr[k], k,")", sep = ""))
  }

  for ( r in rownames(Q))
  {
    Q[r, ] <- length(which(p == r))
  }
  Q['Frq', 1] <- frq
  return(Q)
}
#**********************************************************************
normlize.function <-function(z)
{
  c <- sapply(z, function(x) (x-min(z))/(max(z)-min(z)))
  return(c)
}
#****************************************************************************************

# A function to process data of the control  group
getData <- function(g)
{
  d  <- readRDS("./Data/FullData.Rda")
  d1 <- subset(d, d$Group == g)
  d1 <- d1[!duplicated(d1[,c(6,9, 11)]),]
  y <-  ifelse(d1$Productive.status =='true', 1, 0)
  yp <- y[y==1]
  yn <- y[y==0]
  dp <- d1[d1$Productive.status =='true',]
  dn <- d1[d1$Productive.status =='false', ]
  
  dp$Count <- normalize(dp$Frequency....) # by count
  dn$Count <- normalize(dn$Frequency....)
  
  pmat <- matrix(0, nrow = dim(dp)[1], ncol = 472)
  colnames(pmat) <- rownames(parm())
  
  nmat <- matrix(0, nrow = dim(dn)[1], ncol = 472)
  colnames(nmat) <- rownames(parm())

  for ( i in 1:dim(dp)[1])
  {
    pmat[i, ] <- t(getInputVector(dp[i, 6], d1[i, 9], d1[i, 11], d1[i, 4]))
  }
  
  for ( i in 1:dim(dn)[1])
  {
    nmat[i, ] <- t(getInputVector(dn[i, 6], d1[i, 9], d1[i, 11], d1[i, 4]))
  }
  
  mat <- rbind(pmat, nmat)
  y <- c(yp, yn)
  
  saveRDS(y, paste("./",g, "/", g, "_labels.Rda", sep = ""))
  saveRDS(mat, paste("./",g, "/", g, "_data.Rda", sep = ""))

  print("The data has been processed")
  #return(mat)
}
#*************************************************************
# Function to count the number of missclassified points and the cost of objcetive function
computeCost <- function(mats, y , parm)
{
  Q <- mats %*% parm[, 1]
  h <- 1/(1 + exp(Q))
  cst <- -sum(y*log(h) + (1-y)*log(1-h))/length(y)
  prob <- ifelse( h < 0.5, 0, 1)

  differ <- sum(abs(prob - y))
  return(c(cst, differ))
}

#**********************************************************************************
# Function to measure the performance of the model
computePerformance <- function(mats, y, parm)
{
  Q <-  mats %*% parm[, 1]
  probs <- 1/(1 + exp(Q))
  prob <- ifelse( probs < 0.5, 0, 1)
  results <- data.frame(y, prob, probs, abs(y - prob))
  colnames(results) <- c("Observed", "Predicted", "value", "Diff" )
  #area <- auc(results$"Predicted", results$"EMP")
  print(sum(abs(y - prob)))
  return(results)
}

###############################################################################
# A function to compute the gradient
gradient <- function(j,  mats, y, parm, lam, a)
{
  js <- which(rownames(parm)==j)
  xi <-  mats %*% parm[, 1]
  pi <- 1/(1 + exp(xi))
  wi <- pi*(1 - pi)
  zi <- xi + (pi - y)/wi
  aj <- sum(wi*(mats[,j]^2))/length(y)
  cj <- sum(wi*mats[,j]*(zi - (mats[,-js] %*%parm[-js, 1])))/length(y)

  if (cj < -lam*a)
  {
    thj <- (cj + lam*a)/(aj + lam*(1-a))
  }else if (cj > lam*a)
  {
    thj <- (cj - lam*a)/(aj +lam*(1- a))
  }else
  {
    thj <- 0
  }
  return(thj)
}

#**************************************************************************************
lassoModel <- function(mats, y, n, lam, a)
{

  parm <- parm()
  print(c('Iteration', 'Cost', 'Miss Classification', 'sum of parameters value', 'Value of lambda'))
  print(c(0, computeCost(mats, y, parm), sum(abs(parm)), lam))

  for (k in 1:n)
  {
    for (j in rownames(parm))
    {

      parm[j, 1] <- gradient(j,  mats, y, parm, lam, a)
    }
    if (k %% 2 == 0) {lam <- ifelse( lam > 1/length(y),  lam*0.5, lam) }

    print(c(k, computeCost(mats, y, parm), sum(abs(parm)), lam))
  }

  return(parm)
}
#*****************************************************************************************
# m fold cross-Validation function
model.crossvalidation <- function(g, m,  n, lam, a = 0.5)
{
  set.seed(10)
  ds <- readRDS(paste("./", g, "/", g, "_data.Rda",  sep = ""))
  y <- readRDS(paste("./", g, "/", g, "_labels.Rda",  sep = ""))
  indexs <- sample(1:length(y), length(y), replace = F)
  folds <- split(indexs, seq(1:m) )
  
  dfm <- NULL
  for(k in 1:m)
  {
    z <- unlist(folds[k])
    Q <- lassoModel(ds[-z, ], y[-z], n, lam, a)
    #saveRDS(Q, paste("./",g, "/",g, "Q", k,".Rda", sep = ""))
    dfm <- rbind.data.frame(dfm, computePerformance(ds[z,], y[z], Q))
    print("***********************************************************")
  }
  print("All is done ...")
  #saveRDS(dfm, paste("./",g, "/", g, "_five_fold_Results.Rda", sep = ""))
  return(dfm)
}
#********************************************************************
#**************************************************************************************

myroc1 <- function(ypp, y)
{
  tpr <- c(); fpr <- c()
  lam <- seq(0.0, 1.0, 0.1)
  
  for (j in 1:length(lam))
  {
    yp <- ifelse(ypp < lam[j], 0, 1)
    tp <- 0; tn <- 0; fp <- 0; fn <- 0
    for (i in 1:length(y))
    {
      if ((y[i] == 0 ) & (yp[i] == 0))
      {
        tn <- tn + 1
      }else if ((y[i] == 1 ) & (yp[i] == 1))
      {
        tp <- tp + 1
      }else if ((y[i] == 0 ) & (yp[i] == 1))
      {
        fp <- fp + 1
      }else if ((y[i] == 1 ) & (yp[i] == 0))
      {
        fn <- fn + 1
      }
    }
    #show(c(lam[j], tp, tn, fp, fn))
    tpr <- c(tpr, tp/(length(y[y == 1])))
    fpr <- c(fpr, fp/(length(y[y == 0])))
  }
  
  return(data.frame(fpr, tpr))
}
############################### A function to generate ROC graph
aucgraph1 <- function(dfm)
{
  df <- myroc1(dfm$value, dfm$Observed)
  area <- auc(dfm$Observed, dfm$value)
  g <- ggplot(df, aes(x = fpr, y = tpr)) + geom_line( color = "red", size = 0.3) + geom_point(size=0.2)
  g <- g + scale_x_continuous(name = "False Positive Rate (1 - specificity)", breaks = seq(0,1, 0.1)) + scale_y_continuous(breaks = seq(0,1, 0.1), "True Positive Rate (Sensitivity)")
  g <- g + ggtitle(paste("ROC Curve for the control group (AUC =", round(area[1], 3), ")", sep = " "))
  g <- g +  theme_bw()
  g <- g + theme(plot.title = element_text(size = 10,  vjust = 0.5, lineheight = 0.2)) + coord_equal()
  g <- g +  theme(plot.background = element_rect(fill = "white"))
  g <- g + geom_line(linetype = 2, aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)))
  g <- g + geom_line( aes(x = seq(0, 1, 0.1), y = 0)) + geom_line( aes(y = seq(0, 1, 0.1), x = 0))
  g <- g + geom_line( aes(x = seq(0, 1, 0.1), y = 1)) + geom_line( aes(y = seq(0, 1, 0.1), x = 1))
  g <- g + theme( axis.line = element_line(colour = "darkblue", size = 0.51, linetype = "solid"))
  
  return(g)
}
#***********************************************************************
combined.crossvalidation <- function( m,  n, lam, a = 0.5)
{
  set.seed(2)
  
  dc <- readRDS('./Controls/Full_Controls_data.Rda')
  yc <- readRDS("./Controls/Full_Controls_labels.Rda")
  
  da <- readRDS("./Asymptomatic/Asymptomatic_data.Rda")
  ya  <- readRDS('./Asymptomatic/Asymptomatic_labels.Rda')
  
  dr = readRDS('./Recovery/Recovery_data.Rda')
  yr = readRDS('./Recovery/Recovery_labels.Rda')
  
  ds = readRDS('./Severe/Severe_data.Rda')
  ys = readRDS('./Severe/Severe_labels.Rda')
  
  du = readRDS('./Uncomplicated/Uncomplicated_data.Rda')
  yu = readRDS('./Uncomplicated/Uncomplicated_labels.Rda')
  
  ds = rbind(dc, dr, ds, du, da)
  y = c(yc, yr, ys, yu, ya)

  print("ALL IS COMBINED")
  
  indexs <- sample(1:length(y), length(y), replace = F)
  folds <- split(indexs, seq(1:m) )
  
  dfm <- NULL
  for(k in 1:m)
  {
    z <- unlist(folds[k])
    Q <- lassoModel(ds[-z, ], y[-z], n, lam, a)
    saveRDS(Q, paste("./Combined/combined_Q", k,".Rda", sep = ""))
    dfm <- rbind.data.frame(dfm, computePerformance(ds[z,], y[z], Q))
    print("***********************************************************")
  }
  print("All is done ...")
  saveRDS(dfm, paste("./Combined/combined_full_Results.Rda", sep = ""))
  #return(dfm)
  print("ALL IS DONE")
}