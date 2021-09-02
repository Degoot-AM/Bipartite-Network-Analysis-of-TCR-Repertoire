# Importing required libraries....
library(stringr) # for string processing..
library(pROC)
library(ggplot2)
library(ggpubr)
library(LaplacesDemon)
library(reshape2)
library(latex2exp)
library(corrplot)
library(RColorBrewer)
# This program to process data for the recovery group, which has 55 unique V genes and 22 unique J gense. The amino acid sequences for thier CDR3 regions have maximum lenght of 34

# List of standard amino acids
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

# maximum length of CDR3 of control group = 39
# maximum length of CDR3 of recovery group = 34
# maximum length of CDR3 of severe group = 41
# maximum length of CDR3 of uncomplicated group = 37
# maximum length of CDR3 of asymptomatic group = 38

heatMapAnalysis <- function(q, m)
{
  mt <- matrix(0, ncol = 20, nrow= m)
  colnames(mt) <- AAlist
  
  for (x in AAlist)
  {
    for (i in 1:m)
    {
       mt[i, x] <- q[paste("q(",x ,i,")", sep = ""), 1] 
       
    }
  }
  mt <- melt(mt)
  g <- ggplot(data = mt, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile(color = "black")+ labs(x = "Amino acids", y  = "Residue position") + 
    scale_fill_gradient2(low = "darkred", high = "blue", mid = "white", 
      midpoint = 0, limit = c(min(mt$value), max(mt$value)), name="Values") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 5, vjust = 0.5, size = 8, color = 'black', hjust = 0.5)) +  theme(axis.text.y = element_text(angle = 5, vjust = 0.5, size = 8,color = 'black', hjust = 0.5), axis.title.x = element_text( face="bold", colour="#990000", size=20), axis.title.y = element_text( face="bold", colour="#990000", size=20))
   
  return(g)
}

#function to generate the heat map
generate.heatmap <- function()
{
  qc = readRDS('./Controls/Q_controls_all.Rda')
  qr = readRDS('./Recovery/Q_recovery_all.Rda')
  qu = as.data.frame(readRDS('./Uncomplicated/Uncomplicated_Q_all.Rda')[-1,])
  row.names(qu) = rownames(qc)
  qa = as.data.frame(readRDS('./Asymptomatic/Asymptomatic_Q_all.Rda')[-1, ])
  rownames(qa) = rownames(qc)
  qs = readRDS('./Severe/Q_severe_all.Rda')
  q = (qc + qs + qr + qu + qa)/5
  #q <- readRDS("./Combined/Q_0.3.Rda")
  g = heatMapAnalysis(q,20)
  return(g)
}
#*****************************************************************
# A function to get the  V  genes within a group:"Asymptomatic", "Recovery", "Severe", "Controls", and Uncomplicated".
vParm <- function(grp)
{
  d <- readRDS("./FullData.Rda")
  d1 <- subset(d, d$Group == grp)
  V <- unique(d1$V)
  v <- c(paste("q(", V, ")", sep = ""))
  return(v)
}
#*****************************************************************
# A function to gether all the common V genes among the groups
gether_v <- function()
{
  pcv <-  vParm("Controls")
  puv <- vParm("Uncomplicated")
  prv <- vParm("Recovery")
  psv <- vParm("Severe")
  pav <- vParm("Asymptomatic")
  
  vgenes <-  Reduce(intersect, list(prv, pcv, puv, psv, pav))
  return(vgenes)
}
#*****************************************************************
# A function to plot the selection factors of V gense 
vfactor <- function()
{
  z <- gether_v()
  qc <- readRDS("./Controls/Q_controls_all.Rda")
  qr <- readRDS("./Recovery/Q_recovery_all.Rda")
  qu <- readRDS("./Uncomplicated/Uncomplicated_Q_all.Rda")
  qs <- readRDS("./Severe/Q_severe_all.Rda")
  qa <- readRDS("./Asymptomatic/Asymptomatic_Q_all.Rda")
  df <- data.frame(matrix(0, ncol = 6, nrow = length(z)))
  colnames(df) <- c ("Genes", "Controls", "Recovery", "Asymptomatic", "Severe", "Uncomplicated")
  df$Genes <- z
  df$Controls <- qc[z,]
  df$Recovery <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(TRB", "", z, fixed = T)
  zz <- gsub(")", "", zz, fixed = T)
  df$Genes <- zz
  
  gc <- ggplot(data = df, aes(x = df$Genes, y = df$Controls)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "V Genes", y = "Contribution of selection factor" )
  
  gs <- ggplot(data = df, aes(x = df$Genes, y = df$Severe)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "V Genes", y = "Contribution of selection factor" )
  
  ga <- ggplot(data = df, aes(x = df$Genes, y = df$Asymptomatic)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "V Genes", y = "Contribution of selection factor" )
  
  gr <- ggplot(data = df, aes(x = df$Genes, y = df$Recovery)) + geom_point(size =1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "V Genes", y = "Contribution of selection factor" )
  
  gu <- ggplot(data = df, aes(x = df$Genes, y = df$Uncomplicated)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "V Genes", y = "Contribution of selection factor" )
  
  g <- ggarrange(gc+ rremove("x.title") + rremove( "x.text"),  gr + rremove("x.title") + rremove( "x.text"), ga, gu, gs, labels =  c ("Control", "Recovery", "Asymptomatic",  "Uncomplicated", "Severe"),
    ncol = 3, nrow = 2, hjust = -1.5)
  return(df)
}
#**********************************************************************************
vfactor1 <- function()
{
  z <- gether_v()
  qc <- readRDS("./Controls/Q_controls_all.Rda")
  qr <- readRDS("./Recovery/Q_recovery_all.Rda")
  qu <- readRDS("./Uncomplicated/Uncomplicated_Q_all.Rda")
  qs <- readRDS("./Severe/Q_severe_all.Rda")
  qa <- readRDS("./Asymptomatic/Asymptomatic_Q_all.Rda")
  
  df <- data.frame(matrix(0, ncol = 6, nrow = length(z)))
  colnames(df) <- c ("Genes", "Controls", "Recovery", "Asymptomatic", "Severe", "Uncomplicated")
  df$Genes <- z
  df$Controls <- qc[z,]
  df$Recovery <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(TRB", "", z, fixed = T)
  zz <- gsub(")", "", zz, fixed = T)
  df$Genes <- zz
  dfm<-  melt(df)
  colnames(dfm) <- c("Genes", "Groups", "value" )
  gp <-  ggplot(dfm, aes(x= Genes, y = value, colour = Groups))+ geom_point( size= 1.5)
  gp <-  gp + theme_classic()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10), axis.text.x  = element_text(angle=45, vjust=0.3, size=8)) + labs(x= "V Genes", y = "Contribution of selection factor" )
  return(gp)
}
#**********************************************************************************
vfactor2 <- function()
{
  z <- gether_v()
  qc <- readRDS("./Controls/Q_controls_all.Rda")
  qr <- readRDS("./Recovery/Q_recovery_all.Rda")
  qu <- readRDS("./Uncomplicated/Uncomplicated_Q_all.Rda")
  qs <- readRDS("./Severe/Q_severe_all.Rda")
  qa <- readRDS("./Asymptomatic/Asymptomatic_Q_all.Rda")
  
  df <- data.frame(matrix(0, ncol = 6, nrow = length(z)))
  colnames(df) <- c ("Genes", "Healthy", "Recovered", "Asymptomatic", "Severe", "Uncomplicated")
  df$Genes <- z
  df$Healthy <- qc[z,]
  df$Recovered <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(TRB", "", z, fixed = T)
  zz <- gsub(")", "", zz, fixed = T)
  df$Genes <- zz
  dfm<-  melt(df)
  colnames(dfm) <- c("Genes", "Groups", "value" )
  gp <- ggplot(data=dfm, aes(x=Genes, y=value, fill=Groups)) +  geom_bar(stat="identity",width = 0.9, position=position_dodge(), colour="black")
  gp <-  gp + theme_classic()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10), axis.text.x  = element_text(angle= 90, vjust=0.3, size=6, colour = 'black')) + labs(x= TeX("$V\\beta$ "), y = "Value" ) +  theme(legend.position="top", legend.title = element_blank())
  return(gp)
}
#*****************************************************************
# A function to get the  J  genes within a group:"Asymptomatic", "Recovery", "Severe", "Controls", and Uncomplicated".
jParm <- function(grp)
{
  d <- readRDS("./FullData.Rda")
  d1 <- subset(d, d$Group == grp)
  J <- unique(d1$J)
  j <- c( paste("q(", J, ")", sep = ""))
  return(j)
}
#*****************************************************************
# A function to gether all the common J genes among the groups
gether_j <- function()
{
  pcj <-  jParm("Controls")
  puj <- jParm("Uncomplicated")
  prj <- jParm("Recovery")
  psj <- jParm("Severe")
  paj <- jParm("Asymptomatic")
  
  jgenes <-  Reduce(intersect, list(prj, pcj, puj, psj, paj))
  return(jgenes)
}
#*****************************************************************
#*****************************************************************
# A function to plot the selection factors of J gense 
jfactor <- function()
{
  z <- gether_j()
  qc <- readRDS("../Control/Q_control.rda")
  qr <- readRDS("../Recovery/Q_recovery.Rda")
  qu <- readRDS("../Uncomplicated/Q_uncomplicated.rda")
  qs <- readRDS("../Severe/Q_severe.rda")
  qa <- readRDS("../Asymptomatic/Q_Asymptomatic.rda")
  df <- data.frame(matrix(0, ncol = 6, nrow = length(z)))
  colnames(df) <- c ("Genes", "Controls", "Recovery", "Asymptomatic", "Severe", "Uncomplicated")
  df$Genes <- z
  df$Controls <- qc[z,]
  df$Recovery <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(", "", z, fixed = T)
  zz <- gsub(")", "", zz, fixed = T)
  df$Genes <- zz
  
  gc <- ggplot(data = df, aes(x = df$Genes, y = df$Controls)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "J Genes", y = "Contribution of selection factor" )
  
  gs <- ggplot(data = df, aes(x = df$Genes, y = df$Severe)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "J Genes", y = "Contribution of selection factor" )
  
  ga <- ggplot(data = df, aes(x = df$Genes, y = df$Asymptomatic)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "J Genes", y = "Contribution of selection factor" )
  
  gr <- ggplot(data = df, aes(x = df$Genes, y = df$Recovery)) + geom_point(size =1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "J Genes", y = "Contribution of selection factor" )
  
  gu <- ggplot(data = df, aes(x = df$Genes, y = df$Uncomplicated)) + geom_point(size = 1.2, color = "blue") + theme_minimal()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10),
    axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "J Genes", y = "Contribution of selection factor" )
  
  g <- ggarrange(gc+ rremove("x.title") + rremove( "x.text"),  gr + rremove("x.title") + rremove( "x.text"), ga, gu, gs, labels =  c ("Control", "Recovery", "Asymptomatic",  "Uncomplicated", "Severe"),
    ncol = 3, nrow = 2, hjust = -1.5)
  return(df)
}
#****************************************************************************************
jfactor1 <- function()
{
  z <- gether_j()
  qc <- readRDS("./Controls/Q_controls_all.Rda")
  qr <- readRDS("./Recovery/Q_recovery_all.Rda")
  qu <- readRDS("./Uncomplicated/Uncomplicated_Q_all.Rda")
  qs <- readRDS("./Severe/Q_severe_all.Rda")
  qa <- readRDS("./Asymptomatic/Asymptomatic_Q_all.Rda")
  
  df <- data.frame(matrix(0, ncol = 6, nrow = length(z)))
  colnames(df) <- c ("Genes", "Controls", "Recovery", "Asymptomatic", "Severe", "Uncomplicated")
  df$Genes <- z
  df$Controls <- qc[z,]
  df$Recovery <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(", "", z, fixed = T)
  zz <- gsub(")", "", zz, fixed = T)
  df$Genes <- zz
  dfm<-  melt(df)
  colnames(dfm) <- c("Genes", "Groups", "value" )
  gp <-  ggplot(dfm, aes(x= Genes, y = value, colour = Groups))+ geom_point( size= 1.5)
  gp <-  gp  + theme_light() + theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10), axis.text.x  = element_text(angle=90, vjust=0.2, size=8)) + labs(x= "V Genes", y = "Contribution of selection factor" )
  return(gp)
}
#****************************************************************************************
jfactor2 <- function()
{
  z <- gether_j()
  qc <- readRDS("./Controls/Q_controls_all.Rda")
  qr <- readRDS("./Recovery/Q_recovery_all.Rda")
  qu <- readRDS("./Uncomplicated/Uncomplicated_Q_all.Rda")
  qs <- readRDS("./Severe/Q_severe_all.Rda")
  qa <- readRDS("./Asymptomatic/Asymptomatic_Q_all.Rda")
  df <- data.frame(matrix(0, ncol = 6, nrow = length(z)))
  colnames(df) <- c ("Genes", "Healthy", "Recovered", "Asymptomatic", "Severe", "Uncomplicated")
  df$Genes <- z
  df$Healthy <- qc[z,]
  df$Recovered <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(TRBJ", "", z, fixed = T)
  zz <- gsub(")", "", zz, fixed = T)
  df$Genes <- zz
  dfm<-  melt(df)
  colnames(dfm) <- c("Genes", "Groups", "value" )
  gp <- ggplot(data=dfm, aes(x=Genes, y=value, fill=Groups)) +  geom_bar(stat="identity",width = 0.9, position=position_dodge(), colour="black")
  gp <-  gp + theme_classic()+ theme(axis.title.x = element_text( face="bold", colour="#990000", size=20),  axis.title.y  = element_text( face="bold", colour="#990000", size=10), axis.text.x  = element_text(angle=90, vjust=0.3, size=8, colour = 'black')) + labs(x= TeX("$J\\beta$ "), y = "Value" ) +  theme(legend.position="top", legend.title = element_blank())
  return(gp)
}
#**************************************************************************************
convert.factors <- function(a)
{
  if (a >0) {return("+")}else if (a <0){return('-')} else {return('0')}
}
#**********************************************************************
## mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
#****************************************************************************************
jfactor3 <- function()
{
  z <- gether_j()
  qc <- readRDS("./Controls/Q_controls_all.Rda")
  qr <- readRDS("./Recovery/Q_recovery_all.Rda")
  qu <- readRDS("./Uncomplicated/Uncomplicated_Q_all.Rda")
  qs <- readRDS("./Severe/Q_severe_all.Rda")
  qa <- readRDS("./Asymptomatic/Asymptomatic_Q_all.Rda")
  df <- data.frame(matrix(0, ncol = 5, nrow = length(z)))
  colnames(df) <- c ( "Healthy", "Recovered", "Asymptomatic", "Severe", "Uncomplicated")
  #df$Genes <- z
  df$Healthy <- qc[z,]
  df$Recovered <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(TRB", "", z, fixed = T)
  zz <- gsub("*01)", "", zz, fixed = T)
  zz <- gsub("*02)", "", zz, fixed = T)
  row.names(df) <- zz
  dfm <- df
  for (i in 1:dim(df)[1])
  {
    for (k in 1:dim(df)[2])
    {
     dfm[i, k] <- convert.factors(df[i,k])
    
    }
  }
  
  return(list(df, dfm))
}
#******************************************************************************************
# corrplot(cor(cor(df[[1]])), method = 'number', type="upper", order="hclust", col= brewer.pal(n=8, name="PuOr"), tl.col="black", tl.srt=45, p.mat = cor.mtest(df[[1]]), sig.level = 0.01, insig = "n", diag = F)
vfactor3 <- function()
{
  z <- gether_v()
  qc <- readRDS("./Controls/Q_controls_all.Rda")
  qr <- readRDS("./Recovery/Q_recovery_all.Rda")
  qu <- readRDS("./Uncomplicated/Uncomplicated_Q_all.Rda")
  qs <- readRDS("./Severe/Q_severe_all.Rda")
  qa <- readRDS("./Asymptomatic/Asymptomatic_Q_all.Rda")
  
  df <- data.frame(matrix(0, ncol = 5, nrow = length(z)))
  colnames(df) <- c ( "Healthy", "Recovered", "Asymptomatic", "Severe", "Uncomplicated")
  df$Healthy <- qc[z,]
  df$Recovered <- qr[z,]
  df$Asymptomatic <- qa[z,]
  df$Severe <- qs[z,]
  df$Uncomplicated <- qu[z,]
  
  zz <- gsub("q(TRB", "", z, fixed = T)
  zz <- gsub("*01)", "", zz, fixed = T)
  zz <- gsub("*02)", "", zz, fixed = T)
  row.names(df) <- zz
  dfm <- df
  for (i in 1:dim(df)[1])
  {
    for (k in 1:dim(df)[2])
    {
      dfm[i, k] <- convert.factors(df[i,k])
      
    }
  }
  
  return(list(df, dfm))
}