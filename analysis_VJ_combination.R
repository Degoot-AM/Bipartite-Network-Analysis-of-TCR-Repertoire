# Importing required libraries....
library(stringr) # for string processing..
library(pROC)
library(ggplot2)
library(reshape2)
library(dplyr)
library(igraph)
library(ggpubr)
# This program to process data for the recovery group, which has 55 unique V genes and 22 unique J gense. The amino acid sequences for thier CDR3 regions have maximum lenght of 34

# List of standard amino acids
AAlist <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# List of J segment genes
J <- c("TRBJ1-1*01", "TRBJ1-2*01" , "TRBJ1-3*01" , "TRBJ1-4*01",  "TRBJ1-5*01" , "TRBJ1-6*02",  "TRBJ2-1*01" , "TRBJ2-2*01" , "TRBJ2-2P*01" ,"TRBJ2-3*01" , "TRBJ2-4*01" , "TRBJ2-5*01" , "TRBJ2-6*01" , "TRBJ2-7*01" )
# short name of J genes 
Jx <- c("J1.1", "J1.2", "J1.3", "J1.4", "J1.5", "J1.6", "J2.1", "J2.2", "J2.2P", "J2.3", "J2.4","J2.5", "J2.6", "J2.7")


# List of V Segment genes.
V <- c("TRBV1*01",  "TRBV10-1*01", "TRBV10-2*01" ,"TRBV10-3*01" ,"TRBV11-1*01", "TRBV11-2*01", "TRBV11-3*01", "TRBV12-1*01" ,"TRBV12-2*01" ,"TRBV12-3*01" ,"TRBV12-4*01", "TRBV12-5*01" ,"TRBV13*01" ,  "TRBV14*01" ,  "TRBV15*01" ,  "TRBV16*02",   "TRBV17*01",   "TRBV18*01",  "TRBV19*01" ,  "TRBV2*01",    "TRBV20-1*01", "TRBV23-1*01", "TRBV24-1*01" ,"TRBV25-1*01", "TRBV27*01",   "TRBV28*01",   "TRBV29-1*01" ,"TRBV3-1*01" , "TRBV3-2*01",  "TRBV30*01" ,  "TRBV4-1*01" , "TRBV4-2*01",  "TRBV4-3*01" , "TRBV5-1*01" , "TRBV5-3*01" , "TRBV5-4*01" ,"TRBV5-5*01",  "TRBV5-6*01" , "TRBV5-7*01",  "TRBV5-8*01",  "TRBV6-1*01",  "TRBV6-2*01",  "TRBV6-4*01" , "TRBV6-5*01",  "TRBV6-6*01",  "TRBV6-7*01",  "TRBV6-8*01",  "TRBV6-9*01",  "TRBV7-1*01",  "TRBV7-2*01",  "TRBV7-3*01" , "TRBV7-4*01",  "TRBV7-6*01" , "TRBV7-7*01" , "TRBV7-8*01" ,"TRBV7-9*01" ,"TRBV9*01")
# Short name of V genes 
 Vx <- c("V1.1", "V10.1", "V10.2", "V10.3", "V11.1", "V11.2", "V11.3", "V12.1", "V12.2", "V12.3", "V12.4", "V12.5", "V13.1",
      "V14.1", "V15.1", "V16.1", "V17.1", "V18.1", "V19.1", "V2.1", "V20.1", "V23.1", "V24.1", "V25.1", "V27.1", "V28.1",
     "V29.1", "V3.1", "V3.2", "V30.1", "V4.1", "V4.2", "V4.3", "V5.1", "V5.3", "V5.4", "V5.5", "V5.6", "V5.7",
     "V5.8", "V6.1", "V6.2", "V6.4", "V6.5", "V6.6", "V6.7", "V6.8", "V6.9", "V7.1", "V7.2", "V7.3", "V7.4",
     "V7.6" ,"V7.7", "V7.8", "V7.9", "V9.1")
#******************************************************************************
#* function to check for proline contents
proline_check <- function()
{
  d  <- readRDS("./Data/FullData.Rda")
  d <- d[!duplicated(d[,c(6,9,11)]),]
  pd = d %>%filter(Productive.status=='true') %>% select(CDR3.amino.acid.sequence)
  pd = pd[, 1]
  nd = d %>%filter(Productive.status=='false') %>% select(CDR3.amino.acid.sequence)
  nd = nd[,1]
  print((typeof(pd)))
  #print(pd[1:4 ])
  #print(nd[1:4 ])
  print(c(length(pd), length(nd)))
  pcount = 0 
  ncount = 0 
  for (x in pd)
  {
    #print(nchar(x))
   #pcount = ifelse(str_detect(x, "P"), pcount +1, pcount)
    if (startsWith(x, "CP"))
    {
      pcount = pcount +1
      #print(c("P",pcount, x))
    }
  }
  
  for (x in nd)
  { 
    #print(nchar(x))
    #ncount = ifelse(str_detect(x, "CP"), ncount +1, ncount)
    if (startsWith(x, "CP"))
    {
      ncount = ncount +1
     # print(c("N", ncount, x))
    }
    
  }
  return(c(pcount/length(pd), ncount/length(nd), (pcount+ncount)/(length(pd)+length(nd))))
}
#********************************************************
#function to find average length of cdr3 with V and J segments
VJ_combination <- function(g)
{
   
    d  <- readRDS("./Data/FullData.Rda")
    d <- subset(d, d$Group == g)
    d <- d[!duplicated(d[,c(6,9,11)]),]
    
    pmat <- matrix(0, ncol = length(J), nrow = length(V))
    colnames(pmat) <- J
    rownames(pmat) <- V
    pd <- subset(d, d$Productive.status=='true')
    for (i in 1:dim(pd)[1])
    { 
        v = pd$V[i] 
        j = pd$J[i]
        cdr = pd$CDR3.amino.acid.sequence[i]
        frq <- pd$Frequency....[i]
        pmat[v, j] <- ifelse(pmat[v, j]==0, nchar(cdr), (pmat[v, j]+nchar(cdr))/1)#*****
        #pmat[v, j] <- ifelse(pmat[v, j]==0, frq, mean(pmat[v, j],frq))

    }
    
    nmat <- matrix(0, ncol = length(J), nrow = length(V))
    colnames(nmat) <- J
    rownames(nmat) <- V
    nd <- subset(d, d$Productive.status=='false')
    
    for (i in 1:dim(nd)[1])
    { 
        v = nd$V[i] 
        j = nd$J[i]
        cdr = nd$CDR3.amino.acid.sequence[i]
        frq <- nd$Frequency....[i]
        nmat[v, j] <- ifelse(nmat[v, j]==0, nchar(cdr), (nmat[v, j]+nchar(cdr))/1)#*****
        #nmat[v, j] <- ifelse(nmat[v, j]==0, frq, mean(nmat[v, j],frq))

        
    }
    
    colnames(pmat) <- Jx
    rownames(pmat) <- Vx
    
    colnames(nmat) <- Jx
    rownames(nmat) <- Vx
    return(list(pmat, nmat))
}
#*******************************************************************
# Afunction to convert incidenece matrix to adjecney matrix
make_adjeceny <- function(B)
{
  n  <- dim(B)[1] + dim(B)[2]
  nams <- c(rownames(B), colnames(B))
  A <- matrix(0, ncol = n, nrow = n)
  colnames(A) <- nams
  rownames(A) <- nams
  for (i in rownames(B))
  {
    for (j in colnames(B))
    {
      
        A[i, j] = B[i,j]
        A[j, i] = B[i,j]
      
    }
  }
  
  return(A)
}
#******************************************************************************
#function to find average length of cdr3 with V and J segments from the entire dataset
Full_VJ_combination <- function()
{
  
  d  <- readRDS("./Data/FullData.Rda")
  #d <- subset(d, d$Group == g)
  d <- d[!duplicated(d[,c(6,9, 11)]),]
  
  pmat <- matrix(0, ncol = length(J), nrow = length(V))
  colnames(pmat) <- J
  rownames(pmat) <- V
  pd <- subset(d, d$Productive.status=='true')
  for (i in 1:dim(pd)[1])
  { 
    v = pd$V[i] 
    j = pd$J[i]
    cdr = pd$CDR3.amino.acid.sequence[i]
    frq <- pd$Frequency....[i]
    #pmat[v, j] <- ifelse(pmat[v, j]==0, nchar(cdr), (pmat[v, j]+nchar(cdr))/1)#*******************************************************
    pmat[v, j] <- pmat[v, j]+ nchar(cdr)
    #pmat[v, j] <- ifelse(pmat[v, j]==0, frq, mean(pmat[v, j],frq))
    
  }
  
  nmat <- matrix(0, ncol = length(J), nrow = length(V))
  colnames(nmat) <- J
  rownames(nmat) <- V
  nd <- subset(d, d$Productive.status=='false')
  
  for (i in 1:dim(nd)[1])
  { 
    v = nd$V[i] 
    j = nd$J[i]
    cdr = nd$CDR3.amino.acid.sequence[i]
    frq <- nd$Frequency....[i]
    #nmat[v, j] <- ifelse(nmat[v, j]==0, nchar(cdr), (nmat[v, j]+nchar(cdr))/1)#*******************************************
    nmat[v, j] <- nmat[v, j] + nchar(cdr)
    #nmat[v, j] <- ifelse(nmat[v, j]==0, frq, mean(nmat[v, j],frq))
  }
  
  for (vv in V)
  {
    for (jj in J)
    {
      xx = dim(pd[which(pd$V==vv & pd$J ==jj),])[1]
      yy = dim(nd[which(nd$V==vv & nd$J ==jj),])[1]
      #print(c(xx, yy))
      pmat[vv, jj] = ifelse(xx !=0, pmat[vv,jj]/xx, 0)
      nmat[vv, jj] = ifelse(yy!=0,  nmat[vv,jj]/yy, 0)
    }
  }
  
  
  colnames(pmat) <- Jx
  rownames(pmat) <- Vx
  
  colnames(nmat) <- Jx
  rownames(nmat) <- Vx
  return(list(pmat, nmat))
}
#***********************************************************
# function for network visulaiation
dummy.visualize <- function(A)
{
  #A <- make_adjeceny(A)
  longData <- melt(A)
  longData <- longData[longData$value!=0,]
  
  # creating a graph that represented as incidence matrix with matrix A
  g <- graph.incidence(A, weighted = TRUE, directed = FALSE)
  #g <- graph.adjacency(A, mode = c('undirected'), weighted = TRUE) 
  # clustering the graph g with Louvain algorithm
  lou <- cluster_louvain(g)
  df.lou <- data.frame(lou$names,lou$membership)
  V(g)$community <- lou$membership
  plot(g, edge.arrow.size=.1, vertex.size=13, vertex.label.cex=0.4, layout=layout_randomly, vertex.label.family= 'Helvetica', edge.color="orange")#, vertex.color=colrs[V(g)$community])
}
#************************************************************
dummy.function <- function(A)
{
   #A <- make_adjeceny(A)
   longData <- melt(A)
 
   longData <- longData[longData$value!=0,]
   
    # creating a graph that represented as incidence matrix with matrix A
    g <- graph.incidence(A, weighted = TRUE, directed = FALSE)
  
    #g <- graph.adjacency(A, mode = c('undirected'), weighted = TRUE) 
    # clustering the graph g with Louvain algorithm
    lou <- cluster_louvain(g)
    df.lou <- data.frame(lou$names,lou$membership)
    
    
    
   
    longData <- left_join(longData, df.lou, by=c("Var1"="lou.names"))
    #print(longData[1:5,])
    colnames(longData)[4] <- "Var1_clust"
    longData$Var2 <- as.factor(longData$Var2)
    longData <- left_join(longData, df.lou, by=c("Var2"="lou.names"))
    colnames(longData)[5] <- "Var2_clust"
    longData$colour <- ifelse(longData$Var1_clust==longData$Var2_clust, longData$Var1_clust, 0)
   
    
    
    longData$Var1 <- factor(longData$Var1, levels=unique(arrange(longData, Var1_clust)[,1]))
    longData$Var2 <- factor(longData$Var2, levels=unique(arrange(longData, Var2_clust)[,2]))
    #levels must be names
    longData$colour <- factor(longData$colour)
    longData$value <- factor(longData$value)
    
    #print(longData[1:10, ])
    
    #print(factor(longData$value, levels=unique(arrange(longData, value))))
    #print(unique(longData$colour))
    #print(colnames(longData))
    #for colours variabes must be factors (discrete scale) otherwise ggplot recognize it continous
    
    jcolors <- c("gray80",'blue', "#B40404", "#0B6121", "#FFBF00", "#008080", '#FF00FF', "#00FFFF")

    g <- ggplot(longData, aes(x = Var2, y = Var1, fill=colour)) + 
        geom_raster() + 
        scale_fill_manual(values=jcolors[1:length(unique(longData$colour))]) + 
        labs(x="J segments", y="V segments", fill ="Avg CDR3 length" ) +
        theme_bw() + theme(axis.text.x=element_text(size=8, angle=40, color = 'black', vjust=0.5),
                           axis.text.y=element_text(size=6, angle = 10, color = 'black'),
                           plot.title=element_text(size=15, hjust = 0.5),
                           legend.text=element_text(size=7)) + theme(legend.position = "none")

    return(g)
}
#****************************************************************************************
# function to generate v-j combination cluster for each group.
generate.vj.combination <- function(target.group)
{
    pn <- VJ_combination(target.group)
    gp <- dummy.function(pn[[1]])
    gn <- dummy.function(pn[[2]])
    g <- ggarrange(gp, gn,  labels = c("A", "B"),ncol = 2, nrow = 1)
    return(g)
}
#****************************************************************************************
# function to generate  v-j combination cluster for entire dataset
generate.full.vj.combination <- function()
{
  pn <- Full_VJ_combination()
  gp <- dummy.function(pn[[1]])
  gn <- dummy.function(pn[[2]])
  g <- ggarrange(gp, gn,  labels = c("A", "B"),ncol = 2, nrow = 1)
  return(g)
}
########################################################
plot.all <- function(i)
{
    groups <- c("Asymptomatic", "Recovery","Severe","Controls", "Uncomplicated")
    asm = VJ_combination(groups[1])[[i]]
    asm = dummy.function(asm)
    
    rec = VJ_combination(groups[2])[[i]]
    rec = dummy.function(rec)
    
    sev = VJ_combination(groups[3])[[i]]
    sev = dummy.function(sev)
    
    con = VJ_combination(groups[4])[[i]]
    con = dummy.function(con)
    unc = VJ_combination(groups[5])[[i]]
    unc = dummy.function(unc)
    g = ggarrange(con, rec, asm, unc, sev,  labels = LETTERS[1:5],ncol = 3, nrow = 2)
    
    return(g)
}