library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(ggpubr)
library(ggseqlogo)
library(viridis)


aalist <- c('A', 'C', 'D', 'E', 'F',  'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
#***************************************************************************************************************
#function to get the statistics from the unclean dataset
get_statistics <- function()
{
  groups <- c("Asymptomatic", "Recovery","Severe","Controls", "Uncomplicated")
  dfm <- matrix(0, nrow = 6, ncol = 3, dimnames = list(c(groups, "Total")))
  d  <- readRDS("./Data/FullData.Rda")
  dfm['Total', 1] = dim(subset(d, d$Productive.status=='true'))[1]
  dfm['Total', 2] = dim(subset(d, d$Productive.status=='false'))[1]
  dfm['Total', 3] = dim(d)[1]
  
  for (g in groups)
  {
    d1 <- subset(d, d$Group == g)
    #d1 <- d1[!duplicated(d1[,c(6,9, 11)]),]
    yps <- subset(d1, d1$Productive.status =='true')
    yne <- subset(d1, d1$Productive.status =='false')
    
    dfm[g, 1] = dim(yps)[1]
    dfm[g,2] = dim(yne)[1]
    dfm[g, 3] = dim(d1)[1]
    
    
  }
  dfm <- as.data.frame(dfm)
  colnames(dfm) <- c("Productive", "Non-productive", "Total")
  return(dfm)
}
#***************************************************************************************************************
#function to get the statistics from the clean dataset
get_clean_statistics <- function()
{
  groups <- c("Asymptomatic", "Recovery","Severe","Controls", "Uncomplicated")
  dfm <- matrix(0, nrow = 6, ncol = 3, dimnames = list(c(groups, "Total")))
  d  <- readRDS("./Data/FullData.Rda")
  dfm['Total', 1] = dim(subset(d, d$Productive.status=='true'))[1]
  dfm['Total', 2] = dim(subset(d, d$Productive.status=='false'))[1]
  dfm['Total', 3] = dim(d)[1]
  
  for (g in groups)
  {
    d1 <- subset(d, d$Group == g)
    #d1 <- d1[!duplicated(d1[,c(6,9, 11)]),]
    yps <- subset(d1, d1$Productive.status =='true')
    yne <- subset(d1, d1$Productive.status =='false')
    
    dfm[g, 1] = dim(yps)[1]
    dfm[g,2] = dim(yne)[1]
    dfm[g, 3] = dim(d1)[1]
    
    
  }
  dfm <- as.data.frame(dfm)
  colnames(dfm) <- c("Productive", "Non-productive", "Total")
  return(dfm)
}
#*****************************************************************************
## A function to analyz the CDR3 length of  each  group
getData_density <- function()
{
  groups <- c("Asymptomatic", "Recovery","Severe","Controls", "Uncomplicated")
  #groups <- c('Controls')
  d  <- readRDS("./Data/FullData.Rda")
  pmax_len <- 355177 
  nmax_len <- 38789
  
  pdf <- as.data.frame(matrix(NA, nrow = pmax_len, ncol = 6))
  ndf <- as.data.frame(matrix(NA, nrow = nmax_len, ncol = 6))
  
  colnames(pdf) <- c('Status', groups)
  colnames(ndf) <- c('Status', groups)
  
  pdf$Status <- rep('true', pmax_len)
  ndf$Status <- rep('false', nmax_len)
  
  
  for (i in 1:length(groups))
  {
    group <- groups[i]
    d1 <- subset(d, d$Group == group)
    d1 <- d1[!duplicated(d1[,c(6,9, 11)]),]
    yps <- subset(d1, d1$Productive.status =='true')
    yne <- subset(d1, d1$Productive.status =='false')
    
    yps <- nchar(yps$CDR3.amino.acid.sequence)
    yne <- nchar(yne$CDR3.amino.acid.sequence)
    pdf[group] <- c(yps, rep(NA, pmax_len - length(yps)))
    ndf[group] <- c(yne, rep(NA, nmax_len - length(yne)))
  }
  df <- rbind(pdf, ndf)
  data <- melt(df)
  data <- na.omit(data)
  mu <- ddply(data, "Status", summarise, grp.mean=mean(value))
  labels <- c(flase = 'non-productive', true = 'productive')
  p <- ggplot(data=data, aes(x=value, group=Status, fill=Status)) +geom_histogram(aes(y=..density..),   binwidth=.5, colour="black", fill="white") +geom_density(alpha=.6, adjust =1) + facet_grid(~variable, labeller = labeller(Status = labels)) +theme_minimal()+ theme( legend.position="bottom", axis.text.x = element_text( size = 14 ), axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 20), panel.spacing = unit(0.1, "lines"))+ geom_vline(data=mu, aes(xintercept=grp.mean, color=Status),   linetype="dashed") +   labs(y="Percent", x="Length") + xlim(5,20)
  
  return(p)
}
#********************************************************************************************************

# A function to analyz the CDR3 length of  each  group
getData <- function()
{
  groups <- c("Asymptomatic", "Recovery","Severe","Controls", "Uncomplicated")
  #groups <- c('Controls')
  d  <- readRDS("./Data/FullData.Rda")
  pmax_len <- 355177 
  nmax_len <- 38789
  
  pdf <- as.data.frame(matrix(NA, nrow = pmax_len, ncol = 6))
  ndf <- as.data.frame(matrix(NA, nrow = nmax_len, ncol = 6))
  
  colnames(pdf) <- c('Status', groups)
  colnames(ndf) <- c('Status', groups)
  
  pdf$Status <- rep('true', pmax_len)
  ndf$Status <- rep('flase', nmax_len)
  
  
  for (i in 1:length(groups))
  {
    group <- groups[i]
    d1 <- subset(d, d$Group == group)
    d1 <- d1[!duplicated(d1[,c(6,9, 11)]),]
    yps <- subset(d1, d1$Productive.status =='true')
    yne <- subset(d1, d1$Productive.status =='false')
    
    yps <- nchar(yps$CDR3.amino.acid.sequence)
    yne <- nchar(yne$CDR3.amino.acid.sequence)
    pdf[group] <- c(yps, rep(NA, pmax_len - length(yps)))
    ndf[group] <- c(yne, rep(NA, nmax_len - length(yne)))
  }
  df <- rbind(pdf, ndf)
  data <- melt(df)
  new_order <- with(data, reorder(variable , value, mean , na.rm=T))
  par(mar = c(3,4,3, 1))
  
  myplot <- boxplot(value ~ Status*new_order , data=data, boxwex=0.3 , 
                    ylab="length",
                    main="length of CDR3 AA sequences" ,
                    col=c("slateblue1" , "tomato") ,  xaxt="n")
  # getting the labels
  my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
  my_names <- my_names[seq(1 , length(my_names) , 2)]
  
  # x-axis lables
  axis(1, at = seq(1.5 , 2*length(groups) , 2), labels = my_names , tick=FALSE , cex=0.3)
  # adding vertical lines
  for(i in seq(0.5 , 20 , 2)){ abline(v=i,lty=1, col="grey")}
  
  # Adding  legend
  legend("topleft", legend = c("Productive", "Non-productive"), 
         col=c("slateblue1" , "tomato"),
         pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))
  # ddply(dfm, c("Status", "variable"), meanmarise, mean = mean(value, na.rm=T), sd = sd(value, na.rm=T),sem = sd(value, na.rm=T)/sqrt(length(value)))
  return(data)
}
###############################################################################################################################################
countcontent <- function(cdraa)
{
  c= sapply(aalist, function(x) 100*length(which(unlist(strsplit(cdraa, '')) == x))/nchar(cdraa) )
  return(c)
}
###############################################################################################################################################
aacontent <- function(group)
{
  groups <- c("Asymptomatic", "Recovery","Severe","Controls", "Uncomplicated")
  # groups <- c('Recovery', 'Controls')
  d  <- readRDS("./Data/FullData.Rda")
  
  
  pd <- as.data.frame(matrix(NA,  ncol = 22))
  nd <- as.data.frame(matrix(NA,  ncol = 22))
  
  colnames(pd) <- c('group', 'status', aalist)
  colnames(nd) <- c('group', 'status', aalist)
  
  
  for (i in 1:length(groups))
  {
    group <- groups[i]
    d1 <- subset(d, d$Group == group)
    d1 <- d1[!duplicated(d1[,c(6,9, 11)]),]
    yps <- subset(d1, d1$Productive.status =='true')
    yne <- subset(d1, d1$Productive.status =='false')
    
    pdf <- as.data.frame(matrix(NA, nrow = dim(yps)[1],  ncol = 22))
    ndf <- as.data.frame(matrix(NA, nrow = dim(yne)[1], ncol = 22))
    colnames(pdf) <- c('group', 'status', aalist)
    colnames(ndf) <- c('group', 'status', aalist)
    
    
    
    pdf['group'] <- yps['Group']
    pdf['status'] <- yps['Productive.status']
    
    cp <- yps$CDR3.amino.acid.sequence
    
    xp <- unlist(lapply(cp, countcontent))
    pdf[1:dim(yps)[1], 3:22] <- xp
    
    ndf['group'] <- yne['Group']
    ndf['status'] <- yne['Productive.status']
    cn <- yne$CDR3.amino.acid.sequence
    
    xn <- unlist(lapply(cn, countcontent))
    ndf[1:dim(yne)[1], 3:22] <- xn
    
    pd <- rbind(pd, pdf)
    
    nd <- rbind(nd, ndf)
  }
  
  df <- rbind(pd, nd)
  df <- df[2:dim(df)[1],]
  
  d <- df %>% group_by(group, status) %>% summarize( A = mean(A,  na.rm = TRUE), C = mean(C,  na.rm = TRUE), D= mean(D,  na.rm = TRUE), E = mean(E,  na.rm = TRUE),
                                                     F = mean(F,  na.rm = TRUE), G = mean(G,  na.rm = TRUE), H= mean(H,  na.rm = TRUE), I = mean(I,  na.rm = TRUE),
                                                     K = mean(K,  na.rm = TRUE), L = mean(L,  na.rm = TRUE), M= mean(M,  na.rm = TRUE), N = mean(N,  na.rm = TRUE),
                                                     P= mean(P,  na.rm = TRUE), Q = mean(Q,  na.rm = TRUE), R= mean(R,  na.rm = TRUE), S = mean(S,  na.rm = TRUE),
                                                     T = mean(T,  na.rm = TRUE), V = mean(V,  na.rm = TRUE), W= mean(W,  na.rm = TRUE), Y = mean(Y,  na.rm = TRUE))
  d <- as.data.frame(na.omit(d))
  return(d)
}
#************************************************************************************************************************************
bbox <- function(df)
{
  data <- melt(df)
  p <- ggplot(data, aes(variable, value, colour = status)) +  
    geom_linerange(aes(x = variable, ymin = 0,ymax = value, group = status), size = 0.5,
                   position = position_dodge(0.3))+  geom_point(aes(color = status),
                                                                position = position_dodge(0.3), size = 0.5) + scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))+
    theme_pubclean() +theme(legend.text = element_text(colour="blue", size=10))+ 
    labs(x = "Amino acids", y = " Average aa composition") +   labs(color = "Status\n") +
    scale_color_manual(labels = c("productive", "non-productive"), values = c("#0073C2FF", "#EFC000FF")) 
  return(p)
  
}
#*****************************************************************************************
# fuction to generate position frequency matrix 
pfm <- function(peptide)
{
  peptide <- unlist(strsplit(toupper(peptide), split = ""))
  # n <- length(peptide)
  mat <- matrix(0, ncol = 20, nrow = 20, dimnames = list(aalist))
  for(i in 1:min(length(peptide), 20))
  {
    if (peptide[i] %in% aalist)  
    {
      mat[peptide[i], i] <- 1
    }
  }
  return(mat)
}
#*************************************************************************
ppfm <- function(cdrs)
{
  pmat <- matrix(0, ncol = 20, nrow = 20, dimnames = list(aalist))
  for(x in cdrs)
  {
    pmat = pmat + pfm(x)
  }
  return(pmat)
}
#************************************************************************************************************************************
# function to design positional amino acid composition
composetionAA <- function()
{
  groups <- c("Asymptomatic", "Recovery","Severe","Controls", "Uncomplicated")
  d  <- readRDS("./Data/FullData.Rda")
  d1 <- subset(d, d$Group == groups[])
  d1 <- d1[!duplicated(d1[,c(6,9, 11)]),]
  yps <- subset(d1, d1$Productive.status =='true')
  yne <- subset(d1, d1$Productive.status =='false')
  pcdrs <- yps$CDR3.amino.acid.sequence
  ncdrs <- yne$CDR3.amino.acid.sequence
  pm <- ppfm(pcdrs)
  pm <- t(t(pm)*1/colSums(pm))
  pmm <- pm
  x <- apply(pm, MARGIN=c(2), mean)
  for (a in aalist){for (i in 1:20) if (pmm[a,i] <x[i]) pmm[a,i]=0}
  
  pn <- ppfm(ncdrs)
  pn <- t(t(pn)*1/colSums(pn))
  pnn <- pn
  xx <- apply(pm, MARGIN=c(2), mean)
  for (a in aalist){for (i in 1:20) if (pmm[a,i] <xx[i]) pnn[a,i]=0}
  p1 <- ggseqlogo(pmm,  method='p', seq_type='aa', stack_width = 0.90 ) + theme(axis.text.x = element_blank(), legend.position = 'none')
  p2 <- ggseqlogo(pnn, method ='p', seq_type='aa', stack_width = 0.90)
  gridExtra::grid.arrange(p1, p2)
}

#********************************************************************************************
countFreq <- function(x, txt, d)
{
  D <- subset.data.frame(d, d$Individual==x & d$Productive.status==txt)
  return(sum(D$Frequency....))
  
}
#********************************************************************************************

fre_count <- function()
{
  d  <- readRDS("./Data/FullData.Rda")
  individual <- unique(d$Individual)
  productive <- sapply(unique(d$Individual), countFreq, 'true', d)
  non_productive <- sapply(unique(d$Individual), countFreq, 'false', d)
  df <- data.frame(individual, productive, non_productive)
  df <- df[order(df$individual),]
  dff <- t(df[,2:3])
  barplot(dff, col= c(2,3), cex.names =0.6,
          font.axis=1, beside=T, names.arg= df$individual,
          legend=c('Productive', 'Non-productive'), xlab="individuals", font.lab=1, 
          ylab = 'Frequency (%)', ylim = c(0, 100))
}
#********************************************************************************************
fre_count_group <- function()
{
  d  <- readRDS("./Data/FullData.Rda")
  groups <- unique(d$Group)
  df <- matrix(nrow=5, ncol=2)
  rownames(df) <- groups
  
  for( g in groups)
  {
    dd <- d[d$Group==g, ]
    df[g, 1] <- sum(subset(d$Frequency...., d$Group == g & d$Productive.status=='true'))
    df[g, 2] <- sum(subset(d$Frequency...., d$Group == g & d$Productive.status=='false'))
  }
  colnames(df) <- c('Productive', 'Non productive')
  
  return(df)
}