# Importing required libraries....
library(stringr) # for string processing..
library(pROC)
library(ggplot2)
library(reshape2)
#*****************************************************************************************
# function to combine  auc values and curves from all the group. 
# #**************************************************************************************

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
aucgraph_all <- function()
{
    groups <- c("Con" = 'blue', "Asm" = 'red', "Seve"= 'blueviolet', 'Rec' = 'green')
    dfm_asm = readRDS("./Asymptomatic/Asymptomatic_five_fold_Results.Rda")
    df_asm <- myroc1(dfm_asm$value, dfm_asm$Observed)
    area_asm <- auc(dfm_asm$Observed, dfm_asm$value)
    
    dfm_con = readRDS("./Controls/Controls_five_fold_Results.Rda")
    df_con <- myroc1(dfm_con$value, dfm_con$Observed)
    area_con <- auc(dfm_con$Observed, dfm_con$value)
    
    dfm_sev = readRDS("./Severe/Severe_five_fold_Results.Rda")
    df_sev <- myroc1(dfm_sev$value, dfm_sev$Observed)
    area_sev <- auc(dfm_sev$Observed, dfm_sev$value)
    
    dfm_rec = readRDS("./Recovery/Recovery_five_fold_Results.Rda")
    df_rec <- myroc1(dfm_rec$value, dfm_rec$Observed)
    area_rec <- auc(dfm_rec$Observed, dfm_rec$value)
    
    g <- ggplot() + geom_line(data = df_asm, aes(x = fpr, y = tpr), color = "red", size = 0.3) 
    g <- g + scale_x_continuous(name = "False Positive Rate (1 - specificity)", breaks = seq(0,1, 0.1))
    g <- g + scale_y_continuous(breaks = seq(0,1, 0.1), "True Positive Rate (Sensitivity)")
    g <- g + theme(plot.title = element_text(size = 10,  vjust = 0.5, lineheight = 0.2)) + coord_equal()
    g <- g + geom_line(linetype = 2, aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)))
    g <- g + geom_line( aes(x = seq(0, 1, 0.1), y = 0)) + geom_line( aes(y = seq(0, 1, 0.1), x = 0))
    g <- g + geom_line( aes(x = seq(0, 1, 0.1), y = 1)) + geom_line( aes(y = seq(0, 1, 0.1), x = 1))
    g <- g + theme( axis.line = element_line(colour = "darkblue", size = 0.51, linetype = "solid"))
    g <- g + geom_line( data = df_con, aes(x = fpr, y = tpr), color = "blue", size = 0.3)
    g <- g + geom_line( data = df_sev, aes(x = fpr, y = tpr), color = "blueviolet", size = 0.3)
    g <- g + geom_line( data = df_rec, aes(x = fpr, y = tpr), color = "green", size = 0.3)

    g <- g + theme_bw() 
    
    return(g)
}
#####################################################################
############################### A function to generate ROC graph
aucgraph_all2 <- function()
{
    groups <- c("Con" = 'blue', "Asm" = 'red', "Seve"= 'blueviolet', 'Rec' = 'green')
    dfm_asm = readRDS("./Asymptomatic/Asymptomatic_five_fold_Results.Rda")
    df_asm <- myroc1(dfm_asm$value, dfm_asm$Observed)
    df_asm$Group = "Asymptomatic"
    area_asm <- auc(dfm_asm$Observed, dfm_asm$value)
    
    dfm_con = readRDS("./Controls/Controls_five_fold_Results.Rda")
    df_con <- myroc1(dfm_con$value, dfm_con$Observed)
    df_con$Group = "Control"
    area_con <- auc(dfm_con$Observed, dfm_con$value)
    
    dfm_sev = readRDS("./Severe/Severe_five_fold_Results.Rda")
    df_sev <- myroc1(dfm_sev$value, dfm_sev$Observed)
    df_sev$Group = 'Severe'
    area_sev <- auc(dfm_sev$Observed, dfm_sev$value)
    
    dfm_rec = readRDS("./Recovery/Recovery_five_fold_Results.Rda")
    df_rec <- myroc1(dfm_rec$value, dfm_rec$Observed)
    df_rec$Group = 'Recovery'
    area_rec <- auc(dfm_rec$Observed, dfm_rec$value)
    
    dfm_unc = readRDS("./Uncomplicated/Uncomplicated__full_Results.Rda")
    df_unc <- myroc1(dfm_unc$value, dfm_unc$Observed)
    df_unc$Group = 'Uncomplicated'
    area_unc <- auc(dfm_unc$Observed, dfm_unc$value)
    
    dfm_com = readRDS("./Combined/combined_full_Results.Rda")
    df_com <- myroc1(dfm_com$value, dfm_com$Observed)
    df_com$Group = 'Combined'
    area_com <- auc(dfm_com$Observed, dfm_com$value)
    
    dfm = rbind(df_asm, df_con, df_rec, df_sev, df_unc, df_com)
    dfm = as.data.frame(dfm)
    lbls <- c("Asymptomatic: 0.78", 'Combined: 0.70', "Control: 0.72", "Recovery: 0.79", "Severe: 0.68", "Uncomplicated: 0.70")
    g <- ggplot() + geom_line(data = dfm, aes(x = dfm$fpr, y = dfm$tpr, color = factor(dfm$Group)), size = 0.3) 
     g <- g + scale_color_manual(values=c("blue", "red", "#56B4E9", "#999999", "#E69F00", "green"), 
                                 labels = lbls) + scale_fill_hue(l=45, c=55)
    g <- g + scale_x_continuous(name = "False Positive Rate (1 - specificity)", breaks = seq(0,1, 0.1))
    g <- g + scale_y_continuous(breaks = seq(0,1, 0.1), "True Positive Rate (Sensitivity)")
    g <- g + theme(plot.title = element_text(size = 10,  vjust = 0.5, lineheight = 0.2)) + coord_equal()
    g <- g + geom_line(linetype = 2, aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1))) + theme_bw()
    g <- g + geom_line(aes(x = seq(0, 1, 0.1), y = 0)) + geom_line( aes(y = seq(0, 1, 0.1), x = 0))
    g <- g + geom_line( aes(x = seq(0, 1, 0.1), y = 1)) + geom_line( aes(y = seq(0, 1, 0.1), x = 1))
    g <- g + theme( axis.line = element_line(colour = "darkblue", size = 0.51, linetype = "solid"))
    g <- g + theme(legend.position = c(0.7, 0.3), legend.text = element_text(colour="black", size=7, face="bold"))
    g <- g + theme(legend.title = element_blank())
    print(c(area_asm, area_con, area_rec, area_sev, area_unc, area_com))
    
    return(g)
}