library(PtProcess)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(data.table)
library(tidyverse)

plot_data_initiation <- function(whradjbmi_path, title1, 
                                 level1, level2, level3, level4,
                                 main_title){
  grs_df <- fread(whradjbmi_path)
  
  grs_df$GRS <- title1
  
  #Let's make the nsnps a bit different:
  grs_df$nsnps <- paste("n=", grs_df$nsnps, sep = "")
  
  #Let's make those plots!!
  
  #Now that they are merged, let's add the full effects, but let's be careful.
  #Because some of the 0.10 transform to 0.1
  
  betas_rounded <- ifelse(nchar(round(grs_df$beta, 2)) == 3, paste(round(grs_df$beta, 2), "0", sep=""), round(grs_df$beta, 2))
  lower_ci_rounded <- ifelse(nchar(round(grs_df$lower_ci, 2)) == 3, paste(round(grs_df$lower_ci, 2), "0", sep=""), round(grs_df$lower_ci, 2))
  upper_ci_rounded <- ifelse(nchar(round(grs_df$upper_ci, 2)) == 3, paste(round(grs_df$upper_ci, 2), "0", sep=""), round(grs_df$upper_ci, 2))
  
  grs_df$full_effect <- paste(betas_rounded,
                              " (",
                              lower_ci_rounded,
                              ", ",
                              upper_ci_rounded,
                              ")",
                              sep=""
  )
  
  #STEP 0: let's transform in factors the traits that we columns that we need:
  
  grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3, level4))
  
  grs_df$GRS = factor(grs_df$GRS, levels=c(title1))
  
  plotio <-
    
    #SETTING the info where we are all gonna work
    
    ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=-0.4, ymax=0.55)) +
    
    #Generating geom_point:
    
    geom_point(aes(color = "black", fill="black"), size = 2.5,  position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("black" = "black")) +
    scale_fill_manual(values = c("black" = "black")) +
    scale_shape_manual(name= "GRS", values=c("Smoking initiation-WHRadjBMI GRS" = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    
    #Adding the error bars:
    
    geom_errorbar(aes(ymin=as.numeric(lower_ci), ymax=as.numeric(upper_ci)), width = 0.15, cex= 1, position = position_dodge(width = 0.75), color="black", size=1) +
    
    #Generating a facet to distribute the data beautifully:
    
    #facet_wrap(~type, strip.position="left", nrow=7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype=2) +
    
    #First geom_text:
    
    #geom_text(aes(label=grs_df$nsnps),
    #          color = "black", position = position_dodge(width = 0.75), vjust = 1.6, size=3.5) +
    
    #Second geom_text:
    geom_text(aes(label=as.character(grs_df$full_effect)),
              color = "black", position = position_dodge(width = 0.75), vjust =-1.5, size=4) +
    
    coord_flip() +
    
    #SETTING AXIS OPTIONS:
    
    xlab("") +
    ylab("effect size (sd)") +
    
    #SETTING TITLE OPTIONS:
    
    ggtitle(paste0(main_title)) +
    theme_bw() +
    theme(legend.position="none") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold"), plot.title = element_text(size = 14, face = "bold", hjust=0.5))
  
  return(plotio)
}

plot_data_lifetime <- function(whradjbmi_path, title1, level1, level2, level3, level4, main_title){
  
  grs_df <- fread(whradjbmi_path)
  
  grs_df$GRS <- title1
  
  #Let's make the nsnps a bit different:
  grs_df$nsnps <- paste("n=", grs_df$nsnps, sep = "")
  
  #Now that they are merged, let's add the full effects, but let's be careful.
  #Because some of the 0.10 transform to 0.1
  
  betas_rounded <- ifelse(nchar(round(grs_df$beta, 2)) == 3, paste(round(grs_df$beta, 2), "0", sep=""), round(grs_df$beta, 2))
  lower_ci_rounded <- ifelse(nchar(round(grs_df$lower_ci, 2)) == 3, paste(round(grs_df$lower_ci, 2), "0", sep=""), round(grs_df$lower_ci, 2))
  upper_ci_rounded <- ifelse(nchar(round(grs_df$upper_ci, 2)) == 3, paste(round(grs_df$upper_ci, 2), "0", sep=""), round(grs_df$upper_ci, 2))
  
  grs_df$full_effect <- paste(betas_rounded,
                              " (",
                              lower_ci_rounded,
                              ", ",
                              upper_ci_rounded,
                              ")",
                              sep=""
  )
  
  #STEP 0: let's transform in factors the traits that we columns that we need:
  
  grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3, level4))
  
  grs_df$GRS = factor(grs_df$GRS, levels=c(title1))
  
  plotio <-
    
    #SETTING the info where we are all gonna work
    
    ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=-0.4, ymax=0.55)) +
    # ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=min(grs_df$lower_ci - 0.1), ymax=max(grs_df$upper_ci + 0.1), color = Trait, shape=GRS)) +
    
    #Generating geom_point:
    
    geom_point(aes(color = "black", fill="black"), size = 2.5,  position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("black" = "black")) +
    scale_fill_manual(values = c("black" = "black")) +
    # scale_shape_manual(name= "GRS", values=c("Lifetime smoking-WHRadjBMI GRS" = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    #Adding the error bars:
    
    geom_errorbar(aes(ymin=as.numeric(lower_ci), ymax=as.numeric(upper_ci)), width = 0.15, cex= 1, position = position_dodge(width = 0.75), color="black", size=1) +
    
    #Generating a facet to distribute the data beautifully:
    
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype=2) +
    
    #First geom_text:
    
    #geom_text(aes(label=grs_df$nsnps),
    #          color = "black", position = position_dodge(width = 0.75), vjust = 1.6, size=3.5) +
    
    #Second geom_text:
    geom_text(aes(label=as.character(grs_df$full_effect)),
              color = "black", position = position_dodge(width = 0.75), vjust =-1.5, size=4) +
    
    coord_flip() +
    
    #SETTING AXIS OPTIONS:
    
    xlab("") +
    ylab("effect size (sd)") +
    
    #SETTING TITLE OPTIONS:
    
    ggtitle(paste0(main_title)) +
    theme_bw() +
    theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=12,face="bold"), plot.title = element_text(size = 14, face = "bold", hjust=0.5))
  
  return(plotio)
}

## LIFETIME
p1 <- plot_data_lifetime("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_weighted.txt",
                         "Lifetime smoking-WHRadjBMI GRS", 
                         "Cortisol", "Oestradiol", "Testosterone", "SHBG", 
                         "Lifetime smoking")


# SMOKING INITIATION
p2 <- plot_data_initiation("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_weighted.txt",
                           "Smoking initiation-WHRadjBMI GRS", 
                           "Cortisol", "Oestradiol", "Testosterone", "SHBG",
                           "Smoking initiation")

gt1 <- ggplotGrob(p1)
gt2 <- ggplotGrob(p2)

newWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3])

gt1$widths[2:3] = as.list(newWidth)
gt2$widths[2:3] = as.list(newWidth)

output_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/final_hormones_lifetime_AND_smoking_initiation_whradjbmi_weighted_grs_plot.tiff"
tiff(output_path, width = 4000, height=3000, res = 500)
grid.arrange(gt1, gt2, ncol=2, widths = c(4,3.1))
dev.off()

### SAVING DATA FOR TABLES
## LIFETIME
whradjbmi_path = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_lifetime_whradjbmi_grs_res_weighted.txt"
title1 = "Lifetime smoking-WHRadjBMI GRS"
level1 = "Cortisol"
level2 = "Oestradiol"
level3 = "Testosterone"
level4 = "SHBG"
grs_df <- fread(whradjbmi_path)
grs_df$GRS <- title1
grs_df$nsnps <- paste("n=", grs_df$nsnps, sep = "")
betas_rounded <- ifelse(nchar(round(grs_df$beta, 2)) == 3, paste(round(grs_df$beta, 2), "0", sep=""), round(grs_df$beta, 2))
lower_ci_rounded <- ifelse(nchar(round(grs_df$lower_ci, 2)) == 3, paste(round(grs_df$lower_ci, 2), "0", sep=""), round(grs_df$lower_ci, 2))
upper_ci_rounded <- ifelse(nchar(round(grs_df$upper_ci, 2)) == 3, paste(round(grs_df$upper_ci, 2), "0", sep=""), round(grs_df$upper_ci, 2))
grs_df$full_effect <- paste(betas_rounded," (", lower_ci_rounded,", ",upper_ci_rounded,")", sep="")
grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3, level4))
grs_df$GRS = factor(grs_df$GRS, levels=c(title1))
fwrite(grs_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/code/pepa/output_tables/Hormones_Lifetime.csv")

## INITIATION
whradjbmi_path = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_weighted.txt"
title1 = "Smoking initiation-WHRadjBMI GRS"
level1 = "Cortisol"
level2 = "Oestradiol"
level3 = "Testosterone"
level4 = "SHBG"
grs_df <- fread(whradjbmi_path)
grs_df$GRS <- title1
grs_df$nsnps <- paste("n=", grs_df$nsnps, sep = "")
betas_rounded <- ifelse(nchar(round(grs_df$beta, 2)) == 3, paste(round(grs_df$beta, 2), "0", sep=""), round(grs_df$beta, 2))
lower_ci_rounded <- ifelse(nchar(round(grs_df$lower_ci, 2)) == 3, paste(round(grs_df$lower_ci, 2), "0", sep=""), round(grs_df$lower_ci, 2))
upper_ci_rounded <- ifelse(nchar(round(grs_df$upper_ci, 2)) == 3, paste(round(grs_df$upper_ci, 2), "0", sep=""), round(grs_df$upper_ci, 2))
grs_df$full_effect <- paste(betas_rounded," (", lower_ci_rounded,", ",upper_ci_rounded,")", sep="")
grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3, level4))
grs_df$GRS = factor(grs_df$GRS, levels=c(title1))
fwrite(grs_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/code/pepa/output_tables/Hormones_Initiation.csv")
