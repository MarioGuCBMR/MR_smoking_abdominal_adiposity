##############
#INTRODUCTION#
##############

#This code matches our SNPs data with the traits of interest
#harmonizes them so that they are properly aligned.
#And finally performs GRS, which are saved ina dataframe.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

plot_data <- function(whr_path, whradjbmi_path, title1, title2, level1, level2, level3, level4, level5, level6, main_title){
  whr_w_grs_df <- fread(whr_path)
  whradjbmi_w_grs_df <- fread(whradjbmi_path)
  
  whr_w_grs_df$GRS <- title1
  whradjbmi_w_grs_df$GRS <- title2
  
  #Let's make the nsnps a bit different:
  
  whr_w_grs_df$nsnps <- paste("n=", whr_w_grs_df$nsnps, sep = "")
  whradjbmi_w_grs_df$nsnps <- paste("n=", whradjbmi_w_grs_df$nsnps, sep = "")
  
  #Let's make those plots!!
  
  grs_df <- rbind(whr_w_grs_df, whradjbmi_w_grs_df)
  
  #Now that they are merged, let's do the following:
  
  #Now that they are merged, let's add the full effects, but let's be careful.
  #Because some of the 0.10 transform to 0.1
  
  grs_df$full_effect <- paste(formatC(grs_df$beta, format = "e", digits = 2),
                              " (",
                              formatC(grs_df$lower_ci, format = "e", digits = 2),
                              ",",
                              formatC(grs_df$upper_ci, format = "e", digits = 2),
                              ")",
                              sep="")
  
  #STEP 0: let's transform in factors the traits that we columns that we need:
  
  grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3, level4, level5, level6))
  
  grs_df$type <- c("Fat depot","Fat depot","Fat depot", "Ratio","Ratio","Ratio", "Fat depot","Fat depot","Fat depot", "Ratio","Ratio","Ratio")
  grs_df$type = factor(grs_df$type, levels=c("Fat depot", "Ratio"))
  
  grs_df$GRS = factor(grs_df$GRS, levels=c(title1, title2))
  
  plotio <-
    
    #SETTING the info where we are all gonna work
    
    # ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=-0.2, ymax=0.9, color = Trait, shape=GRS)) +
    ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=min(grs_df$lower_ci - 0.15), ymax=max(grs_df$upper_ci + 0.15), color = Trait, shape=GRS)) +
    
    #Generating geom_point:
    
    geom_point(aes(color = Trait, shape=GRS, fill=Trait), size = 4,  position = position_dodge(width = 0.75)) +
    scale_shape_manual(name= "GRS", values=c("Smoking initiation-WHRadjBMI GRS" = 17, "Smoking initiation-WHR GRS" = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    
    #Adding the error bars:
    
    geom_errorbar(aes(ymin=as.numeric(lower_ci), ymax=as.numeric(upper_ci)), width = 0.15, cex= 1, position = position_dodge(width = 0.75), color="black", size=1) +
    
    #Generating a facet to distribute the data beautifully:
    
    facet_wrap(~type, strip.position="left", nrow=7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype=2) +
    
    #First geom_text:
    
    #geom_text(aes(label=grs_df$nsnps),
    #          color = "black", position = position_dodge(width = 0.75), vjust = 1.6, size=3.5) +
    
    #Second geom_text:
    geom_text(aes(label=as.character(grs_df$full_effect)),
              color = "black", position = position_dodge(width = 0.75), vjust =-1.25, size=3.5) +
    
    # geom_text(aes(label=as.character(grs_df$full_effect)),
    #           color = "black", position = geom_text(aes(label=as.character(grs_df$full_effect)),
    #           color = "black", position = position_dodge(width = 0.75), hjust=-0.70, vjust =0.60, size=3.5) +
    
    coord_flip() +
    
    
    #SETTING AXIS OPTIONS:
    
    xlab("Fat depots") +
    ylab("effect sizes") +
    
    #SETTING TITLE OPTIONS:
    
    ggtitle(paste0(main_title)) +
    theme_bw() +
    theme(legend.position="right") +
    theme(legend.text=element_text(size=11)) +
    theme(strip.background = element_rect(colour="black", fill="white")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text=element_text(size=11), axis.title=element_text(size=10,face="bold"), plot.title = element_text(size = 12, face = "bold", hjust=0.5))
  
  return(plotio)
}

# SMOKING INITIATION
# whr_path = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whr_grs_res_weighted.txt"
# whradjbmi_path = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whradjbmi_grs_res_weighted.txt"
# title1 = "Smoking initiation-WHR GRS"
# title2 = "Smoking initiation-WHRadjBMI GRS"
# level1 = "VAT"
# level2 = "GSAT"
# level3 = "ASAT"
# level4 = "VAT/GSAT"
# level5 = "VAT/ASAT"
# level6 = "ASAT/GSAT"
# main_title = "Smoking initiation-abdominal adiposity genetic scores"
# output_path= "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/fat_lifetime_smoking_whr_whradjbmi_weighted_grs_plot.tiff"

output_path = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/fat_smoking_initiation_whr_whradjbmi_unweighted_grs_plot.tiff"
tiff(output_path, width = 4200, height=4000, res = 500)
plot_data("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whr_grs_res_unweighted.txt",
          "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whradjbmi_grs_res_unweighted.txt",
          "Smoking initiation-WHR GRS", "Smoking initiation-WHRadjBMI GRS", 
          "VAT", "GSAT", "ASAT", "VAT/GSAT", "VAT/ASAT", "ASAT/GSAT",
          "Smoking initiation-abdominal adiposity genetic scores")
dev.off()

output_path = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/fatadjbmi_smoking_initiation_whr_whradjbmi_unweighted_grs_plot.tiff"
tiff(output_path, width = 4200, height=4000, res = 500)
plot_data("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fatadjbmi_smoking_initiation_whr_grs_res_unweighted.txt",
          "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fatadjbmi_smoking_initiation_whradjbmi_grs_res_unweighted.txt",
          "Smoking initiation-WHR GRS", "Smoking initiation-WHRadjBMI GRS", 
          "VATadjBMI", "GSATadjBMI", "ASATadjBMI", "VAT/GSAT", "VAT/ASAT", "ASAT/GSAT",
          "Smoking initiation-abdominal adiposity genetic scores")
dev.off()

