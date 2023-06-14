library(PtProcess)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(data.table)
library(tidyverse)

plot_data_initiation <- function(whr_path, whradjbmi_path, title1, title2, level1, level2, level3, level4, level5, level6, level7, level8, level9, main_title, output_table){ 
  whr_w_grs_df <- fread(whr_path)
  whradjbmi_w_grs_df <- fread(whradjbmi_path)
  
  whr_w_grs_df$GRS <- title1
  whradjbmi_w_grs_df$GRS <- title2
  
  #Let's make the nsnps a bit different:
  
  whr_w_grs_df$nsnps <- paste("n=", whr_w_grs_df$nsnps, sep = "")
  whradjbmi_w_grs_df$nsnps <- paste("n=", whradjbmi_w_grs_df$nsnps, sep = "")
  
  #Let's make those plots!!
  
  grs_df <- rbind(whr_w_grs_df, whradjbmi_w_grs_df)
  
  # Keeping adjBMI and non-adj with respect to the trait
  # If the trait is adjuested, WHRadjusted (otherwise unadjusted)
  subset_dfadjBMI <- grs_df[grepl("adjBMI", grs_df$trait) & grepl("adjBMI", grs_df$GRS), ]
  subset_df <- grs_df[!grepl("adjBMI", grs_df$trait) & !grepl("adjBMI", grs_df$GRS), ]
  grs_df <- rbind(subset_df, subset_dfadjBMI)
  
  #Now that they are merged, let's do the following:
  
  #Now that they are merged, let's add the full effects, but let's be careful.
  #Because some of the 0.10 transform to 0.1
  
  # betas_rounded <- ifelse(nchar(round(grs_df$beta, 2)) == 3, paste(round(grs_df$beta, 2), "0", sep=""), round(grs_df$beta, 2))
  # lower_ci_rounded <- ifelse(nchar(round(grs_df$lower_ci, 2)) == 3, paste(round(grs_df$lower_ci, 2), "0", sep=""), round(grs_df$lower_ci, 2))
  # upper_ci_rounded <- ifelse(nchar(round(grs_df$upper_ci, 2)) == 3, paste(round(grs_df$upper_ci, 2), "0", sep=""), round(grs_df$upper_ci, 2))
  
  grs_df <- grs_df %>%
    mutate(
      betas_rounded = ifelse(round(beta, 2) != 0, round(beta, 2), ifelse(round(beta, 3) != 0, round(beta, 3), sprintf("%.3E", beta))),
      se_rounded = ifelse(round(se, 2) != 0, round(se, 2), round(se, 3)),
      lower_ci_rounded = ifelse(round(lower_ci, 2) != 0, round(lower_ci, 2), round(lower_ci, 3)),
      upper_ci_rounded = ifelse(round(upper_ci, 2) != 0, round(upper_ci, 2), round(upper_ci, 3)),
      pval_rounded = sprintf("%.2E", pval)
    )
  
  grs_df$full_effect <- paste(
    ifelse(round(grs_df$betas_rounded, 2) != 0.00, 
           sprintf("%.2f", as.numeric(grs_df$betas_rounded)),
           sprintf("%.3f", as.numeric(grs_df$betas_rounded))),
    " (",
    sprintf("%.2f", as.numeric(grs_df$lower_ci_rounded)),
    ", ",
    sprintf("%.2f", as.numeric(grs_df$upper_ci_rounded)),
    ")",
    sep = ""
  )
  
  #STEP 0: let's transform in factors the traits that we columns that we need:
  
  grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3, level4, level5, level6, level7, level8, level9))
  
  grs_df$GRS = factor(grs_df$GRS, levels=c(title1, title2))
  
  fwrite(grs_df, output_table)
  
  plotio <-

    #SETTING the info where we are all gonna work

    ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=-0.4, ymax=0.55, color = "black")) +
    # ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=min(grs_df$lower_ci - 0.1), ymax=max(grs_df$upper_ci + 0.1), color = Trait, shape=GRS)) +

    #Generating geom_point:

    geom_point(aes(color="black"), size = 2.5,  position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("black" = "black")) +
    scale_fill_manual(values = c("black" = "black")) +
    # scale_shape_manual(name= "GRS", values=c("Smoking initiation-WHRadjBMI GRS" = 15, "Smoking initiation-WHR GRS" = 15)) +
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
              color = "black", position = position_dodge(width = 0.75), vjust =-1.25, size=4) +

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

plot_data_lifetime <- function(whr_path, whradjbmi_path, title1, title2, level1, level2, level3, level4, level5, level6, level7, level8, level9, main_title, output_table){ 
  whr_w_grs_df <- fread(whr_path)
  whradjbmi_w_grs_df <- fread(whradjbmi_path)
  
  whr_w_grs_df$GRS <- title1
  whradjbmi_w_grs_df$GRS <- title2
  
  #Let's make the nsnps a bit different:
  
  whr_w_grs_df$nsnps <- paste("n=", whr_w_grs_df$nsnps, sep = "")
  whradjbmi_w_grs_df$nsnps <- paste("n=", whradjbmi_w_grs_df$nsnps, sep = "")
  
  #Let's make those plots!!
  
  grs_df <- rbind(whr_w_grs_df, whradjbmi_w_grs_df)
  
  # Keeping adjBMI and non-adj with respect to the trait
  # If the trait is adjuested, WHRadjusted (otherwise unadjusted)
  subset_dfadjBMI <- grs_df[grepl("adjBMI", grs_df$trait) & grepl("adjBMI", grs_df$GRS), ]
  subset_df <- grs_df[!grepl("adjBMI", grs_df$trait) & !grepl("adjBMI", grs_df$GRS), ]
  grs_df <- rbind(subset_df, subset_dfadjBMI)
  
  #Now that they are merged, let's do the following:
  
  #Now that they are merged, let's add the full effects, but let's be careful.
  #Because some of the 0.10 transform to 0.1
  
  # betas_rounded <- ifelse(nchar(round(grs_df$beta, 2)) == 3, paste(round(grs_df$beta, 2), "0", sep=""), round(grs_df$beta, 2))
  # lower_ci_rounded <- ifelse(nchar(round(grs_df$lower_ci, 2)) == 3, paste(round(grs_df$lower_ci, 2), "0", sep=""), round(grs_df$lower_ci, 2))
  # upper_ci_rounded <- ifelse(nchar(round(grs_df$upper_ci, 2)) == 3, paste(round(grs_df$upper_ci, 2), "0", sep=""), round(grs_df$upper_ci, 2))
  
  grs_df <- grs_df %>%
    mutate(
      betas_rounded = ifelse(round(beta, 2) != 0, round(beta, 2), ifelse(round(beta, 3) != 0, round(beta, 3), sprintf("%.3E", beta))),
      se_rounded = ifelse(round(se, 2) != 0, round(se, 2), round(se, 3)),
      lower_ci_rounded = ifelse(round(lower_ci, 2) != 0, round(lower_ci, 2), round(lower_ci, 3)),
      upper_ci_rounded = ifelse(round(upper_ci, 2) != 0, round(upper_ci, 2), round(upper_ci, 3)),
      pval_rounded = sprintf("%.2E", pval)
    )
  
  grs_df$full_effect <- paste(
    ifelse(round(grs_df$betas_rounded, 2) != 0.00, 
           sprintf("%.2f", as.numeric(grs_df$betas_rounded)),
           sprintf("%.3f", as.numeric(grs_df$betas_rounded))),
    " (",
    sprintf("%.2f", as.numeric(grs_df$lower_ci_rounded)),
    ", ",
    sprintf("%.2f", as.numeric(grs_df$upper_ci_rounded)),
    ")",
    sep = ""
  )
  
  #STEP 0: let's transform in factors the traits that we columns that we need:
  
  grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3, level4, level5, level6, level7, level8, level9))
  
  grs_df$GRS = factor(grs_df$GRS, levels=c(title1, title2))
  
  fwrite(grs_df, output_table)
  
  plotio <-

    #SETTING the info where we are all gonna work

    ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=-0.4, ymax=0.55, color = "black")) +
    # ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(beta), ymin=min(grs_df$lower_ci - 0.1), ymax=max(grs_df$upper_ci + 0.1), color = Trait, shape=GRS)) +

    #Generating geom_point:

    geom_point(aes(color = "black", fill = "black"), size = 2.5,  position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("black" = "black")) +
    scale_fill_manual(values = c("black" = "black")) +
    # scale_shape_manual(name= "GRS", values=c("Lifetime smoking-WHRadjBMI GRS" = 15, "Lifetime smoking-WHR GRS" = 15)) +
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
              color = "black", position = position_dodge(width = 0.75), vjust =-1.25, size=4) +

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
p1 <- plot_data_lifetime("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_lifetime_whr_grs_res_weighted.txt",
                         "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_lifetime_whradjbmi_grs_res_weighted.txt",
                         "Lifetime smoking-WHR GRS", "Lifetime smoking-WHRadjBMI GRS", 
                         "VAT", "VATadjBMI", "GSAT", "GSATadjBMI", "ASAT",  "ASATadjBMI", "VAT/GSAT", "VAT/ASAT","ASAT/GSAT", 
                         "Lifetime smoking",
                         "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/MANUSCRIPT/Fat_Lifetime.csv")

# SMOKING INITIATION
p2 <- plot_data_initiation("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whr_grs_res_weighted.txt",
                           "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whradjbmi_grs_res_weighted.txt",
                           "Smoking initiation-WHR GRS", "Smoking initiation-WHRadjBMI GRS", 
                           "VAT", "VATadjBMI", "GSAT", "GSATadjBMI", "ASAT",  "ASATadjBMI", "VAT/GSAT", "VAT/ASAT","ASAT/GSAT", 
                           "Smoking initiation", 
                           "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/MANUSCRIPT/Fat_Initiation.csv")
gt1 <- ggplotGrob(p1)
gt2 <- ggplotGrob(p2)

newWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3])

gt1$widths[2:3] = as.list(newWidth)
gt2$widths[2:3] = as.list(newWidth)

output_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/MANUSCRIPT/final_fat_lifetime_AND_smoking_initiation_whr_whradjbmi_weighted_grs_plot.tiff"
tiff(output_path, width = 4000, height=3000, res = 500)
grid.arrange(gt1, gt2, ncol=2, widths = c(4,3.1))
dev.off()

# ### SAVING DATA FOR TABLES
# out_initiation_fat <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/code/pepa/output_tables/Fat_Initiation.csv"
# out_lifetime_fat <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/code/pepa/output_tables/Fat_Lifetime.csv"
# ## LIFETIME
# p1 <- plot_data_lifetime("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_lifetime_whr_grs_res_weighted.txt",
#                          "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_lifetime_whradjbmi_grs_res_weighted.txt",
#                          "Lifetime smoking-WHR GRS", "Lifetime smoking-WHRadjBMI GRS", 
#                          "VAT", "VATadjBMI", "GSAT", "GSATadjBMI", "ASAT",  "ASATadjBMI", "VAT/GSAT", "VAT/ASAT","ASAT/GSAT", 
#                          "Lifetime smoking", out_lifetime_fat)
# 
# # SMOKING INITIATION
# p2 <- plot_data_initiation("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whr_grs_res_weighted.txt",
#                            "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/fat_smoking_initiation_whradjbmi_grs_res_weighted.txt",
#                            "Smoking initiation-WHR GRS", "Smoking initiation-WHRadjBMI GRS", 
#                            "VAT", "VATadjBMI", "GSAT", "GSATadjBMI", "ASAT",  "ASATadjBMI", "VAT/GSAT", "VAT/ASAT","ASAT/GSAT", 
#                            "Smoking initiation", out_initiation_fat)
