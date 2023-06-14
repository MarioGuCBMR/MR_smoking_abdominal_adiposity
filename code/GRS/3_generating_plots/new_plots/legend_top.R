library(ggplot2)
library(data.table)
library(forecats)
# HORMONES GRS
whr_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whr_grs_res_weighted.txt"
whradjbmi_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/data/grs_output/hormones_smoking_initiation_whradjbmi_grs_res_weighted.txt"
title1 = "Smoking initiation-WHR GRS"
title2 = "Smoking initiation-WHRadjBMI GRS"
level1 = "Cortisol"
level2 = "Testosterone"
level3 = "SHBG"
main_title = "Smoking initiation-abdominal adiposity genetic scores"
output_path= "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/German&Mario&MariaJose/plots/new/hormones_smoking_initiation_whr_whradjbmi_weighted_grs_plot.tiff"

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

betas_rounded <- ifelse(nchar(round(grs_df$beta, 2)) == 3, paste(round(grs_df$beta, 2), "0", sep=""), round(grs_df$beta, 2))
lower_ci_rounded <- ifelse(nchar(round(grs_df$lower_ci, 2)) == 3, paste(round(grs_df$lower_ci, 2), "0", sep=""), round(grs_df$lower_ci, 2))
upper_ci_rounded <- ifelse(nchar(round(grs_df$upper_ci, 2)) == 3, paste(round(grs_df$upper_ci, 2), "0", sep=""), round(grs_df$upper_ci, 2))

grs_df$full_effect <- paste(betas_rounded,
                            " (",
                            lower_ci_rounded,
                            ",",
                            upper_ci_rounded,
                            ")",
                            sep=""
)

#STEP 0: let's transform in factors the traits that we columns that we need:

grs_df$Trait <- factor(grs_df$trait, levels = c(level1, level2, level3))

# grs_df$type <- c("Hormones", "Hormones", "Hormones", "Hormones", "Hormones", "Hormones")
# grs_df$type = factor(grs_df$type, levels=c("Hormones"))

grs_df$GRS = factor(grs_df$GRS, levels=c(title1, title2))

tiff(output_path, width = 4200, height=4000, res = 500)

ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(beta), 
                   ymin = min(grs_df$lower_ci - 0.15), 
                   ymax = max(grs_df$upper_ci + 0.15), 
                   color = Trait, shape = GRS)) +
  geom_point(aes(fill = Trait), size = 4, position = position_dodge(width = 0.75)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  geom_errorbar(aes(ymin = as.numeric(lower_ci), ymax = as.numeric(upper_ci)), 
                width = 0.15, cex = 1, position = position_dodge(width = 0.75), 
                color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
  geom_text(aes(label = as.character(grs_df$full_effect)), 
            color = "black", position = position_dodge(width = 0.75), 
            vjust = -1.25, size = 3.5) +
  coord_flip() +
  xlab("Hormones") +
  ylab("Effect sizes") +
  ggtitle(paste0(main_title)) +
  labs(subtitle = paste0("                              ■ Lifetime smoking-WHRadjBMI GRS   ", 
                         "▲ Lifetime smoking-WHR GRS")) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 11)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 10, face = "bold"), 
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

dev.off()
