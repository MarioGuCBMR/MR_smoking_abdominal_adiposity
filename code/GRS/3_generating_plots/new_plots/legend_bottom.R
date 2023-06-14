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
  theme_bw() +
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0.5, size = 11,
                                    # face = "bold",
                                    margin = margin(t = 10, b = 10))) +
  labs(caption = paste0("■ Lifetime smoking-WHRadjBMI GRS   ", 
                        "▲ Lifetime smoking-WHR GRS")) +
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 10, face = "bold"), 
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
