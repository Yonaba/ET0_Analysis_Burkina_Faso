Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(ggplot2)
library(hrbrthemes)
library(gridExtra)
library(ggpubr)

stations <- read.csv("angstrom_coeff.csv", header = T, sep = ",", dec = ".")
lats <- stations$Latitude

# lm_lat_b <- lm(b ~ Elevation, stations)
# r2_b_lat <- summary(lm_lat_b)$r.squared
# pv_b_lat <- cor.test(stations$Elevation, stations$b, alternative = "two.sided", method = "spearman")$p.value
# 
# lm_lat_a <- lm(a ~ Elevation, stations)
# r2_a_lat <- summary(lm_lat_a)$r.squared
# pv_a_lat <- cor.test(stations$Elevation, stations$a, alternative = "two.sided", method = "spearman")$p.value

lm_lat_ab <- lm(a.b ~ Elevation, stations)
r2_ab_lat <- summary(lm_lat_ab)$r.squared
pv_ab_lat <- cor.test(stations$Elevation, stations$a.b, alternative = "two.sided", method = "spearman")$p.value

p1 <- ggplot(stations, aes(x=Elevation, y=a.b)) +
  annotate("text", x = 300, y = 0.55, label = paste0("y = ",round(lm_lat_ab$coefficients[2], 5),
                                                 "x +", round(lm_lat_ab$coefficients[1], 3)), fontface = 2, size = 5) +
  annotate("text", x = 300, y = 0.54, label = paste0("R²: ",round(r2_ab_lat, 3),
                                                     ", p-value: ",round(pv_ab_lat, 4)), fontface = 4, size = 5) + 
  # annotate("text", x = 275, y = 0.7, label = "a)", fontface = 2, size = 5) +
  geom_point() +
  labs(x = "Station Elevation (masl)", y = "Atmospheric transmissivity (as + bs)") +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  geom_text(aes(label=Name), size=3, vjust = 1) + 
  theme_bw()

# p2 <- ggscatter(stations, x = "Elevation", y = "a.b",
#           color = "black", shape = 21, size = 3, # Points color, shape and size
#           add = "reg.line",  # Add regressin line
#           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#           conf.int = TRUE, # Add confidence interval
#           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
#           cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")) + xlab("Elevation (masl)") + 
#           ylab("Atmospheric transmissivity (as + bs)")

# p2 <- ggplot(stations, aes(x=lats, y=b)) +
#   annotate("text", x = 11, y = 0.35, label = paste0("y = ",round(lm_lat_b$coefficients[2], 3),
#                                                     "x +", round(lm_lat_b$coefficients[1], 3)), fontface = 2) +
#   annotate("text", x = 11, y = 0.34, label = paste0("R²: ",round(r2_b_lat, 3)), fontface = 4) +  
#   annotate("text", x = 11, y = 0.33, label = paste0("p-value: ",round(pv_b_lat, 4)), fontface = 4) + 
#   annotate("text", x = 10.5, y = 0.46, label = "b)", fontface = 2) +
#   geom_point() +
#   labs(x = "Absolute latitude (degrees)", y = "Angstrom coefficient (b)") +
#   geom_smooth(method=lm , color="blue", fill="#69b3a2", se=TRUE) +
#   geom_text(aes(label=stations$Name), size=3, vjust = 2) + 
#   theme_ipsum()
# 
# p2


ggsave(filename = "graphs/angstrom_lat.png", p1, dpi = 500, scale = 2,
       width = 10, height = 10, units = "cm")
