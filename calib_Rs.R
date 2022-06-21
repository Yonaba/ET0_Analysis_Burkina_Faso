Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(lmtest)
library(readxl)
library(writexl)
library(ggpubr)

var_path <- "raw_data/calib_rs_draw.xlsx"
stations <- excel_sheets(var_path)
ss <- read.csv("bf_stations.csv", header = TRUE, sep = ",", dec = ".")
ss$Name <- tolower(ss$Name)
ss <- ss[order(ss$Latitude, decreasing=T),]
ss <- subset(ss, Name != "bogande")

plots <- sapply(1:9, function(x) NULL)
tb <- data.frame(matrix(nrow=0, ncol=8))
#dev.off()
i<- 0

for (s in ss$Name) {
  i <- i + 1
  print(s)
  #s <- "bobo"
  df.s <- data.frame(read_excel(var_path, sheet = s, guess_max = 10000))

  #df.s$Rs.Ra[df.s$Rs.Ra < 0.33] <- 0.33
  #df.s$Rs.Ra[df.s$Rs.Ra > 1] <- 1
  #df.s$n.N[df.s$n.N > 1] <- 1
  df.s <- subset(df.s, (n.N <= 1 & Rs.Ra <= 1 & Rs.Ra >= 0.33))
  lm <- lm(Rs.Ra ~ n.N, data = df.s)
  lmf <- lm$coefficients
  as <- lmf[1]
  bs <- lmf[2]
  plots[[i]] <- ggscatter(df.s, x = "n.N", y = "Rs.Ra",
            color = "black", shape = 1, size = 1, # Points color, shape and size
            add = "reg.line",
            add.params = list(color = "blue", fill = "red", size = 2), # Customize reg. line
            conf.int = T, # Add confidence interval
            cor.coef = F, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"),
            cor.coef.coord = c(0.1,0.9),
            title = paste0(letters[seq( from = 1, to = length(stations) )][i],") ",s)) + xlim(0,1) + ylim(0,1) + xlab("nn/N") + ylab("rs/Ra") +
    annotate(geom="text",x=0.5,y=0.1,label=paste0("R² = ",round(summary(lm)$r.squared,3),", p-value < 2.2e-16"))+
    annotate(geom="text",x=0.5,y=0.2,label=paste0("Rs/Ra = ",round(as,4),"+",round(bs,4),"*(nn/N)"))
  tb[nrow(tb)+1,] <- c(s, ss$Latitude[i], as, bs, as+bs, summary(lm)$r.squared, 
                       as.numeric(summary(lm)$coefficients[,4][1]), nrow(df.s))
}

ggarrange(plotlist = plots)
ggsave("graphs/calib_Rs.png", units="cm", dpi = 400, width = 25, height = 25, scale = 1.1)

colnames(tb) <- c("Station", "Latitude","as", "bs", "as+bs","R²", "p-value","n")
tb$`p-value` <- "< 2.2e-16"

write_xlsx(tb, path = "graphs/calib_rs.xlsx")
