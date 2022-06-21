Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")

library(naniar)
library(readxl)
library(ggpubr)

obs <- read_excel("raw_data/obs.xlsx", sheet = "data", guess_max = 10000)
# obs.f <- read_excel("raw_data/obs_fill_rs.xlsx", sheet = "data", guess_max = 10000)


obs <- obs[,!grepl("_tm", colnames(obs))]
# obs <- obs[,!grepl("_nn", colnames(obs))]
# obs.f <- obs.f[,!grepl("bg_", colnames(obs.f))]
# obs.f <- obs.f[,!grepl("bg_", colnames(obs.f))]

dates <- as.Date(obs$Date)
obs$Date <- NULL
obs$year <- as.numeric(format(dates,"%Y"))
obs <- obs[obs$year >= 1988,]

bc <- colnames(obs)[-length(colnames(obs))]
bc <- matrix(bc,ncol=10, nrow=6,byrow=T)
plot.all <- sapply(1:10,function(x) NULL)

s <- c("Dori", "Ouahigouya","Bogande","Dedougou","Ouagadougou","Fada","Boromo","Bobo","Po","Gaoua")
for (i in 1:10) {
  #i <- 1
  d.obs <- obs[,bc[,i]]
  d.obs$year <- obs$year
  plot.all[[i]] <- gg_miss_fct(x = d.obs, fct = year) + xlab("") + ylab("") + ggtitle(s[i]) +
    scale_y_discrete(labels=c("nn","rh","rs","tn","tx","ws"))
    #theme(plot.margin = margin(rep(0.1,4), "cm"))
}

dev.off()
ggarrange(plotlist = plot.all,nrow = 4, ncol = 3)
ggsave("graphs/miss_var.png", units="cm", dpi = 400, width = 45, height = 30, scale = 0.8)

#obs.f$year <- as.numeric(format(dates,"%Y"))


#obs.f <- obs.f[obs.f$year >= 1988,]

#vis_miss(obs)
#vis_miss(obs.f)
# 
# gg_miss_var(x = obs, show_pct = F)
# gg_miss_var(x = obs.f, show_pct = F)
# 
# obs$year <- as.factor(obs$year)
# obs.f$year <- as.factor(obs.f$year)

# gg_miss_fct(x = obs, fct = year)
# gg_miss_fct(x = obs.f, fct = year)
