Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(openair)
library(readxl)

fitMethod <- "PTFscale"
obs_path <- "biasc/obs_fill_rs_1988_2017.xlsx"
merra2_path <- "biasc/merra2_1988_2017.xlsx"
merra2bc_path <- paste0("biasc/merra2_bc_1988_2017_",fitMethod,".xlsx")

sheets <- excel_sheets(obs_path)
data.all <- data.frame(matrix(nrow=0, ncol = 4))
colnames(data.all) <- c("obsvalues", "modvalues", "model", "var")

sheets <- sheets[c(1:2,4:6)]

for (var in sheets) {
  print(paste("Processing",var))
  #var <- sheets[1]
  obs <- read_excel(obs_path, sheet = var, guess_max = 10000)
  merra2 <- read_excel(merra2_path, sheet = var, guess_max = 10000)
  merra2[,c(2:NCOL(merra2))] <- sapply(merra2[,c(2:NCOL(merra2))], as.numeric)
  merra2bc <- read_excel(merra2bc_path, sheet = var, guess_max = 10000)
  
  dates <- as.Date(obs$Date)
  obs$Date <- merra2$Date <- merra2bc$Date <- NULL
  
  obsvalues <- append(stack(obs)$values,stack(obs)$values)
  modvalues <- append(stack(merra2)$values,stack(merra2bc)$values)
  n <- length(obsvalues)/2
  model <- append(rep("merra2", n),rep("merra2bc", n))
  data <- data.frame(obsvalues, modvalues, model)
  data$var <- var
  data.all <- rbind(data.all, data)
}

res <- 1200
factor <- res/72 
png(paste0("graphs/taylorPlotbc_Merra_",fitMethod,".png"), width=800*factor, height=400 * factor, res=res, pointsize = 1)

TaylorDiagram(data.all, obs = "obsvalues", mod = "modvalues", group= "model", type = "var", 
              normalise = T, 
              cols = c("red", "green1"),
              annotate = "centered\n RMSE", auto.text = F) 
# insert ggplot green
dev.off()
