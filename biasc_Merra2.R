Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(qmap)
library(readxl)
library(writexl)

fitMethod <- "PTF"
doType <- NULL
obs_path <- "biasc/obs_fill_rs_1988_2017.xlsx"
merra2_path <- "biasc/merra2_1988_2017.xlsx"

sheets <- excel_sheets(obs_path)
outxl <- sapply(sheets,function(x) NULL)

for (var in sheets) {
  #var <- sheets[1]
  print(paste("Processing",var))
  obs <- read_excel(obs_path, sheet = var, guess_max = 10000)
  merra2 <- read_excel(merra2_path, sheet = var, guess_max = 10000)
  #merra2[,c(2:NCOL(merra2))] <- sapply(merra2[,c(2:NCOL(merra2))], as.numeric)
  obs$Date <- merra2$Date <- as.Date(obs$Date)
  
  obs.m <- split(obs, format(obs$Date, "%m"))
  merra2.m <- split(merra2, format(merra2$Date, "%m"))
  
  bc <- colnames(obs)[-1]
  bc.df <- sapply(bc,function(x) NULL)
  
  for (mon in 1:12) {
    #mon <- 1  
    print(mon)
    for (s in bc) {
      #s <- bc[1]
      print(s)
      #TAU <- as.numeric(length(obs.m[[mon]][[s]])) * TAU_RATIO
      qfit <- fitQmap(obs.m[[mon]][[s]], merra2.m[[mon]][[s]],method = fitMethod, transfun = "scale", wet.day = F)
      bc.df[[s]][[mon]] <- doQmap(merra2.m[[mon]][[s]], qfit, type = doType)
    }
  }
  
  dates <- c()
  for (mon in 1:12) dates <- append(dates, obs.m[[mon]][["Date"]])
  df.f <- data.frame(matrix(nrow = length(dates), ncol = 0))
  df.f$Date <- dates 
  for (s in names(bc.df)) {
    cc <- c()
    for (mon in 1:12) {
      cc <- append(cc, bc.df[[s]][[mon]])
    }
    df.f[,s] <- cc
  }
  colnames(df.f) <- c("Date",names(bc.df))
  df.f <- df.f[order(df.f$Date),]
  
  outxl[[var]] <- df.f
}

write_xlsx(
  outxl,
  path = paste0("biasc/merra2_bc_1988_2017_",fitMethod,"scale.xlsx"),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

print("All done.")