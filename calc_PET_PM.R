Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
source("PET_Formulas.R")
library(readxl)
library(writexl)

station_data <- function(path, vars, dates, stations, s) {
  n  <- length(dates)
  df <- data.frame(matrix(nrow = n, ncol = 0))
  df$dates <- dates
  for (var in vars) {
    #var <- vars[1]
    dd <- read_excel(path, sheet = var, guess_max = 10000)
    dd$Date <- NULL
    colnames(dd) <- stations
    #print(paste0("var", var))
    df <- cbind(df, dd[,s])
  }
  colnames(df) <- c("Date", vars)
  return (df)
}

ss <- c("bobo", "boromo", "dori", "dedougou", "fada", "gaoua", "ouaga", "ouahigouya", "po")
stations <- read.csv("bf_stations.csv", header = TRUE, sep = ",", dec = ".")
rownames(stations) <- tolower(stations$Name)

obs_path <- "biasc/obs_complete_1988_2017.xlsx"
vars <- excel_sheets(obs_path)

dates <- seq(as.Date("1988-01-01"), as.Date("2017-12-31"), by = "day")
pet.d <- data.frame(matrix(nrow = length(dates), ncol = 0))
pet.d$dates <- dates

for (s in ss) {
  lat <- stations[s,"Latitude"]
  z <- as.numeric(stations[s,"Elevation"])
  print(paste("Processing ",s,"Lat ",lat,"Alt ", z))
  cd <- station_data(obs_path, vars, dates, ss, s)
  
  cd$Date <- as.Date(cd$Date)
  colnames(cd) <- c("Date","tx", "tn","tm","rh","rs","ws")
  pet.d <- cbind(pet.d, PET.PM(cd, lat = lat, z = z)$PET)
  
}

colnames(pet.d) <- c("Date",ss)

write_xlsx(
  pet.d,
  path = "et0/ET0_OBS_1988_2017.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

print("All done.")