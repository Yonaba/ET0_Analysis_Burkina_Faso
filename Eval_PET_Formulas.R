Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(hydroGOF)
library(readxl)
library(writexl)


list.functions.in.source <- function(..., local=NULL) {
  tmp <- new.env(parent=parent.frame())
  source(..., local = tmp)
  funs <- names(tmp)[unlist(eapply(tmp, is.function))]
  return (funs[order(funs, decreasing = T)])
}

PET.funs <- list.functions.in.source("PET_Formulas.R")
PET.funs <- PET.funs[PET.funs != "PET.PM"]
source("PET_Formulas.R")

Gsc <- 0.0820
var_path <- "biasc/obs_complete_1988_2017.xlsx"


station_data <- function(path, vars, dates, stations, s) {
  n  <- length(dates)
  df <- data.frame(matrix(nrow = n, ncol = 0))
  df$dates <- dates
  for (var in vars) {
    dd <- read_excel(path, sheet = var, guess_max = 10000)
    dd$Date <- NULL
    colnames(dd) <- stations
    df <- cbind(df, dd[,s])
  }
  colnames(df) <- c("dates", vars)
  return (df)
}

stations <- read.csv("bf_stations.csv", header = TRUE, sep = ",", dec = ".")
rownames(stations) <- tolower(stations$Name)
locations <- c("bobo", "boromo", "dori", "dedougou", "fada", "gaoua", "ouagadougou", "ouahigouya", "po")
#locations <- "bobo"

vars <- excel_sheets(var_path)
np <- data.frame(read_excel("eval/pet_funs_params.xlsx", guess_max = 10000))
rownames(np) <- np$Formulas
np$Formulas <- NULL

dates <- seq(as.Date("1988-01-01"), as.Date("2017-12-31"), by="day")
metrics <- c("r","R2", "bR2", "d", "md","RMSE", "MAE", "PBIAS %","NSE","rNSE", "KGE")

#df <- data.frame(matrix(nrow=0, ncol=length(metrics)+2))
df <- data.frame(matrix(nrow=0, ncol=2))

for (f in PET.funs) {
  #f <- "PET.MTF.Horton"
  print(paste0("Evaluating: ", f))
  fun.args <- as.list(args(f))
  fun.args[[length(fun.args)]] <- NULL
  fun.args <- names(fun.args)
  k <- np[f,1]
  
  petPMc <- petfc <- c()
  
  for (s in locations) {
    #s <- "bobo"
    print(paste0(" Pooling >> ", s))
    cd <- station_data(var_path, vars, dates, locations, s)
    colnames(cd) <- tolower(colnames(cd))
    colnames(cd)[1] <- "Date"
    lat <- stations[s, "Latitude"]
    z <- as.numeric(stations[s, "Elevation"])
    
    tdiff <- cd$tx - cd$tn
    rem <- which(tdiff < 0)
    if (length(rem)>0) cd <- cd[-c(rem),]
    
    petPM <- PET.PM(cd, lat, z, checkTdiff = T)$PET
    def.args <- list(cd = cd, lat = lat, z = z)

    fargs <- vector("list", length(fun.args))
    names(fargs) <- fun.args
    for (name in names(fargs)) fargs[[name]] <- def.args[[name]]
    petf <- do.call(f, as.list(fargs))$PET
    
    petPMc <- c(petPMc, petPM)
    petfc <- c(petfc, petf)
  }
  
  
  mgof <- gof(petfc, petPMc, na.rm = T, method = "2012")
  mgof <- t(mgof[metrics,])
  
  rsq <- mgof[2]
  n <- length(petPMc)
  adjrsq <- 1-(((1-rsq)*(n-1))/(n-k-1)) 
  mgof <- c(mgof[1:2], adjrsq, mgof[3:length(mgof)])
  
  #df[nrow(df)+1,] <- c(f,mgof)
  df[nrow(df)+1,] <- c(f,n)
}

cnames <- c(metrics[1:2],"adj-R²",metrics[3:length(metrics)])
colnames(df) <- c("Formulas",cnames)

#df[,c(2:NCOL(df))] <- sapply(df[,c(2:NCOL(df))], as.numeric)

write_xlsx(
  df,
  path = "eval/nn.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)
