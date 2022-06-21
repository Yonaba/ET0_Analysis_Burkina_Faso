Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(readxl)
library(nlstools)
library(hydroGOF)
library(transport)
library(Metrics)

LoadFunction <- function(file,...) {
  
  dots <- match.call(expand.dots = FALSE)$...
  dots <- sapply(dots, as.character)
  
  output <- lapply(dots, function(x,file){eval(parse(text=paste(x," <- function(x) {0}",sep="")),envir = .GlobalEnv)
    suppressMessages(insertSource(file, functions=x))
    eval(parse(text=paste(x," <- ",x,"@.Data",sep="")),envir = .GlobalEnv) },file=file)
  
}

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
  colnames(df) <- c("Date", vars)
  return (df)
}

LoadFunction(file="PET_Formulas.R", PET.PM)
LoadFunction(file="PET_Formulas_Opt.R", PET.COMB.Valiantzas2)

Gsc <- 0.0820
var_path <- "biasc/obs_complete_1988_2017.xlsx"
vars <- excel_sheets(var_path)
#vars <- vars[-3]

dates <- seq(as.Date("1988-01-01"), as.Date("2017-12-31"), by="day")
locations <- c("bobo", "boromo", "dori", "dedougou", "fada", "gaoua", "ouagadougou", "ouahigouya", "po")
#locations <- c("ouahigouya")

stations <- read.csv("bf_stations.csv", header = TRUE, sep = ",", dec = ".")
rownames(stations) <- tolower(stations$Name)

petfun <- PET.COMB.Valiantzas2
has_daily_wind <- R

pet0c <- petc <- c()
cdc <- data.frame(matrix(nrow = 0, ncol = 8))

for (s in locations) {
  #s <- "boromo"
  print(paste0("processing ", s))
  cd <- station_data(var_path, vars, dates, locations, s)
  colnames(cd)[-1] <- tolower(colnames(cd)[-1])
  
  lat <- stations[s, "Latitude"]
  z <- as.numeric(stations[s, "Elevation"])
  cd$lat <- lat
  cd$z <- z
  
  tdiff <- cd$tx - cd$tn
  rem <- which(tdiff < 0)
  if (length(rem)>0) cd <- cd[-c(rem),]
  
  cd_opt <- cd
  if (!has_daily_wind) cd_opt$ws <- mean(cd$ws)
  
  pet0 <- PET.PM(cd, lat = lat, z = z)$PET
  #cd$pet0 <- pet0
  pet <- petfun(cd_opt)
  
  pet0c <- c(pet0c, pet0)
  petc <- c(petc, pet)
  cdc <- rbind(cdc, cd_opt)
}

colnames(cdc) <- colnames(cd)
cdc$petx <- pet0c

dev.off()
par(mfrow = c(1,2))
plot(pet0c, petc, xlim = c(0,10), ylim=c(0,10))

nls(petx ~ petfun(cdc,a,b,c,d,e,n1,n2), data = cdc,
    start=list(a=0.051, b=2.4, c=0.00012, d=0.5*0.048, e=0.536*0.048, n1=0.5,n2=2 ),
    # upper = list(a=10),
    # lower = list(a=0),
    trace=TRUE, algorithm = "port", control=list(warnOnly=T))

pet_opt <- c()
for (s in locations) {
  #s <- "bobo"
  print(paste0("processing ", s))
  cd <- station_data(var_path, vars, dates, locations, s)
  colnames(cd)[-1] <- tolower(colnames(cd)[-1])
  
  if (!has_daily_wind) cd$ws <- mean(cd$ws)
  
  lat <- stations[s, "Latitude"]
  z <- as.numeric(stations[s, "Elevation"])  
  cd$lat <- lat
  cd$z <- z  

  tdiff <- cd$tx - cd$tn
  rem <- which(tdiff < 0)
  if (length(rem)>0) cd <- cd[-c(rem),]
  
  #pet <- petfun(cd,a=1.4624 , b=17.0470 , n1=-0.2435  , n2=1.0858, n3=  )
  pet <- petfun(cd, c=0.08284, n1=0.36471 ,n2=1.67277 ,d=6.63834 ,e=0.04250)
  pet_opt <- c(pet_opt, pet)
}
plot(pet0c, pet_opt, xlim = c(0,10), ylim=c(0,10))

metrics <- c("r","R2", "bR2","d", "md", "RMSE", "MAE", "PBIAS %","NSE","rNSE", "KGE")

#data.frame(gof(petc, pet0c, na.rm = T, method = "2012"),gof(pet_opt, pet0c, na.rm = T, method = "2012"))

mgof <- data.frame(gof(petc, pet0c, na.rm = T, method = "2012")[metrics,],
                   gof(pet_opt, pet0c, na.rm = T, method = "2012")[metrics,])

swd <- c(wasserstein1d(pet0c, petc)/sd(c(pet0c, petc)),
        wasserstein1d(pet0c, pet_opt)/sd(c(pet0c, pet_opt)))

mgof[nrow(mgof)+1,] <- swd
rownames(mgof)[nrow(mgof)] <- "swd"

colnames(mgof) <- c("raw", "opt")
t(mgof)
#View(mgof)

