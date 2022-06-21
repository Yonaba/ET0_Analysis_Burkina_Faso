Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(lmtest)
library(readxl)
library(writexl)

LoadFunction <- function(file,...) {
  
  dots <- match.call(expand.dots = FALSE)$...
  dots <- sapply(dots, as.character)
  
  output <- lapply(dots, function(x,file){eval(parse(text=paste(x," <- function(x) {0}",sep="")),envir = .GlobalEnv)
    suppressMessages(insertSource(file, functions=x))
    eval(parse(text=paste(x," <- ",x,"@.Data",sep="")),envir = .GlobalEnv) },file=file)
  
}

station_data <- function(path, vars, pet, s) {
  n  <- length(pet$Date)
  df <- data.frame(matrix(nrow = n, ncol = 0))
  df$dates <- pet$Date
  for (var in vars) {
    #var <- vars[1]
    dd <- read_excel(path, sheet = var, guess_max = 10000)
    dd$Date <- NULL
    colnames(dd) <- colnames(pet)[-1]
    #print(paste0("var", var))
    df <- cbind(df, dd[,s])
  }
  df$pet <- pet[,s]
  colnames(df) <- c("dates", vars, "pet")
  return (df)
}

LoadFunction(file="PET_Formulas.R", PET.PM)

Gsc <- 0.0820
var_path <- "biasc/obs_complete_1988_2017.xlsx"
et0_path <- "et0/ET0_OBS_1988_2017.xlsx"
SENS_MIN <- -0.25
SENS_MAX <- 0.25
SENS_STEP <- 0.025

rsens <- c(seq(SENS_MIN,-SENS_STEP, by=SENS_STEP), seq(SENS_STEP,SENS_MAX, by=SENS_STEP))

stations <- read.csv("bf_stations.csv", header = TRUE, sep = ",", dec = ".")
rownames(stations) <- tolower(stations$Name)

vars <- excel_sheets(var_path)
vars <- vars[-3]
pet.d <- data.frame(read_excel(et0_path, guess_max = 10000))
pet.d$Date <- as.Date(pet.d$Date)

locations <- colnames(pet.d)[-1]

sensitivity <- function(season = F, msg) {
  #season <- F
  print(msg)
  df.sens <- sapply(locations,function(x) NULL)
  df.csens <- data.frame(matrix(nrow = 0, ncol=length(vars)+1))
  
  for (s in locations) {
    #s <- "dori"
    print(paste0("station: ",s))
    lat <- stations[s,"Latitude"]
    z <- as.numeric(stations[s,"Elevation"])
    data.d <- station_data(var_path, vars, pet.d, s)
    
    #season <- "dry"
    if (season != F) {
      data.d$J <- as.numeric(format(data.d$dates, "%j"))
      incl <- (data.d$J < 152 | data.d$J > 304)
      data.d$incl <- if(season == "dry") incl else !incl
      data.d <- data.d[data.d$incl,]
      data.d$J <- data.d$incl <- NULL
    }
    
    dates <- data.d$dates
    data.d$dates <- NULL
    pet0 <- data.d$pet
    data.d$pet <- NULL
    colnames(data.d) <- tolower(colnames(data.d))
    
    cc.df <- data.frame(matrix(nrow = length(rsens), ncol = 0))
    cc.df$rsens <- rsens
    lcsens <- c(s)
    for (var in colnames(data.d)) {
      #var <- "rh"
      #print(var)
      cc <- c()
      for (sens in rsens) {
        #sens <- +0.4
        #print(sens)
        odf <- data.d[,var]
        ndf <- data.frame(data.d[,var])
        colnames(ndf) <- var
        ndf <- ndf*(1+sens)
        
        if (var == "rh") ndf[ndf>100] <- NA
        
        ndf <- cbind(dates, ndf,data.d[,!(colnames(data.d) %in% c(var))])
        colnames(ndf)[1] <- "Date"
        rpet <- PET.PM(ndf, lat = lat, z = z, checkTdiff = T)$PET
        #sensr <- mean(((((rpet-pet0)/pet0)/(1+rsens))*(odf/pet0)), na.rm=T)
        #sensr <- mean(((((rpet-pet0)/pet0)/(1+rsens))*((odf-modf)/pet0)), na.rm=T)
        sensr <- mean(((rpet-pet0)/(ndf[[var]]-odf))*(odf/pet0), na.rm=T)
        cc <- append(cc, sensr)
      }
      cc.df <- cbind(cc.df,cc)
      lmfit <- lm(cc ~ rsens)
      csens <- lmfit$coefficients
      print(paste0(var," (lmfit): ",round(summary(lmfit)$r.squared, digits = 4),
                   " - p-val: ",formatC(summary(lmfit)$coefficients[,4], format = "e", digits = 3)))
      lcsens <- append(lcsens, mean(cc))
    }
    colnames(cc.df) <- c("sens",colnames(data.d)) 
    df.sens[[s]] <- cc.df
    df.csens[nrow(df.csens)+1,] <- lcsens
    
  }
  rownames(df.csens) <- locations
  colnames(df.csens) <- c("station",vars)
  df.csens[,c(2:NCOL(df.csens))] <- sapply(df.csens[,c(2:NCOL(df.csens))], as.numeric)
  
  ll <- (list(df.sens, df.csens))
}

ann_daily_sens <- sensitivity(season = F, "daily annual sensitivity evaluation")
dry_daily_sens <- sensitivity(season = "dry", "daily dry seasonal sensitivity evaluation")
wet_daily_sens <- sensitivity(season = "wet", "daily wet seasonal sensitivity evaluation")

print("Writing varying sensitivities")
write_xlsx(
  ann_daily_sens[[1]],
  path = "sens/sens_var_PET_ann_daily.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

write_xlsx(
  dry_daily_sens[[1]],
  path = "sens/sens_var_PET_dry_daily.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

write_xlsx(
  wet_daily_sens[[1]],
  path = "sens/sens_var_PET_wet_daily.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

write_xlsx(
  ann_daily_sens[[2]],
  path = "sens/sens_var_coeff_ann.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

write_xlsx(
  dry_daily_sens[[2]],
  path = "sens/sens_var_coeff_dry.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

write_xlsx(
  wet_daily_sens[[2]],
  path = "sens/sens_var_coeff_wet.xlsx",
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

print("All done.")
