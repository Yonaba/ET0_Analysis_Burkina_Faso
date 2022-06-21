Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(readxl)
library(lmtest)
library(modifiedmk)

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

var_path <- "biasc/obs_complete_1988_2017.xlsx"
et0_path <- "et0/ET0_OBS_1988_2017.xlsx"

stations <- read.csv("bf_stations.csv", header = TRUE, sep = ",", dec = ".")
rownames(stations) <- tolower(stations$Name)

vars <- excel_sheets(var_path)
pet.d <- data.frame(read_excel(et0_path, guess_max = 10000))
pet.d$Date <- as.Date(pet.d$Date)

stations <- colnames(pet.d)[-1]
pv <- ss <- data.frame(matrix(nrow = 0, ncol=length(vars)+1))

dev.off()

for (s in stations) {
  #s <- "bobo"
  print(s)
  data.d <- station_data(var_path, vars, pet.d, s)
  data.d$Y <- as.numeric(format(data.d$dates,"%Y"))
  
  data.y <- aggregate(x=subset(data.d, select=-c(dates,Y)), by=list(data.d$Y), FUN = mean)
  pet_ann <- data.frame(data.d$Y, data.d$pet)
  pet_ann <- aggregate(x=pet_ann, by=list(pet_ann$data.d.Y), FUN=sum)
  data.y$pet <- pet_ann[,3]
  
  colnames(data.y)[1] <- "year"
  param <- colnames(data.y)[-1]
  pvs <- sss <- c()
  for (var in param) {
    #var <- "Tx"
    print(var)
    dwpval <- dwtest(formula = as.formula(paste0("year ~ ",var)), data = data.y)$p.value
    bgpval <- bgtest(formula = as.formula(paste0("year ~ ",var)), data = data.y)$p.value
    ts.data <- ts(data.y[,var], start = 1988, end = 2017, frequency = 1)
    ts.data <- data.y[,var]
    # plot(ts.data, type="l")
    y <- mkttest(ts.data)
    x <- pwmk(ts.data)  
    pval <- ifelse(((dwpval < 0.05) && (bgpval < 0.05)), x[4], y[5])
    sval <- ifelse(((dwpval < 0.05) && (bgpval < 0.05)), x[3], y[2])
    pvs <- append(pvs, pval)
    sss <- append(sss, sval)
    
  }
  pv[nrow(pv)+1,] <- pvs 
  ss[nrow(ss)+1,] <- sss 
  ts_ann <- ts(data.y$pet, start=1988, end = 2017)
  plot(ts_ann, type="l", main=s)
}
rownames(pv) <- rownames(ss) <- stations
colnames(pv) <- colnames(ss) <- param

