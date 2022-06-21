Sys.setenv(TZ="UTC")
setwd("D:/Recherche/ET0/analysis")
library(zoo)
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
dcal <- 0.3

rsens <- c(seq(SENS_MIN,-SENS_STEP, by=SENS_STEP), seq(SENS_STEP,SENS_MAX, by=SENS_STEP))

stations <- read.csv("bf_stations.csv", header = TRUE, sep = ",", dec = ".")
stations <- subset(stations, stations$Name != "BOGANDE")
rownames(stations) <- stations$Name <- tolower(stations$Name)

vars <- excel_sheets(var_path)
#vars <- vars[-3]
pet.d <- data.frame(read_excel(et0_path, guess_max = 10000))
pet.d$Date <- as.Date(pet.d$Date)
colnames(pet.d)[8] <- "ouagadougou"

locations <- colnames(pet.d)[-1]
stations <- stations[order(stations$Latitude, decreasing = T),]
#sensitivity <- function(season = F, msg) {
  #season <- F
#print(msg)
df.sens <- sapply(locations,function(x) NULL)

for (s in stations$Name) {
  #s <- "dori"
  print(paste0("station: ",s))
  lat <- stations[s,"Latitude"]
  z <- as.numeric(stations[s,"Elevation"])
  data.d <- station_data(var_path, vars, pet.d, s)
  
  #season <- "dry"
  #if (season != F) {
  data.d$J <- as.numeric(format(data.d$dates, "%j"))
  
  dates <- data.d$dates
  data.d$dates <- NULL
  pet0 <- data.d$pet
  data.d$pet <- NULL
  colnames(data.d) <- tolower(colnames(data.d))
  
  cc.df <- data.frame(matrix(nrow = 366, ncol = 0))

  for (var in tolower(vars)) {
    #var <- "tx"
    print(var)
    
    cc <- data.frame(matrix(nrow = 366, ncol = 0))
    cc$j <- 1:366    
    for (sens in rsens) {
      #sens <- +0.2
      #print(sens)
      odf <- data.d[,var]
      ndf <- data.frame(data.d[,var])
      colnames(ndf) <- var
      ndf <- ndf*(1+sens)
      nndf <- ndf[[var]]
      if (var == "rh") ndf[ndf>100] <- NA
      
      ndf <- cbind(dates, ndf,data.d[,!(colnames(data.d) %in% c(var))])
      colnames(ndf)[1] <- "Date"
      rpet <- PET.PM(ndf, lat = lat, z = z, checkTdiff = T)$PET
      sensr <- data.frame(data.d$j, ((rpet-pet0)/(ndf[[var]]-odf))*(odf/pet0))
      colnames(sensr) <- c("j","sens")
      sensr <- aggregate(sensr$sens, list(sensr$j), FUN=mean)
      if (length(which(is.na(sensr$x)))>0) sensr$x <- na.approx(sensr$x)
      cc <- cbind(cc,sensr$x)
    }
    colnames(cc) <- c(colnames(cc)[1],rsens)
    cc$j <- NULL
    cc.df <- cbind(cc.df,rowMeans(cc, na.rm=T))
  }
  
  colnames(cc.df) <- tolower(vars)
  #for (v in colnames(cc.df)) cc.df[,v] <- predict(loess(cc.df[,v] ~ as.numeric(1:366),span=0.04))
  df.sens[[s]] <- cc.df
}

miny <- -0.75
maxy <- 1

colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols <- c("red","red4","blue", "yellow2","green")

plot.loc <- function(df, miny, maxy, cols,s, cexf, n) {
  plot(df$tx, type="l", ylim=c(miny, maxy), 
       xlab=ifelse(n>6,"Julian Day of Year",""),
       ylab=ifelse(((n==1) || (n==4) || (n==7)),"Daily sensitivity coefficient",""),
       col=cols[1], main=sub('^(\\w?)', '\\U\\1', s, perl=T), lwd=2,
       cex.lab = cexf, cex.axis = cexf, cex.main = cexf*1.25)
  abline(v=156, lty=2)
  abline(v=304, lty=2)
  text(70,0.95, "Dry", font=3, cex=cexf*1.1)
  text(225,0.95, "Wet",font=3, cex=cexf*1.1)
  text(340,0.95, "Dry",font=3, cex=cexf*1.1)
  coli <- 1
  for (v in names(df)[-1]) {
    print(v)
    coli <- coli+1
    lines(df[v], col = cols[coli], lwd = 2)
  }
  
}

res <- 400
factor <- (res/72)

png(filename = "graphs/daily_sens.png",
    width = 1200 * factor, height = 800 * factor, res = res)

#dev.off()
layout(matrix(c(1:9,10,10,10), ncol=3, byrow=TRUE), heights=c(3, 3, 3, 1))
par(mar = rep(5, 4))

plotn <- 0
cexf <- 1.5
for (s in stations$Name) {
  print(s)
  plotn <- plotn + 1
  plot.loc(df.sens[[s]] , miny, maxy, cols,s, cexf, plotn)
}

par(mai=rep(0,4))
plot.new()
legend(x="center",ncol = 5, cex = cexf*1.25, legend = tolower(vars),
       lty = 1, col = cols, bty = "n", lwd=2,
       title="Variables",)

dev.off()
