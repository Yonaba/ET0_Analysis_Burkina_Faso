Sys.setenv(TZ="UTC")

Gsc <- 0.0820

#################################### Reference Penman-Monteith FAO ###############

# Allen et al. (1998) - Reference Penman-Monteith
# Source: FAO-56 
PET.PM <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1, checkTdiff = F) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  gamma <- 0.000665*patm

  sigma <- 4.903e-9
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tdiff <- cd$tx - cd$tn
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$year <- as.numeric(format(cd$Date, "%Y"))
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  p$DT <- p$ssv/(p$ssv + gamma*(1+0.34*cd$ws))
  p$PT <- gamma/(p$ssv + gamma*(1+0.34*cd$ws))
  p$TT <- (900/(p$tm + 273.16))*cd$ws
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  p$Rng <- 0.408 * p$Rn
  p$etrad <- p$DT * p$Rng
  p$etws <- p$PT * p$TT * (p$es - p$ea)
  p$PET <- (p$etrad + p$etws)
  if (checkTdiff) p$PET <- with(p, ifelse(p$tdiff < 0, NA, p$PET))
  return (p)
}

#################################### Mass Transfer Equations ####################

#Horton (1919)
# SOurce: https://doi.org/10.1016/j.agwat.2020.106043
PET.MTF.Horton <- function (cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- 0.4*(2-exp(-2*cd$ws))*(p$es-p$ea)
  return (p)
}

#Rohwer (1962)
# SOurce: 10.2166/nh.2021.128
PET.MTF.Rohwer <- function (cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- 0.44*(1+0.27*cd$ws)*(p$es-p$ea)
  return (p)
}

#Mahringer (1970)
# SOurce: 10.1007/s10584-010-9869-7
PET.MTF.Mahringer <- function (cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- 0.15072*sqrt(cd$ws*3.6)*(p$es-p$ea)
  return (p)
}

#Romanenko et al. (1961) modified by Oudin et al. (2005)
# Source : 10.1016/j.jhydrol.2004.08.026
PET.MTF.Romanenko <- function (cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- 4.5*((1+(p$tm/25))^2)*(1-p$ea/p$es)
  return (p)
}

#################################### Temperature based Equations ############

#Baier and Robertson (1965)
# SOurce: 10.1016/j.jhydrol.2015.06.057
PET.TMP.BaierRobertson <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- 0.09*(1.67021*cd$tx+1.68085*(cd$tx-cd$tn)+1.159575*p$Ra-57.3404)
  return (p)
}

# Linacre (1977)
# Source: 10.3390/w9100734 
PET.TMP.Linacre <- function (cd, lat, z) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$td <- (116.91 + 237.3*log(p$ea))/(16.78-log(p$ea))
  p$PET <- ((500*(p$tm + 0.006*z)/(100-lat)) + 15*(p$tm - p$td))/(80-p$tm)
  return (p)
}

#Ahooghalandari et al (2016)
# SOurce: https://doi.org/10.1016/j.agwat.2020.106043
PET.TMP.Ahooghalandari <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws))) 
  
  p$PET <- 0.369*p$Ra + 0.139*cd$tx*(1-cd$rh/100)-1.95
  return (p)
}

#Lobit et al (2017)
# SOurce: https://doi.org/10.1016/j.agwat.2020.106043
PET.TMP.Lobit <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  
  p$PET <- 0.1555*p$Ra*(0.00428*p$tm+0.09967)*(cd$tx-cd$tn)^0.5
  return (p)
}

#Tang et al (2019)
# SOurce: https://doi.org/10.1016/j.agwat.2020.106043
PET.TMP.Tang <- function (cd, lat, z) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  
  p$PET <- 1e-4*p$Ra*(p$tm+36.6)*(7+0.003*z)*(cd$tx-cd$tn)^0.5
  return (p)
}

#McCloud (1955)
# SOurce: http://dx.doi.org/10.1016/j.jhydrol.2015.06.057
PET.TMP.McCloud <- function (cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- 0.254*(1.07^(1.8*p$tm))
  return (p)
}

#Hamon (1961)
# Source: https://ascelibrary.org/doi/abs/10.1061/JYCEAJ.0000599
PET.TMP.Hamon1 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j")) 
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$N <- (24/pi)*p$ws
  
  p$PET <- 13.97*(p$N/12)^2*(p$eT/100)*10
  return (p)
}

#Hamon (1963)
# Source: https://doi.org/10.1016/j.jhydrol.2014.12.006
PET.TMP.Hamon2 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j")) 
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$N <- (24/pi)*p$ws
  
  p$PET <- 0.1651*(p$N/12)*(216.7*p$eT/(p$tm+273.3))*10
  return (p)
}

#Oudin et al. (2005) based on Hamon (1961)
# Source: http://dx.doi.org/10.1016/j.jhydrol.2013.09.005
PET.TMP.Oudin1 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j")) 
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$N <- (24/pi)*p$ws
  
  p$PET <- ((p$N/12)^2)*exp(p$tm/16)
  return (p)
}

#Oudin et al. (2005)
# Source: 10.1016/j.jhydrol.2004.08.026
PET.TMP.Oudin2 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$theta <- 0.4093*sin(((2*pi*p$J)/365)-1.405)
  p$cosGz <- cos(phi-p$theta)
  p$cosGz <- with(p, ifelse(p$cosGz <= 0.001,0.001,p$cosGz))
  p$Gz <- acos(p$cosGz)
  p$cosW <- 1-p$cosGz/(cos(phi)*cos(p$theta))
  p$cosW <- with(p, ifelse(p$cosW > 1, 1, p$cosW))
  p$cosW <- with(p, ifelse(p$cosW < -1, -1, p$cosW))
  p$W <- acos(p$cosW)
  p$eta <- 1+cos(2*pi*p$J/365)/30
  p$cosZ <- p$cosGz + cos(phi)*cos(p$theta)*((sin(p$W)/p$W)-1)
  p$Ge <- 446 * p$cosZ * p$W * p$eta
  p$PET <- p$Ge*((p$tm + 5)/28.5/100)
  p$PET <- with(p, ifelse(p$PET < 0, 0, p$PET))
  return (p)
}

#Hargreaves-Samani(1985)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
PET.TMP.HargreavesSamani <- function (cd, lat) {
  phi <- lat*pi/180
  krs <- 0.17 
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- 0.0135*krs*(p$Ra/2.45)*((cd$tx-cd$tn)^0.5)*(p$tm+17.8)
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p)
}

# Penman-Monteith Temperature
# Source: 10.1007/s00704-016-1996-2 
PET.TMP.PMT1 <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  gamma <- 0.000665*patm
  k1 <- 0.16
  k2 <- 0.5
  lambda <- 2.45
  
  sigma <- 4.903e-9
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$td <- (116.91 + 237.3*log(p$ea))/(16.78-log(p$ea))
  p$eTd <- 0.6108*exp((17.27*p$td)/(p$td + 237.3))
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$denom <- p$ssv+gamma*(1+0.34*cd$ws)

  p$PMTaero <- (gamma*((900*cd$ws)/(p$tm+273.16))*(p$es-p$eTd))/p$denom
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$Rs <- k1*(cd$tx-cd$tn)^k2
  p$Aterm <- 0.77*p$Ra*k1*(cd$tx-cd$tn)^k2
  
  p$Rs0 <- 0.75+2e-5*z
  p$Rs_Rs0 <- p$Rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  
  p$Bterm <- (1.35*p$Rs_Rs0-0.35)
  p$Cterm <- 0.34-0.14*sqrt(p$eTd)
  p$Dterm <- 0.5*((cd$tx+273.15)^4+(cd$tn+273.15)^4)
   
  p$PMTrad <- ((1/lambda)*p$ssv*(p$Aterm-(p$Bterm*sigma*p$Cterm*p$Dterm)))/p$denom
   
  p$PET <- (p$PMTrad + p$PMTaero)
  return (p)
}

# Valiantzas (2013)
# Source: http://dx.doi.org/10.1016/j.jhydrol.2013.09.005
PET.TMP.Valiantzas1 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$ea <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))*(cd$rh/100)
  p$tdew <- (116.91 + 237.3*log(p$ea))/(16.78-log(p$ea))
  p$PET <- 0.00668*p$Ra*sqrt((p$tm+9.5)*(cd$tx-p$tdew))-0.0696*(cd$tx-p$tdew)-0.024*(p$tm+20)*(1-cd$rh/100)-
    0.00455*p$Ra*(cd$tx-p$tdew)^0.5+0.0984*(p$tm+17)*(1.03+0.00055*(cd$tx-cd$tn)^2-cd$rh/100)
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p)
}

#Valiantzas (2018)
# Source: 10.1016/j.agwat.2018.06.028
PET.TMP.Valiantzas2 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- 0.0118*((1-cd$rh/100)^0.2)*((cd$tx-cd$tn)^0.3)*(p$Ra*((p$tm+10)^0.5)-40)+
    0.1*(p$tm+20)*(1-cd$rh/100)
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p)
}

#################################### Radiation based Equations ############

#Ritchie (1972) in Jones and Ritchie (1990)
# Source: 10.1007/s00271-011-0295-z
PET.RAD.Ritchie <- function(cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$alpha <- with(cd,
                  ifelse(cd$tx <= 5,0.01*exp(0.18*(cd$tx+20)),
                         ifelse(cd$tx <= 35,1,
                                1.1+0.05*(cd$tx-35))))
  p$PET <- (0.00387*cd$rs*(0.6*p$tm + 0.4*cd$tn + 29))*p$alpha
  return (p)
}

#Jensen-Haise (1963)
# Source: 10.2166/nh.2021.128
PET.RAD.JensenHaise <- function(cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  lambda <- 2.45
  p$PET <- 0.025 * (p$tm + 3) * cd$rs/lambda
  return (p)
}

#Alexandris et al. (2006)
#Source: 10.1016/j.agwat.2005.08.001
PET.RAD.Alexandris <- function(cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  a1 <- 0.6416
  a2 <- -0.00784
  a3 <- 0.372
  a4 <- -0.00264
  b1 <- -0.0033
  b2 <- 0.00812
  b3 <- 0.101
  b4 <- 0.00584
  m1 <- 0.057
  m2 <- 0.227
  m3 <- 0.643
  m4 <- 0.0124
  p$C1 <- a1+a2*cd$rh + a3*cd$rs + a4*cd$rs*cd$rh
  p$C2 <- b1 + b2*p$tm + b3*cd$rs + b4*cd$rs*p$tm
  p$PET <- m1 + m2*p$C2 + m3*p$C1 + m4*p$C1*p$C2
  return (p)
}

#Doorenbos and Pruit (1977)
# Source: 10.2166/nh.2021.128
PET.RAD.DoorenbosPruitt <- function(cd, z) {
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  gamma <- 0.000665*patm
  sigma <- 4.903e-9
  lambda <- 2.45
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  p$alpha <- 1.066-0.13*1e-2*cd$rh+0.045*cd$ws-0.2e-3*cd$rh*cd$ws-0.135e-4*cd$rh^2-0.11e-2*cd$ws^2
  p$PET <- p$alpha*(p$ssv/(p$ssv+gamma))*(cd$rs/lambda)-0.3
  
  return (p)
}

#Hansen (1984)
# Source: https://doi.org/10.2166/nh.1984.0017
PET.RAD.Hansen <- function(cd, z) {
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  gamma <- 0.000665*patm
  sigma <- 4.903e-9
  lambda <- 2.45
  d <- 0
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  p$PET <- 0.7*(p$ssv/(p$ssv+gamma))*(cd$rs/lambda)-d
  
  return (p)
}

#Trajkovic and Stojnić (2007)
#Source: https://doi.org/10.1016/j.agwat.2020.106043
PET.RAD.TrajkovicStojnic <- function(cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- 0.013*(0.8383-0.0313*cd$ws^2+0.1706*cd$ws)*(23.88*cd$rs+50)*(p$tm/(p$tm+15))
  return (p)
}

# McGuiness-Bordne (1972)
# Source: 10.1111/j.1752-1688.1996.tb04044.x
PET.RAD.McGuinessBordne <- function (cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- (0.00597*p$tm + 0.0838)*cd$rs
  return (p)
}

# Abtew (1996)
# Source: 10.1111/j.1752-1688.1996.tb04044.x
PET.RAD.Abtew1 <- function (cd) {
  lambda <- 2.45
  k <- 0.53
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$PET <- k * cd$rs/lambda
  return (p)
}

# Abtew (1996)
# Source: 10.1111/j.1752-1688.1996.tb04044.x
PET.RAD.Abtew2 <- function (cd) {
  lambda <- 2.45
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$PET <- 0.01786 * cd$rs * cd$tx/lambda
  return (p)
}

# Hargreaves and Allen (2003)
# Source: https://doi.org/10.1061/(ASCE)0733-9437(2003)129:1(53)
PET.RAD.HargreavesAllen <- function (cd) {
  lambda <- 2.45
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$PET <- (0.0135*cd$tn + 0.2403)*(cd$rs/lambda)
  return (p)
}

#Irmak et al. (2003)
#Source: https://doi.org/10.1061/(ASCE)0733-9437(2003)129:5(336)
PET.RAD.Irmak <- function(cd) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- (-0.611 + 0.149*cd$rs +0.079*p$tm)
  return (p)
}

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000502
PET.RAD.Valiantzas1 <- function (cd, lat) {
  phi <- lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- 0.0393*(p$tm+9.5)^0.5*cd$rs-0.19*(cd$rs^0.6)*(phi^0.15)+0.0061*(p$tm+20)*(1.12*p$tm-cd$tn-2)^0.7
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p)
}

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000520
PET.RAD.Valiantzas2 <- function (cd, lat) {
  phi <- lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- 0.0393*(p$tm+9.5)^0.5*cd$rs-0.19*(cd$rs^0.6)*(phi^0.15)+0.0037*(p$tm+20)*(p$tm-cd$tn)^0.7
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p)
}

# Priestley-Taylor (1972)
# Source: 10.1016/j.jhydrol.2004.08.026
PET.RAD.PriestleyTaylor <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  gamma <- 0.000665*patm
  alpha <- 1.26
  sigma <- 4.903e-9
  lambda <- 2.45
  rhoa <- 1
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$PET <- alpha*(p$ssv/(p$ssv+gamma))*(p$Rn/(rhoa*lambda))
  
  return (p)
}


################################## Combinatory methods ###################################

# Penman (1956)
# Source: 10.1016/j.jhydrol.2004.08.026 / https://emf-creaf.github.io/meteolandbook/potentialevapotranspiration.html
PET.COMB.Penman2 <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  rho <- 1
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$aw <- 0.4+1.4*exp(-((p$J-173)/52)^2)
  p$bw <- 0.605+0.345*exp(-((p$J-243)/80)^2)
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/((lambda*rho)*(p$ssv + gamma)))
  p$f2 <- (gamma/((lambda*rho)*(p$ssv + gamma)))
  p$W <- 2.626+1.381*cd$ws
  p$PET <- p$f1*p$Rn + p$f2*p$W*(p$es-p$ea)
  
  return (p)
}

# Kimberly-Penman dans Wright (1982)
# Source: 10.1016/j.jhydrol.2004.08.026
PET.COMB.KimberlyPenman <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  rhoa <- 1
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/((lambda*rhoa)*(p$ssv + gamma)))
  p$f2 <- (gamma/((lambda*rhoa)*(p$ssv + gamma)))
  p$aw <- 0.4+0.14*exp(-((p$J-173)/58)^2)
  p$bw <- 0.605+0.345*exp(-((p$J-243)/80)^2)  
  p$Wf <- p$aw + p$bw*cd$ws
  p$PET <- p$f1*p$Rn + p$f2*p$Wf*(p$es-p$ea)
  
  return (p)
}

# George et al (1985)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
PET.COMB.George <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$aw <- 0.4+1.4*exp(-((p$J-173)/52)^2)
  p$bw <- 0.605+0.345*exp(-((p$J-243)/80)^2)
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/(p$ssv + gamma))
  p$f2 <- (gamma/(p$ssv + gamma))
  p$Wf <- p$aw + p$bw*cd$ws
  p$PET <- p$f1*p$Rn + 0.268*p$f2*p$Wf*(p$es-p$ea)
  
  return (p)
}

# Thom and Oliver (1977)
# Source: 10.1016/j.jhydrol.2004.08.026
PET.COMB.ThomOliver <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  rhoa <- 1
  rs <- 69
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))
  
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  ws <- cd$ws
  ws[ws==0] <- min(min(ws), 0.1)
  p$ra <- 208/ws
  
  p$f1 <- (p$ssv/((lambda*rhoa)*(p$ssv + gamma*(1+rs/p$ra))))
  p$f2 <- (gamma/((lambda*rhoa)*(p$ssv + gamma*(1+rs/p$ra))))
  
  p$Wf <- 2.6*(1+0.536*cd$ws)
  p$PET <- p$f1*p$Rn + 2.5*p$f2*p$Wf*(p$es-p$ea)
  
  return (p)
}

# VanBavel (1966)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
PET.COMB.VanBavel <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  rhoa <- 997
  k <- 0.41
  z1 <- 0.12
  z0 <- 0.123*z1
  d <- 0.67*z1
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))
  
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/(p$ssv + gamma))
  p$f2 <- (gamma/(p$ssv + gamma))
  p$PET <- p$f1*p$Rn + p$f2*(0.622*rhoa*k^2/patm)*((cd$ws*(p$es-p$ea))/((log(z1-d)/z0)^2))
  
  return (p)
}

# Rijtema (1966)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
PET.COMB.Rijtema <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  hc <- 0.12
  r <- 0.123*hc
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$f <- (p$ssv/(p$ssv + gamma))
  p$PET <- p$f*p$Rn + p$f*gamma*r*(cd$ws^0.75)*(p$es-p$ea)
  
  return (p)
}

# Penman (1963)
#Source: http://dx.doi.org/10.1016/j.ejrh.2015.02.002
PET.COMB.Penman3 <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  rhoa <- 1
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$aw <- 0.4+1.4*exp(-((p$J-173)/52)^2)
  p$bw <- 0.605+0.345*exp(-((p$J-243)/80)^2)
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/((lambda*rhoa)*(p$ssv + gamma)))
  p$f2 <- (gamma/((lambda*rhoa)*(p$ssv + gamma)))
  p$W <- 6.43*(1+0.536*cd$ws)
  p$PET <- p$f1*p$Rn + p$f2*p$W*(p$es-p$ea)
  
  return (p)
}

# Matt-Shuttleworth (2006)
#Source: 10.13031/2013.29217
PET.COMB.MattShuttleworth <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  roua <- 1.2
  ca <- 0.001013
  rs <- 70
  hc <- 0.12
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$N <- (24/pi)*p$ws
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  p$rclim <- 86400 * roua * ca * (p$es-p$ea)/(p$ssv * p$Rn)
  p$rclim <- with(p, ifelse(p$rclim == 0, 0.1, p$rclim))
  p$ws <- cd$ws
  p$ws <- with(p, ifelse(p$ws == 0, 0.1, p$ws))
  p$vpd <- (302 * (p$ssv + gamma) + 70 * gamma * p$ws)/(208 *(p$ssv + gamma) + 70 * gamma * p$ws) +
    1/p$rclim * ((302 * (p$ssv + gamma) + 70 * gamma * p$ws)/(208 * (p$ssv + gamma) + 70 * gamma * p$ws) * (208/p$ws) - (302/p$ws))
 
  r_c50 <- 1/((0.41)^2) * log((50 - 0.67 * hc)/(hc)) * 
    log((50 - 0.67 * hc)/(0.0123 * hc)) * log((2 - 0.08)/0.0148)/log((50 - 0.08)/0.0148)
  
  p$PET <-  1/lambda * (p$ssv * p$Rn + (roua * ca * p$ws * (p$es - p$ea))/r_c50 * p$vpd)/(p$ssv + gamma * (1 + rs * p$ws/r_c50))
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p)
}

#Irmak et al. (2003)
#Source: 10.1061/~ASCE!0733-9437~2003!129:5~336!
PET.COMB.Irmak <- function(cd, lat, z, albedo = 0.23, rs_rso_min = 0.33, rs_rso_max = 1) {
  phi <- lat*pi/180
  patm <- 101.3*(((293-0.0065*z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  
  p$Rs0 <- (0.75+0.00002*z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rn <- p$Rns - p$Rnl
  
  p$PET <- 0.489 + 0.289*p$Rn + 0.023*p$tm
  
  return (p)
}

#Valiantzas (2006)
#Source: 10.1016/j.jhydrol.2006.06.012
PET.COMB.Valiantzas1 <- function (cd, lat, z) {
  phi <- lat*pi/180
  alpha <- 0.23
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- 0.051*(1-alpha)*((p$tm+9.5)^0.5)*cd$rs-2.4*((cd$rs/p$Ra)^2)+0.00012*z +
    0.048*(p$tm+20)*(1-cd$rh/100)*(0.5+0.536*cd$ws)
  return (p)
}


#Valiantzas (2006)
#Source: 10.1016/j.jhydrol.2006.06.012
PET.COMB.Valiantzas2 <- function (cd, lat) {
  phi <- lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- 0.038*((p$tm+9.5)^0.5)*cd$rs-2.4*((cd$rs/p$Ra)^2)+0.075*(p$tm+20)*(1-cd$rh/100)
  return (p)
}

#Valiantzas (2013)
# Source: https://doi.org/10.1061/(ASCE)HE.1943-5584.0000590
PET.COMB.Valiantzas3 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  p$Waero <- with(cd, ifelse(cd$rh>65,0.78,1.067))
  p$PET <- 0.0393*(cd$rs)*((p$tm+9.5)^0.5)-0.024*(p$tm+20)*(1-cd$rh/100)-2.4*((cd$rs/p$Ra)^2)+
    0.066*p$Waero*(p$tm+20)*(1-cd$rh/100)*cd$ws^0.6
  return (p)
}

#Valiantzas (2013)
# Source: https://doi.org/10.1061/(ASCE)HE.1943-5584.0000590
PET.COMB.Valiantzas4 <- function (cd, lat) {
  phi <- lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  p$Cu <- with(cd, ifelse(cd$rh > 65, 0.054, 0.083))
  p$PET <- 0.0393*(cd$rs)*((p$tm+9.5)^0.5)-2.4*((cd$rs/p$Ra)^2)+p$Cu*(p$tm+20)*(1-cd$rh/100)
  return (p)
}

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000520
PET.COMB.Valiantzas5 <- function (cd, lat) {
  phi <- lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- 0.0393*(cd$rs)*((p$tm+9.5)^0.5)-0.19*(cd$rs^0.6)*(phi^0.15)+0.048*(p$tm+20)*(1-cd$rh/100)*(cd$ws^0.7)
  return (p)
}

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000520
PET.COMB.Valiantzas6 <- function (cd, lat) {
  phi <- lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- 0.0393*(cd$rs)*((p$tm+9.5)^0.5)-0.19*(cd$rs^0.6)*(phi^0.15)+0.078*(p$tm+20)*(1-cd$rh/100)
  return (p)
}

