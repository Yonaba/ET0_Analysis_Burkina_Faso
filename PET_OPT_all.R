Sys.setenv(TZ="UTC")

#################### Optimization Case 1: Tx, Tn are known ##############

#Tang et al (2019)
# SOurce: https://doi.org/10.1016/j.agwat.2020.106043
#a=2.089308 , b=1.451514 , c=6.312104 , d=0.006747, n1=0.345908
OPET.TMP.Tang <- function (cd, a=2.089308 , b=1.451514 , c=6.312104 , d=0.006747, n1=0.345908) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  
  p$PET <- 1e-4*p$Ra*(a*p$tm+b)*(c+d*cd$z)*(cd$tx-cd$tn)^n1
  return (p$PET)
}

# Penman-Monteith Temperature
# Source: 10.1007/s00704-016-1996-2 
#c = 0.7494 , k1 = 0.1884, k2 = 0.4409 , d = 1.2236, e = 0.3220
#c = 1.2190 , k1 = 0.2249, k2 = 0.2484 , d = 2.7959 , e = 0.2420 (with long-term mean ws)
OPET.TMP.PMT1 <- function(cd, c = 0.7494 , k1 = 0.1884, k2 = 0.4409 , d = 1.2236, e = 0.3220) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  gamma <- 0.000665*patm
  
  lambda <- 2.45
  
  albedo <- 0.23 
  rs_rso_min <- 0.33
  rs_rso_max <- 1
  
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
  
  p$denom <- p$ssv+gamma*(d+e*cd$ws)
  
  p$PMTaero <- (gamma*((900*cd$ws)/(p$tm+273.16))*(p$es-p$eTd))/p$denom
  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  p$Rs <- k1*(cd$tx-cd$tn)^k2
  p$Aterm <- c*p$Ra*k1*(cd$tx-cd$tn)^k2
  
  p$Rs0 <- 0.75+2e-5*cd$z
  p$Rs_Rs0 <- p$Rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  
  p$Bterm <- (1.35*p$Rs_Rs0-0.35)
  p$Cterm <- 0.34-0.14*sqrt(p$eTd)
  p$Dterm <- 0.5*((cd$tx+273.15)^4+(cd$tn+273.15)^4)
  
  p$PMTrad <- ((1/lambda)*p$ssv*(p$Aterm-(p$Bterm*sigma*p$Cterm*p$Dterm)))/p$denom
  
  p$PET <- (p$PMTrad + p$PMTaero)
  return (p$PET)
}

#Oudin et al. (2005)
# Source: 10.1016/j.jhydrol.2004.08.026
#a=-33.497    , b=18.255   , c= -5.474 
OPET.TMP.Oudin2 <- function (cd, a=-33.497    , b=18.255   , c= -5.474 ) {
  phi <- cd$lat*pi/180
  
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
  p$PET <- p$Ge*((p$tm + a)/b/100)-c
  p$PET <- with(p, ifelse(p$PET < 0, 0, p$PET))
  return (p$PET)
}

#Oudin et al. (2005) based on Hamon (1961)
# Source: http://dx.doi.org/10.1016/j.jhydrol.2013.09.005
#a=1.4777   , b=23.4793   , c= -0.2555 , n1= -0.6578
OPET.TMP.Oudin1 <- function (cd,a=1.4777   , b=23.4793   , c= -0.2555 , n1= -0.6578) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j")) 
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$N <- (24/pi)*p$ws
  
  p$PET <- a*((p$N/12)^n1)*exp(p$tm/b)+c
  return (p$PET)
}

#McCloud (1955)
# SOurce: http://dx.doi.org/10.1016/j.jhydrol.2015.06.057
#a=0.0453, b=2.7942
OPET.TMP.McCloud <- function (cd,a=0.0453, b=2.7942) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- a*(1.07^(1.8*p$tm))+b
  return (p$PET)
}

#Lobit et al (2017)
# SOurce: https://doi.org/10.1016/j.agwat.2020.106043
#a= 0.01145, b=0.02544   , n1= 0.32483
OPET.TMP.Lobit <- function (cd, a= 0.01145, b=0.02544   , n1= 0.32483) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  
  p$PET <- 0.1555*p$Ra*(a*p$tm+b)*(cd$tx-cd$tn)^n1
  return (p$PET)
}

# Linacre (1977)
# Source: 10.3390/w9100734 
#a=1.02017   , b=0.02846, c=105.96368, d=2.10185  , e=81.47602
OPET.TMP.Linacre <- function (cd,a=1.02017   , b=0.02846, c=105.96368, d=2.10185  , e=81.47602) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es
  
  p$td <- (116.91 + 237.3*log(p$ea))/(16.78-log(p$ea))
  p$PET <- ((500*(a*p$tm + b*cd$z)/(c-cd$lat)) + d*(p$tm - p$td))/(e-p$tm)
  return (p$PET)
}

#Hargreaves-Samani(1985)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
#krs=0.3231 , n1=0.3248 , b=2.2219 
OPET.TMP.HargreavesSamani <- function (cd, krs=0.3231 , n1=0.3248 , b=2.2219) {
  phi <- cd$lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- 0.0135*krs*(p$Ra/2.45)*((cd$tx-cd$tn)^n1)*(p$tm+b)
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p$PET)
}

#Hamon (1963)
# Source: http://dx.doi.org/10.1016/j.jhydrol.2013.09.005
#c= 0.1171, d=1.0471
OPET.TMP.Hamon2 <- function (cd,c= 0.1171, d=1.0471) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j")) 
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$N <- (24/pi)*p$ws
  
  p$PET <- c*(p$N/12)*(216.7*p$eT/(p$tm+273.3))*10+d
  return (p$PET)
}

#Hamon (1961)
# Source: http://dx.doi.org/10.1016/j.jhydrol.2013.09.005
#a= 9.2625 , n1= -0.2708,b=0.7057
OPET.TMP.Hamon1 <- function (cd,a= 9.2625 , n1= -0.2708,b=0.7057) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j")) 
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$N <- (24/pi)*p$ws
  
  p$PET <- a*(p$N/12)^n1*(p$eT/100)*10+b
  return (p$PET)
}

#Baier and Robertson (1965)
# SOurce: 10.1016/j.jhydrol.2015.06.057
#a=0.19443, b=-0.02723, c= 0.03592, d=3.34133
OPET.TMP.BaierRobertson <- function (cd,a=0.19443, b=-0.02723, c= 0.03592, d=3.34133) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  p$PET <- a*cd$tx+b*(cd$tx-cd$tn)+c*p$Ra-d
  return (p$PET)
}

#################### Optimization Case 2: Tx, Tn, Rh are known ##############

#Romanenko et al. (1961) modified by Oudin et al. (2005)
# Source : 10.1016/j.jhydrol.2004.08.026
#a=-0.7204 , b=21.4283  , c= 3.4682
OPET.MTF.Romanenko <- function (cd,a=-0.7204 , b=21.4283  , c= 3.4682) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- 4.5*((a+(p$tm/b))^2)*(1-p$ea/p$es)+c
  return (p$PET)
}

#Rohwer (1962)
# SOurce: 10.2166/nh.2021.128
#a=0.09061 , b=0.99359   , c= 3.31195 
#a=0.4259  , b=0.7625    , c= 3.2523 (long term ws)
OPET.MTF.Rohwer <- function (cd, a=0.09061 , b=0.99359   , c= 3.31195 ) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- 0.44*(a+b*cd$ws)*(p$es-p$ea)+c
  return (p$PET)
}

#Mahringer (1970)
# SOurce: 10.1007/s10584-010-9869-7
#a=0.3233  , b=3.0914
#a=0.2882  , b=3.2202 (with long term ws)
OPET.MTF.Mahringer <- function (cd,a=0.3233  , b=3.0914) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- a*sqrt(cd$ws*3.6)*(p$es-p$ea)+b
  return (p$PET)
}

#Horton (1919)
# SOurce: https://doi.org/10.1016/j.agwat.2020.106043
#a=4.7179   , b=0.8611  
#a=4.6649    , b=0.9607 
OPET.MTF.Horton <- function (cd,a=4.7179   , b=0.8611 ) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$eTn <- 0.6108*exp((17.27*cd$tn)/(cd$tn + 237.3))
  p$eTx <- 0.6108*exp((17.27*cd$tx)/(cd$tx + 237.3))
  p$es <- 0.5* (p$eTn + p$eTx)
  p$ea <- (cd$rh/100)*p$es  
  p$PET <- 0.4*(a-exp(-b*cd$ws))*(p$es-p$ea)
  return (p$PET)
}

# Valiantzas (2013)
# Source: http://dx.doi.org/10.1016/j.jhydrol.2013.09.005
#a=0.017254  , n1=0.530887 , b=-0.481758  , c=0.260185, d=0.016704  , n2=1.059265 ,f=-0.326022 , g=-0.002027, n3=1.812105
OPET.TMP.Valiantzas1 <- function (cd,a=0.017254  , n1=0.530887 , b=-0.481758  , c=0.260185, d=0.016704  , n2=1.059265 ,f=-0.326022 , g=-0.002027, n3=1.812105) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$ea <- 0.6108*exp((17.27*p$tm)/(p$tm+237.3))*(cd$rh/100)
  p$tdew <- (116.91 + 237.3*log(p$ea))/(16.78-log(p$ea))
  p$PET <- a*p$Ra*((p$tm+9.5)*(cd$tx-p$tdew))^n1-b*(cd$tx-p$tdew)-c*(p$tm+20)*(1-cd$rh/100)-
    d*p$Ra*(cd$tx-p$tdew)^n2+0.0984*(p$tm+17)*(f+g*(cd$tx-cd$tn)^n3-cd$rh/100)
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p$PET)
}

#Ahooghalandari et al (2016)
# Source: https://doi.org/10.1016/j.agwat.2020.10604
#c=0.19083 , d=0.08685 , e=3.85402 
OPET.TMP.Ahooghalandari <- function (cd,c=0.19083 , d=0.08685 , e=3.85402) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws))) 
  
  p$PET <- c*p$Ra + d*cd$tx*(1-cd$rh/100)-e
  return (p$PET)
}

#Valiantzas (2018)
# Source: 10.1016/j.agwat.2018.06.028
#a=0.01474  , n1=0.30133  , n2=0.02954  , n3=0.67394 , b=-0.01046 
OPET.TMP.Valiantzas2 <- function (cd, a=0.01474  , n1=0.30133  , n2=0.02954  , n3=0.67394 , b=-0.01046) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- a*((1-cd$rh/100)^n1)*((cd$tx-cd$tn)^n2)*(p$Ra*((p$tm+10)^n3)-40)+
    b*(p$tm+20)*(1-cd$rh/100)
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p$PET)
}

#################### Optimization Case 3: Tx, Tn, Rs are known ##############

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000520
#c=0.0029685  , n1=1.1859811  , d=0.0007709  , n2=2.5852572  , n3=0.5187935  ,e=0.0150364  ,n4=0.1907575
OPET.RAD.Valiantzas2 <- function (cd,c=0.0029685  , n1=1.1859811  , d=0.0007709  , n2=2.5852572  , n3=0.5187935  ,e=0.0150364  ,n4=0.1907575) {
  phi <- cd$lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- (c*(p$tm+9.5)^n1)*cd$rs-d*(cd$rs^n2)*(phi^n3)+e*(p$tm+20)*(p$tm-cd$tn)^n4
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p$PET)
}

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000502
#c=0.0030255 , n1=1.1738433 , d=0.0007019 , n2=2.6310437 , n3=0.5999451 ,e=0.0140295 ,n4=0.2276241
OPET.RAD.Valiantzas1 <- function (cd,c=0.0030255 , n1=1.1738433 , d=0.0007019 , n2=2.6310437 , n3=0.5999451 ,e=0.0140295 ,n4=0.2276241) {
  phi <- cd$lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- (c*(p$tm+9.5)^n1)*cd$rs-d*(cd$rs^n2)*(phi^n3)+e*(p$tm+20)*(1.12*p$tm-cd$tn-2)^n4
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p$PET)
}

#Trajkovic and Stojnić (2007)
#Source: https://doi.org/10.1016/j.agwat.2020.106043
#a=1.142e+00, b=2.735e-03, c=3.618e-01, n1=3.296e+00, d=3.179e+01,e=2.488e+02,f=7.553e1
#a=1.7559, b=-2.4358, c=-1.0278, n1=0.7295  , d=75.7332 ,e=649.2709 ,f=509.0275 (long term ws)
OPET.RAD.TrajkovicStojnic <- function(cd,a=1.142e+00, b=2.735e-03, c=3.618e-01, n1=3.296e+00, d=3.179e+01,e=2.488e+02,f=7.553e1) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- 0.013*(a-b*cd$ws^n1+c*cd$ws)*(d*cd$rs+e)*(p$tm/(p$tm+f))
  return (p$PET)
}

#Ritchie (1972) in Jones and Ritchie (1990)
# Source: 10.1007/s00271-011-0295-z
#d=-3.082,e=  2.313,f= 93.551
OPET.RAD.Ritchie <- function(cd, d=-3.082,e=  2.313,f= 93.551) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$alpha <- with(cd,
                  ifelse(cd$tx <= 5,0.01*exp(0.18*(cd$tx+20)),
                         ifelse(cd$tx <= 35,1,
                                1.1+0.05*(cd$tx-35))))
  p$PET <- (0.00387*cd$rs*(d*p$tm + e*cd$tn + f))*p$alpha
  return (p$PET)
}

# McGuiness-Bordne (1972)
# Source: 10.1111/j.1752-1688.1996.tb04044.x
#a=0.00740 ,b=-0.05476  ,c=1.57594  
OPET.RAD.McGuinessBordne <- function (cd,a=0.00740 ,b=-0.05476  ,c=1.57594  ) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- (a*p$tm + b)*cd$rs+c
  return (p$PET)
}

#Jensen-Haise (1963)
# Source: 10.2166/nh.2021.128
#a=0.01342  , b= 0.21306
OPET.RAD.JensenHaise <- function(cd, a=0.01342  , b= 0.21306) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  lambda <- 2.45
  p$PET <- (a*p$tm + b) * cd$rs/lambda
  return (p$PET)
}

#Irmak et al. (2003)
#Source: 10.1061/~ASCE!0733-9437~2003!129:5~336!
#a=2.3093 , b=0.1496 , c=0.1416 
OPET.RAD.Irmak <- function(cd,a=2.3093 , b=0.1496 , c=0.1416) {
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- (-a + b*cd$rs +c*p$tm)
  return (p$PET)
}

# Hargreaves and Allen (2003)
# Source: https://doi.org/10.1061/(ASCE)0733-9437(2003)129:1(53)
#a=0.01246 , b=0.17202 ,c=1.14460
OPET.RAD.HargreavesAllen <- function (cd, a=0.01246 , b=0.17202 ,c=1.14460) {
  lambda <- 2.45
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$PET <- (a*cd$tn + b)*(cd$rs/lambda)+c
  return (p$PET)
}

#Hansen (1984)
# Source: https://doi.org/10.2166/nh.1984.0017
#a=0.6433  , b=-0.7850
OPET.RAD.Hansen <- function(cd, a=0.6433  , b=-0.7850) {
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  gamma <- 0.000665*patm
  sigma <- 4.903e-9
  lambda <- 2.45
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  p$PET <- a*(p$ssv/(p$ssv+gamma))*(cd$rs/lambda)-b
  
  return (p$PET)
}

# Abtew (1996)
# Source: 10.1111/j.1752-1688.1996.tb04044.x
#a=0.01191  , b=1.34665 
OPET.RAD.Abtew2 <- function (cd,a=0.01191  , b=1.34665) {
  lambda <- 2.45
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$PET <- a * cd$rs * cd$tx/lambda+b
  return (p$PET)
}

# Abtew (1996)
# Source: 10.1111/j.1752-1688.1996.tb04044.x
#a=0.5048   , b=0.7633
OPET.RAD.Abtew1 <- function (cd, a=0.5048   , b=0.7633) {
  lambda <- 2.45
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$PET <- a * cd$rs/lambda + b
  return (p$PET)
}

# Priestley-Taylor (1972)
# Source: 10.1016/j.jhydrol.2004.08.026
#a=1.1858 , b=-3.7281, c=0.5991 , d= 0.7215 , e=0.2199    , f=0.2342 
OPET.RAD.PriestleyTaylor <- function(cd, a=1.1858 , b=-3.7281, c=0.5991 , d= 0.7215 , e=0.2199    , f=0.2342 ) {
  
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  gamma <- 0.000665*patm
  sigma <- 4.903e-9
  lambda <- 2.45
  rhoa <- 1
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))
  p$Rn <- p$Rns - p$Rnl
  
  p$PET <- a*(p$ssv/(p$ssv+gamma))*(p$Rn/(rhoa*lambda))+b
  
  return (p$PET)
}

# Matt-Shuttleworth (2006)
#Source: 10.13031/2013.29217
#,ca=60.46031  , a=0.87549  ,b=0.97334,c=0.45757  , d= 0.11543, e=0.18257  , f=0.04897
#ca=64.12650  , a=0.67554    ,b=1.48875,c=0.48690  , d= 0.12514  , e=0.18092    , f=0.04599 
OPET.COMB.MattShuttleworth <- function(cd, ca=60.46031  , a=0.87549  ,b=0.97334,c=0.45757  , d= 0.11543, e=0.18257  , f=0.04897) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  roua <- 1.2
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))  
  p$Rn <- p$Rns - p$Rnl
  p$rclim <- 86400 * roua * ca * (p$es-p$ea)/(p$ssv * p$Rn)
  p$rclim <- with(p, ifelse(p$rclim == 0, 0.1, p$rclim))
  p$ws <- cd$ws
  p$ws <- with(p, ifelse(p$ws == 0, 0.1, p$ws))
  p$vpd <- (302 * (p$ssv + gamma) + 70 * gamma * p$ws)/(208 *(p$ssv + gamma) + 70 * gamma * p$ws) +
    1/p$rclim * ((302 * (p$ssv + gamma) + 70 * gamma * p$ws)/(208 * (p$ssv + gamma) + 70 * gamma * p$ws) * (208/p$ws) - (302/p$ws))
  
  r_c50 <- 1/((0.41)^2) * log((50 - 0.67 * hc)/(hc)) * 
    log((50 - 0.67 * hc)/(0.0123 * hc)) * log((2 - 0.08)/0.0148)/log((50 - 0.08)/0.0148)
  
  p$PET <-  1/lambda * (p$ssv * p$Rn + (roua * ca * p$ws * (p$es - p$ea))/r_c50 * p$vpd)/(p$ssv + gamma * (a + b*rs * p$ws/r_c50))
  p$PET <- with(p, ifelse(is.nan(p$PET),0,p$PET))
  return (p$PET)
}

#Irmak et al. (2003)
#Source: 10.1061/~ASCE!0733-9437~2003!129:5~336!
#a=20.6735  , b=0.3586  , g=0.6732  ,c=0.5710 , d= -2.2371, e=0.2141, f=0.2429
OPET.COMB.Irmak <- function(cd, a=20.6735  , b=0.3586  , g=0.6732  ,c=0.5710 , d= -2.2371, e=0.2141, f=0.2429) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))   
  p$Rn <- p$Rns - p$Rnl
  
  p$PET <- a + b*p$Rn + g*p$tm
  
  return (p$PET)
}

#################### Optimization Case 4: Tx, Tn, Rh and Rs are known ##############

#Doorenbos and Pruit (1977)
# Source: 10.2166/nh.2021.128
#a=0.178995, b=-0.00824066, c= 0.424952, d= 0.00446275, e= 5.60509e-05, f= 0.0112964, g= -0.888841
#a=2.187e-01,b=-6.880e-03,c=4.736e-01,d= 4.632e-03,e=4.413e-05,f=3.106e-02,g=-7.649e-01 (with long-term mean ws)
OPET.RAD.DoorenbosPruitt <- function(cd, a=0.178995, b=-0.00824066, c= 0.424952, d= 0.00446275, e= 5.60509e-05, f= 0.0112964, g= -0.888841) {
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  gamma <- 0.000665*patm
  sigma <- 4.903e-9
  lambda <- 2.45
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$eT <- 0.6108*exp((17.27*p$tm)/(p$tm + 237.3))
  p$ssv <- (4098*p$eT)/((p$tm + 237.3)^2)
  p$alpha <- a-b*cd$rh+c*cd$ws-d*cd$rh*cd$ws-e*cd$rh^2-f*cd$ws^2
  p$PET <- p$alpha*(p$ssv/(p$ssv+gamma))*(cd$rs/lambda)-g
  
  return (p$PET)
}

#Alexandris et al. (2006)
#Source: 10.1016/j.agwat.2005.08.001
#m1=0.37373  , m2=0.64075  , m3=0.30084, m4=-0.02272
OPET.RAD.Alexandris <- function(cd,m1=0.37373  , m2=0.64075  , m3=0.30084, m4=-0.02272) {
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
  return (p$PET)
}

# VanBavel (1966)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
#a=5.421e+04, b=-1.430e+00,c=4.392e-01, d= 1.345e-01, e=5.322e-02, f=4.981e-02
#aa=6.121e+04, b=-2.105e+00,c=4.265e-01, d= 1.418e-01, e=3.625e-02, f=3.697e-02  (mean ws)
OPET.COMB.VanBavel <- function(cd, a=5.421e+04, b=-1.430e+00,c=4.392e-01, d= 1.345e-01, e=5.322e-02, f=4.981e-02) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  rhoa <- 997
  k <- 0.41
  z1 <- 0.12
  z0 <- 0.123*z1
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1  
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))  
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/(p$ssv + gamma))
  p$f2 <- (gamma/(p$ssv + gamma))
  p$PET <- p$f1*p$Rn + p$f2*(a*rhoa*k^2/patm)*((cd$ws*(p$es-p$ea))/((log(z1-0.67*z1)/z0)^2))+b
  
  return (p$PET)
}

# Thom and Oliver (1977)
# Source: 10.1016/j.jhydrol.2004.08.026
#a=-0.3300 , b=7.2277  ,c=0.4582, d= 0.1361  , e=0.1898, f=0.0585
#a=0.87607  , b=6.96262   ,c=0.49119 , d= 0.13161   , e=0.18817 , f=0.04697 (with long term mean ws)
OPET.COMB.ThomOliver <- function(cd, a=-0.3300 , b=7.2277  ,c=0.4582, d= 0.1361  , e=0.1898, f=0.0585) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  rhoa <- 1
  rs <- 69
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))   
  p$Rn <- p$Rns - p$Rnl
  
  ws <- cd$ws
  ws[ws==0] <- min(min(ws), 0.1)
  p$ra <- 208/ws
  
  p$f1 <- (p$ssv/((lambda*rhoa)*(p$ssv + gamma*(1+rs/p$ra))))
  p$f2 <- (gamma/((lambda*rhoa)*(p$ssv + gamma*(1+rs/p$ra))))
  
  p$Wf <- a+b*cd$ws
  p$PET <- p$f1*p$Rn + p$f2*p$Wf*(p$es-p$ea)
  
  return (p$PET)
}

# Rijtema (1966)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
#n1=5.27675  ,b=-3.92509, c=0.49733    , d= 0.35144    , e=0.08175    , f=0.11791
# n1=11.78332, b=-3.96879,c=0.47060, d= 0.32928  , e=0.06618  , f=0.10452 (long term mean ws)
OPET.COMB.Rijtema <- function(cd, n1=5.27675  ,b=-3.92509, c=0.49733    , d= 0.35144    , e=0.08175    , f=0.11791) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))    
  p$Rn <- p$Rns - p$Rnl
  
  p$f <- (p$ssv/(p$ssv + gamma))
  p$PET <- p$f*p$Rn + p$f*gamma*r*(cd$ws^n1)*(p$es-p$ea)+b
  
  return (p$PET)
}

# Penman (1963)
#Source: http://dx.doi.org/10.1016/j.ejrh.2015.02.002
#a=1.89789 , b=5.24909 ,c=0.48481 , d= 0.11190   , e=0.19273   , f=0.05117 
#a=3.30041  , b=5.16168  ,c=0.50643  , d= 0.08647    , e=0.18259    , f=0.02827 (with long term mean ws)
OPET.COMB.Penman3 <- function(cd,a=1.89789 , b=5.24909 ,c=0.48481 , d= 0.11190   , e=0.19273   , f=0.05117) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  rhoa <- 1
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  albedo <- 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))   
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/((lambda*rhoa)*(p$ssv + gamma)))
  p$f2 <- (gamma/((lambda*rhoa)*(p$ssv + gamma)))
  p$W <- a+b*cd$ws
  p$PET <- p$f1*p$Rn + p$f2*p$W*(p$es-p$ea)
  
  return (p$PET)
}

# Penman (1956)
# Source: 10.1016/j.jhydrol.2004.08.026 / https://emf-creaf.github.io/meteolandbook/potentialevapotranspiration.html
#a=1.89789   , b=5.24909   ,c=0.48481   , d= 0.11190     , e=0.19273     , f=0.05117 
#a=0.6111   ,b=2.9373,g=5.5723  ,h=-2.7081, c=0.4871, d= 0.1699    , e=0.1913  , f=0.0776 (long term mean ws)
OPET.COMB.Penman2 <- function(cd,a=1.89789   , b=5.24909   ,c=0.48481   , d= 0.11190     , e=0.19273     , f=0.05117) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  lambda <- 2.45
  rho <- 1
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))     
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/((lambda*rho)*(p$ssv + gamma)))
  p$f2 <- (gamma/((lambda*rho)*(p$ssv + gamma)))
  p$W <- a+b*cd$ws
  p$PET <- p$f1*p$Rn + p$f2*p$W*(p$es-p$ea)
  
  return (p$PET)
}

# Kimberly-Penman dans Wright (1982)
# Source: 10.1016/j.jhydrol.2004.08.026
#a=1.84645 ,b=-0.30170,g=5.33543 ,h=-0.56006, c=0.48233, d= 0.10527  , e=0.19259  , f=0.04604
#a=0.6111  ,b=2.9373  ,g=5.5723 ,h=-2.7081, c=0.4871  , d= 0.1699  , e=0.1913  , f=0.0776   (with long term mean ws)
OPET.COMB.KimberlyPenman <- function(cd, a=1.84645 ,b=-0.30170,g=5.33543 ,h=-0.56006, c=0.48233, d= 0.10527  , e=0.19259  , f=0.04604) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  rhoa <- 1
  
  albedo <- 0.23
  rs_rso_min <- 0.33
  rs_rso_max <- 1
  
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))  
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/((lambda*rhoa)*(p$ssv + gamma)))
  p$f2 <- (gamma/((lambda*rhoa)*(p$ssv + gamma)))
  p$aw <- a+b*exp(-((p$J-173)/58)^2)
  p$bw <- g+h*exp(-((p$J-243)/80)^2)  
  p$Wf <- p$aw + p$bw*cd$ws
  p$PET <- p$f1*p$Rn + p$f2*p$Wf*(p$es-p$ea)
  
  return (p$PET)
}

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000520
#a=0.02986 , b=0.05419 , c=0.02630 , n1=0.68916 , n2=1.48945 , n3=0.23140
OPET.COMB.Valiantzas6 <- function (cd,a=0.02986 , b=0.05419 , c=0.02630 , n1=0.68916 , n2=1.48945 , n3=0.23140) {
  phi <- cd$lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- a*(cd$rs)*((p$tm+9.5)^n1)-b*(cd$rs^n2)*(phi^n3)+c*(p$tm+20)*(1-cd$rh/100)
  return (p$PET)
}

#Valiantzas (2013)
#Source: https://doi.org/10.1061/(ASCE)IR.1943-4774.0000520
#a=0.33588 , b=0.51364  , c=0.02632  , n1=0.24423  , n2=1.07988  , n3=0.03573, n4=1.33586 
#a=0.03423  , b=0.07997   , c=0.02316   , n1=0.65623   , n2=1.33236   , n3=0.16385 , n4=1.61825 (with long term ws)
OPET.COMB.Valiantzas5 <- function (cd,a=0.33588 , b=0.51364  , c=0.02632  , n1=0.24423  , n2=1.07988  , n3=0.03573, n4=1.33586 ) {
  phi <- cd$lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$PET <- a*(cd$rs)*((p$tm+9.5)^n1)-b*(cd$rs^n2)*(phi^n3)+c*(p$tm+20)*(1-cd$rh/100)*(cd$ws^n4)
  return (p$PET)
}

#Valiantzas (2013)
# Source: https://doi.org/10.1061/(ASCE)HE.1943-5584.0000590
#a=0.10067  , b=7.01018  , n1=0.32150   , n2=1.69469,c = 0.02849 , d= 0.03925
OPET.COMB.Valiantzas4 <- function (cd,a=0.10067  , b=7.01018  , n1=0.32150   , n2=1.69469,c = 0.02849 , d= 0.03925) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  p$Cu <- with(cd, ifelse(cd$rh > 65, c, d))
  p$PET <- a*(cd$rs)*((p$tm+9.5)^n1)-b*((cd$rs/p$Ra)^n2)+p$Cu*(p$tm+20)*(1-cd$rh/100)
  return (p$PET)
}

#Valiantzas (2013)
# Source: https://doi.org/10.1061/(ASCE)HE.1943-5584.0000590
#c=0.03793 ,n1=0.52220 ,n2=2.33146 ,d=0.02075 ,e=3.52358,n3=0.88773 ,g=0.03854 ,h=0.04622
#c=0.11271  ,n1=0.30412  ,n2=1.41986  ,d=0.01768  ,e=6.46543 ,n3=0.75655  ,g=0.04055  ,h=0.05295 (long term ws)
OPET.COMB.Valiantzas3 <- function (cd,c=0.03793 ,n1=0.52220 ,n2=2.33146 ,d=0.02075 ,e=3.52358,n3=0.88773 ,g=0.03854 ,h=0.04622) {
  phi <- cd$lat*pi/180
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))
  p$Waero <- with(cd, ifelse(cd$rh>65,g,h))
  p$PET <- c*(cd$rs)*((p$tm+9.5)^n1)-d*(p$tm+20)*(1-cd$rh/100)-e*((cd$rs/p$Ra)^n2)+
    p$Waero*(p$tm+20)*(1-cd$rh/100)*cd$ws^n3
  return (p$PET)
}

#Valiantzas (2006)
#Source: 10.1016/j.jhydrol.2006.06.012
#c=0.08284, n1=0.36471 ,n2=1.67277 ,d=6.63834 ,e=0.04250
OPET.COMB.Valiantzas2 <- function (cd,c=0.08284, n1=0.36471 ,n2=1.67277 ,d=6.63834 ,e=0.04250) {
  phi <- cd$lat*pi/180
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- c*((p$tm+9.5)^n1)*cd$rs-d*((cd$rs/p$Ra)^n2)+e*(p$tm+20)*(1-cd$rh/100)
  return (p$PET)
}


#Valiantzas (2006)
#Source: 10.1016/j.jhydrol.2006.06.012
#a=0.0350847, b=3.0567746  , c=0.0001802 , d=-0.0133132, e=0.0399366  , n1=0.6003212,n2=2.3475096
#a=0.1295987  , b=6.2577842   , c=-0.0006359  , d=-0.0018575  , e=0.0408143    , n1=0.3296744  ,n2=1.5151543 (long term ws)
OPET.COMB.Valiantzas1 <- function (cd,a=0.0350847, b=3.0567746  , c=0.0001802 , d=-0.0133132, e=0.0399366  , n1=0.6003212,n2=2.3475096) {
  phi <- cd$lat*pi/180
  alpha <- 0.23
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)
  p$J <- as.numeric(format(cd$Date, "%j"))  
  p$dr <- 1 + 0.033*cos((2*pi*p$J)/365)
  p$ds <- 0.409*sin(((2*pi*p$J)/365)-1.39)
  p$ws <- acos(-tan(phi)*tan(p$ds))
  p$Ra <- (24*60/pi)*Gsc*p$dr*((p$ws*sin(phi)*sin(p$ds))+(cos(phi)*cos(p$ds)*sin(p$ws)))  
  p$PET <- a*(1-alpha)*((p$tm+9.5)^n1)*cd$rs-b*((cd$rs/p$Ra)^n2)+c*cd$z +
    (p$tm+20)*(1-cd$rh/100)*(d+e*cd$ws)
  return (p$PET)
}

# George et al (1985)
# Source: https://doi.org/10.1016/j.agwat.2020.106043
#a=10.59962 , b=6.83140      ,c=0.44164     , d= -0.02144      , e=0.03843       , f=-0.01034
#a=11.78313   , b=7.56283      ,c=0.43977     , d= -0.04012 , e=0.02929 , f=-0.02382 (long term mean ws)
OPET.COMB.George <- function(cd, a=10.59962 , b=6.83140      ,c=0.44164     , d= -0.02144      , e=0.03843       , f=-0.01034) {
  phi <- cd$lat*pi/180
  patm <- 101.3*(((293-0.0065*cd$z)/293)^5.26)
  
  lambda <- 2.45
  gamma <- 0.0016286*patm/lambda
  sigma <- 4.903e-9
  
  albedo = 0.23
  rs_rso_min = 0.33
  rs_rso_max = 1
  
  p <- data.frame(matrix(nrow = length(cd$Date), ncol = 0))
  p$tm <- 0.5*(cd$tx + cd$tn)  
  p$J <- as.numeric(format(cd$Date, "%j"))
  p$aw <- a+1.4*exp(-((p$J-173)/52)^2)
  p$bw <- b+0.345*exp(-((p$J-243)/80)^2)
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
  
  p$Rs0 <- (0.75+0.00002*cd$z)*p$Ra
  p$Rns <- (1-albedo)*cd$rs
  p$Rs_Rs0 <- cd$rs/p$Rs0
  p$Rs_Rs0 <- sapply(p$Rs_Rs0, function(x) max(min(x, rs_rso_max), rs_rso_min))
  #p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(0.34-0.14*sqrt(p$ea))*((1.35*p$Rs_Rs0)-0.35)
  p$Rnl <- sigma * ((((cd$tx+273.16)^4)+((cd$tn+273.16)^4))/2)*(c*p$Rs_Rs0-d-e*p$Rs_Rs0*sqrt(p$ea)+f*sqrt(p$ea))   
  p$Rn <- p$Rns - p$Rnl
  
  p$f1 <- (p$ssv/(p$ssv + gamma))
  p$f2 <- (gamma/(p$ssv + gamma))
  p$Wf <- p$aw + p$bw*cd$ws
  p$PET <- p$f1*p$Rn + 0.268*p$f2*p$Wf*(p$es-p$ea)
  
  return (p$PET)
}

