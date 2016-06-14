## ------------------------------------ Bern ------------------------------------ ##

setwd("Documents/Universität/TUM/CaseStudy/Wetter")
Bern <- read.table("Bern-hour.dat-2", header=FALSE)
Bern <- Bern[,1:5]
Bern$V4 <- NULL
colnames(Bern)[1] <- 'Month'
colnames(Bern)[2] <- 'Day'
colnames(Bern)[3] <- 'Hour'
colnames(Bern)[4] <- 'Radiation'
BernSSW <- Bern[(Bern$Month==6) & (Bern$Day==19),3:4]
#plot(BernSSW$Hour, BernSSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice Bern ", type = "l")
BernWSW <- Bern[(Bern$Month==12) & (Bern$Day==22),3:4]
#plot(BernWSW$Hour, BernWSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Winter Solstice Bern ", type = "l")
################## Test on Bern ###########
BernSummer <- Bern[(Bern$Month==6),2:4]
#plot(BernSummer$Day, BernSummer$Radiation, xlab = "Day", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice Bern ", type = "l")

BernWinter <- Bern[(Bern$Month==12),2:4]
#plot(BernWinter$Day, BernWinter$Radiation, xlab = "Day", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Winter Solstice Bern", type = "l")


##---------------------------------- BRASILIA -----------------------------------------##

setwd("Documents/Universität/TUM/CaseStudy/Wetter")
Brasil <- read.table("BRASILIA_BR-hour.dat", header=FALSE)
Brasil <- Brasil[,1:5]
Brasil$V4 <- NULL
colnames(Brasil)[1] <- 'Month'
colnames(Brasil)[2] <- 'Day'
colnames(Brasil)[3] <- 'Hour'
colnames(Brasil)[4] <- 'Radiation'
BrasilSSW <- Brasil[(Brasil$Month==6) & (Brasil$Day==23),3:4]
#plot(BrasilSSW$Hour, BrasilSSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice Brasilia ", type = "l")
BrasilWSW <- Brasil[(Brasil$Month==12) & (Brasil$Day==18),3:4]
#plot(BrasilWSW$Hour, BrasilWSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2) , main="Winter Solstice Brasilia ", type = "l")
################## Test on Brasil ###########
BrasilSummer <- Brasil[(Brasil$Month==6),2:4]
#plot(BrasilSummer$Day, BrasilSummer$Radiation, xlab = "Day", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice Brasilia ", type = "l")

BrasilWinter <- Brasil[(Brasil$Month==12),2:4]
#plot(BrasilWinter$Day, BrasilWinter$Radiation, xlab = "Day", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Winter Solstice Brasilia ", type = "l")
## ------------------------------------ Johannisburg ------------------------------------ ##

setwd("Documents/Universität/TUM/CaseStudy/Wetter")
Johannisburg <- read.table("Johannesburg_SF-hour.dat", header=FALSE)
Johannisburg <- Johannisburg[,1:5]
Johannisburg$V4 <- NULL
colnames(Johannisburg)[1] <- 'Month'
colnames(Johannisburg)[2] <- 'Day'
colnames(Johannisburg)[3] <- 'Hour'
colnames(Johannisburg)[4] <- 'Radiation'
JohannisburgSSW <- Johannisburg[(Johannisburg$Month==6) & (Johannisburg$Day==22),3:4] # 22 500
#plot(JohannisburgSSW$Hour, JohannisburgSSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice Johannisburg ", type = "l")
JohannisburgWSW <- Johannisburg[(Johannisburg$Month==12) & (Johannisburg$Day==21),3:4] #21 1100, 22 900
#plot(JohannisburgWSW$Hour, JohannisburgWSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Winter Solstice Johannisburg ", type = "l")


## ------------------------------------ Perth ------------------------------------ ##

setwd("Documents/Universität/TUM/CaseStudy/Wetter")
Perth <- read.table("Perth_AS-hour.dat", header=FALSE)
Perth <- Perth[,1:5]
Perth$V4 <- NULL
colnames(Perth)[1] <- 'Month'
colnames(Perth)[2] <- 'Day'
colnames(Perth)[3] <- 'Hour'
colnames(Perth)[4] <- 'Radiation'
PerthSSW <- Perth[(Perth$Month==6) & (Perth$Day==20),3:4] # 20 oder 21 egal
#plot(PerthSSW$Hour, PerthSSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice Perth ", type = "l")
PerthWSW <- Perth[(Perth$Month==12) & (Perth$Day==23),3:4]
#plot(PerthWSW$Hour, PerthWSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Winter Solstice Perth ", type = "l")

## ------------------------------------ SanFrancisco ------------------------------------ ##

setwd("Documents/Universität/TUM/CaseStudy/Wetter")
SanFrancisco <- read.table("San_Francisco_US-hour.dat", header=FALSE)
SanFrancisco <- SanFrancisco[,1:5]
SanFrancisco$V4 <- NULL
colnames(SanFrancisco)[1] <- 'Month'
colnames(SanFrancisco)[2] <- 'Day'
colnames(SanFrancisco)[3] <- 'Hour'
colnames(SanFrancisco)[4] <- 'Radiation'
SanFranciscoSSW <- SanFrancisco[(SanFrancisco$Month==6) & (SanFrancisco$Day==21),3:4]
#plot(SanFranciscoSSW$Hour, SanFranciscoSSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice SanFrancisco ", type = "l")
SanFranciscoWSW <- SanFrancisco[(SanFrancisco$Month==12) & (SanFrancisco$Day==21),3:4]
#plot(SanFranciscoWSW$Hour, SanFranciscoWSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Winter Solstice SanFrancisco ", type = "l")

## -------------------------- Combinded Plots Summer -------------------------------- ##
plot(BernSSW$Hour, BernSSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Summer Solstice", type = "l", col = "green")
lines(SanFranciscoSSW$Hour, SanFranciscoSSW$Radiation, col= "red")
lines(JohannisburgSSW$Hour, JohannisburgSSW$Radiation, col= "blue")
lines(PerthSSW$Hour, PerthSSW$Radiation, col= "orange")
lines(BrasilSSW$Hour, BrasilSSW$Radiation, col= "purple")

legend ("topright",c("Bern","SanFrancisco","Brasil", "Johannisburg","Perth"), lty=c(1,1,1,1,1), lwd=c(1,2,2,2,2),col=c("green","red","purple","blue","orange"), cex=0.6)

## -------------------------- Combinded Plots Winter -------------------------------- ##

plot(JohannisburgWSW$Hour, JohannisburgWSW$Radiation, xlab = "Hour", ylab = expression(plain(Radiation_in_Wh/m)^2), main="Winter Solstice", type = "l", col = "blue")
lines(SanFranciscoWSW$Hour, SanFranciscoWSW$Radiation, col= "red")
lines(PerthWSW$Hour, PerthWSW$Radiation, col= "orange")
lines(BernWSW$Hour, BernWSW$Radiation, col= "green")
lines(BrasilWSW$Hour, BrasilWSW$Radiation, col= "purple")

legend ("topright",c("Bern","SanFrancisco","Brasil", "Johannisburg","Perth"), lty=c(1,1,1,1,1), lwd=c(1,2,2,2,2),col=c("green","red","purple","blue","orange"), cex=0.6)

