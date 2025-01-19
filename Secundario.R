#### Primary and Secondary Model ####
library(readxl)
library(tidyverse)
library(minpack.lm)

CM_n <- function(n,X,Xopt,Xmax,Xmin,muopt){
  z1 <- (X-Xmax)*((X-Xmin)^n)
  z2 <- (Xopt-Xmin)^(n-1)
  aux <- (Xopt-Xmin)*(X-Xopt)-(Xopt-Xmax)*((n-1)*Xopt+Xmin-n*X)
  z3 <- z2*aux
  z <- ifelse(X < Xmax & X > Xmin, (muopt*z1)/z3, 0)
  return (z)
}

tau <- function(X,Xopt,Xmax,Xmin,muopt){
  z1 <- (X-Xmin)*(X-Xmax)
  z2 <- (X-Xmin)*(X-Xmax)
  z3 <- (X-Xopt)^2
  z4 <- z2 - z3
  z <- ifelse(X < Xmax & X > Xmin,(muopt*z1)/z4,0)
  return(z)
}

rho <- function(muopt,pH,pHmin,pH1.2){
  z <- muopt*(1-2^((pHmin-pH)/(pH1.2-pHmin)))/(1-2^((pHmin-7.3)/(pH1.2-pHmin)))
  return(z)
}

sales <- function(X,a,b,c,muopt){
  z <- muopt*(a*X^2+b*X+c)
  return(z)
}

metrics <- function(fit, df){
  df <- na.omit(df)
  residuos <- residuals(fit)
  ajustados <- fitted(fit)
  SST <- sum((df$GR - mean(df$GR))^2)  
  SSE <- sum((df$GR-ajustados)^2)                  
  n <- nrow(df)      
  p <- length(coef(fit))  
  # R^2
  R2 <- 1 - (SSE / SST)
  R2_adj <- 1 - (1 - R2) * ((n - 1) / (n - p - 1))
  MSE <- mean((df$GR-fitted(fit))^2)
  return(c(R2,R2_adj, MSE))
}


hojasaw <- excel_sheets("datos/Completed dataset - NaCl_250723.xls")
dataw <- read_excel("datos/Completed dataset - NaCl_250723.xls", sheet = hojasaw[2])
hojaspH <- excel_sheets("datos/Completed dataset pH_250723_final.xlsx")
datpH <- read_excel("datos/Completed dataset pH_250723_final.xlsx", sheet = hojaspH[2])
nam <- "Meat"
strain <- "PM1"
dfaw <- dataw %>% filter(Strain == strain, Origin == nam)
dfaw$GR <- as.numeric(dfaw$GR)
dfaw$aw <- as.numeric(dfaw$aw)
dfaw$NaCl <- as.numeric(dfaw$NaCl)
dfpH <- datpH %>% filter(Strain == strain, Origin == nam)
dfpH$GR <- as.numeric(dfpH$GR)
dfpH$pH <- as.numeric(dfpH$pH)

png("Figures/PractiseCaseSecondaryModel.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))

# Cardinal aw n=2: ####
fit_aw <- nlsLM(GR ~ CM_n(2, aw, 0.997, 1, awmin, muopt), data = dfaw,
                start = c(awmin = 0.93, muopt = 0.906),
                algorithm = "port",
                control = nls.lm.control(maxiter = 100))
summary(fit_aw)
metrics(fit_aw, dfaw)
caw <- coef(fit_aw)

# Cardinal NaCl n=2: ####
# fit_NaCl <- nlsLM(GR ~ CM_n(2, NaCl, NaClopt, NaClmax, NaClmin, muopt), data = dfaw,
#                 start = c(NaClopt = 4, NaClmax = 8, NaClmin = 1, muopt = 1),
#                 algorithm = "port",
#                 control = nls.lm.control(maxiter = 100))
# a = -0.0124, b = 0.0543, c = 0.5397, muopt = 1.2955
ajustados <- sales(dfaw$NaCl[1:10], a = -0.0124, b = 0.0543, c = 0.5397, muopt = 1.2955)
SST <- sum((dfaw$GR[1:10] - mean(dfaw$GR, na.rm = TRUE))^2)  
SSE <- sum((dfaw$GR[1:10]-ajustados)^2)                  
n <- nrow(dfaw)      
p <- 4
# R^2
R2 <- 1 - (SSE / SST)
R2_adj <- 1 - (1 - R2) * ((n - 1) / (n - p - 1))
MSE <- mean((dfaw$GR[1:10]-ajustados)^2)
cat(c(R2,R2_adj, MSE))

# pH rho (Aryani, 2015)####
fit_pH <- nlsLM(GR ~ rho(muopt,pH,pHmin,pH1.2), data = dfpH,
                 start = c(muopt = 1, pHmin = 4, pH1.2 = 5),
                 algorithm = "port",
                 control = nls.lm.control(maxiter = 100))
summary(fit_pH)
metrics(fit_pH, dfpH)
cpH <- coef(fit_pH)

# pH tau ####
fit_pH2 <- nlsLM(GR ~ tau(pH,pHopt,pHmax,pHmin,muopt), data = dfpH,
                 start = c(pHopt = 5.057667, pHmax = 10, pHmin = 4.576, muopt = 1.310111),
                 algorithm = "port",
                 control = nls.lm.control(maxiter = 100))
summary(fit_pH2)
metrics(fit_pH2, dfpH)
cpH2 <- coef(fit_pH2)

# aw plot ####
awT <- seq(0.92,1,length.out = 100)
muaw <- CM_n(2,awT,0.997,1,caw[1],caw[2])
plot(GR ~ aw, data = dfaw, xlim = c(0.92,1), 
     main = "", ylim = c(0,1),
     xlab = "aw", ylab = expression(mu[max](h^-1)))
points(awT,muaw, type = "l", col = "red")

# NaCl plot
NaClT <- seq(0,8,length.out = 100)
muNaCl <- sales(NaClT, a = -0.0124, b = 0.0543, c = 0.5397, muopt = 1.2955) #CM_n(2,NaClT, cNaCl[1], cNaCl[2], cNaCl[3], cNaCl[3])
plot(GR ~ NaCl, data = dfaw, xlim = c(0,8), 
     main = "", ylim = c(0,1),
     xlab = "NaCl", ylab = expression(mu[max](h^-1)))
points(NaClT,muNaCl, type = "l", col = "red")


# pH rho ####
pHT <- seq(4,7.5,length.out = 100)
mupH <- rho(cpH[1],pHT,cpH[2],cpH[3])
plot(GR ~ pH, data = dfpH, ylim = c(0,1.5), xlim = c(4,7.5),
     main = "",
     xlab = "pH", ylab = expression(mu[max](h^-1)))
points(pHT,mupH, type = "l", col = "red")

# pH tau ####
pHT2 <- seq(4,7.5,length.out = 100)
mupH2 <- tau(pHT2,cpH2[1],cpH2[2],cpH2[3],cpH2[4])
plot(GR ~ pH, data = dfpH, ylim = c(0,1.5), xlim = c(4,7.5),
     main = "",
     xlab = "pH", ylab = expression(mu[max](h^-1)))
points(pHT2,mupH2, type = "l", col = "red")
par(mfrow = c(1,1))


# Mumax ####
mumax1 <- CM_n(2,mean(awT),0.997,1,caw[1],caw[1])*rho(cpH[1],mean(pHT),cpH[2],cpH[3])*sales(mean(NaClT), a = -0.0124, b = 0.0543, c = 0.5397, muopt = 1.2955)
cat("Mumax con rho es:",mumax1,"\n")
mumax2 <- CM_n(2,mean(awT),0.997,1,caw[1],caw[1])*tau(mean(pHT2),cpH2[1],cpH2[2],cpH2[3],cpH2[4])*sales(mean(NaClT), a = -0.0124, b = 0.0543, c = 0.5397, muopt = 1.2955)
cat("Mumax con tau es:",mumax2,"\n")
dev.off()
