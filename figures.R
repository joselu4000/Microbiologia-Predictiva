#### Figuras TFM ####

#### Figura i: Ejemplo de Ratkowsky mejorado ####
g <- function(TT,x){
  return(c*(TT-x))
}
ratkowsky <- function(TT){
  return(
    b*(g(TT,Tmin)/c)*(1-exp(g(TT,Tmax)))
  )
}
# Temperature vector
temp <- seq(from=-12,to=50,length.out =10000)
# Example 1
b <- 0.025
Tmin <- -15
c <- 0.30
Tmax <- 45
par(mar = c(5, 4, 4, 2) + 0.1)
plot(temp, ratkowsky(temp), type = "l",
     xlab = "Temperature (ºC)", 
     ylab = expression(sqrt(mu[max])~(h^{-1/2})),
     xlim = c(-10, 45), ylim = c(0, 1.2),
     xaxt = "n", yaxt = "n",  
     lwd = 2, col = "blue") 
# Example 2
b <- 0.025
Tmin <- -10
c <- 0.30
Tmax <- 40
points(temp, ratkowsky(temp), type = "l", col = "purple") 
z <- optimize(function(TT) -ratkowsky(TT), interval = c(0, 45))
points(rep(z$minimum, 2), c(-10, ratkowsky(z$minimum)), type = "l", lty = 2) 
# Example 3
b <- 0.025
Tmin <- -5
c <- 0.30
Tmax <- 30
points(temp, ratkowsky(temp), type = "l", col = "red")
# Example 4
b <- 0.025
Tmin <- 0
c <- 0.30
Tmax <- 20
points(temp, ratkowsky(temp), type = "l", col = "orange")
# Axis
axis(1, at = seq(-10, 45, by = 5), labels = seq(-10, 45, by = 5)) 
axis(2, at = pretty(c(0, 1.2)), labels = pretty(c(0, 1.2)))  
#### Bisection method ####

x0 <- 20
de <- "l"
z0 <- 35
{
  temp <- seq(from=10,to=40,length.out =1000)
  plot(temp,ratkowsky(temp), type="l", ylim = c(0,1))
  if(de == "l"){
    points(c(-10,40),c(ratkowsky(x0),ratkowsky(x0)),type = "l",col="red")
  }else{
    points(c(-10,40),c(ratkowsky(z0),ratkowsky(z0)),type = "l",col="red")
  }
  i <- 2
  while(e <= 1e-34 && ee <= 1e-34){
    if(de == "l"){
      P <- abs(ratkowsky(temp)-ratkowsky(x0))
      pos.1 <- which(P <= 1e-2)
      z0 <- temp[pos.1][length(pos.1)]
    }else{
      P <- abs(ratkowsky(temp)-ratkowsky(z0))
      pos.1 <- which(P <= 1e-2)
      x0 <- temp[pos.1][length(pos.1)]
    }
    e <- abs(x0-z0)
    ee <- abs(ratkowsky(x0)-ratkowsky(z0))
    p <- (x0+z0)/2
    cat("(x0,p,z0)=",x0,p,z0)
    temp <- seq(from = x0, to = z0, length.out=1000)
    Q <- abs(ratkowsky(temp)-ratkowsky(p))
    pos.2 <- which(Q <= 1e-2)
    q <- temp[pos.2][length(pos.2)]
    if(p <= q && q <= z0){
      x0 <- p
      de <- "l"
    }else{
      z0 <- p
      de <- "r"
    }
    points(c(-10,40),rep(ratkowsky(p),2), type="l", col = i)
    i <- i+1
  }
  topt.2 <- p
  cat("El punto optimo es",p)
}
temp <- seq(from=-10,to=40,length.out =1000)
plot(temp,ratkowsky(temp), type="l", ylim = c(0,1))
points(c(topt.1,topt.1),c(0,ratkowsky(topt.1)), type = "l", col = "red")
points(c(topt.2,topt.2),c(0,ratkowsky(topt.2)), type = "l", col = "blue")

result <- optimize(function(TT) -ratkowsky(TT), interval = c(20, 35))
cat("El máximo valor de la función es:", -result$objective, "\n")
cat("El valor de TT que maximiza la función es:", result$minimum, "\n")

####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+1: Ejemplo de funciones gamma ####

gamT <- function(t,tm,topt){
  return(((t-tm)/(topt-tm))^2)
}
gamaw <- function(aw,awm){
  return((aw-awm)/(1-awm))
}
gampH <- function(pH,pHm,pHopt,pHM){
  return(
    ((pH-pHm)*(pHM-pH))/((pHopt-pHm)*(pHM-pHopt))
  )
}
mumaxG <- function(t,aw,pH){
  return(
    muopt*gamT(t,tm,topt)*gamaw(aw,awm)*gampH(pH,pHm,pHopt,pHM)
  )
}
# Param
tm <- -8
topt <- 31.48915
awm <- 0.93
pHm <- 5.7
pHM <- 9
pHopt <- (pHM+pHm)/2 #7.8
muopt <- 0.3
# Graph
t <- seq(from = tm, to = topt, length.out = 1e4)
aw <- seq(from=awm,to=1,length.out=1e4)
pH <- seq(from=pHm,to=pHM,length.out=1e4)
png("Figures/GammaExample.png", width = 2100, height = 1670, res = 300)
par(mfrow = c(2,2))
plot(t,gamT(t,tm,topt), type = "l", 
     xlab = "Temperature (ºC)", ylab = expression(gamma(T)))
plot(aw,gamaw(aw,awm), type = "l",
     xlab = expression(a[w]), ylab = expression(gamma(a[w])))
plot(pH,gampH(pH,pHm,pHopt,pHM),type="l",
     xlab = expression(pH), ylab = expression(gamma(pH)))
plot(t,mumaxG(t,aw,pH),type = "l",
     xlab = "Temperature (ºC)", ylab = expression(mu[max](T,a[w],pH)))
par(mfrow = c(1,1))
dev.off()
####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+2: Ejemplo de funciones gamma pH, variando pHopt ####

gampH <- function(pH,pHm,pHopt,pHM){
  return(
    ((pH-pHm)*(pHM-pH))/((pHopt-pHm)*(pHM-pHopt))
  )
}
mumaxGpHopt <- function(t,aw,pH,pHopt){
  return(
    muopt*gamT(t,tm,topt)*gamaw(aw,awm)*gampH(pH,pHm,pHopt,pHM)
  )
}
tm <- -8
topt <- 31.48915
awm <- 0.93
pHm <- 5.7
pHM <- 9
pHopt <- (pHM+pHm)/2 #7.8
muopt <- 0.3
pHm <- 5.7
pHM <- 9
pHopt <- (pHM+pHm)/2 #7.8
pH <- seq(from=pHm,to=pHM,length.out=1e4)
z <- seq(from=pHopt-1,to=pHM,length.out=10)
png("Figures/VariandopHoptGamma", width = 2100, height = 1670, res = 300)
par(mfrow = c(1,2))
{
plot(pH,gampH(pH,pHm,pHopt,pHM),type="l", col="red",
     xlab = expression(pH), ylab = expression(gamma(pH)), 
     ylim = c(0,4), lwd=5)
for (i in z){
  points(pH,gampH(pH,pHm,i,pHM),type="l")
}
points(pH,gampH(pH,pHm,pHopt,pHM),type="l", col="red",lwd=5)
}
{
plot(pH,mumaxGpHopt(t,aw,pH,pHopt),type = "l",
     col = "red", lwd = 3, ylim = c(0,0.4),
     xlab = "pH", ylab = expression(mu[max](T,a[w],pH)))
for (i in z){
  points(pH,mumaxGpHopt(t,aw,pH,i),type="l")
}
points(pH,mumaxGpHopt(t,aw,pH,pHopt),type="l", col="red",lwd=5)
}
dev.off()

####------------------------------------------------------------------------####
#### Figura i+3: Ejemplo de funciones CMn variando todo                     ####
####------------------------------------------------------------------------####

CM_n <- function(X, X_min, X_max, X_opt, n) {
  # Función indicatriz
  I <- ifelse(X >= X_min & X <= X_max, 1, 0)
  
  # Definición de p_n(X)
  p_n <- function(X) {
    (X_opt - X_min) * (X - X_opt) - 
      (X_opt - X_max) * ((n - 1) * X_opt + X_min - n * X)
  }
  
  # Cálculo de CM_n(X)
  CM_n_value <- ((X - X_max) * (X - X_min)^n) / 
    ((X_opt - X_min)^(n - 1) * p_n(X)) * I
  
  return(CM_n_value)
}

# Guardar la imagen en proporciones de A4
png("MC_nEvolution_A4.png", width = 2100, height = 2970, res = 300) # A4 con 300 dpi

# Definir los parámetros
X_min <- c(2,3,4)
X_max <- c(8,7,6)
X_opt <- c(4.5,5,5.5)
n <- c(1,2,3)

# Configuración de la cuadrícula 3x3
par(mfrow = c(3,3), mar = c(4.5, 5, 4, 2)) # Márgenes ajustados: abajo, izquierda, arriba, derecha


# Generar las gráficas
for (k in 1:3){
  for (j in 1:3){
    # Configurar la gráfica base
    plot(0, 0, xlim = c(2, 8), ylim = c(0, 1.6), type = "n",
         xlab = "X", ylab = "CM_n(X)",
         main = paste("X_opt =", X_opt[j], ", n =", n[k]), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.2)
    for (i in 1:3){
      X <- seq(from = X_min[i], to = X_max[i], length.out = 1000)
      result <- CM_n(X, X_min[i], X_max[i], X_opt[j], n[k])
      # Graficar las curvas
      lines(X, result, col = i, lwd = 2)
    }
  }
}

# Cerrar el dispositivo gráfico para guardar la imagen
#dev.off()


X2_star <- function(X_max, X_min, X_opt, n) {
  numerator <- X_max^2 * (n - 1) + 2 * X_max * X_min - X_opt * (n * X_max + X_min)
  denominator <- n * (X_max - X_opt) + (X_min - X_opt)
  X2_star <- numerator / denominator
  return(X2_star)
}

png("MC_nXestEvolution.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(3,2))
X_max <- 8
X_min <- 3
X_opt <- 5
n <- 1:6
X <- seq(from=3,to=8,length.out=1000)
for (k in n){
  result <- CM_n(X,X_min, X_max, X_opt, k)
  X2 <- X2_star(X_max, X_min, X_opt, k)
  if (k <= 4){
    L <- 1.5
  }else{ L <- 3.5}
  fXo <- CM_n(X_opt, X_min, X_max, X_opt, k)
  fX2 <- CM_n(X2, X_min, X_max, X_opt, k)
  text1 <- bquote(n == .(k))
  if (fXo > fX2) {
    text2 <- bquote(X[opt] > X[2]^"*" ~ "=" ~ .(round(X2, 4)))
  } else {
    text2 <- bquote(X[opt] < X[2]^"*" ~ "=" ~ .(round(X2, 4)))
  }
  main_title <- bquote(.(text1) ~ "," ~ .(text2))
  plot(
    X, result, type = "l", col = "black", lwd = 2, 
    ylim = c(0, L),
    main = main_title,
    ylab = ""
  )
  points(X2,CM_n(X2,X_min, X_max, X_opt, k), col="red", lwd=2)
  points(X_opt,CM_n(X_opt,X_min, X_max, X_opt, k), col="blue", lwd=2)
  lines(c(3,8),rep(1,2),col="black", lty=2)
  lines(rep(X2,2),c(-2,CM_n(X2,X_min, X_max, X_opt, k)),type = "l", col="red")
  lines(rep(X_opt,2),c(CM_n(X_opt,X_min, X_max, X_opt, k),6),type = "l", col="blue")
}
dev.off()


####------------------------------------------------------------------------####
#### Figura i+4: Ejemplo de funciones polinomios grado 2                    ####
####------------------------------------------------------------------------####
png("polinomialExample.png", width = 2100, height = 2970, res = 300)
{
f <- function(x){
  return(c(1,x,x^2,x^3))
}
pol2 <- function(b,x){
  return(b[1]+b[2]*x+b[3]*x^2+b[4]*x^3)
}
par(mfrow = c(2,1))
layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)

# Aplicar la cuadrícula
layout(layout_matrix)

# pH
x <- c(3,5,7,9)
y <- matrix(c(0,0.3,0.9,0.5),ncol = 1)
plot(x,y, type = "p",ylim = c(0,1), col="red", xlab = "pH", ylab = "h(pH)")
dat <- as.data.frame(cbind(x,y))
X <- matrix(ncol = 4, nrow = 4)
for (i in 1:nrow(X)){
  X[i,] <- f(x[i])
}
A <- t(X)%*%X 
invA <- solve(A)
b <- invA%*%t(X)%*%y
t <- seq(from=2,to=10, length.out=100)
lines(t,pol2(b,t)+rnorm(100,mean = 0, sd = 0.01))
points(x,y, type = "p",ylim = c(0,1), col="red", lwd=3)

# aw
gamaw <- function(aw,awm){
  return((aw-awm)/(1-awm))
}
xw <- seq(from=0.96,to=1, length.out=5)
yw <- matrix(gamaw(xw,0.96),ncol = 1)
plot(xw,yw, type = "p",ylim = c(0,1), col="red", 
     xlab = expression(a[w]), ylab = expression(h(a[w])))
dat <- as.data.frame(cbind(xw,yw))
Xw <- matrix(ncol = 4, nrow = length(xw))
for (i in 1:nrow(Xw)){
  Xw[i,] <- f(xw[i])
}
Aw <- t(Xw)%*%Xw  
invAw <- solve(Aw)
bw <- invAw%*%t(Xw)%*%yw
tw <- seq(from=0.96,to=1, length.out=100)
lines(tw,pol2(bw,tw)+rnorm(100,mean = 0, sd = 0.01))
points(xw,yw, type = "p",ylim = c(0,1), col="red", lwd=3)

# mu_max
gamaw <- function(aw,awm){
  return((aw-awm)/(1-awm))
}
gampH <- function(pH,pHm,pHopt,pHM){
  return(
    ((pH-pHm)*(pHM-pH))/((pHopt-pHm)*(pHM-pHopt))
  )
}
mumaxG <- function(aw,pH){
  return(
    muopt*gamaw(aw,awm)*gampH(pH,pHm,pHopt,pHM)
  )
}
# Param
awm <- 0.93
pHm <- 5.7
pHM <- 9
pHopt <- (pHM+pHm)/2
muopt <- 0.3
g <- function(x,y){
  return(c(1,x,y,x*y,x^2,y^2,x^3,y^3,x*y^2,x^2 * y))
}
pol.T <- function(b,x,y){
  return(
    b[1] +b[2]*x +b[3]*y +b[4]*x*y +b[5]*x^2 +b[6]*y^2 +b[7]*x^3 +b[8]*y^3 + b[9]*x*y^2 + b[10]*x^2 * y
  )
}
m <- 10
aw <- seq(from=awm,to=1,length.out=m)
pH <- seq(from=pHm,to=pHM,length.out=m)
mumax <- mumaxG(aw,pH)
plot(pH,mumax, type = "p", col="red", xlab = "pH", ylab = expression(h(pH,a[w])),
     ylim = c(0,0.22))
dat <- as.data.frame(cbind(pH,mumax))
X <- matrix(ncol = m, nrow = length(pH))
for (i in 1:nrow(X)){
  X[i,] <- g(pH[i],aw[i])
}
A <- t(X)%*%X +diag(1e-8, ncol = 10, nrow = 10)
invA <- solve(A)
b.T <- invA%*%t(X)%*%mumax
t1 <- seq(from=pHm,to=pHM, length.out=100)
t2 <- seq(from=awm,to=1, length.out=100)
lines(t1,pol.T(b.T,t1,t2)+rnorm(100,mean = 0, sd = 0.01))
points(pH,mumax, type = "p", col="red", lwd=3)
}
dev.off()


####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+5: Ejemplo de funciones gamma ####

gamT <- function(t,tm,topt){
  return(((t-tm)/(topt-tm))^2)
}
gamaw <- function(aw,awm){
  return((aw-awm)/(1-awm))
}
gampH <- function(pH,pHm,pHopt,pHM){
  return(
    ((pH-pHm)*(pHM-pH))/((pHopt-pHm)*(pHM-pHopt))
  )
}
mumaxG <- function(t,aw,pH){
  return(
    muopt*gamT(t,tm,topt)*gamaw(aw,awm)*gampH(pH,pHm,pHopt,pHM)
  )
}
# Param
tm <- -8
topt <- 31.48915
awm <- 0.93
pHm <- 5.7
pHM <- 9
pHopt <- (pHM+pHm)/2
muopt <- 0.3
# Graph
noisy <- 0.01
{
l <- 10
tx <- seq(from = tm, to = topt, length.out = l)
aw <- seq(from=awm,to=1,length.out=l)
pH <- seq(from=pHm,to=pHM,length.out=l)
mumaxx <- mumaxG(tx,aw,pH)
p <- 5
MUMU <- matrix(ncol = p, nrow = l)
for (i in 1:p){
  MUMU[,i] <- mumaxx+rnorm(l,0,noisy)
}
l <- 20
t <- seq(from = tm, to = topt, length.out = l)
aw <- seq(from=awm,to=1,length.out=l)
pH <- seq(from=pHm,to=pHM,length.out=l)
mumax <- mumaxG(t,aw,pH)
# IC
ICmumax <- function(mumax){
  meanmumax <- mean(mumax)
  sc <- sum((mumax-meanmumax)^2)/(length(mumax)-1)
  se_mu_max <- sqrt(sc)
  t_critical <- qt(0.975, df = length(t) - 2)
  return(t_critical * se_mu_max)
}
critical <- unlist(apply(MUMU,1,function(x) ICmumax(x)))
critical <- ifelse(is.na(critical),1e-7,critical)
upper_bound <- mumaxx+critical
lower_bound <- mumaxx-critical
par(mfrow = c(1,1))
plot(rep(tx,p), MUMU, type = "p",
     xlab = "Temperature (ºC)", ylab = expression(mu[max](T,pH,a[w])),
     ylim = c(-0.025,0.12))
lines(t, mumax, col = "black", lty = 1)
lines(tx, upper_bound, col = "red", lty = 2)
lines(tx, lower_bound, col = "red", lty = 2)
}

####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+6: Trifase dat2 ####
#png("Figures/TrifaseComparationNumerical.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))
source("Trifase Primario/TrifaseANRD1.R")
source("Trifase Primario/TrifaseANR.R")
source("Trifase Primario/TrifaseGNRD1True.R")
source("Trifase Primario/TrifaseGNRDTrue.R")
#dev.off()

#png("Figures/TrifaseComparationNumerical2.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))
source("Trifase Primario/TrifaseANAdamD1.R")
source("Trifase Primario/TrifaseANAdamD2.R")
source("Trifase Primario/TrifaseAdamD1.R")
source("Trifase Primario/TrifaseAdamD2.R")
#dev.off()

par(mfrow = c(1,1))

####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+6: Logistico dat2 ####
#png("Figures/LogisticoComparationNumerical.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))
source("Logistico Primario/GomperztANRD1.R")
source("Logistico Primario/GomperztANR.R")
source("Logistico Primario/GomperztGNRD1.R")
source("Logistico Primario/GomperztGNR.R")
#dev.off()

#png("Figures/LogisticoComparationNumerical2.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))
source("Logistico Primario/GomperztANAdamD1.R")
source("Logistico Primario/GomperztANAdamD2.R")
source("Logistico Primario/GomperztAdamD1.R")
source("Logistico Primario/GomperztAdamD2.R")
#dev.off()

par(mfrow = c(1,1))


# ####------------------------------------------------------------------------####
# ####------------------------------------------------------------------------####
# #### Figura i+7: Logistico Modificado dat2 ####
# #png("Figures/ModificadoComparationNumerical.png", width = 2100, height = 2970, res = 300)
# par(mfrow = c(2,2))
# source("Modificado Primario/ModificadoANRD1.R")
# source("Modificado Primario/ModificadoANR.R")
# source("Modificado Primario/ModificadoANAdamD1.R")
# source("Modificado Primario/ModificadoANAdamD2.R")
# #dev.off()
# 
# par(mfrow = c(1,1))

####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+8: Baranyi dat2 ####
#png("Figures/BaranyiComparationNumerical1.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))
source("Baranyi Primario/BaranyiANRD1.R")
source("Baranyi Primario/BaranyiANR.R")
source("Baranyi Primario/BaranyiGND1.R")
source("Baranyi Primario/BaranyiGNR.R")
#dev.off()

# Adam
#png("Figures/BaranyiComparationNumerical2.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))
source("Baranyi Primario/BaranyiANAdamD1.R")
source("Baranyi Primario/BaranyiANAdamD2.R")
source("Baranyi Primario/BaranyiAdamD1.R")
source("Baranyi Primario/BaranyiAdamD2.R")
#dev.off()

par(mfrow = c(1,1))

####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+8: LagExpo dat2 ####
#png("Figures/LagExpoComparationNumerical.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,2))
source("LagExpo Primario/LagExpoGND1R.R")
source("LagExpo Primario/LagExpoGNR.R")
source("LagExpo Primario/LagExpoAdamD1.R")
source("LagExpo Primario/LagExpoAdamD2.R")
#dev.off()

par(mfrow = c(1,1))


####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+8: LagExpo dat2 ####
library(readxl)
datos <- read_excel("datos/ComBaseExport.xlsx")
datos <- datos[1:151,]
A <- function(x, q, nu) {
  x + (1/nu) * log((exp(-nu * x) + q) / (1 + q))
}

f <- function(q, mumax, nu, m, y0, ymax, x) {
  a <- A(x, q, nu)
  y0 + mumax * a - (1/m) * log(1 + (exp(m * mumax * a) - 1) / exp(m * (ymax - y0)))
}
# Recta ideal 
rect <- function(t,mumax,y0,tlag){
  z <- y0 + mumax*(t-tlag)
  return(z)
}
# Recta de Baranyi
termA <- function(t,mumax,y0,nu,q){
  a <- A(t, q, nu) 
  z <- y0 + mumax*t
  return(z)
}
# Recta general
recta <- function(t,m,n){
  return(m*t+n)
}
dat0 <- datos$`Time (Hours)`
y_vals <- datos$`Log10 Concentration`
q0 = 499
mumax=0.39
nu=1
m=1.5
y0=min(y_vals)
ymax=max(y_vals)
tlag = 4.45

png("Figures/Estudioh.png", width = 2100, height = 1670, res = 300)
par(mfrow = c(1,1), # grid
    mgp = c(0, 1, 1), # ejes
    mar = c(4, 4, 2, 2) + 0.5 # ventana de textos
    )
plot(0,0,type="n",xlim=c(0,max(dat0)), 
     lwd =2, axes = FALSE,
     xlab = "Tiempo (h)",
     ylim = c(0, max(y_vals) + 1),
     ylab = expression(log(ufc/g)),
)
# Curva
lines(dat0,y_vals,col="black",lwd=3)
# Recta de crecimiento
lines(tempo,rect(tempo,mumax,y0,tlag),col="red",lwd=3,lty=2)
# Limitadores tlag
points(c(0,tlag),c(y0,y0),col="blue",lwd=2,lty=2, type="l")
points(c(0),c(y0),col="blue",lwd=2,pch=16)
text(0,y0,labels = expression(y[0]),pos=2)
points(0,rect(0,mumax,y0,tlag),col="red",lwd=2,pch=16)
text(0,rect(0,mumax,y0,tlag),labels = expression(y[id]),pos=2)
points(c(tlag,tlag),c(0,y0),col="blue",lwd=2,lty=2, type="l")
# h
h <- y0-rect(0,mumax,y0,tlag)
text(tlag,y0/2,labels = paste("h =",h),pos=4)
text(tlag,0,labels = expression(t[lag]),pos=1)
# Ejes
points(c(0,40),rep(0,2),col="black",type="l")
points(rep(0,2),c(0,9),col="black",type="l")
dev.off()
# h(t) = log(1+1/q(t)) => q = 1/(exp(h)-1)
q <- 1/(exp(h)-1)
cat("q es",q)

# Resolución del sistema tiempos maximos
# H <- h^{-1}
# A <- matrix(c(1,-H*tlag,-mumax,1), ncol=2, byrow = TRUE)+diag(1e-2,2,2)
# b <- matrix(c(tlag*(1-y0),y0-mumax*tlag),ncol=1, byrow = TRUE)
# sol <- solve(A)%*%b # t_max, y_max
# sol

####------------------------------------------------------------------------####
####------------------------------------------------------------------------####
#### Figura i+9: ####
f  <- function(x){
  return(x^2 + x + 1)
}

P <- function(x,l,u){
  lapply(x, function(x) max(l,min(x,u)))
}

library(readxl)
datos <- read_excel("datos/ComBaseExport.xlsx")
datos <- datos[1:151,]
A <- function(x, q, nu) {
  x + (1/nu) * log((exp(-nu * x) + q) / (1 + q))
}

f <- function(q, mumax, nu, m, y0, ymax, x) {
  a <- A(x, q, nu)
  y0 + mumax * a - (1/m) * log(1 + (exp(m * mumax * a) - 1) / exp(m * (ymax - y0)))
}

dat0 <- datos$`Time (Hours)`
y_vals <- datos$`Log10 Concentration`
q0 = 499
mumax=0.39
nu=1
m=1.5
y0=min(y_vals)
ymax=max(y_vals)
tlag = 4.45

#png("Figures/Pproyection.png", width = 2100, height = 2970, res = 300)
par(mfrow = c(2,1), # grid
    mgp = c(0, 1, 1), # ejes
    mar = c(5, 4, 4, 2) + 0.5 # ventana de textos
)

# 1st
t <- seq(-0.5-5,-0.5+5,length.out = 100)
plot(t,f(t), type = "l", lwd = 2, axes = FALSE,
     xlab = "", ylab = expression(pol[2](Tiempo)))
lines(t, P(f(t), 0, 15), col="red", lwd=3, lty=2)
lines(t, P(f(t), 5, 15), col="blue", lwd=3, lty=3)

# 2nd
plot(0,0,type="n",xlim=c(0,max(dat0)), 
     lwd =2, axes = FALSE,
     xlab = "Tiempo (h)",
     ylim = c(0, max(y_vals) + 1),
     ylab = expression(log[10](ufc/g)),
)
# Curva
lines(dat0, y_vals,col="black",lwd=3)
lines(dat0, P(y_vals,4,6),col="red", lwd=3, lty=2)
lines(dat0, P(y_vals,0,5),col="blue", lwd=3, lty=3)
#dev.off()

par(mfrow = c(1,1))




































