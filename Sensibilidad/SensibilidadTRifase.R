library(readxl)
datos <- read_excel("datos/datos.xlsx")

trifase <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, y_0,
         ifelse(t<=t_max, y_0+mu*(t-t_lag), y_max))
  return(z)
}

tempo <- c(0:552)
param <- c(2.161,0,307.761,7.812)
TT <- trifase(tempo,param)
plot(tempo, TT, type = "l")
points(datos$x,datos$y, col = "red")

#### Analisis de sensibilidad####

# Por filas es el tempo
# Por columnas derivadas parametros
dy_0 <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 1-(t-t_lag)/(t_max-t_lag),
              ifelse(t<=t_max, 1, 0))
  return(z)
}

dt_lag <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 0,
              ifelse(t<=t_max, mu*(t-t_max), 
                     0))
  return(z)
}

dt_max <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 0,
              ifelse(t<=t_max, mu*(t-t_lag)/(t_max-t_lag), 0))
  return(z)
}

dy_max <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- 1/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 0,
              ifelse(t<=t_max, mu*(t-t_lag), 1))
  return(z)
}

N <- 1e3
tempo <- seq(min(datos$x), max(datos$x), length.out = N)
S <- matrix(nrow = N, ncol = 4)
for (i in 1:N){
  S[i,] <- c(dt_lag(tempo[i],param),
            dt_max(tempo[i],param),
            dy_0(tempo[i],param),
            dy_max(tempo[i],param))
}
colnames(S) <- c("t_lag","t_max","y_0","y_max")

####
plot(c(), c(), xlim = c(0, 552), ylim = c(-6, 2), 
     xlab = "Tiempo (min)", ylab = expression(f[theta[i]]))
for (i in 1:4) {
  lines(tempo, S[, i], col = i, lty = i, lwd = 3)
}
legend("right", legend = colnames(S), col = 1:4, lty = 1:4, lwd = 3)

#### Estudio correlacion ####

MS <- t(S)%*%S
heatmap(cor(MS))

delta <- seq(0,1,length.out=100)
cov_tot <- function(MS,ci, delta){
  s <- 0
  CC <- cor(MS)
  for (i in 1:ncol(MS)){
    if (i != ci){
      p <- CC[i,ci]
      s <- s + p*(p >= 1-delta)
    }
  }
  return(s)
}

EV <- matrix(nrow=length(delta),ncol = ncol(MS))
colnames(EV) <- colnames(MS)
#png("Figures/Pproyection.png", width = 2100, height = 2970, res = 300)
plot(c(),c(),ylim=c(0,2),xlim=c(0,1),ylab="",xlab=expression(delta))
for (j in 1:ncol(MS)){
  ev <- c()
  for (k in 1:length(delta)){
    ev <- c(ev,cov_tot(MS,j,delta[k]))
  }
  EV[,j] <- ev
  lines(delta, ev, col=j, type = "l", lwd=4, lty=j)
}
legend("right", legend = colnames(EV), col = 1:ncol(EV), lty = 1:4)
# Parece que el parametro t_lag no tiene correlacion con el resto

for (i in 1:4){
  cat(paste(colnames(EV)[i],":", cov_tot(MS,i,0.25)),"\n")
}

#### Norma columnas ####
for (i in 1:4){
  cat(paste(colnames(EV)[i],":", sum(MS[,i]^2),"\n"))
}

#### PCA ####

pca <- prcomp(MS, scale. = TRUE)
# Scree plot
plot(pca, type = "l", main = "Scree Plot")
# Biplot de los primeros dos componentes principales
biplot(pca, main = "Biplot de PCA")

library(ggplot2)
# Extraer los puntajes de los primeros dos componentes
pca_data <- as.data.frame(pca$x)

# Extraer las cargas (loadings)
loadings <- as.data.frame(pca$rotation)
# Gráfico de las cargas de las variables
ggplot(loadings, aes(x = PC1, y = PC2, label = rownames(loadings))) +
  geom_point(color = "red", size = 3) +
  geom_text(vjust = 1, hjust = 1) +
  ggtitle("Gráfico de Cargas: Variables en el espacio PC1 vs PC2") +
  theme_minimal()

# Se podria estimar que existe una relación fuerte entre las variables y_0, y_max y t_max
# de manera que se podria reducir alguna de ellas.

eigen(MS)

#### Datos 2 ####

datos <- read_excel("datos/Datos2.xlsx")

trifase <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, y_0,
              ifelse(t<=t_max, y_0+mu*(t-t_lag), y_max))
  return(z)
}

tempo <- c(0:552)
param <- c(2.161,0,307.761,7.812)
TT <- trifase(tempo,param)
plot(tempo, TT, type = "l")
points(datos$x,datos$y, col = "red")

#### Analisis de sensibilidad####

# Por filas es el tempo
# Por columnas derivadas parametros
dy_0 <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 1-(t-t_lag)/(t_max-t_lag),
              ifelse(t<=t_max, 1, 0))
  return(z)
}

dt_lag <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 0,
              ifelse(t<=t_max, mu*(t-t_max), 
                     0))
  return(z)
}

dt_max <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 0,
              ifelse(t<=t_max, mu*(t-t_lag)/(t_max-t_lag), 0))
  return(z)
}

dy_max <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- 1/(t_max-t_lag)
  z <- ifelse(t <= t_lag, 0,
              ifelse(t<=t_max, mu*(t-t_lag), 1))
  return(z)
}

N <- 1e3
tempo <- seq(min(datos$x), max(datos$x), length.out = N)
S <- matrix(nrow = N, ncol = 4)
for (i in 1:N){
  S[i,] <- c(dt_lag(tempo[i],param),
             dt_max(tempo[i],param),
             dy_0(tempo[i],param),
             dy_max(tempo[i],param))
}
colnames(S) <- c("t_lag","t_max","y_0","y_max")

####
plot(c(), c(), xlim = c(0, 552), ylim = c(-6, 2), 
     xlab = "Tiempo (min)", ylab = expression(f[theta[i]]))
for (i in 1:4) {
  lines(tempo, S[, i], col = i, lty = i, lwd = 3)
}
legend("right", legend = colnames(S), col = 1:4, lty = 1:4, lwd = 3)

#### Estudio correlacion ####

MS <- t(S)%*%S
heatmap(cor(MS))

delta <- seq(0,1,length.out=100)
cov_tot <- function(MS,ci, delta){
  s <- 0
  CC <- cor(MS)
  for (i in 1:ncol(MS)){
    if (i != ci){
      p <- CC[i,ci]
      s <- s + p*(p >= 1-delta)
    }
  }
  return(s)
}

EV <- matrix(nrow=length(delta),ncol = ncol(MS))
colnames(EV) <- colnames(MS)
#png("Figures/Pproyection.png", width = 2100, height = 2970, res = 300)
plot(c(),c(),ylim=c(0,2),xlim=c(0,1),ylab="",xlab=expression(delta))
for (j in 1:ncol(MS)){
  ev <- c()
  for (k in 1:length(delta)){
    ev <- c(ev,cov_tot(MS,j,delta[k]))
  }
  EV[,j] <- ev
  lines(delta, ev, col=j, type = "l", lwd=4, lty=j)
}
legend("right", legend = colnames(EV), col = 1:ncol(EV), lty = 1:4)
# Parece que el parametro t_lag no tiene correlacion con el resto

for (i in 1:4){
  cat(paste(colnames(EV)[i],":", cov_tot(MS,i,0.25)),"\n")
}

#### Norma columnas ####
for (i in 1:4){
  cat(paste(colnames(EV)[i],":", sum(MS[,i]^2),"\n"))
}

#### PCA ####

pca <- prcomp(MS, scale. = TRUE)
# Scree plot
plot(pca, type = "l", main = "Scree Plot")
# Biplot de los primeros dos componentes principales
biplot(pca, main = "Biplot de PCA")

library(ggplot2)
# Extraer los puntajes de los primeros dos componentes
pca_data <- as.data.frame(pca$x)

# Extraer las cargas (loadings)
loadings <- as.data.frame(pca$rotation)
# Gráfico de las cargas de las variables
ggplot(loadings, aes(x = PC1, y = PC2, label = rownames(loadings))) +
  geom_point(color = "red", size = 3) +
  geom_text(vjust = 1, hjust = 1) +
  ggtitle("Gráfico de Cargas: Variables en el espacio PC1 vs PC2") +
  theme_minimal()

# Se podria estimar que existe una relación fuerte entre las variables y_0, y_max y t_max
# de manera que se podria reducir alguna de ellas.

eigen(MS)


