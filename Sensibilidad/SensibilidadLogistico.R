library(readxl)
library(ggplot2)

##### Datos 1 ####
{
datos <- read_excel("datos/datos.xlsx")

logistico <- function(t, param){
  A <- param[1]
  B <- param[2]
  C <- param[3]
  M <- param[4]
  z <- A + C/(1+exp(-B*(t-M)))
  return(z)
}

par(mfrow = c(1,1))
tempo <- c(0:552)
param <- c(1.1592,0.0112,6.9267,148.5352)
TT <- logistico(tempo,param)
plot(tempo, TT, type = "l", ylim=c(0,max(datos$y)+1))
points(datos$x,datos$y, col = "red")

#### Analisis de sensibilidad####

dA <- function(t, param){
  A <- param[1]
  B <- param[2]
  C <- param[3]
  M <- param[4]
  z <- rep(1, length(t))
  return(z)
}

dB <- function(t, param){
  A <- param[1]
  B <- param[2]
  C <- param[3]
  M <- param[4]
  z <- C*(t-M)*exp(-B*(t-M))/((1+exp(-B*(t-M)))^2)
  return(z)
}

dC <- function(t, param){
  A <- param[1]
  B <- param[2]
  C <- param[3]
  M <- param[4]
  z <- 1/(1+exp(-B*(t-M)))
  return(z)
}

dM <- function(t, param){
  A <- param[1]
  B <- param[2]
  C <- param[3]
  M <- param[4]
  z <- -B*C*exp(-B*(t-M))/((1+exp(-B*(t-M)))^2)
  return(z)
}

N <- 1e3
tempo <- seq(min(datos$x), max(datos$x), length.out = N)
S <- matrix(nrow = N, ncol = length(param))
for (i in 1:N){
  S[i,] <- c(dA(tempo[i],param),
             dB(tempo[i],param),
             dC(tempo[i],param),
             dM(tempo[i],param))
}
colnames(S) <- c("A","B","C","M")

####
par(mar = c(4, 4, 2, 1))
plot(c(), c(), xlim = c(0, 552), ylim = c(0, 150), 
     xlab = "Tiempo (min)", ylab = expression(f[theta[i]]))
for (i in 1:4) {
  lines(tempo, S[, i], col = i, lty = i, lwd = 3)
}
legend("left", legend = colnames(S), col = 1:4, lty = 1:4, lwd = 3)

# Ajustar la ventana gráfica para la subgráfica
par(fig = c(0.55, 1, 0.25, 1), new = TRUE)

# Subgráfica
plot(c(), c(), xlim = c(0, 552), ylim = c(0, 1),
     xlab = "", ylab = "", axes = FALSE)
rect(0, 0, 552, 1, col = "white", border = NA)
box()
axis(1)
axis(2)
for (i in 1:4) {
  lines(tempo, S[, i], col = i, lty = i, lwd = 3)
}

#### Estudio correlacion ####

MS <- t(S)%*%S
heatmap(cor(MS))
cor(MS)

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

#### Autovalores ####
eigen(MS)

#### PCA ####
pca <- prcomp(MS, scale. = TRUE)
# Scree plot
plot(pca, type = "l", main = "Scree Plot")
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
}

##### Datos 2 ####
{
  datos <- read_excel("datos/Datos2.xlsx")
  
  logistico <- function(t, param){
    A <- param[1]
    B <- param[2]
    C <- param[3]
    M <- param[4]
    z <- A + C/(1+exp(-B*(t-M)))
    return(z)
  }
  
  par(mfrow = c(1,1))
  tempo <- c(0:max(datos$x))
  param <- c(1.8342,0.015,6.0921,168.1591)
  TT <- logistico(tempo,param)
  plot(tempo, TT, type = "l", ylim=c(0,max(datos$y)+1))
  points(datos$x,datos$y, col = "red")
  
  #### Analisis de sensibilidad####
  
  dA <- function(t, param){
    A <- param[1]
    B <- param[2]
    C <- param[3]
    M <- param[4]
    z <- rep(1, length(t))
    return(z)
  }
  
  dB <- function(t, param){
    A <- param[1]
    B <- param[2]
    C <- param[3]
    M <- param[4]
    z <- C*(t-M)*exp(-B*(t-M))/((1+exp(-B*(t-M)))^2)
    return(z)
  }
  
  dC <- function(t, param){
    A <- param[1]
    B <- param[2]
    C <- param[3]
    M <- param[4]
    z <- 1/(1+exp(-B*(t-M)))
    return(z)
  }
  
  dM <- function(t, param){
    A <- param[1]
    B <- param[2]
    C <- param[3]
    M <- param[4]
    z <- -B*C*exp(-B*(t-M))/((1+exp(-B*(t-M)))^2)
    return(z)
  }
  
  N <- 1e3
  tempo <- seq(min(datos$x), max(datos$x), length.out = N)
  S <- matrix(nrow = N, ncol = length(param))
  for (i in 1:N){
    S[i,] <- c(dA(tempo[i],param),
               dB(tempo[i],param),
               dC(tempo[i],param),
               dM(tempo[i],param))
  }
  colnames(S) <- c("A","B","C","M")
  ####
  par(mar = c(4, 4, 2, 1))
  plot(c(), c(), xlim = c(0, 552), ylim = c(0, 100), 
       xlab = "Tiempo (min)", ylab = expression(f[theta[i]]))
  for (i in 1:4) {
    lines(tempo, S[, i], col = i, lty = i, lwd = 3)
  }
  legend("left", legend = colnames(S), col = 1:4, lty = 1:4, lwd = 3)
  
  # Ajustar la ventana gráfica para la subgráfica
  par(fig = c(0.55, 1, 0.25, 1), new = TRUE)
  
  # Subgráfica
  plot(c(), c(), xlim = c(0, 552), ylim = c(0, 1),
       xlab = "", ylab = "", axes = FALSE)
  rect(0, 0, 552, 1, col = "white", border = NA)
  box()
  axis(1)
  axis(2)
  for (i in 1:4) {
    lines(tempo, S[, i], col = i, lty = i, lwd = 3)
  }
  
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
  #### PCA ####
  pca <- prcomp(MS, scale. = TRUE)
  # Scree plot
  plot(pca, type = "l", main = "Scree Plot")
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
  
  print(MS)
  print(cor(MS))
  print(sort(colMeans(EV)))
  eigen(MS)
}




