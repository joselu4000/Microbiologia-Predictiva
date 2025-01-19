library(readxl)
library(pracma) # Para funciones matemáticas avanzadas
library(Deriv) # Para derivadas simbólicas

##### Datos 1 ####
{
  # Leer datos desde el archivo Excel
  library(readxl)
  datos <- read_excel("datos/datos.xlsx")
  
  # Función logística de Baranyi
  A <- function(x, q, nu) {
    x + (1/nu) * log((exp(-nu * x) + q) / (1 + q))
  }
  
  f <- function(q, mumax, nu, m, y0, ymax, x) {
    a <- A(x, q, nu)
    y0 + mumax * a - (1/m) * log(1 + (exp(m * mumax * a) - 1) / exp(m * (ymax - y0)))
  }
  
  q = 4.7371
  mumax=0.0171 
  nu=4.6965 
  m=1.7274 
  y0=2.0891 
  ymax=7.8764
  
  par(mfrow = c(1,1))
  tempo <- seq(0,552,length.out=100)
  TT <- f(q = 4.7371, mumax=0.0171, nu=4.6965, m=1.7274, y0=2.0891, ymax=7.8764, tempo)
  plot(tempo, TT, type = "l", ylim=c(0,max(datos$y)+1))
  points(datos$x,datos$y, col = "red")
  
  #### Analisis de sensibilidad####
  
  f_q <- Deriv(f, "q")
  f_mumax <- Deriv(f, "mumax")
  f_nu <- Deriv(f, "nu")
  f_m <- Deriv(f, "m")
  f_y0 <- Deriv(f, "y0")
  f_ymax <- Deriv(f, "ymax")
  
  N <- 1e3
  tempo <- seq(min(datos$x), max(datos$x), length.out = N)
  S <- matrix(nrow = N, ncol = 6)
  for (i in 1:N){
    S[i,] <- c(f_q(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_mumax(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_nu(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_m(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_y0(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_ymax(q, mumax, nu, m, y0, ymax, tempo[i]))
  }
  colnames(S) <- c("q0","mumax","nu","m","y0","ymax")
  
  ####
  par(mar = c(4, 4, 2, 1))
  plot(c(), c(), xlim = c(0, 552), ylim = c(0, 250),
       xlab = "Tiempo (min)", ylab = expression(f[theta[i]]))
  for (i in 1:6) {
    lines(tempo, S[, i], col = i, lty = i, lwd = 3)
  }
  legend("left", legend = colnames(S), col = 1:6, lty = 1:6, lwd = 3)
  
  # Ajustar la ventana gráfica para la subgráfica
  par(fig = c(0.55, 1, 0.25, 1), new = TRUE)
  
  # Subgráfica
  plot(c(), c(), xlim = c(0, 552), ylim = c(0, 1),
       xlab = "", ylab = "", axes = FALSE)
  rect(0, 0, 552, 1, col = "white", border = NA)
  box()
  axis(1)
  axis(2)
  for (i in 1:6) {
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
  plot(c(),c(),ylim=c(0,5),xlim=c(0,1),ylab="",xlab=expression(delta))
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
  
  for (i in 1:6){
    cat(paste(colnames(EV)[i],":", cov_tot(MS,i,0.25)),"\n")
  }
  
  #### Norma columnas ####
  for (i in 1:6){
    cat(paste(colnames(EV)[i],":", sum(MS[,i]^2),"\n"))
  }
  
  #### Autovectores ####
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
    ggtitle("") +
    theme_minimal()
}

##### Datos 2 ####
{
  # Leer datos desde el archivo Excel
  library(readxl)
  datos <- read_excel("datos/Datos2.xlsx")
  datos <- datos[16:nrow(datos),c(1,2)]
  
  # Función logística de Baranyi
  A <- function(x, q, nu) {
    x + (1/nu) * log((exp(-nu * x) + q) / (1 + q))
  }
  
  f <- function(q, mumax, nu, m, y0, ymax, x) {
    a <- A(x, q, nu)
    y0 + mumax * a - (1/m) * log(1 + (exp(m * mumax * a) - 1) / exp(m * (ymax - y0)))
  }
  
  q = 598.6903
  mumax=0.0175 
  nu=0.1 
  m=2.8585 
  y0=2.1209 
  ymax=7.8235
  
  par(mfrow = c(1,1))
  tempo <- seq(0,max(datos$x),length.out=100)
  TT <- f(q = 4.7371, mumax=0.0171, nu=4.6965, m=1.7274, y0=2.0891, ymax=7.8764, tempo)
  plot(tempo, TT, type = "l", ylim=c(0,max(datos$y)+1))
  points(datos$x,datos$y, col = "red")
  
  #### Analisis de sensibilidad####
  
  f_q <- Deriv(f, "q")
  f_mumax <- Deriv(f, "mumax")
  f_nu <- Deriv(f, "nu")
  f_m <- Deriv(f, "m")
  f_y0 <- Deriv(f, "y0")
  f_ymax <- Deriv(f, "ymax")
  
  N <- 1e3
  tempo <- seq(min(datos$x), max(datos$x), length.out = N)
  S <- matrix(nrow = N, ncol = 6)
  for (i in 1:N){
    S[i,] <- c(f_q(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_mumax(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_nu(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_m(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_y0(q, mumax, nu, m, y0, ymax, tempo[i]),
               f_ymax(q, mumax, nu, m, y0, ymax, tempo[i]))
  }
  colnames(S) <- c("q0","mumax","nu","m","y0","ymax")  
  
  ####
  par(mar = c(4, 4, 2, 1))
  plot(c(), c(), xlim = c(0, 552), ylim = c(0, 250),
       xlab = "Tiempo (min)", ylab = expression(f[theta[i]]))
  for (i in 1:6) {
    lines(tempo, S[, i], col = i, lty = i, lwd = 3)
  }
  legend("left", legend = colnames(S), col = 1:6, lty = 1:6, lwd = 3)
  
  # Ajustar la ventana gráfica para la subgráfica
  par(fig = c(0.55, 1, 0.25, 1), new = TRUE)
  
  # Subgráfica
  plot(c(), c(), xlim = c(0, 552), ylim = c(0, 1),
       xlab = "", ylab = "", axes = FALSE)
  rect(0, 0, 552, 1, col = "white", border = NA)
  box()
  axis(1)
  axis(2)
  for (i in 1:6) {
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
  plot(c(),c(),ylim=c(0,5),xlim=c(0,1),ylab="",xlab=expression(delta))
  for (j in 1:ncol(MS)){
    ev <- c()
    for (k in 1:length(delta)){
      ev <- c(ev,cov_tot(MS,j,delta[k]))
    }
    EV[,j] <- ev
    lines(delta, ev, col=j, type = "l", lwd=4, lty=j)
  }
  legend("right", legend = colnames(EV), col = 1:ncol(EV), lty = 1:4)
  
  for (i in 1:6){
    cat(paste(colnames(EV)[i],":", cov_tot(MS,i,0.25)),"\n")
  }
  
  #### Norma columnas ####
  for (i in 1:6){
    cat(paste(colnames(EV)[i],":", sum(MS[,i]^2),"\n"))
  }
  
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
    ggtitle("") +
    theme_minimal()
  
  print(MS)
  print(cor(MS))
  print(sort(colMeans(EV)))
  eigen(MS)
}




