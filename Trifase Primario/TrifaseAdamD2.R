library(readxl)

# Definir la función por tramos
f <- function(params,x) {
  t_lag <- params[1]
  t_max <- params[2]
  y0 <- params[3]
  y_max <- params[4]
  
  ifelse(
    x <= t_lag, y0,
    ifelse(
      x > t_lag & x <= t_max, y0 + ((y_max - y0) / (t_max - t_lag)) * (x - t_lag),
      y_max
    )
  )
}

# Función de pérdida
error_function <- function(params, x_data, y_data) {
  qA <- params[1]
  B <- params[2]
  C <- params[3]
  M <- params[4]
  
  y_pred <- f(params,x_data)
  loss <- mean((y_pred - y_data)^2)
  return(loss)
}

# Implementación de Adam
adam_optimization <- function(initial_params, x_data, y_data, lr = 0.01, 
                              beta1 = 0.9, beta2 = 0.999, 
                              eps = 1e-8, epochs = 5000) {
  params <- as.numeric(initial_params)
  m <- rep(0, length(params))
  v <- rep(0, length(params))
  t <- 0
  
  for (epoch in 1:epochs) {
    t <- t + 1
    
    # Calcular gradiente numérico
    grad <- rep(0, length(params))
    h <- 1e-5
    for (i in seq_along(params)) {
      params_plus <- params
      params_minus <- params
      params_plus[i] <- params_plus[i] + h
      params_minus[i] <- params_minus[i] - h
      
      grad[i] <- (error_function(params_plus, x_data, y_data) - error_function(params_minus, x_data, y_data)) / (2 * h)
    }
    
    # Actualizar momentos
    m <- beta1 * m + (1 - beta1) * grad
    v <- beta2 * v + (1 - beta2) * (grad^2)
    
    # Corregir sesgos
    m_hat <- m / (1 - beta1^t)
    v_hat <- v / (1 - beta2^t)
    
    # Actualizar parámetros
    params <- params - lr * m_hat / (sqrt(v_hat) + eps)
    if(params[1] < 0){
      params[1] <- 0
    }
    
    
    # Mostrar progreso cada 500 épocas
    if (epoch %% 500 == 0) {
      loss <- error_function(params, x_data, y_data)
      cat(sprintf("Epoch %d: Loss = %.6f\n", epoch, loss))
    }
  }
  
  return(params)
}

file_path <- "datos/Datos2.xlsx"
sdf <- read_excel(file_path)
df <- sdf[16:nrow(sdf), c(1, 2)]  # Filtrar filas y columnas necesarias
colnames(df) <- c("x", "y")  # Renombrar columnas
df <- na.omit(df)  # Eliminar valores faltantes

# Convertir columnas a numéricas
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)

# Extraer datos
x_data <- df$x
y_data <- df$y

# Parámetros iniciales
initial_params <- c(0, 300, 2, 7)

# Optimización
optimized_params <- adam_optimization(initial_params, x_data, y_data)

# Generar puntos para la curva ajustada
x_fit <- seq(min(x_data), max(x_data), length.out = 100)
y_fit <- f(optimized_params,x_fit)
# Crear la gráfica
plot(x_data, y_data, col = "red", lwd=2,
     xlab = "Tiempo (min)",
     ylim = c(0, max(y_data) + 1),
     ylab = expression(log(ufc/g)),
     main = "Trifase Adam Datos 1")
lines(x_fit, y_fit, col = "black", lwd = 2)
lines(rep(optimized_params[1],2),c(0,optimized_params[3]),lty = 2, col = "blue")
lines(rep(optimized_params[2],2),c(0,optimized_params[4]),lty = 2, col = "blue")


# Calcular R^2 ajustado
y_pred <- f(optimized_params,x_data)
SSE <- sum((y_data - y_pred)^2)
SST <- sum((y_data - mean(y_data))^2)
R2 <- 1 - (SSE / SST)
n <- length(y_data)
k <- length(optimized_params)
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)

cat("Parámetros optimizados: (tlag,tmax,y0,ymax)\n", optimized_params)

# Calcular el MSE
MSE <- mean((y_pred - y_data)^2)
cat(sprintf("Error Cuadrático Medio (MSE): %.4f\n", MSE))
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))


