library(readxl)
library(deSolve)

# Definir el modelo para trabajar directamente con los valores originales de y
model <- function(t, state, parameters) {
  y <- state[1]
  y0 <- parameters["y0"]
  t_lag <- parameters["t_lag"]
  t_max <- parameters["t_max"]
  y_max <- parameters["y_max"]
  
  # Calcular mu
  mu <- (y_max - y0) / (t_max - t_lag)
  
  # Definir las fases del modelo
  if (t <= t_lag) {
    dy <- 0  # Fase 1: Constante
  } else if (t > t_lag && t <= t_max) {
    dy <- mu  # Fase 2: Crecimiento lineal
  } else {
    dy <- 0  # Fase 3: Constante en y_max
  }
  
  list(c(dy))
}

# Función de pérdida para ajuste de parámetros
error_function <- function(parameters, times, y_data) {
  state <- c(y = y_data[1])  # Condición inicial
  out <- ode(y = state, times = times, func = model, parms = parameters)
  
  y_pred <- out[, "y"]
  loss <- mean((y_pred - y_data)^2)
  return(loss)
}

# Implementación de Adam
adam_optimization <- function(initial_params, times, y_data, lr = 0.01, 
                              beta1 = 0.9, beta2 = 0.999, 
                              eps = 1e-8, epochs = 5000) {
  params <- initial_params
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
      
      grad[i] <- (error_function(params_plus, times, y_data) - error_function(params_minus, times, y_data)) / (2 * h)
    }
    
    # Actualizar momentos
    m <- beta1 * m + (1 - beta1) * grad
    v <- beta2 * v + (1 - beta2) * (grad^2)
    
    # Corregir sesgos
    m_hat <- m / (1 - beta1^t)
    v_hat <- v / (1 - beta2^t)
    
    # Actualizar parámetros
    params <- params - lr * m_hat / (sqrt(v_hat) + eps)
    
    # Mostrar progreso cada 500 épocas
    if (epoch %% 500 == 0) {
      loss <- error_function(params, times, y_data)
      cat(sprintf("Epoch %d: Loss = %.6f\n", epoch, loss))
    }
  }
  
  return(params)
}

# Lectura de datos
file_path <- "datos/Datos2.xlsx"
sdf <- read_excel(file_path)
df <- sdf[16:nrow(sdf), c(1, 2)]  # Filtrar filas y columnas necesarias
colnames(df) <- c("x", "y")  # Renombrar columnas
df <- na.omit(df)  # Eliminar valores faltantes

# Convertir columnas a numéricas
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)

# Extraer datos
times <- df$x
y_data <- df$y

# Parámetros iniciales nombrados
initial_params <- c(y0 = 2, t_lag = 1, t_max = 320, y_max = 7.825)

# Optimización
optimized_params <- adam_optimization(initial_params, times, y_data)

# Resolver la EDO con parámetros optimizados
state <- c(y = y_data[1])
out <- ode(y = state, times = times, func = model, parms = optimized_params)

# Crear la gráfica
plot(x, y, col = "red", lwd = 2, 
     xlab = "Tiempo (min)",
     ylim = c(0, max(y_data) + 1),
     ylab = expression(log(ufc/g)),
     main = "Trifase ode-Adam Datos 2")
#lines(out[, "time"], out[, "y"], col = "black", lwd = 2)
segments(0, y_data[1], optimized_params["t_lag"], y_data[1], col = "black", lwd = 2)  # Fase 1
segments(optimized_params["t_lag"], y_data[1], optimized_params["t_max"], optimized_params["y_max"], col = "black", lwd = 2)  # Fase 2
segments(optimized_params["t_max"], optimized_params["y_max"], max(times), optimized_params["y_max"], col = "black", lwd = 2)  # Fase 3
lines(rep(optimized_params["t_lag"], 2), c(0, y_data[1]), lty = 2, col = "blue")
lines(rep(optimized_params["t_max"], 2), c(0, optimized_params["y_max"]), lty = 2, col = "blue")


# Calcular R^2 ajustado
y_pred <- out[, "y"]
SSE <- sum((y_data - y_pred)^2)
SST <- sum((y_data - mean(y_data))^2)
R2 <- 1 - (SSE / SST)
n <- length(y_data)
k <- length(optimized_params)
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
# Calcular el MSE
MSE <- mean((y_pred - y_data)^2)

# Mostrar resultados optimizados
cat("Parámetros optimizados:\n")
print(optimized_params)
cat(sprintf("Error Cuadrático Medio (MSE): %.4f\n", MSE))
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))

