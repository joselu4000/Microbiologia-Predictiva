library(readxl)
library(deSolve)

# Definir el modelo diferencial
model <- function(t, state, parameters) {
  y <- state[1]
  mu <- parameters["mu"] # mu 
  ymax <- parameters["ymax"] # xmax
  m <- parameters["m"] # m
  nu <- parameters["nu"] # nu
  q0 <- parameters["q0"] # q0
  
  a <- q0 / (q0 + exp(-nu * t))
  term1 <- max(1 - (y / ymax)^m, 0)
  
  dy <- mu * a * y * term1 
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
file_path <- "datos/datos.xlsx"  # Cambia este path según sea necesario
data <- read_excel(file_path)
data <- na.omit(data)
data <- aggregate(y ~ x, data, mean)  # Filtrar y agregar datos

# Extraer datos
times <- data$x
y_data <- data$y

# Parámetros iniciales
initial_params <- c(mu = 0.02, ymax = 7.8764, m = 1.7274, nu = 4.6965, q0 = 4.7371)

# Optimización
optimized_params <- adam_optimization(initial_params, times, y_data)

# Resolver la EDO con parámetros optimizados
state <- c(y = y_data[1])
out <- ode(y = state, times = times, func = model, parms = optimized_params)

# Crear la gráfica
plot(times, y_data, col = "red", lwd = 2, 
     xlab = "Tiempo (min)",
     ylim = c(0, max(y_data) + 1),
     ylab = expression(log(ufc/g)),
     main = "Baranyi ode-Adam Datos 1")
lines(out[, "time"], out[, "y"], col = "black", lwd = 2)

# Calcular R^2 ajustado
y_pred <- out[, "y"]
SSE <- sum((y_data - y_pred)^2)
SST <- sum((y_data - mean(y_data))^2)
R2 <- 1 - (SSE / SST)
n <- length(y_data)
k <- 6
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
# Calcular el MSE
MSE <- mean((y_pred - y_data)^2)

# Mostrar resultados optimizados
cat("Parámetros optimizados:\n")
print(optimized_params)
cat(sprintf("Error Cuadrático Medio (MSE): %.4f\n", MSE))
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))

