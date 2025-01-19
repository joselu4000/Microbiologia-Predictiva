library(readxl)

# Definir las funciones A y f
A <- function(x, q, nu) {
  x + (1 / nu) * log((exp(-nu * x) + q) / (1 + q))
}

f <- function(x, q, mumax, nu, m, y0, ymax) {
  a <- A(x, q, nu)
  y0 + mumax * a - (1 / m) * log(1 + (exp(m * mumax * a) - 1) / exp(m * (ymax - y0)))
}

# Función de pérdida
error_function <- function(params, x_data, y_data) {
  q <- params[1]
  mumax <- params[2]
  nu <- params[3]
  m <- params[4]
  y0 <- params[5]
  ymax <- params[6]
  
  y_pred <- f(x_data, q, mumax, nu, m, y0, ymax)
  loss <- mean((y_pred - y_data)^2)
  return(loss)
}

# Implementación de Adam
adam_optimization <- function(initial_params, x_data, y_data, lr = 0.01, beta1 = 0.9, beta2 = 0.999, eps = 1e-8, epochs = 5000) {
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
    
    # Mostrar progreso cada 500 épocas
    if (epoch %% 500 == 0) {
      loss <- error_function(params, x_data, y_data)
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
x_data <- data$x
y_data <- data$y

# Parámetros iniciales
initial_params <- c(4.7371, 0.0171, 4.6965, 1.7274, 2.0891, 7.8764)

# Optimización
optimized_params <- adam_optimization(initial_params, x_data, y_data)

# Mostrar resultados optimizados
q_opt <- optimized_params[1]
mumax_opt <- optimized_params[2]
nu_opt <- optimized_params[3]
m_opt <- optimized_params[4]
y0_opt <- optimized_params[5]
ymax_opt <- optimized_params[6]

# Generar puntos para la curva ajustada
x_fit <- seq(min(x_data), max(x_data), length.out = 100)
y_fit <- f(x_fit, q_opt, mumax_opt, nu_opt, m_opt, y0_opt, ymax_opt)
# Crear la gráfica
plot(x_data, y_data, col = "red", lwd=2,
     xlab = "Tiempo (min)",
     ylim = c(0, max(y_data) + 1),
     ylab = expression(log(ufc/g)),
     main = "Baranyi Adam Datos 1")
lines(x_fit, y_fit, col = "black", lwd = 2)


# Calcular R^2 ajustado
y_pred <- f(x_data, q_opt, mumax_opt, nu_opt, m_opt, y0_opt, ymax_opt)
SSE <- sum((y_data - y_pred)^2)
SST <- sum((y_data - mean(y_data))^2)
R2 <- 1 - (SSE / SST)
n <- length(y_data)
k <- length(optimized_params)
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)

cat(sprintf("Parámetros optimizados:\n q = %.4f, mumax = %.4f, nu = %.4f, m = %.4f, y0 = %.4f, ymax = %.4f\n", 
            q_opt, mumax_opt, nu_opt, m_opt, y0_opt, ymax_opt))

# Calcular el MSE
MSE <- mean((y_pred - y_data)^2)
cat(sprintf("Error Cuadrático Medio (MSE): %.4f\n", MSE))
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))


