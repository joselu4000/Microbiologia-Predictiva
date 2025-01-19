# Gomperzt

# Cargar librerías necesarias
library(readxl)
library(deSolve)
library(stats)
library(ggplot2)

# Leer los datos desde el archivo Excel
file_path <- "datos/Datos2.xlsx"
sdf <- read_excel(file_path)
df <- sdf[16:nrow(sdf), c(1, 2)]  # Filtrar filas y columnas necesarias
colnames(df) <- c("x", "y")  # Renombrar columnas
df <- na.omit(df)  # Eliminar valores faltantes

# Convertir columnas a numéricas
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)

# Deslogaritmizar las concentraciones
t_data <- df$x
y_data <- df$y

# Definir el modelo diferencial
model <- function(t, state, parameters) {
  y <- state[1]
  b1 <- parameters["b1"]
  b2 <- parameters["b2"]
  
  term1 <- max(1 - y / b2, 0)
  dy <- b1 * y * term1
  list(c(dy))
}

# Definir la función de error con escala logarítmica
error_function <- function(params, t_data, y_data) {
  y0 <- params[1]
  b1 <- params[2]
  b2 <- params[3]
  
  parameters <- c(b1 = b1, b2 = b2)
  state <- c(y = y0)
  
  # Resolver el sistema diferencial
  times <- seq(min(t_data), max(t_data), length.out = length(t_data))
  sol <- ode(y = state, times = times, func = model, parms = parameters, 
             method = "ode45")
  
  # Interpolar para coincidir con los datos
  interpolated <- approx(sol[, "time"], sol[, "y"], xout = t_data)$y
  error <- mean((interpolated - y_data)^2)
  return(error)
}

# Parámetros iniciales para la optimización
initial_params <- c(min(y_data), 0.0219, max(y_data))

# Ajustar los límites de los parámetros
lower_bounds <- c(0, 0.001, 7)
upper_bounds <- c(3, 0.5, 9)

# Ejecutar la optimización con L-BFGS-B
result <- optim(
  par = initial_params,
  fn = error_function,
  t_data = t_data,
  y_data = y_data,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(maxit = 5000)
)

# Obtener los parámetros optimizados
y0_opt <- result$par[1]
b1_opt <- result$par[2]
b2_opt <- result$par[3]

# Generar solución con mayor resolución
parameters <- c(b1 = b1_opt, b2 = b2_opt)
state <- c(y = y0_opt)
times <- seq(min(t_data), max(t_data), length.out = 100)
sol <- ode(y = state, times = times, func = model, parms = parameters, 
           method = "ode45")
sol <- as.data.frame(sol)

# Graficar los datos y el modelo ajustado
plot(t_data,y_data, type="p", col="red",lwd=2,
     ylim = c(0, max(y_data) + 1),
     xlab = "Tiempo (min)",
     ylab = expression(log(ufc/g)),
     main = "Logistico ode-L-BFGS-B Datos 2")
lines(sol$time,sol$y, col = "black", lwd = 2)

# Calcular el error cuadrático medio (MSE)
interpolated <- approx(sol$time, sol$y, xout = t_data)$y
MSE <- mean((interpolated - y_data)^2)

# Calcular R^2 ajustado
y_pred <- approx(sol$time, sol$y, xout = t_data)$y
SSE <- sum((y_data - y_pred)^2)
SST <- sum((y_data - mean(y_data))^2)
R2 <- 1 - (SSE / SST)
n <- length(y_data)
k <- 2
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))

# Imprimir los parámetros optimizados
cat(sprintf("Parámetros optimizados:\n y0 = %.3f, b1 = %.3f, b2 = %.3e, MSE =%.3e\n", 
            y0_opt, b1_opt, b2_opt, MSE))
