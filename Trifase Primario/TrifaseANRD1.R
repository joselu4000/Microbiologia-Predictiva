# Trifase AN

# Cargar librerías necesarias
library(readxl)
library(deSolve)
library(stats)
library(ggplot2)

# Leer los datos desde el archivo Excel
file_path <- "datos/datos.xlsx"
sdf <- read_excel(file_path)
df <- sdf[1:nrow(sdf), c(1, 2)]  # Filtrar filas y columnas necesarias
colnames(df) <- c("x", "y")  # Renombrar columnas
df <- na.omit(df)  # Eliminar valores faltantes

# Convertir columnas a numéricas
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)

# Deslogaritmizar las concentraciones
t_data <- df$x
y_data <- df$y

# Definir el modelo para trabajar directamente con los valores originales de y
model_direct_y <- function(t, state, parameters) {
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

# Definir la función por tramos
model_piecewise <- function(x, t_lag, t_max, y0, y_max) {
  ifelse(
    x <= t_lag, y0,
    ifelse(
      x > t_lag & x <= t_max, y0 + ((y_max - y0) / (t_max - t_lag)) * (x - t_lag),
      y_max
    )
  )
}

# Definir la función de error
error_function_direct_y <- function(params, t_data, y_data) {
  y0 <- params[1]
  t_lag <- params[2]
  t_max <- params[3]
  y_max <- params[4]
  
  parameters <- c(y0 = y0, t_lag = t_lag, t_max = t_max, y_max = y_max)
  state <- c(y = y0)
  
  # Resolver el sistema
  times <- seq(min(t_data), max(t_data), length.out = 100)
  sol <- ode(y = state, times = times, func = model_direct_y, parms = parameters, 
             method = "ode45")
  
  # Interpolar para coincidir con los datos
  interpolated <- approx(sol[, "time"], sol[, "y"], xout = t_data)$y
  
  # Calcular el error cuadrático medio
  error <- mean((interpolated - y_data)^2)
  return(error)
}

# Parámetros iniciales para la optimización
initial_params_with_y0 <- c(2, 0, 358, 7.825)

# Ajustar los límites de los parámetros
lower_bounds <- c(0,0,100,3)
upper_bounds <- c(3, 100, 500, 8)

# Ejecutar la optimización
result_direct_y <- optim(
  par = initial_params_with_y0,
  fn = error_function_direct_y,
  t_data = t_data,
  y_data = y_data,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(maxit = 5000)
)

# Obtener los parámetros optimizados
y0_opt <- result_direct_y$par[1]
t_lag_opt <- result_direct_y$par[2]
t_max_opt <- result_direct_y$par[3]
y_max_opt <- result_direct_y$par[4]

# Generar solución con mayor resolución
parameters <- c(y0 = y0_opt, t_lag = t_lag_opt, t_max = t_max_opt, y_max = y_max_opt)
state <- c(y = y0_opt)
times <- seq(min(t_data), max(t_data), length.out = 500)
sol_final_direct_y <- ode(y = state, times = times, func = model_direct_y, parms = parameters, method = "ode45")

# Graficar resultados

x <- df$x
y <- df$y
t <- seq(min(x),max(x),length.out=100)
trifase <- model_piecewise(t, t_lag_opt, t_max_opt, y0_opt, y_max_opt)
plot(x,y, type="p", col="red",lwd=2,
     xlab = "Tiempo (min)",
     ylab = expression(log(ufc/g)),
     main = "Trifase ode-L-FBGS-B Datos 1",
     ylim = c(0,max(y_data)+1))
lines(t,trifase, col = "black", lwd = 2)
lines(rep(t_lag_opt,2),c(0,y0_opt),lty = 2, col = "blue")
lines(rep(t_max_opt,2),c(0,y_max_opt),lty = 2, col = "blue")

# Calcular el error cuadrático medio (MSE)
interpolated <- approx(sol$time, sol$y, xout = t_data)$y
MSE <- mean((interpolated - y_data)^2)

# Calcular R^2 ajustado
y_pred <- approx(sol$time, sol$y, xout = t_data)$y
SSE <- sum((y - y_pred)^2)
SST <- sum((y - mean(y))^2)
R2 <- 1 - (SSE / SST)
n <- length(y)
k <- 4
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))

# Imprimir los parámetros optimizados
cat("Parámetros optimizados:\n")
cat(sprintf("y0 = %.3f, t_lag = %.3f, t_max = %.3f, y_max = %.3f, MSE = %3.f\n", y0_opt, t_lag_opt, t_max_opt, y_max_opt, MSE))
cat(MSE)
cat(sprintf("Tasa de crecimiento calculada: mu = %.6f\n", (y_max_opt - y0_opt) / (t_max_opt - t_lag_opt)))


