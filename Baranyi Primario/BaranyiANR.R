# Logistico Modificado Dataos 2

# Cargar librerías necesarias
library(readxl)
library(deSolve)
library(stats)

# Leer los datos desde el archivo Excel
file_path <- "datos/Datos2.xlsx"
sdf <- read_excel(file_path)
df <- sdf[16:nrow(sdf), c(1, 2)]  # Filtrar filas y columnas necesarias
colnames(df) <- c("x", "y")  # Renombrar columnas
df <- na.omit(df)  # Eliminar valores faltantes
x <- df$x
y <- df$y
df <- aggregate(y ~ x, data = df, FUN = mean)

# Convertir columnas a numéricas
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)

# Deslogaritmizar las concentraciones
t_data <- df$x
y_data <- df$y

# Definir el modelo diferencial
model <- function(t, state, parameters) {
  y <- state[1]
  mu <- parameters["mu"] # mu 
  ymax <- parameters["ymax"] # xmax
  m <- parameters["m"] # m
  nu <- parameters["nu"] # nu
  q0 <- parameters["q0"] # q0
  
  a <- q0/(q0+exp(-nu*t))
  term1 <- max(1 - (y / ymax)^m, 0)
  
  dy <- mu *  a *  y * term1 
  list(c(dy))
}

# Definir la función de error con escala logarítmica
error_function <- function(params, t_data, y_data) {
  y0 <- params[1]
  mu <- params[2]
  ymax <- params[3]
  m <- params[4]
  nu <- params[5]
  q0 <- params[6]
  
  parameters <- c(mu = mu, ymax = ymax, m = m, nu = nu, q0 = q0)
  state <- c(y = y0)
  
  # Resolver el sistema diferencial
  times <- seq(min(t_data), max(t_data), length.out = max(100,length(t_data)))
  sol <- ode(y = state, times = times, func = model, parms = parameters, 
             method = "ode45")
  
  # Interpolar para coincidir con los datos
  interpolated <- approx(sol[, "time"], sol[, "y"], xout = t_data)$y
  error <- mean((interpolated - y_data)^2)
  return(error)
}

# Parámetros iniciales para la optimización
initial_params <- c(min(y_data), 0.01, max(y_data), 1.5, 0.5, 2)
# q0 no aporta nada

# Ajustar los límites de los parámetros
lower_bounds <- c(1e-4, 0.001, 7, 1, 1e-4, 1)
upper_bounds <- c(4, 1, 9, 3, 3, 500)

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
# mu = mu, ymax = ymax, m = m, nu = nu, q0 = q0
y0_opt <- result$par[1]
mu_opt <- result$par[2]
ymax_opt <- result$par[3]
m_opt <- result$par[4]
nu_opt <- result$par[5]
q0_opt <- result$par[6]

# Generar solución con mayor resolución
parameters <- c(mu = mu_opt, ymax = ymax_opt, m = m_opt, nu = nu_opt, q0 = q0_opt)
state <- c(y = y0_opt)
times <- seq(min(t_data), max(t_data), length.out = 100)
sol <- ode(y = state, times = times, func = model, parms = parameters, 
           method = "ode45")
sol <- as.data.frame(sol)

# Graficar los datos y el modelo ajustado
plot(x,y, type="p", col="red",lwd=2,
     ylim = c(0, max(y_data) + 1),
     xlab = "Tiempo (min)",
     ylab = expression(log(ufc/g)),
     main = "Baranyi ode-L-BFGS-B Datos 2")
lines(sol$time,sol$y, col = "black", lwd = 2)

# Calcular el error cuadrático medio (MSE)
interpolated <- approx(sol$time, sol$y, xout = t_data)$y
MSE <- mean((interpolated - y_data)^2)

# Mostrar resultados
# Calcular R^2 ajustado
y_predT <- interpolated
SSE <- sum((df$y - y_predT)^2)
SST <- sum((df$y - mean(df$y))^2)
R2 <- 1 - (SSE / SST)
n <- length(sdf$y) # el original
k <- 6
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)

# Imprimir los parámetros optimizados
cat(sprintf("Parámetros optimizados:\n y0 = %.3f, mu = %.3f, ymax = %.3f, m = %.3f, nu = %.3f, q0 = %.3f, MSE = %.3f\n", 
            y0_opt, mu_opt, ymax_opt, m_opt, nu_opt, q0_opt, MSE))
cat("R^2:",R2,"R^2adj:",R2_adj)
