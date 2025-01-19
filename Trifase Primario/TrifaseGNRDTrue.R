# Cargar las librerías necesarias
library(readxl)

# Leer los datos desde el archivo Excel
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

# Definir el modelo trifásico
f <- function(x, t_lag, t_max, y0, y_max) {
  ifelse(
    x <= t_lag, y0,
    ifelse(
      x > t_lag & x <= t_max, y0 + ((y_max - y0) / (t_max - t_lag)) * (x - t_lag),
      y_max
    )
  )
}

# Derivadas parciales del modelo trifásico calculadas explícitamente
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

# Definir el Jacobiano del modelo trifásico usando las derivadas explícitas
trifase_jacobian <- function(params, x) {
  t_lag <- params[1]
  t_max <- params[2]
  y0 <- params[3]
  y_max <- params[4]
  
  jacobian <- matrix(0, nrow = length(x), ncol = 4)
  
  for (i in seq_along(x)) {
    jacobian[i, ] <- c(
      dy_0(x[i], params),
      dt_lag(x[i], params),
      dt_max(x[i], params),
      dy_max(x[i], params)
      
    )
  }
  
  return(jacobian)
}

# Inicialización de parámetros
params <- c(t_lag = 0.5, t_max = 307.761, y0 = 2.270, y_max = 7.812)
tolerancia <- 1e-6  # Tolerancia para el cambio de parámetros
max_iter <- 1000  # Número máximo de iteraciones

# Ajuste usando Gauss-Newton
iteracion <- 0
convergencia <- FALSE

while (iteracion < max_iter && !convergencia) {
  # Calcular las predicciones actuales
  pred <- sapply(x_data, function(x) f(x, params[1], params[2], params[3], params[4]))
  
  # Calcular la matriz Jacobiana
  jacobian <- trifase_jacobian(params, x_data)
  
  # Calcular los residuos
  residuals <- y_data - pred
  
  # Calcular el paso de actualización usando Gauss-Newton
  jacobian_t <- t(jacobian)
  hessian_approx <- jacobian_t %*% jacobian + diag(1e-12, 4, 4)  # Regularización
  gradient <- jacobian_t %*% residuals
  delta <- solve(hessian_approx, gradient)
  
  # Actualizar los parámetros
  params <- params + as.numeric(delta)
  
  # Condición de convergencia
  if (sqrt(sum(delta^2)) < tolerancia) {
    convergencia <- TRUE
  }
  if(params[1] < 0){
    params[1] <- 0
  }
  
  iteracion <- iteracion + 1
}

# Generar predicciones finales
pred_line <- sapply(seq(min(x_data), max(x_data), length.out = 500), 
                    function(x) f(x, params[1], params[2], params[3], params[4]))

# Crear la gráfica final
plot(x_data, y_data, col = "red", ylim = c(0, max(y_data) + 1), lwd = 2,
     xlab = "Tiempo (min)",
     ylab = expression(log(ufc/g)),
     main = "Trifase GN Datos 2")
lines(seq(min(x_data), max(x_data), length.out = 500), pred_line, col = "black", lwd = 2)
lines(rep(params[1],2),c(0,params[3]),lty = 2, col = "blue")
lines(rep(params[2],2),c(0,params[4]),lty = 2, col = "blue")

# Calcular el error cuadrático medio (MSE)
MSE <- mean((y_data - sapply(x_data, function(x) f(x, params[1], params[2], params[3], params[4])))^2)

# Calcular R^2 ajustado
SSE <- sum((y_data - pred)^2)
SST <- sum((y_data - mean(y_data))^2)
R2 <- 1 - (SSE / SST)
n <- length(y_data)
k <- length(params)
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)

# Mostrar resultados finales
cat(sprintf("Convergencia alcanzada en la iteración %d\n", iteracion))
cat(sprintf("Parámetros optimizados: t_lag = %.4f, t_max = %.4f, y0 = %.4f, y_max = %.4f\n",
            params[1], params[2], params[3], params[4]))
cat(sprintf("Error Cuadrático Medio (MSE): %.4f\n", MSE))
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))
