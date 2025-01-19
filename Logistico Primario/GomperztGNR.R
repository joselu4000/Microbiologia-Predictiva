# Gomperzt Datos 2

# Cargar librerías necesarias
library(dplyr)
library(ggplot2)
library(pracma)

# Leer los datos desde el archivo Excel
library(readxl)
file_path <- "datos/Datos2.xlsx"
sdf <- read_excel(file_path)
df <- sdf[16:nrow(sdf), c(1, 2)]  # Filtrar filas y columnas necesarias
colnames(df) <- c("x", "y")  # Renombrar columnas
df <- na.omit(df)  # Eliminar valores faltantes

# Convertir columnas a numéricas
df <- df %>% mutate(x = as.numeric(x), y = as.numeric(y))

# Asignar los valores de las columnas 'x' (tiempo) y 'y' (log de las concentraciones)
dat0 <- df$x  
y_vals <- df$y  

# Definir la función logística
logistic_function <- function(params, x) {
  A <- params[1]
  B <- params[2]
  C <- params[3]
  M <- params[4]
  A + C / (1 + exp(-B * (x - M)))
}

# Derivadas parciales de la función logística
logistic_jacobian <- function(params, x) {
  A <- params[1]
  B <- params[2]
  C <- params[3]
  M <- params[4]
  
  exp_term <- exp(-B * (x - M))
  denom <- (1 + exp_term)^2
  
  dA <- rep(1, length(x))  # Derivada con respecto a A
  dB <- C * (x - M) * exp_term / denom  # Derivada con respecto a B
  dC <- 1 / (1 + exp_term)  # Derivada con respecto a C
  dM <- -C * B * exp_term / denom  # Derivada con respecto a M
  
  cbind(dA, dB, dC, dM)
}

# Inicialización de parámetros
params <- c(A = 1, B = 0.011, C = 10, M = 50)  # Valores iniciales
tolerancia <- 1e-6  # Tolerancia para el cambio de parámetros
max_iter <- 1000  # Número máximo de iteraciones

# Gráfica inicial con los datos originales
plot(dat0, y_vals, col = "red", ylim = c(0, max(y_vals) + 1), lwd = 2,
     xlab = "Tiempo (min)",
     ylab = expression(log(ufc/g)),
     main = "Gompertz GN Datos 2")

iteracion <- 0
convergencia <- FALSE

while (iteracion < max_iter && !convergencia) {
  # Calcular las predicciones actuales
  pred <- logistic_function(params, dat0)
  
  # Calcular la matriz Jacobiana
  jacobian <- logistic_jacobian(params, dat0)
  
  # Calcular la diferencia entre los valores reales y las predicciones
  residuals <- y_vals - pred
  
  # Calcular el paso de actualización usando Gauss-Newton
  jacobian_t <- t(jacobian)
  hessian_approx <- jacobian_t %*% jacobian + diag(1e-8,4,4)
  gradient <- jacobian_t %*% residuals
  delta <- solve(hessian_approx, gradient)
  
  # Actualizar los parámetros
  params <- params + as.numeric(delta)
  
  # Añadir la curva ajustada a la gráfica
  pred_line <- logistic_function(params, seq(min(dat0), max(dat0), length.out = 500))
  
  
  # Condición de convergencia
  if (sqrt(sum(delta^2)) < tolerancia) {
    convergencia <- TRUE
  }
  
  iteracion <- iteracion + 1
}
lines(seq(min(dat0), max(dat0), length.out = 500), pred_line, col = "black", lwd = 2)

# Calcular el error cuadrático medio (MSE)
interpolated <- logistic_function(params, dat0)
MSE <- mean((interpolated - y_vals)^2)

# Calcular R^2 ajustado
y_pred <- logistic_function(params, dat0)
SSE <- sum((y_vals - y_pred)^2)
SST <- sum((y_vals - mean(y_vals))^2)
R2 <- 1 - (SSE / SST)
n <- length(y_vals)
k <- 2
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
cat(sprintf("Coeficiente de determinación (R^2): %.4f\n", R2))
cat(sprintf("Coeficiente de determinación ajustado (R^2 ajustado): %.4f\n", R2_adj))

# Mostrar resultados finales
cat(sprintf("Convergencia alcanzada en la iteración %d\n", iteracion))
cat(sprintf("Parámetros optimizados: A = %.4f, B = %.4f, C = %.4f, M = %.4f, MSE =%.4f\n",
            params[1], params[2], params[3], params[4],MSE))
