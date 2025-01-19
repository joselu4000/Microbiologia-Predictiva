library(pracma) # Para funciones matemáticas avanzadas
library(Deriv) # Para derivadas simbólicas
library(readxl) # Para leer archivos Excel

# Leer los datos desde el archivo Excel
file_path <- "datos/Datos2.xlsx"
sdf <- read_excel(file_path)
dff <- sdf[16:nrow(sdf), 1:2] # Filtrar filas y columnas necesarias
colnames(dff) <- c("x", "y")  # Renombrar columnas
df <- aggregate(y ~ x, data = dff, FUN = mean)

# Asignar los valores de las columnas 'x' (tiempo) y 'y' (log de las concentraciones)
dat0 <- df$x  # Tiempo (x)
y_vals <- df$y  # Log de las concentraciones (y)

# Inicialización de parámetros
mumax0 <- 0.01696
y00 <- 2.0928
ymax0 <- 7.8374
tlambda0 <- 5  # Tiempo lambda

# Definir límites razonables para los parámetros
param_limits <- list(
  mumax0 = c(0.001, 5),
  y00 = c(0, 10),
  ymax0 = c(0, 10),
  tlambda0 = c(0,10)
)

# Tolerancia y número máximo de iteraciones
tolerancia <- 1e-6
max_iter <- 100

# Función modificada para la ecuación
f <- function(mumax, y0, ymax, tlambda, t) {
  ifelse(t <= tlambda,
         y0,
         ymax + mumax * (t - tlambda) - log(
           exp(ymax - y0) - 1 + exp(mumax * (t - tlambda))
         )
  )
}

# Derivadas de f respecto a cada parámetro
f_mumax <- Deriv(f, "mumax")
f_y0 <- Deriv(f, "y0")
f_ymax <- Deriv(f, "ymax")
f_tlambda <- Deriv(f,"tlambda")

# Iterar para ajustar los parámetros usando Gauss-Newton
iteracion <- 0
while (iteracion < max_iter) {
  # Evaluar la función f en los datos con los parámetros actuales
  f_vals <- sapply(dat0, function(t) f(mumax0, y00, ymax0, tlambda0, t))
  f_vals[is.na(f_vals) | is.infinite(f_vals)] <- 0
  
  # Evaluar derivadas en los datos
  f_mumax_vals <- sapply(dat0, function(t) f_mumax(mumax0, y00, ymax0, tlambda0, t))
  f_y0_vals <- sapply(dat0, function(t) f_y0(mumax0, y00, ymax0,tlambda0, t))
  f_ymax_vals <- sapply(dat0, function(t) f_ymax(mumax0, y00, ymax0,tlambda0, t))
  f_ymax_vals <- sapply(dat0, function(t) f_ymax(mumax0, y00, ymax0, tlambda0, t))
  f_tlambda_vals <- sapply(dat0, function(t) f_ymax(mumax0, y00, ymax0, tlambda0, t))
  
  # Matriz de derivadas (Jacobiana)
  F0 <- cbind(f_mumax_vals, f_y0_vals, f_ymax_vals, f_tlambda_vals)
  
  # Calcular (F0' * F0)
  F0_t <- t(F0)
  F0_p <- F0_t %*% F0
  
  # Regularización dinámica
  regularization_factor <- 1e-6
  while (kappa(F0_p) > 1e12) {
    regularization_factor <- regularization_factor * 10
    F0_p <- F0_p + regularization_factor * diag(ncol(F0_p))
  }
  
  # Inversa regularizada
  F0_inv <- solve(F0_p)
  
  # Calcular la diferencia entre valores reales y actuales
  Y0 <- y_vals - f_vals
  
  # Delta: vector de actualización de parámetros
  delta <- F0_inv %*% (F0_t %*% Y0)
  
  # Actualizar parámetros con restricciones
  mumax0 <- max(min(mumax0 + delta[1], param_limits$mumax0[2]), param_limits$mumax0[1])
  y00 <- max(min(y00 + delta[2], param_limits$y00[2]), param_limits$y00[1])
  ymax0 <- max(min(ymax0 + delta[3], param_limits$ymax0[2]), param_limits$ymax0[1])
  tlambda0 <- max(min(tlambda0 + delta[3], param_limits$tlambda0[2]), param_limits$tlambda0[1])
  
  # Verificar condición de parada
  if (norm(delta, type = "2") < tolerancia) {
    message("Convergencia alcanzada en la iteración ", iteracion)
    break
  }
  
  iteracion <- iteracion + 1
}

# Graficar el ajuste final
plot(dff$x, dff$y, col = "red", lwd = 2, 
     xlab = "Tiempo (min)",
     ylim = c(0, max(y_vals) + 1),
     ylab = expression(log(ufc/g)),
     main = "LagExpo GN Datos 2")

tempo <- seq(from=min(dat0),to=max(dat0),length.out = 100)
ff <- sapply(tempo, function(t) f(mumax0, y00, ymax0,tlambda0, t))
lines(tempo, ff, col = "black", lwd = 2)

ff <- sapply(dff$x, function(x) f(mumax0, y00, ymax0, tlambda0, x))
MSE <- mean((dff$y-ff)^2)

# Mostrar resultados
# Calcular R^2 ajustado
y_predT <- ff
SSE <- sum((dff$y - y_predT)^2)
SST <- sum((dff$y - mean(dff$y))^2)
R2 <- 1 - (SSE / SST)
n <- length(dff$y)
k <- 4
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
cat("R^2:",R2,"R^2adj:",R2_adj)

# Mostrar resultados
cat(sprintf("Parámetros ajustados: mumax0=%.4f, y00=%.4f, ymax0=%.4f, lambda0=%.4f, MSE=%.4f \n", 
            mumax0, y00, ymax0, tlambda0, MSE))
