library(pracma) # Para funciones matemáticas avanzadas
library(Deriv) # Para derivadas simbólicas

# Leer datos desde el archivo Excel
library(readxl)
df <- read_excel("datos/datos.xlsx")

# Asignar los valores de las columnas 'x' (tiempo) y 'y' (log de las concentraciones)
dat0 <- df$x # Tiempo
y_vals <- df$y # Log de las concentraciones

# Inicialización de parámetros
q0 <- 2
mumax0 <- 0.0170
nu0 <- 0.5
m0 <- 2
y00 <- min(y_vals)
ymax0 <- max(y_vals)

# Definir epsilon para evitar valores de log negativos
epsilon <- 1e-6

# Función logística de Baranyi
A <- function(x, q, nu) {
  x + (1/nu) * log((exp(-nu * x) + q) / (1 + q))
}

f <- function(q, mumax, nu, m, y0, ymax, x) {
  a <- A(x, q, nu)
  y0 + mumax * a - (1/m) * log(1 + (exp(m * mumax * a) - 1) / exp(m * (ymax - y0)))
}

# Derivadas de f respecto a cada parámetro
f_q <- Deriv(f, "q")
f_mumax <- Deriv(f, "mumax")
f_nu <- Deriv(f, "nu")
f_m <- Deriv(f, "m")
f_y0 <- Deriv(f, "y0")
f_ymax <- Deriv(f, "ymax")

# Tolerancia y número máximo de iteraciones
tolerancia <- 1e-6
max_iter <- 100

# Graficar datos originales
plot(dat0, y_vals, col = "red", lwd = 2, 
     xlab = "Tiempo (min)",
     ylim = c(0, max(y_data) + 1),
     ylab = expression(log(ufc/g)),
     main = "Baranyi GN Datos 1")

# Iterar para ajustar los parámetros
iteracion <- 0

while (iteracion < max_iter) {
  # Evaluar la función f en los datos con los parámetros actuales
  f_vals <- sapply(dat0, function(x) f(q0, mumax0, nu0, m0, y00, ymax0, x))
  
  # Reemplazar NaN e infinitos
  f_vals[is.na(f_vals) | is.infinite(f_vals)] <- 0
  
  # Evaluar derivadas en los datos
  f_q_vals <- sapply(dat0, function(x) f_q(q0, mumax0, nu0, m0, y00, ymax0, x))
  f_mumax_vals <- sapply(dat0, function(x) f_mumax(q0, mumax0, nu0, m0, y00, ymax0, x))
  f_nu_vals <- sapply(dat0, function(x) f_nu(q0, mumax0, nu0, m0, y00, ymax0, x))
  f_m_vals <- sapply(dat0, function(x) f_m(q0, mumax0, nu0, m0, y00, ymax0, x))
  f_y0_vals <- sapply(dat0, function(x) f_y0(q0, mumax0, nu0, m0, y00, ymax0, x))
  f_ymax_vals <- sapply(dat0, function(x) f_ymax(q0, mumax0, nu0, m0, y00, ymax0, x))
  
  # Matriz de derivadas (Jacobiana)
  F0 <- cbind(f_q_vals, f_mumax_vals, f_nu_vals, f_m_vals, f_y0_vals, f_ymax_vals)
  
  # Calcular (F0' * F0)
  F0_t <- t(F0)
  F0_p <- F0_t %*% F0
  
  # Regularización dinámica
  regularization_factor <- 1e-6
  while (kappa(F0_p) > 1e12) { # Condición numérica para asegurar invertibilidad
    regularization_factor <- regularization_factor * 10
    F0_p <- F0_p + regularization_factor * diag(ncol(F0_p))
  }
  
  # Inversa regularizada
  F0_inv <- solve(F0_p)
  
  # Calcular la diferencia entre valores reales y actuales
  Y0 <- y_vals - f_vals
  
  # Delta: vector de actualización de parámetros
  delta <- F0_inv %*% (F0_t %*% Y0)
  
  # Actualizar parámetros
  q0 <- q0 + delta[1]
  mumax0 <- mumax0 + delta[2]
  nu0 <- nu0 + delta[3]
  m0 <- m0 + delta[4]
  y00 <- y00 + delta[5]
  ymax0 <- ymax0 + delta[6]
  
  # Verificar condición de parada
  if (norm(delta, type = "2") < tolerancia) {
    message("Convergencia alcanzada en la iteración ", iteracion)
    break
  }
  
  iteracion <- iteracion + 1
}

# Añadir curva ajustada
lines(dat0, f_vals, col = "black", lwd = 2)

MSE <- mean((y_vals-f_vals)^2)

# Mostrar resultados
# Calcular R^2 ajustado
y_predT <- f_vals
SSE <- sum((df$y - y_predT)^2)
SST <- sum((df$y - mean(df$y))^2)
R2 <- 1 - (SSE / SST)
n <- length(df$y)
k <- 6
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)

cat(sprintf("Parámetros ajustados: q0=%.4f, mumax0=%.4f, nu0=%.4f, m0=%.4f, y00=%.4f, ymax0=%.4f, MSE=%.4f \n", 
            q0, mumax0, nu0, m0, y00, ymax0, MSE))
cat("R^2:",R2,"R^2adj:",R2_adj)


