library(shiny)
library(bslib)
library(readxl)

# Cargar datos
datos <- read_excel("datos.xlsx")

# Definir la función logística
Baranyi <- function(t, param, datos) {
  q <- param[1]
  mumax <- param[2]
  nu <- param[3]
  m <- param[4]
  A <- t + (1/nu) * log((exp(-nu * t) + q) / (1 + q))
  z <- min(datos$y) + mumax * A - (1/m) * log(1 + (exp(m * mumax * A) - 1) / (exp(m * (max(datos$y) - min(datos$y)))))
  return(list(A = A, z = z))
}

# Definir la interfaz de usuario
ui <- fluidPage(
  titlePanel("Baranyi-Robert Exploring"),
  sidebarLayout(
    sidebarPanel(
      # Input: q
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "q",
                 label = "Initial media concentration:",
                 min = 0,
                 max = 1000,
                 value = 40,
                 step = 0.01,
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "q_numeric",
                 label = NULL,
                 value = 40,
                 min = 0,
                 max = 1000,
                 width = "200%"
               ))
      ),
      # Input: Mu_max
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "Mu_max",
                 label = "Maximum slope:",
                 min = 0,
                 max = 1,
                 value = 0.01,
                 step = 0.0001,
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "Mu_max_numeric",
                 label = NULL,
                 value = 0.01,
                 min = 0,
                 max = 1,
                 width = "200%"
               ))
      ),
      # Input: Nu
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "nu",  # Cambiar 'Nu' a 'nu'
                 label = "Influence value:",
                 min = 0,
                 max = 1,
                 value = 0.5,
                 step = 0.0001,
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "nu_numeric",  # Cambiar 'Nu_numeric' a 'nu_numeric'
                 label = NULL,
                 value = 0.5,
                 min = 0,
                 max = 1,
                 width = "200%"
               ))
      ),
      # Input: m
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "m",
                 label = "Scale term:",
                 min = 0,
                 max = 10,
                 value = 1,
                 step = 1,
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "m_numeric",
                 label = NULL,
                 value = 1,
                 min = 0,
                 max = 10,
                 width = "200%"
               ))
      )
    ),
    mainPanel(
      plotOutput(outputId = "plot", height = "600px")
    )
  )
)

# Definir la lógica del servidor
server <- function(input, output, session) {
  # Sincronizar q
  observe({
    updateSliderInput(session, "q", value = input$q_numeric)
  })
  observe({
    updateNumericInput(session, "q_numeric", value = input$q)
  })
  
  # Sincronizar Mu_max
  observe({
    updateSliderInput(session, "Mu_max", value = input$Mu_max_numeric)
  })
  observe({
    updateNumericInput(session, "Mu_max_numeric", value = input$Mu_max)
  })
  
  # Sincronizar nu
  observe({
    updateSliderInput(session, "nu", value = input$nu_numeric)
  })
  observe({
    updateNumericInput(session, "nu_numeric", value = input$nu)
  })
  
  # Sincronizar m
  observe({
    updateSliderInput(session, "m", value = input$m_numeric)
  })
  observe({
    updateNumericInput(session, "m_numeric", value = input$m)
  })
  
  # Graficar
  output$plot <- renderPlot({
    par(mfrow = c(2, 1))
    
    plot(datos$x, datos$y, col = "red", type = "p", pch = 2,
         xlab = "Time", ylab = "log(CFU/g)",
         ylim = c(0, 8))
    
    param <- c(input$q, input$Mu_max, input$nu, input$m)
    TT <- Baranyi(datos$x, param, datos)
    points(datos$x, TT$z, col = "blue", type = "l")
    
    plot(datos$x, TT$A, type = "p", col = "blue",
         xlab = "Time", ylab = "A(t)")
    points(datos$x,datos$x , col = "black", type = "l")
  })
}

# Definir la aplicación Shiny
shinyApp(ui = ui, server = server)
