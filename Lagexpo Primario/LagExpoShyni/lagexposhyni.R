library(shiny)
library(bslib)
library(readxl)

# Cargar datos
datos <- read_excel("datos.xlsx")

# Definir la función logística
lagexpo <- function(t, param, datos) {
  lambda <- param[1]
  mumax <- param[2]
  z <- ifelse(t <= lambda, min(datos$y), 
              max(datos$y) + mumax*(t-lambda) - log(exp(max(datos$y)-min(datos$y))-1+exp(mumax*(t-lambda))))
  return(z)
}

# Definir la interfaz de usuario
ui <- fluidPage(
  titlePanel("Logistic Exploring"),
  sidebarLayout(
    sidebarPanel(
      # Input: Lambda
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "Lambda",
                 label = "Lag Term:",
                 min = min(datos$x),
                 max = max(datos$x),
                 value = 1.5,
                 step = 0.001,  
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "Lambda_numeric",
                 label = NULL,
                 value = 1.5,
                 min = min(datos$x),
                 max = max(datos$x),
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
                 value = 0.1,
                 step = 0.0001,  
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "Mu_max_numeric",
                 label = NULL,
                 value = 0.1,
                 min = 0,
                 max = 1,
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
  # Sincronizar A
  observe({
    updateSliderInput(session, "Lambda", value = input$Lambda_numeric)
  })
  observe({
    updateNumericInput(session, "Lambda", value = input$Lambda)
  })
  
  # Sincronizar B
  observe({
    updateSliderInput(session, "Mu_max", value = input$Mu_max_numeric)
  })
  observe({
    updateNumericInput(session, "Mu_max_numeric", value = input$Mu_max)
  })
  
  
  # Graficar
  output$plot <- renderPlot({
    plot(datos$x, datos$y, col = "red", type = "p", pch = 2,
         xlab = "Time", ylab = "log(CFU/g)",
         ylim = c(0, 12))
    
    param <- c(input$Lambda, input$Mu_max)
    z <- lagexpo(datos$x, param, datos)
    points(datos$x, z, col = "blue", type = "l")
  })
}

# Definir la aplicación Shiny
shinyApp(ui = ui, server = server)
