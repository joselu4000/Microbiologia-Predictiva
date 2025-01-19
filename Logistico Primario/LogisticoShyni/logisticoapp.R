library(shiny)
library(bslib)
library(readxl)

# Cargar datos
datos <- read_excel("datos.xlsx")

# Definir la función logística
logistico <- function(t, param) {
  A <- param[1]
  B <- param[2]
  C <- param[3]
  M <- param[4]
  z <- A + C / (1 + exp(-B * (t - M)))
  return(z)
}

# Definir la interfaz de usuario
ui <- fluidPage(
  titlePanel("Logistic Exploring"),
  sidebarLayout(
    sidebarPanel(
      # Input: A
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "A",
                 label = "Independent Term:",
                 min = 0,
                 max = 10,
                 value = 1.5,
                 step = 0.001,  # Cambiado a 0.001
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "A_numeric",
                 label = NULL,
                 value = 1.5,
                 min = 0,
                 max = 10,
                 width = "200%"
               ))
      ),
      # Input: B
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "B",
                 label = "Scale term:",
                 min = 0,
                 max = 2,
                 value = 1,
                 step = 0.001,  # Cambiado a 0.001
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "B_numeric",
                 label = NULL,
                 value = 1,
                 min = 0,
                 max = 2,
                 width = "200%"
               ))
      ),
      # Input: C
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "C",
                 label = "Dilatation Term:",
                 min = 0,
                 max = 10,
                 value = 7,
                 step = 0.001,  # Cambiado a 0.001
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "C_numeric",
                 label = NULL,
                 value = 7,
                 min = 0,
                 max = 10,
                 width = "200%"
               ))
      ),
      # Input: M
      fluidRow(
        column(8,
               sliderInput(
                 inputId = "M",
                 label = "Lag Term:",
                 min = 0,
                 max = max(datos$x),
                 value = 150,
                 step = 0.001,  # Cambiado a 0.001
                 width = "100%"
               )),
        column(4,
               numericInput(
                 inputId = "M_numeric",
                 label = NULL,
                 value = 150,
                 min = 0,
                 max = max(datos$x),
                 width = "200%"
               ))
      )
    ),
    mainPanel(
      plotOutput(outputId = "plot", height = "600px")  # Ajustar la altura de la gráfica
    )
  )
)

# Definir la lógica del servidor
server <- function(input, output, session) {
  # Sincronizar A
  observe({
    updateSliderInput(session, "A", value = input$A_numeric)
  })
  observe({
    updateNumericInput(session, "A_numeric", value = input$A)
  })
  
  # Sincronizar B
  observe({
    updateSliderInput(session, "B", value = input$B_numeric)
  })
  observe({
    updateNumericInput(session, "B_numeric", value = input$B)
  })
  
  # Sincronizar C
  observe({
    updateSliderInput(session, "C", value = input$C_numeric)
  })
  observe({
    updateNumericInput(session, "C_numeric", value = input$C)
  })
  
  # Sincronizar M
  observe({
    updateSliderInput(session, "M", value = input$M_numeric)
  })
  observe({
    updateNumericInput(session, "M_numeric", value = input$M)
  })
  
  # Graficar
  output$plot <- renderPlot({
    plot(datos$x, datos$y, col = "red", type = "p", pch = 2,
         xlab = "Time", ylab = "log(CFU/g)",
         ylim = c(0, 12))
    
    param <- c(input$A, input$B, input$C, input$M)
    z <- logistico(datos$x, param)
    points(datos$x, z, col = "blue", type = "l")
  })
}

# Definir la aplicación Shiny
shinyApp(ui = ui, server = server)
