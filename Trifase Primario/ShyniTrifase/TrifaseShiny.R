library(shiny)
library(bslib)
library(readxl)
datos <- read_excel("datos.xlsx")

trifase <- function(t, param){
  y_0 <- param[1]
  t_lag <- param[2]
  t_max <- param[3]
  y_max <- param[4]
  mu <- (y_max-y_0)/(t_max-t_lag)
  z <- ifelse(t <= t_lag, y_0,
              ifelse(t<=t_max, y_0+mu*(t-t_lag), y_max))
  return(z)
}

# Define UI for app that draws a histogram ----
ui <- page_sidebar(
  # App title ----
  title = "Trifase Exploring",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    # Input: y_0
    sliderInput(
      inputId = "y_0",
      label = "Initial value:",
      min = min(datos$y),
      max = max(datos$y),
      value = mean(datos$y)
    ),
    # Input: t_lag
    sliderInput(
      inputId = "t_lag",
      label = "Lag time:",
      min = min(datos$x),
      max = max(datos$x),
      value = 0
    ),
    # Input: t_max
    sliderInput(
      inputId = "t_max",
      label = "Maximum Time:",
      min = min(datos$x),
      max = max(datos$x),
      value = max(datos$x)
    ),
    # Input: y_max
    sliderInput(
      inputId = "y_max",
      label = "Maximum value:",
      min = min(datos$y),
      max = max(datos$y),
      value = max(datos$y)
    )
  ),
  # Output: Histogram ----
  plotOutput(outputId = "plot")
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$plot <- renderPlot({
    
    plot(datos$x,datos$y, col = "red", type = "p", pch = 2,
         xlab = "Time", ylab = "log(CFU/g")
    
    param <- c(input$y_0, input$t_lag, input$t_max, input$y_max)
    z <- trifase(datos$x, param)
    points(datos$x, z, col = "blue", type = "l")
  })
  
}

# Define ShinyApp
shinyApp(ui = ui, server = server)



