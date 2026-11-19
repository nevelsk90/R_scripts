library(shiny)

counterButton <- function(id, label = "Counter") {
  ns <- NS(id)
  tagList(
    actionButton(ns("button"), label = label),
    verbatimTextOutput(ns("out"))
  )
}

counterServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      count <- reactiveVal(0)
      observeEvent(input$button, {
        count(count() + 1)
      })
      output$out <- renderText({
        count()
      })
      count
    }
  )
}

ui <- fluidPage(
  counterButton(id= "counter1", label="Counter #1"),
  counterButton(id= "counter2", label="Counter #2")
  
)

server <- function(input, output, session) {
  counterServer( id= "counter1")
  counterServer( id= "counter2")
  
}

shinyApp(ui, server)

