library(shiny)

source(server.R)
source(ui.R)
source(neofox_ui.R)
source(custom_ui.R)

options(shiny.host = '127.0.0.1')
options(shiny.port = 3333)

shinyApp(ui, server)
