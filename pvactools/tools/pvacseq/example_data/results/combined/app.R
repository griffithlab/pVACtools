library(shiny)

source(server.R)
source(ui.R)
source(neofox_ui.R)
source(custom_ui.R)

options(shiny.host = '0.0.0.0')
options(shiny.port = 3333)

shinyApp(ui, server)
