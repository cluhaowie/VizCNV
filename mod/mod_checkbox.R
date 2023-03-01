
mod_checkbox_UI <- function(id, value = T) {
  ns <- NS(id)
  tagList(
    checkboxInput(ns("checkbox"), paste0(id), value = value)
  )
}

mod_checkbox_Server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      return(
        list(
          box_state = reactive({input$checkbox})
        )
      )      
    }
  )
}