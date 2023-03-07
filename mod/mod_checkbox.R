#' Checkbox UI
#'
#' Creates a single checkboxInput element
#'
#' @param value method to change default checkbox value, default = T
#'
#
#' @return tagList of Shiny UI elements
#'
#' @examples
#' mod_checkbox_UI("plot", value = F)
#'
#' @export
mod_checkbox_UI <- function(id, value = T) {
  ns <- NS(id)
  tagList(
    checkboxInput(ns("checkbox"), paste0(id), value = value)
  )
}


#' Checkbox Server
#'
#' returns checkbox state value
#'
#' @param 
#
#' @return list with box_state value
#'
#' @examples
#' mod_checkbox_Server("plot")
#'
#' @export
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