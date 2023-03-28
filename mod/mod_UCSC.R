

mod_UCSC_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12, HTML(as.character(uiOutput(ns("ucsc_link")))))
    )
  )
}

mod_UCSC_Server <- function(id, gv, chr, ranges) {
  moduleServer(
    id,
    function(input, output, session) {
      req(!is.null(ranges$x))
      # Generate the UCSC Genome Browser URL link
      output$ucsc_link <- renderUI({
        ucsc_link <- paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=", gv, "&position=",
                            chr, ":", round(ranges$x[1]), "-", round(ranges$x[2]))
        a("UCSC Genome Browser", tags$i(class = "fa fa-external-link", style = "padding-right: 5px;"),
          href = ucsc_link, target = "_blank")
      })
    }
  )
}