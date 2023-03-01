
source("./mod/mod_checkbox.R")


anno_track_UI <- function(id, height = 100) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("plot"),
               brush = brushOpts(id = ns("brush"),
                                 resetOnNew = T,
                                 direction = "x"),
               dblclick = ns("dblclick"),
               height = height)
  )
}

anno_track_Server <- function(id, p, ranges, zoom = T){
  moduleServer(
    id,
    function(input, output, session) {
      if (isTRUE(zoom)){
        output$plot <- renderPlot({
          p +
            coord_cartesian(xlim= ranges$x, ylim = ranges$y, expand = F)
        })
      } else {
        output$plot <- renderPlot({
          p +
            coord_cartesian(expand = F)
        })
      }
      
      observeEvent(input$dblclick, {
        brush <- input$brush
        if (!is.null(brush)) {
          ranges$x <- c(brush$xmin, brush$xmax)
        } else {
          ranges$x <- NULL
        }
      })
      return (ranges) 
    }
  )
}

anno_table_UI <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("download"), "Download (.tsv)"),
    dataTableOutput(ns("table"))
  )
}

anno_table_Server <- function(id, df, ranges, chrn) {
  moduleServer(
    id,
    function(input, output, session) {
      
      gr <- reactive({
        GRanges(seqnames = chrn, ranges = IRanges(ranges$x[1], ranges$x[2]))
      })
      
      
      df_table <- reactive({
        as.data.frame(subsetByOverlaps(as(df,"GRanges"),gr()))
      })
      
      output$table <- renderDataTable({
        df_table()
      })
      
      output$download <- downloadHandler(
        filename = function() {
          paste0(id,".tsv")
        },
        content = function(file) {
          write.table(df, file, sep = "\t", row.names = F, quote = F)
        }
      )
    }
  )
}

mod_plot_switch_UI <- function(id) {
  ns <- NS(id)
  shinyjs::useShinyjs()
  shinyjs::hidden(fluidRow(
    id =ns("panel"),
    tagList(
      anno_track_UI(ns("plot"))
    )
  )
  )
}

mod_plot_switch_Server <- function(id, cbox, p, ranges) {
  moduleServer(
    id,
    function(input, output, session) {
      ranges <- anno_track_Server("plot", p, ranges)
      observe({
        if (isTRUE(cbox())) {
          shinyjs::show(id = "panel")
        } else {
          shinyjs::hide(id = "panel")
        }
      })
      return(ranges)
    }
  )
}

