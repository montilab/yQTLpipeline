table_preview_ui <- tabPanel(
  "Preview",
  titlePanel("Preview Uploaded Result"),
  fluidRow(column(12, paste(rep("-", 100), collapse = ""))),
  reactableOutput("QTL_res_preview")
)

table_preview_server <- function(input, output, session, QTLres_var) {
  output$QTL_res_preview <- renderReactable({
    updateSelectizeInput(
      session = session,
      inputId = "select_phenotype",
      choices = c("[All phenotypes]", unique(QTLres_var$phenotype)), selected = "[All phenotypes]"
    )
    return(reactable(head(QTLres_var, 10)))
  })
}
