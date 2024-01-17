library(shiny)
library(dqshiny)
library(reactable)
library(igraph)
library(qqman)
library(dplyr)
options(shiny.maxRequestSize = 30 * 1024^2)
source("phenotype.R", local = TRUE)
source("network.R", local = TRUE)

ui <- fluidPage(
  titlePanel("QTL Visualization"),
  fileInput(
    inputId = "QTLres", label = "Choose QTL Result",
    multiple = FALSE
  ),
  HTML(c(
    "<p>The input needs to be a dataframe with columns: chr, pos, snpID or variant.id, phenotype, pvalue or Score.pval, beta or Est, mac or MAC.</p>",
    "<p>Preview available after upload...</p>"
  )),
  reactableOutput("QTL_df_preview"),
  tabsetPanel(
    type = "tabs",
    phenotype_ui,
    network_ui
  )
) # end of FluidPage

server <- function(input, output, session) {
  QTLres_var <- reactive({
    ## Read input QTL result RDS and format columns
    req(input$QTLres)
    res <- readRDS(input$QTLres$datapath)
    if ("Score.pval" %in% colnames(res)) res <- res %>% rename(pvalue = Score.pval)
    if ("Est" %in% colnames(res)) res <- res %>% rename(beta = Est)
    if ("mac" %in% colnames(res) & !"MAC" %in% colnames(res)) res <- res %>% rename(MAC = mac)
    if (!"snpID" %in% colnames(res)) res$snpID <- paste0("snp", res$variant.id)
    res$snpID[which(is.na(res$snpID))] <- paste0("snp", res$variant.id[which(is.na(res$snpID))])
    return(res)
  })
  output$QTL_df_preview <- renderReactable({
    updateSelectizeInput(
      session = session,
      inputId = "select_phenotype",
      choices = c("[All phenotypes]", unique(QTLres_var()$phenotype)), selected = "[All phenotypes]"
    )
    return(reactable(head(QTLres_var(), 3)))
  })

  phenotype_server(input, output, session, QTLres_var = QTLres_var())
  network_server(input, output, session, QTLres_var = QTLres_var())
}

# Run the application
shinyApp(ui = ui, server = server)
