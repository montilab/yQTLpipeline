library(shiny)
library(reactable)
library(igraph)
library(SeqArray)
library(dplyr)
library(ggplot2)

options(shiny.maxRequestSize = 30 * 1024^2)
source("table_preview.R", local = TRUE)
source("phenotype.R", local = TRUE)
source("genopheno.R", local = TRUE)
source("network.R", local = TRUE)

ui <- fluidPage(
  tagList(
    tags$head(
      tags$link(type = "text/css", rel = "stylesheet", href = "styling.css")
    )
  ),
  div(titlePanel("QTL Visualization"),
    fileInput(
      inputId = "QTLres", label = "Choose QTL Result",
      multiple = FALSE
    ),
    id = "title_section"
  ),
  HTML("<p>The input needs to be a dataframe with columns: chr, pos, snpID or variant.id, phenotype, pvalue or Score.pval, beta or Est, mac or MAC.</p>"),
  navbarPage(
    title = "Visualization",
    position = "fixed-top",
    table_preview_ui,
    phenotype_ui,
    genopheno_ui,
    network_ui
  )
)

server <- function(input, output, session) {
  QTLres_var <- reactive({
    ## Read input QTL result RDS and format columns
    req(input$QTLres)
    res <- readRDS(input$QTLres$datapath)
    if ("Score.pval" %in% colnames(res)) res <- res %>% rename(pvalue = Score.pval)
    if ("Est" %in% colnames(res)) res <- res %>% rename(beta = Est)
    if ("Est.ca" %in% colnames(res)) res <- res %>% rename(beta = Est.ca)
    if ("beta.ca" %in% colnames(res)) res <- res %>% rename(beta = beta.ca)
    if ("mac" %in% colnames(res) & !"MAC" %in% colnames(res)) res <- res %>% rename(MAC = mac)
    if (!"snpID" %in% colnames(res)) res$snpID <- paste0("snp", res$variant.id)
    res$snpID[which(is.na(res$snpID))] <- paste0("snp", res$variant.id[which(is.na(res$snpID))])
    return(res)
  })

  table_preview_server(input, output, session, QTLres_var = QTLres_var())
  phenotype_server(input, output, session, QTLres_var = QTLres_var())
  genopheno_server(input, output, session)
  network_server(input, output, session, QTLres_var = QTLres_var())
}

# Run the application
shinyApp(ui = ui, server = server)
