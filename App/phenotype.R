phenotype_ui <- tabPanel(
  "phenotype",
  titlePanel("Analyze phenotype"),
  fluidRow(column(12, paste(rep("-", 100), collapse = ""))),
  selectizeInput(inputId = "select_phenotype", label = "Select phenotype", choices = NULL, multiple = FALSE),
  textAreaInput(inputId = "select_snp", label = "Select SNP (seperate by ',')", value = "", placeholder = "", rows = 4),
  sliderInput(
    inputId = "pval_cutoff",
    label = "p-value cutoff 1e-?",
    value = 7.5, min = 1, max = 15, step = 0.1
  ),
  numericInput(
    inputId = "mac_cutoff",
    label = "MAC cutoff >= ?",
    value = 3, min = 1, max = NA, step = 1
  ),
  actionButton("pheno_anal_GO", label = "Analyze phenotype"),
  dq_space(),
  htmlOutput("pheno_anal_text"),
  conditionalPanel(
    condition = "output.pheno_anal_text",
    reactableOutput("pheno_anal_res_df"),
    # downloadButton("pheno_anal_mhtplot_download", "Download plot"),
    plotOutput("pheno_anal_mhtplot")
  ),
  fluidRow(column(12, paste(rep("-", 100), collapse = "")))
)

phenotype_server <- function(input, output, session, QTLres_var) {
  pheno_anal_text_var <- eventReactive(input$pheno_anal_GO, {
    res <- QTLres_var
    if (input$select_phenotype == "[All phenotypes]") {
      pheno_name <- "All phenotypes"
    } else {
      pheno_name <- input$select_phenotype
      res <- res %>% filter(phenotype == input$select_phenotype)
    }
    res <- res %>%
      filter(pvalue <= 10^(-input$pval_cutoff)) %>%
      arrange(pvalue)

    return(list(
      "nrow_res" = nrow(res),
      "select_phenotype" = input$select_phenotype, "QTL_num" = length(unique(res$snpID)),
      "pval_cutoff" = input$pval_cutoff,
      "top_QTL" = res$snpID[1], "top_QTL_chr" = res$chr[1], "top_QTL_pos" = res$pos[1]
    ))
  })

  output$pheno_anal_text <- renderText({
    if (pheno_anal_text_var()$nrow_res == 0) {
      return(c(
        "<p>Analyze phenotype: ", pheno_anal_text_var()$select_phenotype, "</p>",
        "<p>No QTL found exceeding the specified p-value cutoff.</p>"
      ))
    } else {
      return(c(
        "<p>Analyze phenotype: ", pheno_anal_text_var()$select_phenotype, "</p>",
        "<p>Total number of QTLs exceed specified p-value cutoff:", pheno_anal_text_var()$QTL_num, "</p>",
        "<p>Top QTL: ", pheno_anal_text_var()$top_QTL, " , on chromosome ", pheno_anal_text_var()$top_QTL_chr, " pos ", pheno_anal_text_var()$top_QTL_pos, "</p>"
      ))
    }
  })

  pheno_anal_res_df_var <- eventReactive(eventExpr = input$pheno_anal_GO, {
    res <- QTLres_var
    if (input$select_phenotype != "[All phenotypes]") {
      res <- res %>% filter(phenotype == input$select_phenotype)
    }
    if (input$select_snp != "") {
      select_snpid <- input$select_snp
      select_snpid <- gsub(" ", "", select_snpid, fixed = TRUE)
      select_snpid <- gsub("\n", "", select_snpid, fixed = TRUE)
      select_snpid <- unlist(strsplit(select_snpid, ","))
      res <- res %>% filter(snpID %in% select_snpid)
    }
    res <- res %>%
      filter(pvalue <= 1 / 10^(input$pval_cutoff), MAC >= input$mac_cutoff) %>%
      arrange(pvalue)
    return(reactable(res))
  })

  output$pheno_anal_res_df <- renderReactable({
    pheno_anal_res_df_var()
  })

  pheno_anal_mhtplot_var <- eventReactive(eventExpr = input$pheno_anal_GO, {
    res <- QTLres_var
    pheno_name <- input$select_phenotype
    if (pheno_name != "[All phenotypes]") {
      res <- res %>% filter(phenotype == pheno_name)
    }
    res <- res %>% filter(MAC >= input$mac_cutoff)
    return(manhattan(
      x = res, chr = "chr", bp = "pos", snp = "snpID", p = "pvalue",
      genomewideline = pheno_anal_text_var()$pval_cutoff,
      main = paste(pheno_name, "( mac", input$mac_cutoff, ")")
    ))
  })

  output$pheno_anal_mhtplot <- renderPlot({
    pheno_anal_mhtplot_var()
  })

  output$pheno_anal_mhtplot_download <- downloadHandler(
    filename = paste0("mhtplot", pheno_anal_text_var()$select_phenotype, ".png"),
    content = function(file) {
      png(file)
      pheno_anal_mhtplot_var()
      dev.off()
    }
  )
}
