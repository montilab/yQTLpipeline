source("manhattan_plot.R")

phenotype_ui <- tabPanel(
  "Analyze Phenotype",
  titlePanel("Manhattan & Miami Plot of a Phenotype"),
  fluidRow(column(12, paste(rep("-", 100), collapse = ""))),
  fluidRow(
    column(
      3,
      selectizeInput(inputId = "pheno_select_phenotype", label = "Select phenotype", choices = NULL, multiple = FALSE),
      textAreaInput(inputId = "pheno_select_snp", label = "Select a subset of SNPs in snpID column (seperate by ',') to color in the plots", value = "", placeholder = "", rows = 4),
      sliderInput(
        inputId = "pheno_pval_cutoff",
        label = "p-value cutoff 1e-?",
        value = 7.5, min = 1, max = 15, step = 0.1
      ),
      numericInput(
        inputId = "pheno_mac_cutoff",
        label = "MAC cutoff >= ?",
        value = 3, min = 1, max = NA, step = 1
      ),
      actionButton("pheno_GO", label = "Plot")
    ),
    column(
      9,
      htmlOutput("pheno_anal_text"),
      conditionalPanel(
        condition = "output.pheno_anal_text",
        HTML("<p><b>Manhattan plot</b>:</p> "),
        plotOutput("pheno_anal_mhtplot"),
        HTML("<p><b>Miami plot</b>: (note: you will receive error if all SNPs are in the same direction)</p>"),
        plotOutput("pheno_anal_miamiplot"),
        HTML("<p><b>Top SNPs: </b></p>"),
        reactableOutput("pheno_anal_res_df")
      )
    )
  )
)

phenotype_server <- function(input, output, session, QTLres_var) {
  pheno_anal_text_var <- eventReactive(input$pheno_GO, {
    res <- QTLres_var
    if (input$pheno_select_phenotype == "[All phenotypes]") {
      pheno_name <- "All phenotypes"
    } else {
      pheno_name <- input$pheno_select_phenotype
      res <- res %>% filter(phenotype == input$pheno_select_phenotype)
    }
    res <- res %>%
      filter(pvalue <= 10^(-input$pheno_pval_cutoff)) %>%
      arrange(pvalue)

    return(list(
      "nrow_res" = nrow(res),
      "pheno_select_phenotype" = input$pheno_select_phenotype, "QTL_num" = length(unique(res$snpID)),
      "pheno_pval_cutoff" = input$pheno_pval_cutoff,
      "top_QTL" = res$snpID[1], "top_QTL_chr" = res$chr[1], "top_QTL_pos" = res$pos[1]
    ))
  })

  output$pheno_anal_text <- renderText({
    if (pheno_anal_text_var()$nrow_res == 0) {
      return(c(
        "<p>Analyze phenotype: ", pheno_anal_text_var()$pheno_select_phenotype, "</p>",
        "<p>No QTL found exceeding the specified p-value cutoff.</p>"
      ))
    } else {
      return(c(
        "<p>Analyze phenotype: ", pheno_anal_text_var()$pheno_select_phenotype, "</p>",
        "<p>Total number of QTLs exceed specified p-value cutoff:", pheno_anal_text_var()$QTL_num, "</p>",
        "<p>Top QTL: ", pheno_anal_text_var()$top_QTL, " , on chromosome ", pheno_anal_text_var()$top_QTL_chr, " pos ", pheno_anal_text_var()$top_QTL_pos, "</p>"
      ))
    }
  })

  SNP_subet <- eventReactive(eventExpr = input$pheno_GO, {
    if (input$pheno_select_snp != "") {
      pheno_select_snpid <- input$pheno_select_snp
      pheno_select_snpid <- gsub(" ", "", pheno_select_snpid, fixed = TRUE)
      pheno_select_snpid <- gsub("\n", "", pheno_select_snpid, fixed = TRUE)
      pheno_select_snpid <- unlist(strsplit(pheno_select_snpid, ","))
    }
    return(pheno_select_snpid)
  })

  pheno_anal_res_df_var <- eventReactive(eventExpr = input$pheno_GO, {
    res <- QTLres_var
    if (input$pheno_select_phenotype != "[All phenotypes]") {
      res <- res %>% filter(phenotype == input$pheno_select_phenotype)
    }
    if (input$pheno_select_snp != "") {
      res <- res %>% filter(snpID %in% SNP_subet())
    }
    res <- res %>%
      filter(pvalue <= 1 / 10^(input$pheno_pval_cutoff), MAC >= input$pheno_mac_cutoff) %>%
      arrange(pvalue) %>%
      head(10)
    return(reactable(res))
  })

  output$pheno_anal_res_df <- renderReactable({
    pheno_anal_res_df_var()
  })

  pheno_anal_mhtplot_var <- eventReactive(eventExpr = input$pheno_GO, {
    res <- QTLres_var
    pheno_name <- input$pheno_select_phenotype
    if (pheno_name != "[All phenotypes]") {
      res <- res %>% filter(phenotype == pheno_name)
    }
    res <- res %>% filter(MAC >= input$pheno_mac_cutoff)
    highlight_SNPs <- NULL
    if (input$pheno_select_snp != "") {
      highlight_SNPs <- SNP_subet()
    }
    mht_plot <- manhattan(
      x = res, chr = "chr", bp = "pos", snp = "snpID", p = "pvalue",
      highlight = highlight_SNPs,
      genomewideline = pheno_anal_text_var()$pheno_pval_cutoff,
      main = paste(pheno_name, "( mac", input$pheno_mac_cutoff, ")")
    )
    return(mht_plot)
  })

  pheno_anal_miamiplot_var <- eventReactive(eventExpr = input$pheno_GO, {
    res <- QTLres_var
    pheno_name <- input$pheno_select_phenotype
    if (pheno_name != "[All phenotypes]") {
      res <- res %>% filter(phenotype == pheno_name)
    }
    res <- res %>% filter(MAC >= input$pheno_mac_cutoff)
    highlight_SNPs <- NULL
    if (input$pheno_select_snp != "") {
      highlight_SNPs <- SNP_subet()
    }
    par(mfrow = c(2, 1))
    par(mar = c(1.3, 3, 3, 3))
    ylim <- max(-log10(res$pvalue)) + 0.2

    if (all(res$beta > 0) | all(res$beta < 0)) {
      stop("Unable to draw Miami plot. Estimates for all SNPs are in one direction. \n")
    }

    miami_plot <- manhattan(
      res %>% filter(beta > 0),
      ylim = c(0, ylim),
      chr = "chr", bp = "pos", snp = "snpID", p = "pvalue",
      highlight = highlight_SNPs,
      genomewideline = pheno_anal_text_var()$pheno_pval_cutoff,
      main = paste(pheno_name, "( mac", input$pheno_mac_cutoff, ")")
    )
    par(mar = c(3, 3, 1.3, 3))
    miami_plot <- c(miami_plot, manhattan(
      res %>% filter(beta < 0),
      ylim = c(ylim, 0), xlab = "", xaxt = "n",
      chr = "chr", bp = "pos", snp = "snpID", p = "pvalue",
      highlight = highlight_SNPs,
      genomewideline = pheno_anal_text_var()$pheno_pval_cutoff
    ))
    return(miami_plot)
  })

  output$pheno_anal_mhtplot <- renderPlot({
    pheno_anal_mhtplot_var()
  })

  output$pheno_anal_miamiplot <- renderPlot({
    pheno_anal_miamiplot_var()
  })

  output$pheno_anal_mhtplot_download <- downloadHandler(
    filename = paste0("mhtplot", pheno_anal_text_var()$pheno_select_phenotype, ".png"),
    content = function(file) {
      png(file)
      pheno_anal_mhtplot_var()
      dev.off()
    }
  )
}
