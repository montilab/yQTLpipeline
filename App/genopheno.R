genopheno_ui <- tabPanel(
  title = "Genotype-Phenotype",
  div(
    class = "navbar_content",
    titlePanel("Boxplot of a phenotype vs the genotype of a SNP"),
    HTML("<p>(1) Upload genotype data GDS in RDS format and Phenotype data in RDS or CSV format. Type in a SNP ID that matches an entry in annotation/id node in the GDS file.</p>
       <p>(2) Click \"Retrieve Data\" to retrieve the data. Two [Success] messages will print out if both genotype and phenotype data are uploaded successfully. </p>
       <p>(3) Select a phenotype in the new drop down menu. Click \"Plot Boxplot\" to draw a genotype-phenotype boxplot.</p>"),
    fluidRow(column(12, paste(rep("-", 100), collapse = ""))),
    fluidRow(
      column(
        4,
        fileInput(inputId = "genotype_gds", label = "Genotype Data (GDS)", multiple = FALSE),
        fileInput(inputId = "phenotype_rds", label = "Phenotype Data (RDS,CSV)", multiple = FALSE),
        textInput(inputId = "genopheno_select_snp", label = "Input a SNP in 'annotation/id' entry in GDS", value = "", placeholder = ""),
        actionButton("genopheno_READ", label = "Retrieve Data"),
        HTML("<p><i>You need to retrieve data again after selecting another SNP. </i></p>"),
        br(),
        conditionalPanel(
          condition = "output.genopheno_readpheno_text",
          selectizeInput(inputId = "genopheno_select_phenotype", label = "Select phenotype", choices = NULL, selected = NULL, multiple = FALSE),
          actionButton("genopheno_GO", label = "Plot Boxplot")
        ),
        br(), br()
      ),
      column(
        8,
        htmlOutput("genopheno_readgeno_text"),
        htmlOutput("genopheno_readpheno_text"),
        htmlOutput("genopheno_plot_text"),
        conditionalPanel(
          condition = "output.genopheno_plot_text",
          plotOutput("genopheno_boxplot")
        )
      )
    )
  )
)

genopheno_server <- function(input, output, session) {
  genopheno_genodat_var <- eventReactive(input$genopheno_READ, {
    gds <- seqOpen(input$genotype_gds$datapath)

    ## obtain SNP information
    SNP_info <- data.frame(
      "variant.id" = seqGetData(gds, "variant.id"),
      "chr" = seqGetData(gds, "chromosome"),
      "pos" = seqGetData(gds, "position"),
      "REF" = seqGetData(gds, "$ref"),
      "ALT" = seqGetData(gds, "$alt"),
      "snpID" = seqGetData(gds, "annotation/id")
    )

    if (!input$genopheno_select_snp %in% SNP_info$snpID) {
      seqClose(gds)
      stop("Error: Input SNP is not available in the uploaded GDS. Please retry. ")
    }

    SNP_info <- SNP_info %>% filter(snpID == input$genopheno_select_snp)
    seqSetFilter(gds, variant.sel = SNP_info$variant.id)

    ## obtain genotype information for plotting
    geno_dat <- seqGetData(gds, "$dosage")
    rownames(geno_dat) <- seqGetData(gds, "sample.id")
    colnames(geno_dat) <- "SNP"

    geno_dat_formated <- as.data.frame(geno_dat) %>%
      mutate(across(everything(), ~ factor(case_when(
        . == 0 ~ "Alt/Alt",
        . == 1 ~ "Ref/Alt",
        . == 2 ~ "Ref/Ref",
        TRUE ~ as.character(.)
      ), levels = c("Ref/Ref", "Ref/Alt", "Alt/Alt"))))
    geno_dat_formated$sample.id <- rownames(geno_dat_formated)

    seqResetFilter(gds)
    seqClose(gds)

    snp_name <- paste0(SNP_info$snpID, ", ", SNP_info$chr, "_", SNP_info$pos, "_", SNP_info$REF, "_", SNP_info$ALT)
    return(list(
      "message" = paste0("[Success] Genotype data retrieving success. Selected SNP:", snp_name, "."),
      "geno_dat" = geno_dat_formated,
      "snp_name" = snp_name
    ))
  })

  genopheno_phenodat_var <- eventReactive(input$genopheno_READ, {
    if (any(endsWith(input$phenotype_rds$datapath, c(".rds", ".RDS")))) {
      pheno_dat <- readRDS(input$phenotype_rds$datapath)
    } else if (any(endsWith(input$phenotype_rds$datapath, c(".csv", ".CSV")))) {
      pheno_dat <- read.csv(input$phenotype_rds$datapath)
    } else {
      stop("Error: Input phenotype data must be in RDS or CSV format.")
    }
    if (!"sample.id" %in% colnames(pheno_dat)) {
      stop("Error: Input phenotype data must have a column \"sample.id\".")
    }

    updateSelectizeInput(
      session = session,
      inputId = "genopheno_select_phenotype",
      choices = c("[Select phenotype...]", colnames(pheno_dat)), selected = "[Select phenotype...]"
    )

    return(list(
      "message" = "[Success] Phenotype data retrieving success.",
      "pheno_dat" = pheno_dat
    ))
  })

  # Note: must ask genopheno_phenodat_var() to give output, otherwise it will not run. i.e. will not updateSelectizeInput.
  output$genopheno_readgeno_text <- renderText({
    genopheno_genodat_var()$message
  })
  output$genopheno_readpheno_text <- renderText({
    genopheno_phenodat_var()$message
  })


  genopheno_plot_var <- eventReactive(input$genopheno_GO, {
    if (input$genopheno_select_phenotype == "[Select phenotype...]") {
      stop("Please select a phenotype first.")
    }
    phenodat <- genopheno_phenodat_var()$pheno_dat %>% select(sample.id, !!as.name(input$genopheno_select_phenotype))
    plot_dat <- merge(
      x = genopheno_genodat_var()$geno_dat,
      y = phenodat,
      by = "sample.id", all = F
    )

    plot <- ggplot(plot_dat, aes(x = SNP, y = !!as.name(input$genopheno_select_phenotype), fill = SNP)) +
      geom_boxplot() +
      theme_minimal() +
      xlab(genopheno_genodat_var()$snp_name) +
      theme(legend.position = "none") +
      scale_fill_manual(values = c("#E5A1C3", "#A1A1E5", "#A1E5C3")) +
      ggtitle(paste(input$genopheno_select_phenotype, "x", genopheno_genodat_var()$snp_name))

    return(list(
      "message" = "Finish plotting.",
      "plot" = plot
    ))
  })

  output$genopheno_plot_text <- renderText({
    genopheno_plot_var()$message
  })

  output$genopheno_boxplot <- renderPlot(
    genopheno_plot_var()$plot,
    width = 600, height = 400, res = 120
  )
}
