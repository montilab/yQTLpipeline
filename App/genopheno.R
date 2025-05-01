genopheno_ui <- tabPanel(
  title = "Genotype-Phenotype",
  div(
    class = "navbar_content",
    titlePanel("Boxplot of a phenotype vs the genotype of a SNP"),
    HTML("<p>(1) Upload the phenotype data in RDS or CSV format. </p>
       <p>(2) Upload genotype data GDS format. Upon successful uploading, a searchable table will appear.
       <p>(3) Type in a SNP from the \"snpID\" column, then click \"Retrieve Data\". Two [Success] messages will print out. </p>
       <p>(4) Select a phenotype in the new drop down menu, and click \"Plot Boxplot\" to draw a genotype-phenotype boxplot.</p>"),
    fluidRow(column(12, paste(rep("-", 100), collapse = ""))),
    fluidRow(
      column(
        4,
        fileInput(inputId = "phenotype_rds", label = "Phenotype Data (RDS,CSV)", multiple = FALSE),
        fileInput(inputId = "genotype_gds", label = "Genotype Data (GDS)", multiple = FALSE),
        textInput(inputId = "genopheno_select_snp", label = "Input a SNP in 'snpID'", value = "", placeholder = ""),
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
        ),
        reactableOutput("genopheno_gds_snp_preview")
      )
    )
  )
)

genopheno_server <- function(input, output, session) {
  output$genopheno_gds_snp_preview <- renderReactable({
    req(input$genotype_gds$datapath)

    gds <- try(seqOpen(input$genotype_gds$datapath), silent = T)
    if (class(gds) == "try-error") {
      stop("Input genotype file is not a valid GDS file. ")
    }

    ## obtain SNP information
    SNP_info <- try(data.frame(
      "snpID" = seqGetData(gds, "annotation/id"),
      "chr" = seqGetData(gds, "chromosome"),
      "pos" = seqGetData(gds, "position"),
      # "allele" = seqGetData(gds, "allele"),
      "REF" = seqGetData(gds, "$ref"),
      "ALT" = seqGetData(gds, "$alt")
    ), silent = T)

    if (class(SNP_info) == "try-error") {
      stop("Input genotype file does not contain at least one of the following node: chromosome, position, annotation/id, $ref, $alt.")
    }

    seqClose(gds)

    return(reactable(SNP_info, searchable = TRUE, defaultPageSize = 5))
  })


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
      stop("Input SNP ID is not available in the input GDS file. Please double check. ")
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
      "message" = paste0(
        "<p><font color=\"lightseagreen\">[Success]</font> Genotype data retrieving success. </p>",
        "<p>Selected SNP:<font color=\"slateblue\">", snp_name, "</font>.</p>"
      ),
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
      stop("Input phenotype data must be in RDS or CSV format.")
    }
    if (!"sample.id" %in% colnames(pheno_dat)) {
      stop("Input phenotype data must have a column \"sample.id\".")
    }

    updateSelectizeInput(
      session = session,
      inputId = "genopheno_select_phenotype",
      choices = c("[Select phenotype...]", colnames(pheno_dat)), selected = "[Select phenotype...]"
    )

    return(list(
      "message" = "<font color=\"lightseagreen\">[Success]</font> Phenotype data retrieving success.",
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
      "message" = "<font color=\"lightseagreen\">[Success]</font> Finish plotting.",
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
