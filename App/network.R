network_ui <- tabPanel(
  title = "Network",
  div(
    class = "navbar_content",
    titlePanel("Network of SNPs vs Phenotypes"),
    fluidRow(column(12, paste(rep("-", 100), collapse = ""))),
    sliderInput(
      inputId = "nw_pval_cutoff",
      label = "significant p-value cutoff 1/10^?",
      value = 7.5, min = 1, max = 30, step = 0.1
    ),
    numericInput(
      inputId = "nw_mac_cutoff",
      label = "MAC cutoff >= ?",
      value = 3, min = 1, max = NA, step = 1
    ),
    actionButton("nw_anal_GO", label = "Plot Association"),
    br(), br(),
    htmlOutput("nw_anal_text"),
    conditionalPanel(
      condition = "output.nw_anal_text",
      # downloadButton("nw_anal_plot_download", "Download plot"),
      plotOutput("nw_anal_plot")
    )
  )
)

network_server <- function(input, output, session, QTLres_var) {
  nw_anal_plot_var <- eventReactive(input$nw_anal_GO, {
    sub <- QTLres_var %>%
      filter(pvalue < 1 / 10^(input$nw_pval_cutoff), MAC >= input$nw_mac_cutoff) %>%
      arrange(pvalue)
    pheno <- unique(sub$phenotype)
    res <- matrix(0, nrow = 0, ncol = ncol(sub))
    colnames(res) <- colnames(sub)
    for (i in c(1:length(pheno))) {
      pheno0 <- pheno[i]
      dat0 <- sub[which(sub$phenotype == pheno0), ]
      dat0 <- dat0[order(dat0$pvalue, decreasing = F), ]
      chr0 <- unique(dat0$chr)
      for (chridx in chr0) {
        a <- which(dat0$chr == chridx)[1]
        res <- rbind(res, dat0[a, ])
      }
    }
    gdat_wchr <- data.frame(
      "from" = res$snpID,
      "to" = res$phenotype,
      "pval" = res$pvalue
    )
    gdat_wchr <- rbind(
      gdat_wchr,
      distinct(data.frame(
        "from" = res$snpID,
        "to" = paste0("chr", res$chr),
        "pval" = 1
      ))
    )
    g_wchr <- graph_from_data_frame(gdat_wchr, directed = F)
    SNP_length <- length(unique(res$snpID))
    pheno_length <- length(unique(res$phenotype))
    chr_length <- length(unique(res$chr))
    vertex_attr(g_wchr, name = "color") <- c(rep("lightsteelblue1", SNP_length), rep("mistyrose", pheno_length), rep("darkseagreen1", chr_length))

    ## saving igraph plot in a variable does not work. use recordPlot().
    ## plot it in R device
    plot(g_wchr, vertex.label.cex = 1.3, vertex.size = 8, vertex.label.family = "Roboto", vertex.label.color = "grey20")
    ## record current plot
    nw_plot <- recordPlot()
    ## clean up R plot device
    plot.new()

    return(list(
      "message" = "[Success]",
      "plot" = nw_plot
    ))
  })

  output$nw_anal_text <- renderText({
    nw_anal_plot_var()$message
  })

  output$nw_anal_plot <- renderPlot(
    {
      nw_anal_plot_var()$plot
    },
    width = 800,
    height = 800
  )

  # output$nw_anal_plot_download <- downloadHandler(
  #   filename = paste0("plot.png"),
  #   content = function(file) {
  #     png(file)
  #     nw_anal_plot_var()$plot
  #     dev.off()
  #   }
  # )
}
