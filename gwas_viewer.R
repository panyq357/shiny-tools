library(tidyverse)
library(shiny)
library(Gviz)


## TODOs
## - chr, peak plot brush y axis recognize
## - chr, peak plot position connect with input frame
## - snp info hover
## - 

source("gwas_viewer_prep_3k.R")

plot_manhattan = function(result_df) {
    manhattan <- ggplot(result_df) +
        geom_point(aes(x = pos, y = -log10(p))) +
        facet_wrap(~chr, nrow = 1, scales = "free_x") +
        labs(x = NULL, y = expression(-log[10](p.value))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
        theme_bw() +
        theme(panel.spacing.x = unit(0,"line"))
    return(manhattan)
}

plot_global_manhattan = function(result_df) {
    result_df |>
        plot_manhattan() +
        theme(axis.text.x=element_blank())
}

plot_chromosome_manhattan = function(result_df, chromosome) {
    result_df |>
        subset(chr == chromosome) |>
        plot_manhattan() +
        scale_x_continuous(breaks = 0:10 * 10 * 1e6, labels = paste0(0:10 * 10, "Mb"))
}

plot_peak_manhattan = function(result_df, chromosome, start, end) {
    result_df |>
        subset(chr == chromosome & pos >= start & pos <= end) |>
        plot_manhattan()
}



plot_zoom_in_tracks = function(chromosome, start, end, result_df) {

    result_df = subset(result_df, chr==chromosome & pos >= start & pos <= end)
    gwas_gr = with(result_df, data.frame(chr, start=pos, end=pos, neg_log10_p=-log10(p))) |>
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    manhattan = DataTrack(gwas_gr, name = "-log10(p.value)")

    gene_region <- GeneRegionTrack(config$txdb)

    axis = GenomeAxisTrack()

    plotTracks(
        list(manhattan, gene_region, axis),
        transcriptAnnotation="transcript_id",
        chromosome=chromosome,
        from=start,
        to=end
    )
}

get_gene_annotation = function(chromosome, pos_start, pos_end) {

    gene_gr <- genes(config$txdb)
    gene_gr <- subset(gene_gr, start(gene_gr) > pos_start & end(gene_gr) < pos_end & seqnames(gene_gr) == chromosome)

    gene_annotation = data.frame(
        GeneID=gene_gr$gene_id,
        Location=sprintf("%s:%s-%s(%s)", seqnames(gene_gr), start(gene_gr), end(gene_gr), strand(gene_gr))
    )

    gene_annotation = gene_annotation[order(gene_annotation$Location),]

    gene_annotation = cbind(
        gene_annotation, config$gene_anno[match(gene_annotation$GeneID, config$gene_anno[[1]]),-1]
    )

    return(gene_annotation)
}


### shiny app ###

ui <- fluidPage(
    titlePanel(config$title),
    fluidRow(
        column(width = 2,
            selectInput("result_name", label = "Result", choices = names(config$result_list)),
        ),
        column(width = 2,
            sliderInput("min_neg_log10_p_value", "Minimun -log10(p.value)", value = 2, min = 0, max = 10)
        ),
        column(width = 2,
            downloadLink("downloadData", "Download Annotation")
        )
    ),
    fluidRow(
        column(width = 5,
            plotOutput("global_manhattan", click = "global_manhattan_click")
        ),
        column(width = 3,
            plotOutput("chromosome_manhattan", brush = "chromosome_manhattan_brush")
        ),
        column(width = 3,
            plotOutput("peak_manhattan", brush = "peak_manhattan_brush")
        ),
    ),
    fluidRow(
        column(width = 11,
            plotOutput("zoom_in_tracks")
        )
    ),
    fluidRow(
        column(width = 11,
            dataTableOutput("gene_annotation")
        )
    )
)

server <- function(input, output, session) {
    max_p_value <- reactive(10^(-input$min_neg_log10_p_value))

    gwas_result <- reactive(config$result_list[[input$result_name]] |> subset(p <= max_p_value()))

    output$global_manhattan <- renderPlot({
        plot_global_manhattan(gwas_result())
    })

    chromosome <- reactive(input$global_manhattan_click$panelvar1)

    output$chromosome_manhattan <- renderPlot({
        plot_chromosome_manhattan(gwas_result(), chromosome())
    })

    chr_xmin <- reactive(input$chromosome_manhattan_brush$xmin)
    chr_xmax <- reactive(input$chromosome_manhattan_brush$xmax)

    output$peak_manhattan <- renderPlot({
        plot_peak_manhattan(gwas_result(), chromosome(), chr_xmin(), chr_xmax())
    })

    peak_xmin <- reactive(as.integer(input$peak_manhattan_brush$xmin))
    peak_xmax <- reactive(as.integer(input$peak_manhattan_brush$xmax))

    output$zoom_in_tracks <- renderPlot({
         plot_zoom_in_tracks(chromosome(), peak_xmin(), peak_xmax(), gwas_result())
    })

    output$gene_annotation <- renderDataTable({
        get_gene_annotation(chromosome(), peak_xmin(), peak_xmax())
    })

    output$downloadData <- downloadHandler(
        filename = function() {
          paste("chromosome", chromosome(), "-", peak_xmin(), "-", peak_xmax(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(
            get_gene_annotation(chromosome(), peak_xmin(), peak_xmax()), file)
        }
    )
}

shinyApp(ui, server)

