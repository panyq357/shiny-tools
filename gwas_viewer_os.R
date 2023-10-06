library(tidyverse)
library(shiny)
library(biomaRt)
library(Gviz)
library(vroom)               # For fast table reading.


## TODOs
## - chr, peak plot brush y axis recognize
## - chr, peak plot position connect with input frame
## - snp info hover
## - 

config = list(
    "title" = "3K-RG-2K-samples-GWAS",
    "biomart" = useEnsemblGenomes(biomart="plants_mart", dataset="osativa_eg_gene"),
    "col_names" = list(
        "chr" = "Chromosome",
        "pos" = "Position",
        "p" = "P.value"
    ),
    "results" = list(
        "GL_Blink" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/GL/GAPIT.Blink.GL.GWAS.Results.csv",
        "GL_CMLM" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/GL/GAPIT.CMLM.GL.GWAS.Results.csv",
        "GL_GLM" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/GL/GAPIT.GLM.GL.GWAS.Results.csv",
        "GW_Blink" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/GW/GAPIT.Blink.GW.GWAS.Results.csv",
        "GW_CMLM" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/GW/GAPIT.CMLM.GW.GWAS.Results.csv",
        "GW_GLM" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/GW/GAPIT.GLM.GW.GWAS.Results.csv",
        "TGW_Blink" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/TGW/GAPIT.Blink.TGW.GWAS.Results.csv",
        "TGW_CMLM" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/TGW/GAPIT.CMLM.TGW.GWAS.Results.csv",
        "TGW_GLM" = "/mnt/d/Analysis/2022-04-01_LYH_3K-RG-2K-samples-GWAS/results/TGW/GAPIT.GLM.TGW.GWAS.Results.csv"
    )
)

get_gwas_result = function(result_name) {
    result_file_path = config[["results"]][[result_name]]
    raw_result = vroom::vroom(result_file_path)
    clean_result = data.frame(
        chr = raw_result[[config[["col_names"]][["chr"]]]],
        pos = raw_result[[config[["col_names"]][["pos"]]]],
        p = raw_result[[config[["col_names"]][["p"]]]]
    )
    return(clean_result)
}

plot_manhattan = function(gwas_result) {
    ggplot(gwas_result) +
        geom_point(aes(x = pos, y = -log10(p))) +
        facet_wrap(~chr, nrow = 1, scales = "free_x") +
        labs(x = NULL, y = expression(-log[10](p.value))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
        theme_bw() +
        theme(panel.spacing.x = unit(0,"line"))
}

plot_global_manhattan = function(gwas_result) {
    gwas_result |>
        plot_manhattan() +
        theme(axis.text.x=element_blank())
}

plot_chromosome_manhattan = function(gwas_result, chromosome) {
    gwas_result |>
        subset(chr == chromosome) |>
        plot_manhattan() +
        scale_x_continuous(breaks = 0:10 * 10 * 1e6, labels = paste0(0:10 * 10, "Mb"))
}

plot_peak_manhattan = function(gwas_result, chromosome, start, end) {
    gwas_result |>
        subset(chr == chromosome & pos >= start & pos <= end) |>
        plot_manhattan()
}

plot_zoom_in_tracks = function(chromosome, start, end, gwas_result) {

    gwas_result = subset(gwas_result, chr==chromosome & pos >= start & pos <= end)
    gwas_gr = with(gwas_result, data.frame(chr, start=pos, end=pos, neg_log10_p=-log10(p))) |>
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    manhattan_track = DataTrack(gwas_gr, name = "-log10(p.value)")

    gene_track = BiomartGeneRegionTrack(
        start=start, end=end, chromosome=chromosome,
        biomart=config[["biomart"]], name="Gene Models", transcriptAnnotation = "transcript"
    )
    axis_track = GenomeAxisTrack()
    plotTracks(list(manhattan_track, gene_track, axis_track), from=start, to=end)
}

get_gene_annotation = function(chromosome, start, end) {
    query_result = getBM(
        attributes=c(
            "ensembl_gene_id",
            "external_gene_name",
            "chromosome_name",
            "start_position",
            "end_position",
            "strand",
            "description"
        ),
        filters=list("chromosome_name"=chromosome, "start"=start, "end"=end),
        mart=config[["biomart"]]
    )

    query_result$strand[query_result$strand == -1] = "-"
    query_result$strand[query_result$strand ==  1] = "+"

    gene_annotation = with(query_result, data.frame(
        GeneID=ensembl_gene_id,
        Name=external_gene_name,
        Location=sprintf("%s:%s-%s(%s)", chromosome_name, start_position, end_position, strand),
        Description=description
    ))

    gene_annotation = gene_annotation[order(gene_annotation$Location),]

    return(gene_annotation)
}


### shiny app ###

ui <- fluidPage(
    titlePanel(config[["title"]]),
    fluidRow(
        column(width = 2,
            selectInput("result_name", label = "Result", choices = names(config[["results"]])),
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

    gwas_result <- reactive(get_gwas_result(input$result_name) |> subset(p <= max_p_value()))

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

