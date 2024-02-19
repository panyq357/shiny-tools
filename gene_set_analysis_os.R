library(readxl)
library(readr)
library(clusterProfiler)
library(shiny)

config = list(
    title = "gene_set_analysis_os.R",
    oryzabase_ontologies = "/home/panyq/Tools/onto-scripts/results/oryzabase-ontologies.xlsx"
)

sheets <- excel_sheets(config$oryzabase_ontologies)
onto_list <- lapply(
    setNames(sheets, sheets),
    function(x) read_excel(config$oryzabase_ontologies, sheet=x)
)

enricher_os <- function(gene, universe, onto_name) {
    res <- enricher(
        gene = gene,
        universe = universe,
        TERM2GENE = onto_list[[onto_name]][c("OntoID", "GeneID")],
        TERM2NAME = onto_list[[onto_name]][c("OntoID", "Description")],
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
    out <- list(
        dotplot = dotplot(res),
        table = as.data.frame(res)
    )
    return(out)
}


ui <- fluidPage(
    titlePanel(config[["title"]]),
    fluidRow(
        column(width = 2, selectInput("onto_name", label = "Ontology", choices = excel_sheets(config$oryzabase_ontologies)))
    ),
    fluidRow(
        column(width = 3, textAreaInput("gene", label = "Gene")),
        column(width = 3, textAreaInput("universe", label = "Universe")),
        column(width = 1, fluidRow(actionButton("go", label = "Go"))),
        column(width = 1, fluidRow(downloadButton("download", label = "Download")))
    ),
    fluidRow(
        column(width = 8, plotOutput("dotplot")),
    )
)

server <- function(input, output, session) {
    res <- eventReactive(input$go, {
        gene <- strsplit(input$gene, "\n")[[1]]
        if (nchar(input$universe) > 0) {
            universe <- strsplit(input$universe, "\n")[[1]]
        } else {
            universe <- NULL
        }
        return(enricher_os(gene, universe, input$onto_name))
    })
    output$dotplot = renderPlot(res()$dotplot)

    output$download <- downloadHandler(
        filename = function() {
            sprintf("enrich_res.%s.csv", input$onto_name)
        },
        content = function(file) {
            write_excel_csv(res()$table, file)
        }
    )
}

shinyApp(ui, server)

