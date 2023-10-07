library(readxl)
library(clusterProfiler)
library(shiny)

config = list(
    title = "gene_set_analysis_os.R",
    oryzabase_ontologies = "/home/panyq/Tools/oryzabase-ontologies/oryzabase-ontologies-2023-05-27.xlsx"
)

enricher_os = function(gene, universe, onto) {
    enrich_result = enricher(
        gene=gene,
        universe=universe,
        TERM2GENE=onto[c("OntoID", "GeneID")],
        TERM2NAME=onto[c("OntoID", "Description")],
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
    out = list(
        dotplot = dotplot(enrich_result),
        table = as.data.frame(enrich_result)
    )
    return(out)
}


ui <- fluidPage(
    titlePanel(config[["title"]]),
    fluidRow(
        column(width = 2, selectInput("onto_name", label = "Ontology", choices = excel_sheets(config$oryzabase_ontologies)))
    ),
    fluidRow(
        column(width = 2, textAreaInput("enrich_gene", label = "Gene")),
        column(width = 2, textAreaInput("enrich_universe", label = "Universe")),
        column(width = 2,
            fluidRow(actionButton("enrich_go", label = "Go")),
            fluidRow(downloadButton("enrich_download", label = "Download Annotation"))
        ),
        column(width = 4, plotOutput("enrich_dotplot")),
    ),
    fluidRow(
        column(width = 11, dataTableOutput("temp"))
    )
)

server <- function(input, output, session) {
    onto <- eventReactive(input$enrich_go, {
        read_excel(config$oryzabase_ontologies, sheet=input$onto_name)
    })

    output$temp = renderDataTable({
        onto()
    })
    output$enrich_dotplot = renderPlot(plot(iris))
}

shinyApp(ui, server)

