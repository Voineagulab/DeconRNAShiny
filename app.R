library("shiny")
library("shinyjs")
library("ipc")
library("future")
library("ggplot2")
library("promises")

if(Sys.info()['sysname'] == "Windows") {
    plan(multisession) # forks process so uses > 1GB memory
} else {
    plan(multicore) # not available on windows
}

#Provides "sigsBrain"
load("Signatures - Brain.rda")

allSignatures <- c("F5", "IP", "DM", "MM", "VL", "NG", "CA", "TS", "LK")
allTissues <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelia")
steps <- c(0.1, 0.8, 0.1)
steps_start <- c(0.0, 0.1, 0.9)
empty <- list(data.frame(), data.frame(), data.frame(algorithm=c("CIBERSORT", "CIBERSORT", "CIBERSORT", "DeconRNASeq", "DeconRNASeq", "DeconRNASeq"), r=c(0, 0, 0, 0, 0, 0)))

# Runs deconvolution pipeline with simple progress and interruption callbacks
source("core.R")
accessibleAnalysisFunction <- function(mixture, signature, interruptCallback, progressSet) {
    load.deps()

    # Step 1/3: Running DeconRNASeq
    interruptCallback()
    progressSet(value=steps_start[1], message="Step 1/3: Running DeconRNASeq")
    x <- run.DRS(mixture,signature)

    # Step 2/3: Running CIBERSORT
    y <- run.CIB(mixture,signature,interruptCallback,progressSet,steps_start[2],steps[2])

    # Step 3/3 Calculating GoFs
    interruptCallback()
    progressSet(value=steps_start[3], message="Step 3/3 Calculating GoFs")
    x2 <- write.gof(mixture, x, signature)
    y2 <- write.gof(mixture, y, signature)

    progressSet(value=1.0, message="Done")
    results <- list(x, y, data.frame(algorithm=append(rep("DeconRNASeq", ncol(mixture)), rep("CIBERSORT", ncol(mixture))), r=append(x2$r, y2$r)))

    return(results)
}

# Shiny UI
ui <- fluidPage(
    useShinyjs(),
    titlePanel("DeconRNAShiny"),
    sidebarLayout(
        sidebarPanel(
            selectInput("signature", "Signature:",
                    allSignatures),
            selectInput("tissues", "Tissues:",
                allTissues, 
                selected=allTissues,
                multiple=TRUE),
            fileInput("file", "Mixture:", accept = c(".csv", ".tsv")),
            
            p(strong("Run Deconvolution:")),

            fluidRow(
                column(width = 12, align="center",
                    actionButton('run', 'run', width='45%'),
                    actionButton('cancel', 'cancel', width='45%')
                ),
            ),
        ),
        mainPanel(
            fluidRow(
                column(width = 12, align="right",
                    downloadButton('getX', 'x', icon=icon("file-text")),
                    downloadButton('getY', 'y', icon=icon("file-text")),
                    downloadButton('getPlot', 'plot', icon = icon("download"))
                )
            ),
            plotOutput("violin"),
        )
    )
)

# Shiny server
server <- function(input, output, session) {
    inter <- AsyncInterruptor$new()
    
    result_val <- reactiveVal(empty)
    is_running <- reactiveVal(FALSE)

    # React to changes in is_running
    is_running_observer <- observe({
        toggleState("run", !is_running())
        toggleState("cancel", is_running())
    })

    sigUpdate <- reactive({
        #TODO: ensure there are at least 2 tissues
        mask <- intersect(input$tissues, colnames(sigsBrain[[input$signature]]))
        sigsBrain[[input$signature]][mask]
    })

    # Handle cancel button click
    observeEvent(input$cancel,{
        if(!is_running()) return(NULL)

        disable("cancel")
        showNotification("cancel request sent to server")
        inter$interrupt("cancelled")
    })

    # Handle run button click
    observeEvent(input$run,{
        if(is_running() || !mixUpdate()) return(NULL)
        is_running(TRUE)

        progress <- AsyncProgress$new(message="Initializing", min=0.0, max=1.0, value=0.0)
        sig <- sigUpdate()
        mix <- mixUpdate()

        promises::future_promise({
            accessibleAnalysisFunction(mix, sig, inter$execInterrupts, progress$set)
        }) %>% then(
            onFulfilled = function(value) {
                result_val(value)
            },
            onRejected = function(err) {
                showNotification(err$message)
            }
        ) %>% finally(function() {
            # Hide loading bar whether it finished or not
            progress$close()
            is_running(FALSE)
        })

        # Return something other than the future so we don't block the UI
        return(NULL)
    })

    #File updater to prevent rereading of file every time "signal" is changed
    mixUpdate <- reactive({
        validate(need(input$file, 'No file; Please upload a .csv or.tsv file'))
        ext <- tools::file_ext(input$file$name)
        switch(ext,
            csv = read.table(input$file$datapath, sep=",", header=TRUE, row.names=1),
            tsv = read.table(input$file$datapath, sep="\t", header=TRUE, row.names=1),
            validate("Invalid file; Please upload a .csv or .tsv file")
        )
    })

    # Update plot for rendering and downloading
    plotOutput <- reactive({
        req(result_val())
        ggplot(result_val()[[3]], aes(x=algorithm, y=r, fill=algorithm)) + geom_violin()
    })

    # Render output
    output$violin <- renderPlot({
        plotOutput()
    }, res = 96)

    # Downloading files
    output$getX <- downloadHandler(
        filename = function() {"x.csv"},
        content = function(file) {write.csv(result_val()[[1]], file)}
    )
    output$getY <- downloadHandler(
        filename = function() {"y.csv"},
        content = function(file) {write.csv(result_val()[[2]], file)}
    )

    # Downloading plot
    output$getPlot <- downloadHandler(
        filename = function() {"plot.png"},
        content = function(file) {ggsave(file, plot = plotOutput(), device = "png")}
    )

    # Clean up ipc interruptor
    session$onSessionEnded(function() {
        inter$destroy()
    })
}

# Run the application
shinyApp(ui = ui, server = server)