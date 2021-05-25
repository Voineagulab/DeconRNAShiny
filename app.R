library(shiny)
library(shinyWidgets)
library(shinyjs)
library(ipc)
library(future)
library(ggplot2)
library(promises)

if(Sys.info()['sysname'] == "Windows") {
    plan(multisession) # forks process so uses > 1GB memory
} else {
    plan(multicore) # not available on windows
}

#Provides "sigsBrain"
load("Signatures - Brain.rda")

allSig <- c("F5", "IP", "DM", "MM", "VL", "NG", "CA", "TS", "LK")
allCT <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelia")
allAlg <- c("DeconRNASeq", "CIBERSORT")
steps <- c(0.1, 0.8, 0.1)
stepsStart <- c(0.0, 0.1, 0.9)
empty <- list(NULL, NULL, data.frame(Algorithm=append(rep(allAlg[1], 3), rep(allAlg[2], 3)), r=rep(0, 6)))

# Runs deconvolution pipeline with simple progress and interruption callbacks
source("core.R")
accessibleAnalysisFunction <- function(mixture, signature, algorithms, interruptCallback, progressSet) {
    load.deps()

    totalSteps = length(algorithms) + 1
    currStep = 1

    x <- y <- x2 <- y2 <- NULL

    if(is.element("DeconRNASeq", algorithms)) {
        # Step 1/3: Running DeconRNASeq
        interruptCallback()
        progressSet(value=stepsStart[1], message=sprintf("Step %d/%d: Running DeconRNASeq", currStep, totalSteps))
        x <- run.DRS(mixture,signature)
        currStep <- currStep + 1
    }

    if(is.element("CIBERSORT", algorithms)) {
         # Step 2/3: Running CIBERSORT
        y <- run.CIB(mixture,signature,interruptCallback,progressSet,stepsStart[2],steps[2], currStep, totalSteps)
        currStep <- currStep + 1
    }

    # Step 3/3 Calculating GoFs
    interruptCallback()
    progressSet(value=stepsStart[3], message=sprintf("Step %d/%d Calculating GoFs", currStep, totalSteps))
    if(isTruthy(x)) x2 <- write.gof(mixture, x, signature)
    if(isTruthy(y)) y2 <- write.gof(mixture, y, signature)

    progressSet(value=1.0, message="Done")
    results <- list(x, y, data.frame(Algorithm=append(rep("DeconRNASeq", length(x2$r)), rep("CIBERSORT", length(y2$r))), r=append(x2$r, y2$r)))

    return(results)
}

# Shiny UI
ui <- fluidPage(
    useShinyjs(),

    #Fixes fading out of placeholder when "All" selected
    tags$style(".bs-placeholder {color: #000000 !important;}"),

    titlePanel("BrainDeconvShiny"),
    sidebarLayout(
        sidebarPanel(
            selectInput("signature", "Signature:",
                    allSig),
            pickerInput("celltypes", "Cell types:",
                allCT, 
                multiple=TRUE,
                selected=character(0),
                options=pickerOptions(noneSelectedText = "All")
            ),
            pickerInput("algorithms", "Algorithms:",
                allAlg,
                multiple=TRUE,
                selected=character(0),
                options=pickerOptions(noneSelectedText = "All")
            ),
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
                    downloadButton('getDeconRNASeq', 'DeconRNASeq', icon=icon("file-text")),
                    downloadButton('getCIBERSORT', 'CIBERSORT', icon=icon("file-text")),
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
        cts <- input$celltypes
        if(!isTruthy(cts)) cts <- allCT
        
        mask <- intersect(cts, colnames(sigsBrain[[input$signature]]))
        return(sigsBrain[[input$signature]][mask])
    })

    algUpdate <- reactive({
        algs <- input$algorithms
        if(!isTruthy(algs)) algs <- allAlg
        return(algs)
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
        algs <- algUpdate()

        promises::future_promise({
            accessibleAnalysisFunction(mix, sig, algs, inter$execInterrupts, progress$set)
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
        ggplot(result_val()[[3]], aes(x=Algorithm, y=r, fill=Algorithm)) + geom_violin() + labs(y = "Goodness of fit (r)")
    })

    # Render output
    output$violin <- renderPlot({
        plotOutput()
    }, res = 96)

    # Showing download buttons
    result_val_observer <- observe({
        toggleState("getDeconRNASeq", result_val()[[1]])
        toggleState("getCIBERSORT", result_val()[[2]])
        toggleState("getPlot", result_val()[[1]] != NULL || result_val()[[2]] != NULL)
    })

    # Downloading files
    output$getDeconRNASeq <- downloadHandler(
        filename = function() {"DeconRNASeq.csv"},
        content = function(file) {write.csv(result_val()[[1]], file)}
    )
    output$getCIBERSORT <- downloadHandler(
        filename = function() {"CIBERSORT.csv"},
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