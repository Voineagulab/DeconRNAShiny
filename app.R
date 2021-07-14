library(shiny)
library(shinyWidgets)
library(shinyjs)
library(ipc)
library(future)
library(ggplot2)
library(promises)
library(shinythemes)

if(Sys.info()['sysname'] == "Windows") {
    plan(multisession) # forks process so uses > 1GB memory
} else {
    plan(multicore) # not available on windows
}

load("sigsBrain.rda")

allSig <- c("F5", "IP", "DM", "MM", "VL", "NG", "CA", "TS", "LK", "MB")
defaultCT <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelia")
otherCT <- c("OPCs", "Excitatory", "Inhibatory")
allCT <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelia", "OPCs", "Excitatory", "Inhibatory")
choicesCT <- list(Default=defaultCT, Other=otherCT)
allAlg <- c("dtangle", "CIBERSORT")
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

    if(is.element("dtangle", algorithms)) {
        # Step 1/3: Running dtangle
        interruptCallback()
        progressSet(value=stepsStart[1], message=sprintf("Step %d/%d: Running dtangle", currStep, totalSteps))
        x <- run.DTA(mixture,signature)
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
    if(isTruthy(x)) x2 <- write.gof.v2(mixture, x, signature)
    if(isTruthy(y)) y2 <- write.gof.v2(mixture, y, signature)

    progressSet(value=1.0, message="Done")
    results <- list(x, y, data.frame(Algorithm=append(rep("dtangle", length(x2$r)), rep("CIBERSORT", length(y2$r))), r=append(x2$r, y2$r)))

    return(results)
}

# Shiny UI
ui <- navbarPage(theme = shinytheme("paper"), "BrainDeconvShiny", 
    tabPanel("About",
        tags$head(
        # Note the wrapping of the string in HTML()
        tags$style(HTML("
        dl {
            display: grid;
            grid-template-columns: max-content auto;
            column-gap: 1rem;
        }
        dt { grid-column-start: 1;}
        dd { grid-column-start: 2;}
        .well { background-color: #f9f9f9;}
        "))),
        verticalLayout(
            fluidRow(
                column(width = 12, align="center",
                    wellPanel(
                        h3("Welcome to BrainDeconvShiny"),
                        HTML('
                        <span style="text-align:center;">
                            <h5>Developed by VoineaguLab, view the source code on <a href="https://github.com/Voineagulab/DeconRNAShiny">GitHub</a>.</h5>
                        </span>
                        ')
                    ),
                ),
            ),
            wellPanel(
                fluidRow(
                    column(width = 12, align="center",
                        div(img(src="ShinyAppFig.png",width="80%"), style="text-align: center; max-width: 800px;")
                    )
                )
            ),
            wellPanel(
                fluidRow(
                    column(width = 12, align="left",
                        HTML('
                            <span style="text-align:left; color:black;">
                            <h5>BrainDeconvShiny implements the following deconvolution methods:</h5>
                            <dl>
                                <dt>CIBERSORT v1.04 </dt>
                                <dd>obtained from https://cibersort.stanford.edu with default parameters.</dd>
                                <dt>dtangle v2.09 </dt>
                                <dd>using the dtangle Cran R package with default parameters.</dd>
                            </dl>
                            <h5>And the following cell-type signatures:</h5>
                            <dl>
                                <dt>F5 (FANTOM5) </dt>
                                <dd>Cap Analysis of Gene Expression (CAGE) data for robust CAGE peaks was obtained from the FANTOM5 consortium: http://fantom.gsc.riken.jp/5/data/44.</dd>
                                <dt>IP (immuno-purified) </dt>
                                <dd>RNA-seq data from cells immunopurified from human adult brain tissue extracted during surgery were obtained from Zhang et al. 201641.</dd>
                                <dt>MM (Mus musculus) </dt>
                                <dd>RNA-seq data from immunopurified mouse brain tissue was obtained from Zhang et al. 201442.</dd>
                                <dt>DM (Darmanis) </dt>
                                <dd>Human brain single-cell gene expression data from the middle temporal gyrus generated by Darmanis et al. (2015)13.</dd>
                                <dt>LK (Lake)</dt>
                                <dd>Gene expression data for 10,319 human adult frontal cortex nuclei were accessed from Lake et al. 201840.</dd>
                                <dt>VL (Velmeshev) </dt>
                                <dd>10X Chromium for single-nucleus data from the post-mortem adult human brain were accessed Velmeshev et al.36.</dd>
                                <dt>CA (Cell Atlas) </dt>
                                <dd>Count-level exon expression data for NeuN+ sorted adult nuclei from the middle temporal gyrus were acquired from the Human Cell Atlas37.</dd>
                                <dt>NG (Nagy) </dt>
                                <dd>10X Chromium single-nucleus expression data from the adult human post-mortem human prefrontal cortex were accessed from Nagy et al.39.</dd>
                                <dt>TS (Tasic) </dt>
                                <dd>Exon-level SmartSeq2 single-cell expression data from the adult mouse cortex were accessed from Tasic et al.43.</dd>
                                <dt>MB (Multibrain) </dt>
                                <dd>this composite signature was generated by quantile-normalising and averaging the rpkm-level expression of the CA, IP, DM, NG, and VL signatures for five cell-types (neurons, astrocytes, oligodendrocytes, microglia, and endothelia).</dd>
                            </dl>
                            <p>All signatures are cortical in origin but represent a range of purification protocols (scRNA-seq by SmartSeq (DM), snRNA-seq by 10X (VL, NG), snRNA-seq by SmartSeq (CA), and immuno-panning (IP)). Detailed information on data processing and normalisation for cell-ype signatures is available in the Methods section, Sutton et al. 2021 (https://www.biorxiv.org/content/10.1101/2020.06.01.126839v1).</p>

                            <p>Goodness of fit is evaluated using a Pearson correlation coefficient calculated for each sample between the bulk gene expression data and the gene expression data reconstructed with a given cell-type signature and the corresponding estimated cell type proportions.</p>
                            </span>'),
                    )
                )
            )
        )
    ),
    tabPanel("Run",
        # Enables button showing/hiding
        useShinyjs(),

        # Fixes fading out of placeholder when "All" selected
        tags$style(".bs-placeholder {color: #000000 !important;}"),
        tags$head(tags$style("#violin{height:80vh !important;}")),

        # Fixes file select input text overlap
        tags$head(tags$style("input.form-control{margin-left:10px !important;}")),

        sidebarLayout(
            sidebarPanel(
                pickerInput("signature", "Signature:",
                        allSig),
                pickerInput("celltypes", "Cell types:",
                    choicesCT, 
                    multiple=TRUE,
                    selected=character(0),
                    options=pickerOptions(noneSelectedText = "All Default")
                ),
                pickerInput("algorithms", "Algorithms:",
                    allAlg,
                    multiple=TRUE,
                    selected=character(0),
                    options=pickerOptions(noneSelectedText = "All")
                ),
                a(href="mixture.csv", "Mixture: (see example)", download=NA, target="_blank"),
                fileInput("file", NULL, accept = c(".csv", ".tsv")),
                
                p(strong("Run Deconvolution:")),

                fluidRow(
                    column(width = 12, align="center",
                        actionButton('run', 'Run', width='45%'),
                        actionButton('cancel', 'Cancel', width='45%')
                    ),
                ),
            ),
            mainPanel(
                fluidRow(
                    column(width = 12, align="right",
                        downloadButton('getdtangle', 'dtangle', icon=icon("file-text")),
                        downloadButton('getCIBERSORT', 'CIBERSORT', icon=icon("file-text")),
                        downloadButton('getPlot', 'Plot', icon = icon("download"))
                    )
                ),
                plotOutput("violin"),
            )
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

    userSelectCT <- FALSE
    observeEvent(input$celltypes, {
        userSelectCT <<- (length(input$celltypes) > 0)
    }, ignoreNULL = FALSE)

    observeEvent(input$signature, {
        # Get vector of selected (character(0) if all)
        selectedCT <- input$celltypes
        enabledCT <- colnames(sigsBrain[[input$signature]])

        # Convert character(0) to all
        if(!isTruthy(selectedCT)) selectedCT <- defaultCT

        disabledCT <- !allCT %in% enabledCT
        maskCT <- intersect(selectedCT, enabledCT)

        # Default to all cell types if all user selections are disabled
        if(!length(maskCT)) maskCT <- enabledCT

        newSelectedCT <- if(!userSelectCT && setequal(maskCT, intersect(defaultCT, enabledCT)) && !length(setdiff(maskCT, defaultCT))) character(0) else maskCT

        #Update cell type dropdown based on signature
        updatePickerInput(
            session=session, inputId="celltypes",
            choices=choicesCT,
            selected=newSelectedCT,
            choicesOpt=list(disabled = disabledCT)
        )
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
        showNotification("Cancel request sent to server")
        inter$interrupt("Cancelled")
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
    violinOutput <- reactive({
        req(result_val())
        ggplot(result_val()[[3]], aes(x=Algorithm, y=r, fill=Algorithm)) + geom_violin() + labs(y = "Goodness of fit (r)")
    })

    # Render output
    output$violin <- renderPlot({
        violinOutput()
    }, res = 96)

    # Showing download buttons
    result_val_observer <- observe({
        toggleState("getdtangle", result_val()[[1]])
        toggleState("getCIBERSORT", result_val()[[2]])
        toggleState("getPlot", result_val()[[1]] != NULL || result_val()[[2]] != NULL)
    })

    # Downloading files
    output$getdtangle <- downloadHandler(
        filename = function() {"dtangle.csv"},
        content = function(file) {write.csv(result_val()[[1]], file)}
    )
    output$getCIBERSORT <- downloadHandler(
        filename = function() {"CIBERSORT.csv"},
        content = function(file) {write.csv(result_val()[[2]], file)}
    )

    # Downloading plot
    output$getPlot <- downloadHandler(
        filename = function() {"plot.png"},
        content = function(file) {ggsave(file, plot = violinOutput(), device = "png")}
    )

    # Clean up ipc interruptor
    session$onSessionEnded(function() {
        inter$destroy()
    })
}

# Run the application
shinyApp(ui = ui, server = server)