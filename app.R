library(shiny)
library(vroom)

shinyApp(
    ui = shinyUI(fluidPage(
        titlePanel("DeconRNAShiny"),
        sidebarLayout(
            #Sidebar with inputs
            sidebarPanel(
                selectInput("signature", "Signature:",
                    allSignatures, 
                    multiple=FALSE),
                selectInput("tissues", "Tissues:",
                    allTissues, 
                    selected=allTissues,
                    multiple=TRUE),
                fileInput("file", NULL, accept = c(".csv", ".tsv")),

                numericInput("m", "Number of samples:", 0, min = 0, max = 100),
                actionButton("do", "Start")
            ),

            #Show a plot of the generated distribution
            mainPanel(
                plotOutput("hist")
            )
        )
    )),

    server = function(input, output, session) {
        #File updater to prevent rereading of file every time "signal" is changed
        fUpdate <- reactive({
            validate(need(input$file, 'No file; Please upload a .csv or.tsv file'))
            ext <- tools::file_ext(input$file$name)
            switch(ext,
                csv = vroom::vroom(input$file$datapath, delim = ","),
                tsv = vroom::vroom(input$file$datapath, delim = "\t"),
                validate("Invalid file; Please upload a .csv or .tsv file")
            )
        })

        #Mean updater on button press
        mUpdate <- eventReactive(input$do, {
            validate(need(input$m != 0, 'Invalid sample; Please choose a sample > 0'))
            replicate(1e4, mean(runif(input$m)))
        })

        output$hist <- renderPlot({
            file <- fUpdate()
            means <- mUpdate()
            hist(means, breaks = 20, main="Placeholder Histogram")
        }, res = 96)
    }
)

# TODO:
# source(core.R)
# x <- run.DRS(mixture,signature)
# y <- run.CIB(mixture,signature)
# x2 <- write.gof(mixture, x, signature)
# y2 <- write.gof(mixture, y, signature)
# For long computation see https://shiny.rstudio.com/articles/bookmarking-state.html 