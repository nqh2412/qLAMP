library(shiny)
library(shinyBS)
library(imager)
library(colocr)
library(dplyr)
library(qpcR)

source('img_func.R')
source('qpcr_func.R')

# Define UI for application that draws a histogram
ui <-  fluidPage(
  

  
  h3("qLAMP", span("", style = "font-weight: 300"), 
     style = "color: #000000; text-align: center;
background-color:#00AD87;

     padding: 10px"),
  column(6, 
         p("A tool for realtime colorimetric quantitative LAMP analysis from timelapse images of reaction chambers."),
  br()),
  ## Main
  fluidRow( ),
                 # Tabs
                 tabsetPanel(
                   
                   ## Select ROI
                   tabPanel('Image Processing',
                            sidebarLayout(
                              sidebarPanel(
                            tags$br(),
                            
                          
                            fileInput('image1', 'UPLOAD IMAGES', multiple = TRUE),
                            bsTooltip('image1',
                                      'Select and upload the chamber images.',
                                      'right', options = list(container = "body")),
                            tags$hr(),
                            
                        
                            
                           tags$p('Adjust parameters to remove the background,'),
                            sliderInput('threshold', 'Threshold', 1, 99, 29, 1),
                            bsTooltip('threshold',
                                      'Choose a threshold for excluding the background,',
                                      'right', options = list(container = "body")),
                            
                           
                                   sliderInput('shrink', 'Shrink', 1, 10, 1, 1),
                            bsTooltip('shrink',
                                      'Shrink the selected area by eroding the bounderies around it.',
                                      'right', options = list(container = "body")),
                            
                                   sliderInput('grow', 'Grow', 1, 10, 1, 1),
                            bsTooltip('grow',
                                      'Grow the selected area by dilating the bounderies around it.',
                                      'right', options = list(container = "body"))),
                            mainPanel(
                            plotOutput("image_plot", height = "3000px"),
                            tags$h4('Hue value table'),
                            tableOutput('table')
                            
                   ))),
                   
                   ## Realamp 
                   
                   tabPanel('Data Processing',
                            sidebarLayout(
                              sidebarPanel(
                                tags$br(),
                                radioButtons('curve', 'MODE OF ANALYSIS', choices = c("Standard Calibration Curve","Quantification of Real Sample")),
                                tags$hr(),
                        
                                numericInput('dil', 'Lowest concentration', 1,min = 0, max = NA, step = NA ),
                                numericInput('dil_factor', 'Dilution factor', 10,min = 0, max = NA, step = NA),
                                tags$hr(),
                                
                                
                                fileInput('StdCurve', 'Upload the Calibration Curve',
                                          multiple = FALSE)),
                              
                              mainPanel(
                                fluidRow(
                                plotOutput('graph_plot', height = "700px"),
                                tableOutput('table2'),
                                tags$br(),
                                
                                downloadButton("qLAMP.csv", "Save as Calibration Curve")
                              ))))
                            
                   
                   
                   
           ))


# Define server
server <- function(input, output) {
  # intiate interactive values
  values <- reactiveValues()
  
  # load images
  img1 <- reactive({
    image_load(input$image1$datapath[order(input$image1$name, decreasing = FALSE)])
    
  })
  
  ## calculate the pixset
  px <- reactive({
    roi_select(img1(),
               threshold = input$threshold,
               shrink = input$shrink,
               grow = input$grow)
  })
  
  ## calculate correlation
  b1 <- reactive({
    roi_check(px())  
  })
  
  # Output Views
  ## Select ROI
  
  # plots
  output$image_plot <- renderPlot({
    req(input$image1)
    n <- length(input$image1$name)/2
    par(mfrow=c(n+1,4), mar = rep(1, 4))
    roi_show(px())
  })
  
  s1 <- reactive({
    read.csv(input$StdCurve$datapath)
  })
  dil1 <- reactive({
    c(input$dil,
      input$dil*input$dil_factor^1,
      input$dil*input$dil_factor^2,
      input$dil*input$dil_factor^3,
      input$dil*input$dil_factor^4,
      input$dil*input$dil_factor^5)
  })
 
  ## sigmoid plot
  output$graph_plot <- renderPlot({
    req(input$image1)
    par(mfrow=c(2,1))
    rlamp(df = b1(), standardcurve = input$curve, df2 = s1(), dil2 = dil1(), thr_2 = input$thr_1, thr_01 = input$thr_0)
    
})
  ## G/R ratio table
  output$table <- renderTable({
    req(input$image1)
    b1()
  })
  output$table2 <- renderTable({
    req(input$image1)
    rlamp(df = b1(), standardcurve = input$curve, df2 = s1(), dil2 = dil1(), thr_2 = input$thr_1, thr_01 = input$thr_0)
  })
  # download plot

  # download table
  output$RLAMP_data.csv <- downloadHandler(
    #specify the filename
    filename = function() {
      paste("qLAMP", ".csv", sep = "")
    },
    
    content = function(file) {
      req(input$image1)
      write.csv(b1(), file, row.names = FALSE)
    })
  
}

# Run the application
shinyApp(ui = ui, server = server)
