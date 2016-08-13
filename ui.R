
shinyUI(pageWithSidebar(
  
  # Header:
  headerPanel("Angular regression model"),
  
  # Input in sidepanel:
  sidebarPanel(
    tags$style(type='text/css', ".well { max-width: 20em; }"),
    # Tags:
    tags$head(
      tags$style(type="text/css", "select[multiple] { width: 100%; height:10em}"),
      tags$style(type="text/css", "select { width: 100%}"),
      tags$style(type="text/css", "input { width: 19em; max-width:100%}")
    ),
    
    # Select filetype:
    selectInput("readFunction", "Function to read data:", c(
      # Base R:
      "read.table",
      "read.csv",
      "read.csv2",
      "read.delim",
      "read.delim2"

      # # foreign functions:
      # "read.spss",
      # "read.arff",
      # "read.dbf",
      # "read.dta",
      # "read.epiiinfo",
      # "read.mtp",
      # "read.octave",
      # "read.ssd",
      # "read.systat",
      # "read.xport",
      # 
      # # Advanced functions:
      # "scan",
      # "readLines"
      )),
    
    # Argument selecter:
    htmlOutput("ArgSelect"),
    
    # Argument field:
    htmlOutput("ArgText"),
    
    # Upload data:
    fileInput("file", "Upload data-file:"),
    
    # Variable selection:
    htmlOutput("varselectDirection"),
    htmlOutput("varselectDistance"),
    
    br(),
    # Only show this panel if the plot type is a histogram
    
  
    
    
    textInput("formula.reg","Formula:","y~.") #,
    #radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
     #            inline = TRUE),
    #downloadButton('downloadReport')

  
  ),
  
  # Main:
  mainPanel(tabsetPanel(

    tabPanel("Data and trajectory",tableOutput("table"),
             plotOutput("plotTrajectoire")),
    tabPanel("Angular regression model", 
             verbatimTextOutput("regressionAng"),
             plotOutput("goodnessOfFitAng")
             ),
    tabPanel("Consensus regression model", 
             verbatimTextOutput("regressionCons"),
             plotOutput("goodnessOfFitCons")
    )
    
  )
    
  )
))