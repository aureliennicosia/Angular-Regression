shinyServer(function(input, output) {
  
  ### Argument names:
  ArgNames <- reactive({
    Names <- names(formals(input$readFunction)[-1])
    Names <- Names[Names!="..."]
    return(Names)
  })
  
  # Argument selector:
  output$ArgSelect <- renderUI({
    if (length(ArgNames())==0) return(NULL)
    
    selectInput("arg","Argument:",ArgNames())
  })
  
  ## Arg text field:
  output$ArgText <- renderUI({
    fun__arg <- paste0(input$readFunction,"__",input$arg)
    
    if (is.null(input$arg)) return(NULL)
    
    Defaults <- formals(input$readFunction)
    
    if (is.null(input[[fun__arg]]))
    {
      textInput(fun__arg, label = "Enter value:", value = deparse(Defaults[[input$arg]])) 
    } else {
      textInput(fun__arg, label = "Enter value:", value = input[[fun__arg]]) 
    }
  })
  
  
  ### Data import:
  Dataset <- reactive({
    if (is.null(input$file)) {
      # User has not uploaded a file yet
      return(data.frame())
    }
    
    args <- grep(paste0("^",input$readFunction,"__"), names(input), value = TRUE)
    
    argList <- list()
    for (i in seq_along(args))
    {
      argList[[i]] <- eval(parse(text=input[[args[i]]]))
    }
    names(argList) <- gsub(paste0("^",input$readFunction,"__"),"",args)
    
    argList <- argList[names(argList) %in% ArgNames()]
    
    Dataset <- as.data.frame(do.call(input$readFunction,c(list(input$file$datapath),argList)))
    return(Dataset)
  })
  fitCons  <- reactive({
    ifelse( (is.null(input$file) ||is.null(input$formula.reg) || length(input$formula.reg)==0 ),
      return(NULL),
    return(consensus(input$formula.reg, data = Dataset())))
  })
  
  fitAng  <- reactive({
    ifelse( (is.null(input$file) ||is.null(input$formula.reg) || length(input$formula.reg)==0 ),
            return(NULL),
            return(angular(input$formula.reg, data = Dataset())))
  })
  
  trajectoire <- reactive({
    if( (is.null(input$file) || is.null(input$varsDir) ||is.null(input$varsDist))){
            return(NULL)}
    y <- Dataset()[,input$varsDir]
    d <- Dataset()[,input$varsDist]
    n <- length(y)
    position <- matrix(0,n+1,2)
    position[1,] <- c(0,0)
    for (i in (1:n)){
      position[i+1,] <-  position[i,] +d[i]*c(cos(y[i]),sin(y[i]))
    }   
    return(position)
  })
  
  
  # Select variables:
  
  output$plotTrajectoire <- renderPlot({
    if((is.null(input$varsDir) || is.null(input$varsDist) ||
              length(input$varsDir) !=1 ||length(input$varsDist) !=1)){ return(NULL)}
    traj <- as.data.frame(trajectoire())
    colnames(traj) <- c("Xpos","Ypos")
    p <- ggplot(data= traj, aes(x=Xpos, y = Ypos)) + geom_path()
    return(p)
    
  })
  
  output$varselect <- renderUI({
    
    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)
    
    # Variable selection:    
    selectInput("vars", "Variables to use:",
                names(Dataset()), names(Dataset()), multiple =TRUE)            
  })
  output$varselectDirection <- renderUI({
    
    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)
    
    # Variable selection:    
    selectInput("varsDir", "Direction:",
                names(Dataset()), names(Dataset()), multiple =TRUE)            
  })
  
  output$varselectDistance <- renderUI({
    
    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)
    
    # Variable selection:    
    selectInput("varsDist", "Distance:",
                names(Dataset()), names(Dataset()), multiple =TRUE)            
  })
  
  # Summary statistic
  output$polarTable <- renderDataTable({
    if( (is.null(input$file) || is.null(input$varsDir) ||is.null(input$varsDist))){
      return(NULL)}
    y <- Dataset()[,input$varsDir]
    d <- Dataset()[,input$varsDist]
    n <- length(y)
    position <- matrix(0,n+1,2)
    position[1,] <- c(0,0)
    for (i in (1:n)){
      position[i+1,] <-  position[i,] +d[i]*c(cos(y[i]),sin(y[i]))
    }   
    return(as.data.frame(position))
 
    
    
  })
  # Show table:
  output$table <- renderTable({
    
    if (is.null(input$file) ) return(NULL)
    
    return(head(Dataset(),10))
  })
  
  output$regressionCons <- renderPrint({
    if(is.null(input$file) ||is.null(input$formula.reg) || length(input$formula.reg)==0 )
      return(NULL);
    print(fitCons())
    cat("\nBeta parameters:\n")
    print(fitCons()$parambeta)
    
    
  })
  
  output$goodnessOfFitCons <- renderPlot({
    if(is.null(input$file) ||is.null(input$formula.reg) || length(input$formula.reg)==0 )
      return(NULL);
    gg_qq(fitted = fitCons()$mui, resid = sqrt(fitCons()$long)*(fitCons()$y - fitCons()$mui))
    
    
  })
  
  output$regressionAng <- renderPrint({
    if(is.null(input$file) ||is.null(input$formula.reg) || length(input$formula.reg)==0 )
      return(NULL);
    print(fitAng())
    
    
  })
  
  output$goodnessOfFitAng <- renderPlot({
    if(is.null(input$file) ||is.null(input$formula.reg) || length(input$formula.reg)==0 )
      return(NULL);
    gg_qq(fitted = fitAng()$mui, resid =(fitAng()$y - fitAng()$mui))
    
    
  })
  
  
  ### Download dump:

  output$downloadDump <- downloadHandler(
    filename = "Rdata.R",
    content = function(con) {
    
    assign(input$name, Dataset()[,input$vars,drop=FALSE])
      
    dump(input$name, con)
    }
  )
  
  ### Download save:

  output$downloadSave <- downloadHandler(
    filename = "Rdata.RData",
    content = function(con) {
      
      assign(input$name, Dataset()[,input$vars,drop=FALSE])
      
      save(list=input$name, file=con)
    }
  )
  
  
  #### Download report
  output$downloadReport <- downloadHandler(
    filename = 'myreport.pdf',
    
    content = function(file) {
      out = knit2pdf('input.Rnw', clean = TRUE)
      file.rename(out, file) # move pdf to file for downloading
    },
    
    contentType = 'application/pdf'
  )
  
})


