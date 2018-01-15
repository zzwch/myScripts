#setwd("K:/307Hospital/Workspace/EC-RNA-Seq/sprint/20171211/fromHouSiyuan/")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, shinydashboard, tiff, EBImage, SDMTools)
if(!all(p_isloaded(shiny, shinydashboard, tiff, EBImage, SDMTools))) 
  stop("some package is not loaded! You may contact lizc07@vip.qq.com")

## app.R ##
library(shiny)
library(shinydashboard)
library(tiff)
library(EBImage)
require(SDMTools)
options(shiny.maxRequestSize=1000*1024^2, stringsAsFactors = F)

## custom functions
tiffBinary <- function(red, green, blue, red.cutoff = 0, green.cutoff = 0, blue.cutoff = 0){
  require(tiff)
  red[red < red.cutoff] <- 0
  green[green < green.cutoff] <- 0
  blue[blue < blue.cutoff] <- 0
  rgb <- red
  rgb[,,2] <- green[,,2]
  rgb[,,3] <- blue[,,3]
  return(list(rgb = rgb, red = red, green = green, blue = blue))
}
triplet <- function(mat, cutoff = 0.5){
  ind <- which(mat > cutoff, arr.ind = T)
  return(data.frame(x = ind[,2], y = -ind[,1]))
}
overlapStat <- function(list){
  nn <- length(list)
  sapply(list, function(x){
    sapply(list, function(y) {
      length(intersect(x, y))
    })
  })
}
######
ui <- dashboardPage(title = "Widget for Confocal Image",
  dashboardHeader(title = "Widget for Confocal Image", titleWidth = "100%", disable = F),
  dashboardSidebar(width = 200,
    sidebarMenu(id = NULL,
      menuItem("Confocal Images",
               fileInput(inputId = "tiff.red", label = NULL, multiple = F, accept = c(".tiff",".tif"), 
                                width = "100%", buttonLabel = "Red", placeholder = "No file selected"),
               
                      fileInput(inputId = "tiff.green", label = NULL, multiple = F, accept = c(".tiff",".tif"), 
                                width = "100%", buttonLabel = "Green", placeholder = "No file selected"),
              
                      fileInput(inputId = "tiff.blue", label = NULL, multiple = F, accept = c(".tiff",".tif"), 
                                width = "100%", buttonLabel = "Blue", placeholder = "No file selected"),
               
               
                      sliderInput(inputId = "cutoff.red", label = NULL, min = 0, max = 1, value = 0.5, 
                                  step = 0.01, round = -2, ticks = T, width = "100%", pre = "Red "),
            
                      sliderInput(inputId = "cutoff.green", label = NULL, min = 0, max = 1, value = 0.5, 
                                  step = 0.01, round = -2, ticks = T, width = "100%", pre = "Green "),
              
                      sliderInput(inputId = "cutoff.blue", label = NULL, min = 0, max = 1, value = 0.5, 
                                  step = 0.01, round = -2, ticks = T, width = "100%", pre = "Blue ")
               ),
      menuItem("Display Image",
              
                      sliderInput(inputId = "resize", label = "Resize Image", min = 1, max = 10, value = 10, 
                                  step = 1, round = T, ticks = T, width = "100%", pre = "Ratio "),
               
              
                      sliderInput(inputId = "plot.width", label = "Width", min = 600, max = 2000, value = 600, 
                                  step = 1, round = T, ticks = T, width = "100%", pre = NULL),
             
                      sliderInput(inputId = "plot.height", label = "Height", min = 400, max = 1600, value = 500, 
                                  step = 1, round = T, ticks = T, width = "100%", pre = NULL)
               ),
      menuItem("Region Statistics",
               actionButton(inputId = "action.new", label = "Start to Select", width = 160),
               actionButton(inputId = "action.submit", label = "Complete selection", width = 160),
               actionButton(inputId = "action.stat", label = "Statistics", width = 160)
               )
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                               .box { margin-bottom: 1; } 
                              [class*="col-lg-"],[class*="col-md-"],
                              [class*="col-sm-"],[class*="col-xs-"]{
                              padding-right:1 !important;
                              padding-left:0 !important;
                              }
                              '))),
    column(width = 8,
           tabBox(id = "rasterImage", selected = "rasterImage.rgb", title = "Raster Images", width = "100%", height = NULL, side = "left",
                  #
                  tabPanel(title = "RGB", value = "rasterImage.rgb",uiOutput("ui.rgb")),
                  tabPanel(title = "Red", value = "rasterImage.red",uiOutput("ui.red")),
                  tabPanel(title = "Green", value = "rasterImage.green",uiOutput("ui.green")),
                  tabPanel(title = "Blue", value = "rasterImage.blue",uiOutput("ui.blue"))
                  )
           ),
    column(width = 4,
           box(width = "100%",
               h3("Statistics Results:"),
               tableOutput("stat.table"),
               verbatimTextOutput("stat.text"),
               h4("Copy the following  to clipboard"),
               verbatimTextOutput("stat.copy")
               )
           )

    )
)

server <- function(input, output) {
  tiff <- reactive({
    # input$tiff.* will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    req(input$tiff.red, input$tiff.green, input$tiff.blue)
    
    tiff1 <- readTIFF(source = input$tiff.red$datapath, native = F)
    tiff2 <- readTIFF(source = input$tiff.green$datapath, native = F)
    tiff3 <- readTIFF(source = input$tiff.blue$datapath, native = F)
    
    tiffBinary(tiff1, tiff2, tiff3, input$cutoff.red, input$cutoff.green, input$cutoff.blue)
  })
  # rasterImage
  action.submit.plot.rgb <- eventReactive(eventExpr = c(input$action.new, input$action.submit, input$resize, tiff()), {
    raster <- tiff()$rgb
    dims <- dim(raster)
    raster.resize <- EBImage::resize(EBImage::as.Image(raster), dims[1]/input$resize, dims[2]/input$resize)@.Data 
    
    op <- par(bg = "black", mai = c(0,0,0,0), mar = c(0,0,0,0))
    plot(c(0, dims[2]), c(0, -dims[1]), type = "n", xlab = "", ylab = "", 
         xaxs = "i", yaxs = "i", xlim = c(0, dims[2]), ylim = c(-dims[1],0))
    rasterImage(raster.resize, 0, -dims[1], dims[2], 0, interpolate = F)
    points(rbind(click_data, click_data[1,]), col = "orange", type = "l", new = F)
    par(op)
  })
  output$plot.rgb <- renderPlot({
    action.submit.plot.rgb()
  },res = 150)
  output$ui.rgb <- renderUI(plotOutput(outputId = "plot.rgb", click = "plot.rgb.click", width = input$plot.width, height = input$plot.height))
  
  output$plot.red <- renderPlot({
    raster <- tiff()$red
    dims <- dim(raster)
    raster.resize <- EBImage::resize(EBImage::as.Image(raster), dims[1]/input$resize, dims[2]/input$resize)@.Data 
    
    op <- par(bg = "black", mai = c(0,0,0,0), mar = c(0,0,0,0))
    plot(c(0, dims[2]), c(0, -dims[1]), type = "n", xlab = "", ylab = "",
         xaxs = "i", yaxs = "i", xlim = c(0, dims[2]), ylim = c(0,-dims[1]))
    rasterImage(raster.resize, 0, -dims[1], dims[2], 0, interpolate = F)
    par(op)
  },res = 150)
  output$ui.red <- renderUI(plotOutput(outputId = "plot.red", click = "plot.red.click", width = input$plot.width, height = input$plot.height))
  
  output$plot.green <- renderPlot({
    raster <- tiff()$green
    dims <- dim(raster)
    raster.resize <- EBImage::resize(EBImage::as.Image(raster), dims[1]/input$resize, dims[2]/input$resize)@.Data 
    
    op <- par(bg = "black", mai = c(0,0,0,0), mar = c(0,0,0,0))
    plot(c(0, dims[2]), c(0, -dims[1]), type = "n", xlab = "", ylab = "",
         xaxs = "i", yaxs = "i", xlim = c(0, dims[2]), ylim = c(0,-dims[1]))
    rasterImage(raster.resize, 0, -dims[1], dims[2], 0, interpolate = F)
    par(op)
  },res = 150)
  output$ui.green <- renderUI(plotOutput(outputId = "plot.green", click = "plot.green.click", width = input$plot.width, height = input$plot.height))
  
  output$plot.blue <- renderPlot({
    raster <- tiff()$blue
    dims <- dim(raster)
    raster.resize <- EBImage::resize(EBImage::as.Image(raster), dims[1]/input$resize, dims[2]/input$resize)@.Data 
    
    op <- par(bg = "black", mai = c(0,0,0,0), mar = c(0,0,0,0))
    plot(c(0, dims[2]), c(0, -dims[1]), type = "n", xlab = "", ylab = "",
         xaxs = "i", yaxs = "i", xlim = c(0, dims[2]), ylim = c(0,-dims[1]))
    rasterImage(raster.resize, 0, -dims[1], dims[2], 0, interpolate = F)
    par(op)
  },res = 150)
  output$ui.blue <- renderUI(plotOutput(outputId = "plot.blue", click = "plot.blue.click", width = input$plot.width, height = input$plot.height))
  
  # click events
  click_data <- data.frame(x = as.numeric(NULL), y = as.numeric(NULL))
  makeReactiveBinding("click_data")
  
  observeEvent(input$action.new, {
    click_data <<- data.frame(x = as.numeric(NULL), y = as.numeric(NULL))
  })
  
  observeEvent(input$plot.rgb.click,{
    #print(str(input$plot.rgb.click))
    click_data <<- rbind(click_data, data.frame(x = input$plot.rgb.click$x, y = input$plot.rgb.click$y))
  })
  
  observeEvent(input$action.stat, {
    ## select points
    print(click_data)
    
    if(dim(click_data)[1] < 3){
      points.red <- triplet(tiff()$rgb[,,1], input$cutoff.red)
      points.green <- triplet(tiff()$rgb[,,2], input$cutoff.green)
      points.blue <- triplet(tiff()$rgb[,,3], input$cutoff.blue)
      points.red$pip <- points.green$pip <- points.blue$pip <- 1
    }else{
      points.red <- SDMTools::pnt.in.poly(pnts = triplet(tiff()$rgb[,,1], input$cutoff.red), poly.pnts = click_data)
      points.green <- SDMTools::pnt.in.poly(pnts = triplet(tiff()$rgb[,,2], input$cutoff.green), poly.pnts = click_data)
      points.blue <- SDMTools::pnt.in.poly(pnts = triplet(tiff()$rgb[,,3], input$cutoff.blue), poly.pnts = click_data)
    }
    
    points.red.str <- paste(points.red[which(points.red$pip == 1),1], points.red[which(points.red$pip == 1),2], "_")
    points.green.str <- paste(points.green[which(points.green$pip == 1),1], points.green[which(points.green$pip == 1),2], "_")
    points.blue.str <- paste(points.blue[which(points.blue$pip == 1),1], points.blue[which(points.blue$pip == 1),2], "_")
    statRes <- matrix(NA, 3, 3, dimnames = list(c("Red", "Green", "Blue"), c("Red", "Green", "Blue")))
    statRes[1, ] <- c(length(points.red.str), 
                      length(intersect(points.red.str, points.green.str)), 
                      length(intersect(points.red.str, points.blue.str)))
    statRes[2,2:3] <- c(length(points.green.str), 
                        length(intersect(points.green.str, points.blue.str)))
    statRes[3,3] <- c(length(points.blue.str))
    output$stat.table <- renderTable(statRes)
    
    ep <- statRes[1,1]
    ep235 <- statRes[1,2]
    ep6 <- statRes[1,3]
    ep7 <- ep - ep235 - ep6
    textRes <- paste0("(Red)              All EP: ", ep, "\n",
                      "(Red & Green)       EP235: ", ep235, "\n",
                      "(Red & Blue)          EP6: ", ep6, "\n",
                      "(Red & ^(Green|Blue)) EP7: ", ep7)
    output$stat.text <- renderText(textRes)
    copyRes <- paste(ep,ep235,ep6,ep7, sep = "\t")
    output$stat.copy <- renderText(copyRes)
  })
}

shinyApp(ui, server)