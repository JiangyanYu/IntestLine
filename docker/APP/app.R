# Load libraries ----
library(shiny)

rm(list = ls())
library(igraph)
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(ggtext)
library(data.table)
library(readr)

source(file = "functions.R")

####
# initialize global variable to record selected (clicked) rows
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)
columns = c("x", "y", "pos", "order")
selected_backbone = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(selected_backbone) = columns
GLOBAL_SEED = 12345
set.seed(GLOBAL_SEED)

vv_DAPI = data.frame()


rad2deg <- function(rad) {
  (rad * 180) / (pi)
}


# UI ----

ui = shinyUI(fluidPage(
  
  titlePanel(title = div(img(src = "IntestLine-Logo.png", height = 100, width = 100),
                         "IntestLine - Digitally unroll your intestinal images")),
  
  ## panel upper left:  upload data and select regions ---- 
  column(
    wellPanel( width = 6,
               style = "background:#b3b3b3" ,
               tags$head(tags$style(
                 HTML("hr {border-top: 1px solid #000000;}")
               )),
               tags$head(
                 tags$style(
                   type = "text/css",
                   "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          "
                 )
               ),
               # Colors: #CCFF66
               p(),
               #br(),
               tags$b(tags$em(p("Report errors ", 
                                tags$a(href ="https://www.github.com/JiangyanYu/IntestLine",
                                       tags$b(tags$span(style = "color: #cc6699", icon("github"), "GitHub")),
                                       target = "_blank")
                                # tags$a(href ="mailto:jyuu@uni-bonn.de",
                                #        tags$b(tags$span(style = "color: #666699", icon("envelope-open"), "Email")),
                                #        target = "_blank")
               ))),
               p(),
               #br(),
               tags$b(tags$em(p("Cite IntestLine:"))),
               
               conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                                tags$div("Loading...", id = "loadmessage")),
               
               #br(),
               p(),
               tags$b(tags$em(p("Download demo dataset:", tags$a(href = "https://github.com/JiangyanYu/IntestLine",
                                                                 tags$b(tags$span(style = "color: #6699ff", icon("github"), "GitHub")),
                                                                 target = "_blank")
               ))),
               p(),
               
               #br(),
               hr(),
               fileInput(
                 'codex_upload',
                 'Upload CODEX .csv file (comma separated, max. 1Gb)',
                 accept = c('text/csv',
                            'text/comma-separated-values,text/plain',
                            '.csv')
               ),
               
               fileInput(
                 'outline_upload',
                 'Upload saved outline file if available (comma separated)',
                 accept = c('text/csv',
                            'text/comma-separated-values',
                            '.csv')
               ),
               
               fluidRow(
                 column( width = 3,
                         offset = 1,
                         actionButton(
                           "uploadData",
                           "1. Upload Data",
                           icon = icon("upload"),
                           style = "background:#d279a6;color:white;"
                         )
                 ),
                 column(
                   width = 6,
                   offset = 1,
                   numericInput(
                     "pointSize",
                     "Select plotting point size (1 to 100):",
                     value = 50,
                     min = 1,
                     max = 100,
                   )
                 )
               ),
               fluidRow(
                 plotOutput("plotOverview",
                            height = "655px",
                            width = "100%"),
                 hr(),
                 fluidRow(
                   column(
                     width = 6,
                     numericInput("chooseX0",
                                  "Center point X0:",
                                  min = 0,
                                  value = 15975,)
                   ),
                   column(
                     width = 6,
                     numericInput("chooseY0",
                                  "Center point Y0:",
                                  min = 0,
                                  value = 12710,)
                   )
                 ),
                 fluidRow(
                   column(width = 6,
                          numericInput(
                            "chooseA",
                            "Ellipse radius a:",
                            min = 1,
                            value = 8400
                          )),
                   column(width = 6,
                          numericInput(
                            "chooseB",
                            "Ellipse radius b:",
                            min = 1,
                            value = 8850
                          ))
                 ),
                 column(width = 12,
                        offset = 4,
                        actionButton(
                          "selectArea", "1.5. Select Area",
                          icon = icon("chart-area"),
                          style = "background:#ff6666;color:white;"
                        )),
                 br(),
                 br(),
                 p(),
                 p(),
                 p(),
                 plotOutput("plotROI",
                            height = "655px",
                            width = "100%")
               )
    ),
    width = 6),
  
  ## panel upper right: select backbone ----
  column( width = 6,
          wellPanel(
            style = "background:#99ccff",
            sliderInput(
              "percentCells",
              "Select percent of cells to visualize in Shiny:",
              min = 1,
              max = 100,
              value = 40,
              step = 1
            ),
            fluidRow(column(
              width = 6,
              offset = 3,
              actionButton(
                "uploadROI",
                "2. View ROI and then select backbone",
                icon = icon("hand-point-down"),
                style = "background:#6699ff;color:white;"
              )
            )),
            br(),
            h5(HTML(
              paste(
                tags$b(span("Please selecte backbone points from proximal to
                 distal in order. The order you select your points will be used
                 to order your base layer. You can skip this step if you have 
                 uploaded one outline file and go to step 2.5.", 
                            style = "color:#993333;font-size:18;")),
                #'A minimum number of outline points must be selected for
                #' effective analysis.',
                sep = "<br/>"
              )
            ), align = "center"),
            #hr(),
            plotOutput("plotSelectBackbone", click = "clicked",
                       height = "550px",width = "100%"),
            br(),
            fluidRow(column(
              width = 4,
              offset = 4,
              downloadButton('downloadData', 'Download Backbone',
                             style = "background:#cccc00;color:white;")
            )),
            hr(),
            fluidRow(
              column(width = 4,
                     offset = 4,
                     actionButton("backbonePlotButton",
                                  "2.5. Plot backbone on ROI",
                                  icon = icon("diagram-project"),
                                  style = "background:#669900;color:white;")
              )
            ),
            br(),
            plotOutput("backbonePlot",height = "550px",width = "100%"),
            hr(),
            fluidRow(
              column(
                width =4,
                offset = 4,
                actionButton("qualityControl", "QC of backbone selection",  
                             icon = icon("clipboard-check"),
                             style = "background:#ff9966;color:white;")
              )
            ),
            br(),
            plotlyOutput("quality_plot")
          )
  ),
  ## panel middle: linear structure ----
  column(width = 12,
         wellPanel(
           style = "background:#99ffcc",
           fluidRow(column(
             width = 2,
             offset = 2,
             actionButton(
               "unrollData",
               "3. Unroll image",
               icon = icon("digital-ocean"),
               style = "background:#cc6699;color:white;"
             )
           ),
           column(
             width = 4,
             offset = 3,
             downloadButton(
               "downloadStretch",
               "Download converted linear structure",
               style = "background:#999966;color:white;"
             ))
           ),
           # hr(),
           # h5(HTML(
           #   paste(
           #     "All Markers are visualized according to their strength, overlayed over all (100%) of DAPI2+ (Marker>0.0) cells.",
           #     sep = "<br/>"
           #     
           #   )
           # ), align = "center"),
           # hr(),
           p(),
           fluidRow(
             column( width = 6,
                     numericInput(
                       "selectDistance",
                       "Distance treshold (Upload Data! To remove points behind backbone.):",
                       value = Inf,
                     )
             ),
             column( width = 6,
                     numericInput(
                       "selectZscore",
                       "Z-score treshold (Upload Data! To remove projection outliers.):",
                       value = Inf,
                     )
             )
           ),
           fluidRow(
             column(
               width = 3,
               plotOutput("distance_loss_hist")
             ),
             column(
               width = 3,
               plotOutput("distance_loss")
             ),
             # ),
             # p(),
             # fluidRow(
             column(
               width = 3,
               plotOutput("zscore_loss_hist")
             ),
             column(
               width = 3,
               plotOutput("zscore_loss")
             )
           ),
           # hr(),
           p(),
           fluidRow(
             column(width = 12,
                    plotOutput("line_plot",
                               height = "300px",width = "100%"
                    )
             )
             # column(width = 6,
             #        plotOutput("lost_cells")
             # )
           ),
           ## panel bottom: linear structure ----
           hr(),
           fluidRow(column(
             width = 2,
             offset = 2,
             actionButton(
               "unrollData",
               "3.5. Overlay markers",
               icon = icon("layer-group"),
               style = "background:#8585ad;color:white;"
             )
           ),
           column( width = 4,
                   offset = 2,
                   selectizeInput("selectMarker", "Choose marker:", choices = "Upload Data!")
           )
           ),
           fluidRow(
             column(width = 12,
                    plotOutput("marker_plot",
                               height = "300px",width = "100%"
                    )
             )
           )
           # hr(),
           # h5(HTML(
           #   paste(
           #     "Developed by: Dr. Jiangyan Yu and Altay Yuzeir",
           #     'Department of Quantitative Systems biology at LIMES Institute in Bonn, Germany',
           #     "Under the leadership of Prof. Dr. Andreas Schlitzer, April 2022",
           #     sep = "<br/>"
           #     
           #   )
           # ), align = "center"),
           # hr()
         ))
))

# Server ----
server = shinyServer(function(input, output, session) {
  data_upload <- reactive({
    inFile <- input$codex_upload
    if (is.null(inFile))
      return(NULL)
    df <-
      fread(
        file = inFile$datapath,
        header = TRUE,
        sep = ",",
        data.table = F,
        check.names = T
      )
    return(df)
  })
  outline_upload <- reactive({
    inFile <- input$outline_upload
    if (is.null(inFile))
      return(NULL)
    df <-
      fread(
        file = inFile$datapath,
        header = TRUE,
        sep = ",",
        data.table = F,
        check.names = T
      )
    return(df)
  })
  ## Upload main data ----
  observeEvent(input$uploadData, {
    input_image = data_upload()
    output$plotOverview = renderPlot({
      plot.new()
    })
    output$plotROI = renderPlot({
      plot.new()
    })
    output$plotSelectBackbone = renderPlot({
      plot.new()
    })
    output$backbonePlotButton <- renderPlot({
      plot.new()
    })
    output$line_plot = renderPlot({
      plot.new()
    })
    output$lost_cells = renderPlot({
      plot.new()
    })
    output$backbonePlot = renderPlot({
      plot.new()
    })
    output$distance_loss = renderPlot({
      plot.new()
    })
    output$zscore_loss = renderPlot({
      plot.new()
    })
    output$distance_loss_hist = renderPlot({
      plot.new()
    })
    output$zscore_loss_hist = renderPlot({
      plot.new()
    })
    output$marker_plot = renderPlot({
      plot.new()
    })
    output$plotOverview = renderPlot({
      plot(
        #asp = 1,
        input_image$x,
        input_image$y,
        cex = input$pointSize/100,
        pch = 16,
        col = "black",
        ylab = "CODEX coordinate y",
        xlab = "CODEX coordinate x",
        main = "Overview of the uploaded image"
      )
    })
    observeEvent(input$selectArea,
                 {## define center of the whole image
                   input_x0 = input$chooseX0
                   input_y0 = input$chooseY0
                   input_A = input$chooseA
                   input_B = input$chooseB
                   output$plotROI = renderPlot({
                     plot(
                       #asp = 1,
                       input_image$x,
                       input_image$y,
                       cex = input$pointSize/100,
                       pch = 16,
                       ylab = "CODEX coordinate y",
                       xlab = "CODEX coordinate x",
                       main = "Region of interest (ROI)"
                     )
                     plotrix::draw.ellipse(
                       x = input_x0,
                       y = input_y0,
                       a = input_A,
                       b = input_B,
                       lwd = 2,
                       border = "#ff6666"
                     )
                     points(input_x0,input_y0,cex = 2,pch=16,col="#ff6666")
                   })
                 },
                 ignoreInit = T,
                 priority = -1)
  })
  ## Select backbone ----
  observeEvent(input$uploadROI, {
    image_ROI = data_upload()
    x0 = input$chooseX0
    y0 = input$chooseY0
    ex = input$chooseEx
    ey = input$chooseEy
    a = input$chooseA
    b = input$chooseB
    ### select ROI ----
    image_ROI$distance2center = (image_ROI$x - x0) ^ 2 / a ^ 2 + (image_ROI$y - y0) ^ 2 / b ^ 2
    image_ROI = filter(image_ROI, distance2center < 1)
    image_ROI$pos = paste0(image_ROI$x,"_",image_ROI$y)
    output$plotSelectBackbone <- renderPlot({
      plot.new()
    })
    output$line_plot <- renderPlot({
      plot.new()
    })
    output$lost_cells <- renderPlot({
      plot.new()
    })
    output$backbonePlotButton <- renderPlot({
      plot.new()
    })
    output$distance_loss = renderPlot({
      plot.new()
    })
    output$zscore_loss = renderPlot({
      plot.new()
    })
    output$distance_loss_hist = renderPlot({
      plot.new()
    })
    output$zscore_loss_hist = renderPlot({
      plot.new()
    })
    ### select backbone ----
    
    upload_backbone = outline_upload()
    selected_backbone = selected_backbone[0, ]
    percent_cell = input$percentCells
    set.seed(GLOBAL_SEED)
    shiny_plot = sample_n(image_ROI, size = round(percent_cell / 100 * nrow(image_ROI)))
    shiny_plot$order = c(1:nrow(shiny_plot))
    selected <- reactive({
      # add clicked
      selected_backbone <<-
        rbind(
          selected_backbone,
          nearPoints(
            shiny_plot,
            input$clicked,
            threshold = 10,
            maxpoints = 1
          )
        )
      # remove _all_ duplicates if any (toggle mode)
      # http://stackoverflow.com/a/13763299/3817004
      selected_backbone <<-
        selected_backbone[!(duplicated(selected_backbone) |
                              duplicated(selected_backbone, fromLast = TRUE)),]
      #str(selected_backbone)
      return(selected_backbone)
    })
    
    ### download backbone ----
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(
          'selected_backbone-',
          format(Sys.time(), "%d-%b-%Y_%Hh%Mmin"),
          '.csv',
          sep = ''
        )
      },
      content = function(file) {
        tmp = selected()
        tmp = select(tmp, "x", "y","pos", "order")
        write.table(x = tmp,
                    file,
                    sep = ",",
                    row.names = F)
      }
    )
    output$plotSelectBackbone <- renderPlot({
      plot.new()
    })
    output$line_plot <- renderPlot({
      plot.new()
    })
    output$lost_cells <- renderPlot({
      plot.new()
    })
    output$plotSelectBackbone <- renderPlot({
      ggplot(shiny_plot, aes(x = x, y = y)) +
        geom_point() +
        geom_point(data = selected(),
                   colour = "red",
                   size = 0.5) +
        xlab("CODEX coordinate x") + ylab("CODEX coordinate y")
    })
    ### View backbone plus image ----
    observeEvent(input$backbonePlotButton,{
      backbone_points = data.frame()
      if (isTruthy(upload_backbone)) {
        backbone_points <- upload_backbone
        if (nrow(backbone_points) == 0 | nrow(backbone_points) == 1) {
          showNotification("Please select backbone points or upload a 
          backbone file !", type = "error")
        } else {
          output$backbonePlot <- renderPlot({
            plot(
              image_ROI$x,
              image_ROI$y,
              cex = input$pointSize/100,
              pch = 16,
              col = "black",
              ylab = "CODEX coordinate y",
              xlab = "CODEX coordinate x",
              main = "Overview of the uploaded image and selected backbone"
            )
            points(backbone_points$x,backbone_points$y,pch=16,col="red",cex = 0.5)
            showNotification("Use uploaded backbone points.", type = "message")
          })
        }
      } else {
        backbone_points <- selected()
        if (nrow(backbone_points) == 0 | nrow(backbone_points) == 1) {
          showNotification("Please select backbone points or upload a backbone file !", type = "error")
        } else {
          output$backbonePlot <- renderPlot({
            plot(
              image_ROI$x,
              image_ROI$y,
              cex = input$pointSize/100,
              pch = 16,
              col = "black",
              ylab = "CODEX coordinate y",
              xlab = "CODEX coordinate x",
              main = "Overview of the uploaded image and selected backbone"
            )
            points(backbone_points$x, backbone_points$y, pch=16,
                   col="red", cex = 0.5)
            showNotification("Use currently selected backbone points.",
                             type = "warning")
          })
        }
      }
    })
    output$backbonePlot = renderPlot({
      plot.new()
    })
    ### Quality control of backbone selection ----
    observeEvent(input$qualityControl,
                 {
                   output$quality_plot <- renderPlotly({
                     if (isTruthy(upload_backbone)) {
                       base_points1 <- upload_backbone
                     } else{
                       base_points1 <- selected()
                     }
                     # FIXME: Is the line below correct really?
                     base_points1$base_order = c(1:nrow(base_points1))
                     base_points1$Outline_points = base_points1$base_order
                     fig <- plot_ly(
                       base_points1,
                       x = ~ x,
                       y = ~ y,
                       z = ~ Outline_points,
                       marker = list(size = 3)
                     )
                     fig = fig %>% layout(xaxis = list(title = "CODEX coordinate x"),
                                          yaxis = list(title = "CODEX coordinate y"))
                     fig <-
                       fig %>% add_markers(hovertemplate = 'x=%{x}\ny=%{y}\nz=%{z}<extra></extra>')
                     fig
                   })
                 },
                 ignoreInit = T,
                 priority = -2)
    ## Unrolling ----
    observeEvent(input$unrollData,{
      input_image = data_upload()
      x0 = input$chooseX0
      y0 = input$chooseY0
      ex = input$chooseEx
      ey = input$chooseEy
      a = input$chooseA
      b = input$chooseB
      if (isTruthy(outline_upload())) {
        base_points = outline_upload()
      } else{
        base_points = selected()
      }
      ### Order base and calculate length 
      ordered_base = order_base(backbone_points = base_points)
      ordered_base$pos = paste0(ordered_base$x, "_", ordered_base$y)
      ordered_base$distance2center = (ordered_base$x - x0) ^ 2 / a ^ 2 + (ordered_base$y - y0) ^ 2 / b ^ 2
      ### Project all other points
      input_image$pos = paste0(input_image$x, "_", input_image$y)
      input_image$distance2center = (input_image$x - x0) ^ 2 / a ^ 2 + (input_image$y - y0) ^ 2 / b ^ 2
      input_image = subset(input_image, distance2center < 1)
      linear_image = project_points2base(backbone_points = ordered_base, query_points = input_image)
      
      ### Filtering steps  (distance plus outlier)
      ordered_base = zscore_per_backbone_point(converted_image = linear_image,backbone_points = ordered_base)
      linear_image = qc_zscore_outlier(converted_image = linear_image,x0=x0,y0=y0,backbone_points = ordered_base)
      ##############################
      ##
      #
      # We need to add "readr" package in install.R file for the FG Docker 
      linear_image_clean = subset(linear_image,note=="Successfully projected" )
      #linear_image_clean$nn_dist = parse_number(linear_image_clean$nn_dist)

      linear_image_clean$zscore = parse_number(linear_image_clean$zscore)
      
      download_data = merge(linear_image, input_image, by = "pos", no.dups = T, suffixes = c("","1"))
      download_data = select(download_data, -any_of(c("V1", "x1", "y1", "z1", "distance2center1",
                                                      "shortest_path_order", "nn_index")))
      download_data = rename(download_data, Thickness = nn_dist, Length = shortest_path_length)
      write.csv(download_data,"linear_image_clean.csv", row.names = F)
      updateNumericInput(session,
                         inputId = "selectDistance", 
                         label = paste0("Select distance filter (",
                                        round(min(linear_image_clean$nn_dist)),
                                        " to ",
                                        round(max(linear_image_clean$nn_dist))+1,
                                        ")"),
                         min = round(min(linear_image_clean$nn_dist)),
                         max = round(max(linear_image_clean$nn_dist))+1,
                         value = round(0.30*max(linear_image_clean$nn_dist))
      )
      
      updateNumericInput(session,
                         inputId = "selectZscore", 
                         label = paste0("Select z-score to filter (",
                                        round(min(linear_image_clean$zscore, na.rm = T)),
                                        " to ",
                                        round(max(linear_image_clean$zscore, na.rm = T))+1,
                                        ")"),
                         min = round(min(linear_image_clean$zscore, na.rm = T)),
                         max = round(max(linear_image_clean$zscore, na.rm = T))+1,
                         value = round(0.75*max(linear_image_clean$zscore, na.rm = T))
      )
      all_names = colnames(input_image)
      wrong_names = c("pos", "x", "y", "z", "distance2center", "nn_index",
                      "nn_dist", "nearest_Row", "nearest_Column",
                      "shortest_path_order", "shortest_path_length", "note",
                      "nn_dist", "zscore", "V1", "cell_id", "region",
                      "tile_num", "x1", "y1", "z1", "x_tile", "y_tile","size",
                      "tsne_x", "tsne_y", "homogeneity", "distance2center1")
      names <- all_names[! all_names %in% wrong_names]
      updateSelectizeInput(session,
                           "selectMarker",
                           choices = names,
                           selected = names[1],
                           server = TRUE)
      output$distance_loss = renderPlot({
        plot_distance_loss = subset(linear_image_clean, nn_dist > input$selectDistance )
        title_distance_loss = paste(nrow(plot_distance_loss),
                                 "out of",
                                 nrow(linear_image_clean),
                                 "cells are removed \ndue to projection distance filtering")
        plot(
          linear_image_clean$x,
          linear_image_clean$y,
          cex = input$pointSize/100,
          pch = 16,
          ylab = "CODEX coordinate y",
          xlab = "CODEX coordinate x",
          main = title_distance_loss,
          col.main = "red"
        )
        points(
          plot_distance_loss$x,
          plot_distance_loss$y,
          cex = input$pointSize/100,
          pch = 16,
          col = "red"
        )
      })
      output$zscore_loss = renderPlot({
        plot_zscore_loss = subset(linear_image_clean, zscore > input$selectZscore )
        title_zscore_loss = paste(nrow(plot_zscore_loss),
                                  "out of",
                                  nrow(linear_image_clean),
                                  "cells are removed \ndue to Z-score filtering")
        plot(
          linear_image_clean$x,
          linear_image_clean$y,
          cex = input$pointSize/100,
          pch = 16,
          ylab = "CODEX coordinate y",
          xlab = "CODEX coordinate x",
          main = title_zscore_loss,
          col.main = "orange"
        )
        points(
          plot_zscore_loss$x,
          plot_zscore_loss$y,
          cex = input$pointSize/100,
          pch = 16,
          col = "orange"
        )
      })
      ### Histograms ----
      # output$distance_loss_hist = renderPlot({ 
      #   distance = seq(from = round(min(linear_image_clean$nn_dist)),
      #               to = round(max(linear_image_clean$nn_dist))+1,
      #               length.out = 20)
      #   percent_cells = vector()
      #   for(i in 1:20){
      #     tmp = subset(linear_image_clean, nn_dist < distance[i])
      #     percent_cells[i] = (nrow(tmp)/nrow(linear_image_clean))*100
      #   }
      #   for_plot_hist = as.data.frame(cbind(distance, percent_cells))
      #   plot(x = for_plot_hist$distance, y = for_plot_hist$percent_cells,
      #        xlab = "Projection distance", ylab = "Percent of cells",
      #        main = "Cumulative histogram \nof projection distances")
      # })
      
      output$distance_loss_hist = renderPlot({ 
        plot_distance_loss = subset(linear_image_clean, nn_dist > input$selectDistance )
        
        plot(x = linear_image_clean$shortest_path_order,
             y = linear_image_clean$nn_dist,
             cex = input$pointSize/100,
             pch = 16,
             xlab = "Length", 
             ylab = "Distance",
             main = "Mis-projected points")
        points(
          x = plot_distance_loss$shortest_path_order,
          y = plot_distance_loss$nn_dist,
          cex = input$pointSize/100,
          pch = 16,
          col = "red"
        )
      })
      
      output$zscore_loss_hist = renderPlot({ 
        z_score = seq(from = round(min(linear_image_clean$zscore, na.rm = T)),
                      to = round(max(linear_image_clean$zscore, na.rm = T))+1,
                      length.out = 20)
        percent_cells = vector()
        for (i in 1:20) {
          tmp = subset(linear_image_clean, zscore < z_score[i])
          percent_cells[i] = (nrow(tmp)/nrow(linear_image_clean))*100
        }
        for_plot_hist = as.data.frame(cbind(z_score, percent_cells))
        plot(x = for_plot_hist$z_score, y = for_plot_hist$percent_cells,
             xlab = "Z-score", ylab = "Percent of cells",
             main = "Cumulative histogram \nof Z-score")
      })
      ### plot linear structure ----
      output$line_plot <- renderPlot({
        ### 
        plot_linear <- subset(linear_image_clean, nn_dist <= input$selectDistance &
                                zscore <= input$selectZscore)
        plot(
          plot_linear$shortest_path_length,
          plot_linear$nn_dist,
          cex = input$pointSize/100,
          pch = 16,
          col = "black",
          ylab = "Thickness (Serosa-luminal axis)",
          xlab = "Length (Proximal/outer-distal/inner axis)",
          main = "Linear structure of input image"
        )
      })
      ## Overlay markers ----
      output$marker_plot <- renderPlot({
        ### 
        plot_linear <- subset(linear_image_clean, nn_dist <= input$selectDistance &
                                zscore <= input$selectZscore )
        marker_overlay = merge(plot_linear, input_image, by = "pos",
                               no.dups = T, suffixes = c("","1"))
        
        selected_marker = input$selectMarker
        marker_overlay = marker_overlay[,c("x","y","z","pos","shortest_path_length","nn_dist",selected_marker)]
        marker_overlay[, selected_marker] = scale_marker(marker_overlay[, selected_marker])
        ## overlay large values on top
        marker_overlay = marker_overlay[order(marker_overlay[, selected_marker],decreasing = FALSE),]
        min_scale = min(marker_overlay[, selected_marker])
        median_scale = median(marker_overlay[, selected_marker])
        max_scale = max(marker_overlay[, selected_marker])
        ggplot(NULL) +
          geom_point(data = marker_overlay, aes(
            x = shortest_path_length,
            y = nn_dist,
            color = eval(parse(text = selected_marker))),
            size = input$pointSize/100
          ) +
          labs(color = selected_marker) +
          
          scale_colour_gradientn(colors = c("#e6e6e6",
                                            "#e60000",
                                            "#e60000",
                                            "#e60000", 
                                            "#e60000"),
                                 breaks = c(0,
                                            0.5,
                                            1)
          ) +
          ylab("Thickness (Serosa-luminal axis)") + 
          xlab("Length (Proximal/outer-distal/inner axis)") +
          labs(title = paste0("Digitally unrolled for ", selected_marker," marker")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(
            axis.line = element_line(color = 'black'),
            panel.background = element_rect(fill = 'white', color = 'black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank()
          )
      })
      output$downloadStretch <- downloadHandler(
        filename = function() {
          paste(
            "IntestLine-",
            "linear_structure",
            "_",
            #chosen_distance,
            #"_degree_distance-",
            format(Sys.time(), "%d-%b-%Y_%Hh%Mmin"),
            ".csv",
            sep = ""
          )
        },
        content = function(file) {
          write.table(x = download_data,
                      file,
                      sep = ",",
                      row.names = FALSE)
        }
      )},
      ignoreInit = TRUE,
      priority = -1)
  })
})

shinyApp(ui, server)
