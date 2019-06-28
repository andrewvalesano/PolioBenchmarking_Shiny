
### Project: Polio Mixing Experiment for Bangladesh Sample Sequencing
### Purpose: Empirically derive the optimal filtering criteria for variant calling

library(tidyverse)
library(shiny)
library(ggplot2)
library(pROC)
library(wesanderson)
library(gtable)
library(gridExtra)
library(plyr)
library(reshape2)

source("./supportFunctions.R")

ui <- fluidPage(
  titlePanel("Benchmarking Experiment for Poliovirus Variant Calling with Nextera DNA Flex"),
  
  fluidRow(
    column(6, h4("Quality Distribution"), plotOutput("QualDist")),
    column(6, h4("Distribution on Read"), plotOutput("ReadDist")),
    column(6, h4("Frequency Distribution"), plotOutput("FreqDist")),
    column(6, h4("ROC Curve"), plotOutput("ROC")),
    column(6, h4("Variant Frequency: Expected vs. Observed"), plotOutput("ObservedExpected")),
    column(6, h4("Output Table"), tableOutput('table'))
  ),
  
  hr(),
  
  fluidRow(
    column(3,
           sliderInput(inputId = "MapQ",
                       label = "Mean MapQ Cutoff",
                       min = 30, max = 44, value = 30),
           sliderInput(inputId = "Phred",
                       label = "Mean Phred Cutoff",
                       min = 35, max = 39, value = 30),
           
           sliderInput(inputId = "freq.var",
                        label = "Frequency Cutoff",
                       min = 0, max = 0.1, value = 0),
           sliderInput(inputId = "p.val",
                        label = "p-value Cutoff",
                       min = 0, max = 0.1, value = 0.01)
    ),
    column(3,
           radioButtons("dups",
                        label = "Sequencing Duplicates",
                        choices = list("Variants found in both replicates" = "collapsed", 
                                       "Using only the first replicate" = "first",
                                       "Using only the second replicate" = "second"),
                        selected = "collapsed"),
           radioButtons("disp",
                        label = "deepSNV Dispersion Model",
                        choices = list("Binomial" = "binomial", 
                                       "Betabinominal one-sided" = "onesided"),
                                       #"Betabinominal two-sided" = "twosided"),
                        selected = "onesided")
           
    ),
    
    column(3,
           radioButtons("inputLevel",
                        label = "Genome Copy Input",
                        choices = list("10^5" = "E5", "10^4" = "E4", "10^3" = "E3", "10^2" = "E2", "10^1" = "E1"),
                        selected = "E5")
    ),
    column(3,
           sliderInput("pos",
                       label="Read Position Cutoff",
                       min = 0, max = 250, value = c(0, 250))
    )
  )           
)

server <- function(input, output)
{
  dataInput <- reactive({
    
    replicate = ifelse(input$dups == "first" | input$dups == "second", "notcollapsed", "collapsed")
    file <- paste0("./data/raw/shiny.", input$disp, ".", replicate, ".", "variants.csv")
    variants_raw <- read_csv(file)
    
    if(input$dups == "first") {
      data <- filter(variants_raw, Rep == 1)
    } else if(input$dups == "second") {
      data <- filter(variants_raw, Rep == 2)
    } else {
      data <- variants_raw
    }
    
    data <- filter(data, InputLevel == input$inputLevel)
    
    primer_fwd_outer <- 95 # 95 to 115
    primer_rev_outer <- 7440 # 7415 to 7440
    possible_vars <- (primer_rev_outer - primer_fwd_outer - 1)*3 # Positions in primer range, times 3.
    data <- mutate(data, exp.freq = PercWT/100)
    data <- mutate(data, Id = as.factor(as.character(PercWT)))
    data <- filter(data, pos > primer_fwd_outer & pos < primer_rev_outer)
    
    range_factor <- 5 # this filter prevents there from being observed true positives that are due to lack of Sabin1 amplification in low copy number samples
    data <- mutate(data, category = ifelse(freq.var > exp.freq*range_factor | freq.var < exp.freq*(1/range_factor), FALSE, category))
    
    return(data)
  })
  
  # Plot the ROC curves.
  output$ROC <- renderPlot({
    
    data <- dataInput()
    
    primer_fwd_outer <- 95 # 95 to 115
    primer_rev_outer <- 7440 # 7415 to 7440
    expectedVariants <- read_csv("./data/reference/MixingStudyExpectedVariants.csv")
    expectedTruePos <- nrow(expectedVariants) # Expect to see 66 variants.
    possible_vars <- (primer_rev_outer - primer_fwd_outer - 1)*3 # Positions in primer range, times 3.

    filtered_data <- filter(data, p.val < input$p.val, MapQ > input$MapQ & Phred > input$Phred & freq.var > input$freq.var & Read_pos <= input$pos[2] & Read_pos >= input$pos[1])

    palette <- wes_palette("Darjeeling1")
    roc.df <- miseq.roc(filtered_data, expectedTruePos, possible_vars, ">")
    roc.plot <- ggplot(roc.df, aes(x = 1-adj.specificity, y = adj.sensitivity)) + geom_step(aes(color = samp), size = 0.7) + xlab(bquote("1-Specificity (x"~10^-3~")")) + ylab("Sensitivity") + scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,0.005),breaks=c(0,0.001,0.002,0.003,0.004,0.005),labels=c(0:5)) + theme_minimal() + theme(legend.position="right") + scale_color_manual(values = palette, name="Frequency (%)")
    roc.plot
  })
  
  # Plot MapQ and Phred, color-coded by TP/FP
  output$QualDist <- renderPlot({
    data <- dataInput()
    palette <- wes_palette("Darjeeling1")
    ggplot(data, aes(x = Phred, y = MapQ, color = category)) + geom_point() + xlab("Phred") + ylab("MapQ") + theme_minimal() + scale_color_manual(values = palette[c(4,5)]) + ylim(c(30,44)) + xlim(c(35,39)) + geom_vline(xintercept = input$Phred, linetype = "dotted", color = "black", size = 1) + geom_hline(yintercept = input$MapQ, linetype = "dotted", color = "black", size = 1)   
  })
  
  # Plot average read position, color-coded by TP/FP
  output$ReadDist <- renderPlot({
    data <- dataInput()
    palette <- wes_palette("Darjeeling1")
    ggplot(data, aes(x = Read_pos, fill = category)) + geom_histogram(position = "dodge") + xlab("Read Position") + ylab("Count") + scale_fill_manual(values = palette[c(4,5)]) + theme_minimal() + geom_vline(xintercept = input$pos, linetype = "dotted", color = "black", size = 1)
  })
  
  # Plot frequency histogram, color-coded by TP/FP
  output$FreqDist <- renderPlot({
    data <- dataInput()
    palette <- wes_palette("Darjeeling1")
    ggplot(data, aes(x = freq.var, fill = category)) + geom_histogram(binwidth = 0.001, position = "dodge") + xlab("Frequency") + ylab("Count") + scale_fill_manual(values = palette[c(4,5)]) + theme_minimal() + geom_vline(xintercept = input$freq.var, linetype = "dotted", color = "black", size = 1) + xlim(c(0, 0.1))
  })
  
  # Plot observed vs. expected frequency
  output$ObservedExpected <- renderPlot({
    data <- dataInput()
    data <- mutate(data, exp.freq = PercWT/100)
    palette <- wes_palette("Darjeeling1")
    
    ggplot(data, aes(x = exp.freq, y = freq.var, color = category)) + 
      geom_point(size = 1) +
      scale_color_manual(values = palette[c(4,5)]) +
      theme_minimal() + 
      geom_abline(intercept = 0, slope = 1,linetype = 2, size = 1) + 
      xlab("Expected Frequency") + 
      ylab("Observed Frequency") + 
      scale_y_log10(limits=c(0.001,0.1),breaks=c(0.001,0.002,0.005,0.01,0.02,0.05,0.1)) +
      scale_x_log10(limits=c(0.001,0.1),breaks=c(0.001,0.002,0.005,0.01,0.02,0.05,0.1))
  })
  
  # Make the table
  output$table <- renderTable({
    data <- dataInput()
    
    primer_fwd_outer <- 95 # 95 to 115
    primer_rev_outer <- 7440 # 7415 to 7440
    expectedVariants <- read_csv("./data/reference/MixingStudyExpectedVariants.csv")
    expectedTruePos <- nrow(expectedVariants) # Expect to see 66 variants.
    possible_vars <- (primer_rev_outer - primer_fwd_outer - 1)*3 # Positions in primer range, times 3.
    
    filtered_data <- filter(data, p.val < input$p.val, MapQ > input$MapQ & Phred > input$Phred & freq.var > input$freq.var & Read_pos <= input$pos[2] & Read_pos >= input$pos[1])
    
    dd = 4
    m.roc.table <- miseq.roc.table(filtered_data, 1, expectedTruePos, possible_vars, ">")
    m.roc.table <- rename(m.roc.table, c("exp.freq"="Frequency","adj.sensitivity"="Sensitivity","TP"="True\nPositives","adj.specificity"="Specificity","FP"="False\nPositives"))
    m.roc.table$Frequency <- c("100%", "10%", "5%", "2%", "1%")
    m.roc.table$Sensitivity <- round(m.roc.table$Sensitivity, digits = dd)
    m.roc.table$Specificity <- round(m.roc.table$Specificity, digits = dd)
    m.roc.table
    
  }, digits = 4)
  
}

shinyApp(ui = ui, server = server)

