
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
    column(6, h4("Coverage"), plotOutput("Coverage")),
    column(6, h4("Output Table"), tableOutput('table')),
    column(6, h4("Expected Variants by Position"), plotOutput("Variants"))
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
                       min = 0, max = 0.1, value = 0)
    ),
    column(3,
           sliderInput("pos",
                       label="Read Position Cutoff",
                       min = 0, max = 250, value = c(0, 250)),
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
                        selected = "onesided"),
           radioButtons("tna",
                        label = "Spiked into stool total nucleic acid",
                        choices = list("Spike-in" = "yes", 
                                       "No spike-in" = "no"),
                        selected = "yes")
    ),
    
    column(3,
           radioButtons("inputLevel",
                        label = "Genome Copy Input",
                        choices = list("4x10^4" = "4_E4", "9x10^3" = "9_E3", "9x10^2" = "9_E2", "9x10^1" = "9_E1"),
                        selected = "4_E4"),
           radioButtons("subset",
                        label = "Subset of Reads",
                        choices = list("100%" = "100", 
                                       "50%" = "50",
                                       "25%" = "25",
                                       "10%" = "10"),
                        selected = "100")
    )
  )           
)

server <- function(input, output)
{
  dataInput <- reactive({
    
    replicate = ifelse(input$dups == "first" | input$dups == "second", "notcollapsed", "collapsed")
    file <- paste0("./data/processed/shiny2.", input$disp, ".", input$subset, ".", replicate, ".variants.csv")
    variants_raw <- read_csv(file)
    
    if(input$dups == "first") {
      data <- filter(variants_raw, Rep == 1)
    } else if(input$dups == "second") {
      data <- filter(variants_raw, Rep == 2)
    } else {
      data <- variants_raw
    }
    
    data <- filter(data, InputLevel == input$inputLevel, inTNA == input$tna)
    data <- mutate(data, SampleNumber = Id)
    
    primer_fwd_outer <- 95 # 95 to 115
    primer_rev_outer <- 7440 # 7415 to 7440
    possible_vars <- (primer_rev_outer - primer_fwd_outer - 1)*3 # Positions in primer range, times 3.
    data <- mutate(data, exp.freq = PercWT/100)
    data <- mutate(data, Id = as.factor(as.character(PercWT)))
    data <- filter(data, pos > primer_fwd_outer & pos < primer_rev_outer)
    
    range_factor <- 5 # this filter prevents there from being observed true positives that are due to lack of Sabin1 amplification in low copy number samples
    #data <- mutate(data, category = ifelse(freq.var > exp.freq*range_factor | freq.var < exp.freq*(1/range_factor), FALSE, category))
    
    return(data)
  })
  
  covInput <- reactive({
    
    file <- paste0("./data/processed/shiny2.", input$disp, ".", input$subset, ".coverage.csv")
    coverage_raw <- read.csv(file)
    
    return(coverage_raw)
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
      scale_y_log10(limits=c(0.001,0.5),breaks=c(0.001,0.002,0.005,0.01,0.02,0.05,0.1, 0.25, 0.5)) +
      scale_x_log10(limits=c(0.001,0.5),breaks=c(0.001,0.002,0.005,0.01,0.02,0.05,0.1))
  })
  
  # Plot coverage for the given samples
  output$Coverage <- renderPlot({
    data <- dataInput()
    coverage_raw <- covInput()
    coverage <- filter(coverage_raw, Id %in% data$SampleNumber)
    
    base = 12000
    interval = 300
    linesize = 1.2
    ggplot(coverage, aes(x = chr.pos, y = coverage, group = Id)) + geom_line(aes(color = Id)) + 
      theme_bw() + 
      xlab("Genome Position") + 
      ylab("Coverage (Read Depth)") + 
      theme(legend.position = "none") +
      xlim(c(95, 7440)) +
      geom_segment(data = coverage, size = linesize, color = "orange1", aes(x = 98, y = base + interval*0, xend = 2436, yend = base + interval*0)) +
      geom_segment(data = coverage, size = linesize, color = "mediumseagreen", aes(x = 1472, y = base + interval*1, xend = 4151, yend = base + interval*1)) +
      geom_segment(data = coverage, size = linesize, color = "indianred4", aes(x = 3214, y = base + interval*2, xend = 5861, yend = base + interval*2)) +
      geom_segment(data = coverage, size = linesize, color = "slateblue3", aes(x = 4966, y = base + interval*3, xend = 7400, yend = base + interval*3))
  })
  
  # Plot position of true positives and false negatives
  output$Variants <- renderPlot({
    data <- dataInput()
    data <- filter(data, p.val < input$p.val, MapQ > input$MapQ & Phred > input$Phred & freq.var > input$freq.var & Read_pos <= input$pos[2] & Read_pos >= input$pos[1])
    data <- mutate(data, level = ifelse(PercWT == 1, "1_percent", 
                                        ifelse(PercWT == 2, "2_percent", 
                                               ifelse(PercWT == 5, "5_percent", 
                                                      ifelse(PercWT == 10, "10_percent", 
                                                             ifelse(PercWT == 100, "100_percent", NA))))))
    data <- mutate(data, mutation_level = paste0(mutation, "_", level))
    
    expectedVariants <- read_csv("./data/reference/MixingStudyExpectedVariants_ByLevel.csv")
    expectedVariants <- mutate(expectedVariants, mutation_level = paste0(mutation, "_", level))
    expectedVariants <- mutate(expectedVariants, found = ifelse(mutation_level %in% data$mutation_level, "True Positive", "False Negative"))
    
    levels <- c("100_percent", "10_percent", "5_percent", "2_percent", "1_percent")
    chrs <- data.frame("level" = levels)
    chrs <- mutate(chrs, start = 0, stop = 7440)
    palette <- wes_palette("Darjeeling1")
    expectedVariants$level <- factor(expectedVariants$level, levels = rev(c("100_percent","10_percent","5_percent","2_percent","1_percent")))
    chrs$level <- factor(chrs$level, levels = levels(expectedVariants$level))
    
    ggplot(expectedVariants, aes(x = Position, y = level)) +
      geom_point(aes(color = found), size = 5, shape = 108) +
      geom_segment(data = chrs, aes(x = start, y = level, xend = stop, yend = level)) +
      ylab("Expected Frequency Group") +
      xlab("") +
      theme_minimal() +
      scale_color_manual(name = "", values = palette[c(1,2)]) +
      theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
      scale_x_continuous(breaks = c()) +
      theme(panel.grid.major = element_line(colour = "gray96"), panel.grid.minor = element_line(colour = "white")) + theme(legend.position = "bottom") +
      scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 7500))
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

