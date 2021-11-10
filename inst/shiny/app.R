#This is a shiny app to explore the data from the 2021 CLL Combinatorial Screen. 
#There are 5 tabs; one home page and  4 with plots to explore the screening data and associated meta data. 
#Shiny app written by Holly Giles. 

#Set up ----------------------
##Load Libraries--------------
library(shiny)
library(shinydashboard)
library(plyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(magrittr)
library(RColorBrewer)
library(gtable)
library(grid)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)
library(survival)
library(survminer)

##Define functions ---------------------------------
HGmakelegends <- function (legendFor, colors) 
{
  x = NULL
  y = NULL
  colors = colors[names(colors) %in% legendFor]
  nleg = length(colors)
  
  #edit these widths to change alignment with lasso plots in patchwork code
  
  if ("M" %in% legendFor) {
    wdths = c(0.01, 1.5, 1.5, 1.5, 0.01)
    hghts = c(2)
  } else {
    wdths = c(0.01, 2, 2, 0.01)
    hghts = c(2)
  }
  

  
  
  gtl = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
  n = 2
  if ("M" %in% names(colors)) {
    Mgg = ggplot(data = data.frame(x = 1, 
                                   y = factor(c("LP", "IP", "HP"), 
                                              levels = c("LP", "IP", "HP"))), 
                 aes(x = x, y = y, fill = y)) + 
      geom_tile() + 
      scale_fill_manual(name = "Methylation cluster", 
                        values = setNames(colors[["M"]], 
                                          nm = c("LP", "IP","HP"))) + 
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12))
    
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Mgg), "guide-box"), 1, n)
    n = n + 1
  }
  
  if ("I" %in% names(colors)) {
    Igg = ggplot(data = data.frame(x = 1, y = factor(c("Unmutated", 
                                                       "Mutated"), levels = c("Unmutated", "Mutated"))), 
                 aes(x = x, y = y, fill = y)) + geom_tile() + scale_fill_manual(name = "IGHV", 
                                                                                values = setNames(colors[["I"]], nm = c("Unmutated", "Mutated"))) + 
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Igg), 
                                             "guide-box"), 1, n)
    n = n + 1
  }
  
  if ("G" %in% names(colors)) {
    Ggg = ggplot(data = data.frame(x = 1, y = factor(c("Wild Type", 
                                                       "Mutated"), levels = c("Wild Type", "Mutated"))), 
                 aes(x = x, y = y, fill = y)) + geom_tile() + scale_fill_manual(name = "Gene", 
                                                                                values = setNames(colors[["G"]], nm = c("Wild Type", "Mutated"))) + 
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12))
    gtl = gtable_add_grob(gtl, gtable_filter(ggplotGrob(Ggg), 
                                             "guide-box"), 1, n)
    n = n + 1
  }
  
  return(list(plot = gtl, width = sum(wdths), height = sum(hghts)))
}

## Load data -----------------------------------------
load("./shinyAppData.Rdata") 

thedrugs <- unique(df_patmeta$Drug) %>% setdiff("DMSO")
thecytokines <- unique(df_patmeta$Cytokine) %>% setdiff("No Cytokine")

geneticData <- patMeta[c(5, 7:ncol(patMeta))]
geneticData <- geneticData[, colSums(geneticData == 1|geneticData=="M", na.rm = TRUE) > 3]  
geneticData <- geneticData[,names(sort(colSums(geneticData == 1|geneticData=="M", na.rm = TRUE), decreasing = TRUE))]
thegenes <-   colnames(geneticData)

##Set themes ----------------------------------------
t1 <-
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_line(),
    panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.line = element_line(size=.4),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.text.x  = element_text(angle=90, size=16, 
                                face="bold", hjust = 1, vjust = 0.4),
    axis.text.y = element_text(size = 18),
    axis.ticks.x = element_line(linetype = "dotted"),
    axis.ticks.length = unit(0.3,"cm"),
    axis.title.x = element_text(face="bold", size=18), 
    axis.title.y = element_text(face="bold", size=18),
    plot.title = element_text(face="bold", size=18, hjust = 0.5)
)

t2 <- 
  t1+
  theme( axis.text.x  = element_text(angle=0, size=16, face="bold", hjust = 0.5, vjust = 1))


colors <- c("#A1BE1F", #green
            "#F4C61F", #yellow
            "#734595", #purple
            "#D41645", #red
            "#3B6FB6", #blue
            "#F49E17", #orange
            "#E2E868", #light green
            "#FDD757", #light yellow
            "#CBA3D8", #light purple
            "#E58F9E", #light purple
            "#8BB8E8", #light blue
            "#EFC06E",
            "#303030") #light orange

offwhite <- "#f8f8ff"
#for negatives only: 
palblues <- c("#003DA5", "#2055B0", "#406EBC", "#6086C7", "#809ED2", "#9FB6DD", "#BFCFE9", "#DFE7F4")

#for positives only:
palreds <- c("#F4E0E7", "#E9C2CF", "#DEA3B6", "#D3849E", "#C76586", "#BC476E", "#B12855", "#A6093D")

# Dashboard Page Setup ----------------------------------------------------
header <-dashboardHeader(title = "Drug - Cytokine Combinatorial Screen", 
                         titleWidth  = "400")
        
        # Dashboard Sidebar -------------------------------------------------------
sidebar <- dashboardSidebar(width = "400px",
            sidebarMenu(
              #TAB 1 : Home Page
              menuItem("Home Page", tabName = "homePage", icon = icon("home")
              ), 
              
              #TAB 2 : Viability plots
              menuItem("Drug and stimulus responses", tabName = "viabPlots", icon = icon("chart-line")
              ),  
              
              #TAB 3 : Drug and Cytokine Beeswarms
              menuItem("Effects of mutations on drug & stimulus responses", tabName = "beeswarms", icon = icon("dna")
                       ),
              
              #TAB 4 : Genetic predictors of response
              menuItem("Genetic predictors of drug & stimulus responses", tabName = "lassoPlots", icon = icon("project-diagram")
                       ), 
              #TAB 5 : Drug and Cytokine Beeswarms
              menuItem("Genetic predictors of drug & stimulus interactions", tabName = "lassoIntPlots", icon = icon("signal")
                       )
            )
)


#Dashboard Body------------------------------------      
body <- dashboardBody(
  
    tabItems(
      #TAB 1 --------------------------------
      tabItem(tabName = "homePage",
              fluidRow(
                box(#header box
                width = 12,
                status = "primary",
                title = "Welcome to the Shiny App to Explore the Drug-Cytokine Combinatorial Screen Dataset",
                helpText("Please click on the tabs on the left to interactively explore the data. 
                This app accompanies the 2021 paper: Combinatorial drug-microenvironment interaction mapping reveals cell-extrinsic drug resistance mechanisms and clinically relevant patient subgroups in CLL. If you use this app to support published research, please cite the paper. 
                The data can also be explored manually using the online vignette."))), 
               
              fluidRow(#boxes for tabs 2 and 3
                box(title = "Drug and stimulus responses",
                    helpText("Explore log-transformed viabilites for single and combinatorial treatments and annnotate samples by genetic features"), 
                    width = 6),
                
                box(title = "Effects of mutations on drug & stimulus responses",
                    helpText("View log-transformed viability data stratified by genetic features"),
                  width = 6)),
                
                fluidRow(#boxes for tabs 4 and 5
                  box(title = "Genetic predictors of drug & stimulus responses", 
                      helpText("View predictor profiles for viabilities after treatment with drugs and stimuli"),
                
                      width = 6),
                  box(title = "Genetic predictors of drug & stimulus interactions",
                      helpText("View predictor profiles for the size of interactions between drugs and stimuli"),
                      width = 6)),
              
              fluidRow(#boxes for info
                box(title = "Maintenance",
                    helpText("This shiny app is maintained by Holly Giles. If you have questions, please email me at holly.giles@embl.de."),
                    
                    width = 6),
                box(title = "Security and Licence",
                    helpText("Copyright 2021"),
                    width = 6)
                )),
              
              
      #TAB 2 --------------------------------
        tabItem(tabName = "viabPlots",
                fluidRow(
                    infoBoxOutput("drugTarget2", width = 6),
                    infoBoxOutput("cytTarget2", width = 6)),
                fluidRow(
                  box(
                    width = 12,
                    helpText("Select a drug : stimulus combination to view viabilities for all patient samples, coloured by chosen genetic feature."), 
                    #Drug
                    selectInput("drugP2", "Drug:", 
                                choices = thedrugs),
                    #Cytokine
                    selectInput("cytP2", "Stimulus:",
                                choices = thecytokines),
                    #Omic 1
                    selectInput("omicP2", "Genomic Feature:",
                                choices = thegenes))),
                fluidRow(
                  box(plotOutput("linePlot"), title = "Viabilities with single and combinatorial treatments:", status = "primary",
                        width = 12)), 
                fluidRow(
                    box(plotOutput("linePlot_omic1"), 
                        title = "For mutated samples only:", 
                        status = "primary", 
                        width = 6), 
                    box(plotOutput("linePlot_omic2"),  
                        title = "For unmutated samples only:", 
                        status = "primary", 
                        width = 6)
                    )
                ), 
        
        #TAB 3 -----------------------------------------
        tabItem(tabName = "beeswarms",
                fluidRow(
                  infoBoxOutput("drugTarget1", width = 6),
                  infoBoxOutput("cytTarget1", width = 6)),
                
                fluidRow(
                  box(plotOutput("beeswarm", width = "100%"), title = "Beeswarm plot of log transformed viability values, stratified by chosen feature(s)", status = "primary", width = 9),
                  box("Select inputs", 
                      helpText("Select the drug and/or stimulus to visualise log transformed viability values stratified by the selected genetic feature. Optional: Select a second genetic feature to visualise both simultaneously."),
                      width = 3, 
                      #Drug
                      selectInput("drugP1", "Drug:", 
                                  c(thedrugs, "DMSO"), selected = "DMSO"),
                      #Cytokine
                      selectInput("cytP1", "Stimulus:",
                                  c(thecytokines, "No Cytokine"), selected = "Resiquimod"),
                      #Omic 1
                      selectInput("omic1P1", "Genetic Feature 1:",
                                  thegenes, selected = "IGHV.status"),
                      #Omic 2
                      selectInput("omic2P1", "Genetic Feature 2 (optional):",
                                  c("NA",thegenes), 
                                  selected = "NA")
                      
                  )
                )
                
                
        ), 
        #TAB 4 --------------------------
        tabItem(tabName = "lassoPlots",
                fluidRow(
                  infoBoxOutput("drugTarget3", width = 6),
                  infoBoxOutput("cytTarget3", width = 6)
                ),
                fluidRow(
                  
                  box(width = 6,
                      height = 250, #Input: Select drug and cytokine combination 
                      
                      helpText("Select drug and/or stimulus to view genetic predictors of viability with selected treatment."),
                      #Drug
                      selectInput("drugP3", "Drug:", 
                                  c(thedrugs, "DMSO"),
                                  selected = "DMSO" ),
                      #Cytokine
                      selectInput("cytP3", "Stimulus:",
                                  c(thecytokines, "No Cytokine"), 
                                  selected = "Resiquimod")),
                  #Omic 1
                  box(width = 6, 
                      height = 250,
                      helpText("Adjust Frequency Threshold (the proportion of cases in which a predictor must be significant with each bootstapped repeat of the regression) and Coefficient Threshold (minimum value of predictors)."),
                      numericInput("freqCutP3", "Frequency Threshold:", 0.75,  0.1, 1.0, 0.1),
                      numericInput("coefCutP3", "Coefficient Threshold:", 0,0, 0.3, 0.01))
                ),
                
                fluidRow(
                    box(plotOutput("lassoPlot", height = "600px"), 
                        title = "Lasso Regularised Regression to identify genetic features that predict viability with selected treatment",
                        helpText("The predictor profile plot depicts genetic features that predict response to selected treatment. 
                        Bar plot on left indicates size and sign of coefficients for the named predictors. 
                        Positive coefficients indicate higher viability after stimulation, if the feature is present.
                        Scatter plot and heatmap indicate how each genetic feature relates to 
                        patient sample viabilities: Scatter plot indicates log(viability) values, in order of magnitude, for each individual sample.
                        Heatmap shows patient mutation status for each of genetic predictors for corresponding sample in scatter plot."),
                        width = 12))
                
                    
              
                    
                ),
                 
        #TAB 5 ------------------------
        tabItem(tabName = "lassoIntPlots", 
                fluidRow(
                        infoBoxOutput("drugTarget4", width = 6),
                        infoBoxOutput("cytTarget4", width = 6)
                    ),
                fluidRow(
                  
                  box(width = 6,#Input: Select drug and cytokine combination 
                      
                      helpText("Select drug and stimulus to view genetic predictors for the size of interaction (beta int) between chosen combination."),
                      #Drug
                      selectInput("drugP4", "Drug:", 
                                  thedrugs, selected = "Ibrutinib"),
                      #Cytokine
                      selectInput("cytP4", "Stimulus:",
                                  thecytokines, selected = "IL-4")),
                  #Omic 1
                  box(width = 6, 
                      helpText("Adjust Frequency Threshold (the proportion of cases in which a predictor must be significant with each bootstrapped repeat of the regression) and Coefficient Threshold (minimum value of coefficients)."),
                      numericInput("freqCutP4", "Frequency Threshold:", 0.9, 0.1, 1.0, 0.1),
                      
                      numericInput("coefCutP4", "Coefficient Threshold:", 0.0,0.0, 0.5, 0.01))
                ),
                
                fluidRow(
                  box(plotOutput("lassoIntPlot", height = "600px"), 
                      title = "Lasso Regularised Regression to identify genetic features that predict the size of drug - stimulus interactions.",
                      helpText("Predictor profile plot depicts genetic features that modulate the size of interaction between chosen drug and stimulus.
                      To generate predictor profile, a linear model was fitted in a sample - specific manner, to calculate drug- stimulus interaction coefficients (beta int) for each patient sample. 
                      Ranked patient-specific beta int values are shown in lower scatter plot.  
                      Associations between the size of beta int and genetic features were identified using multivariate regression with L1 (lasso) regularisation, 
                      with gene mutations and IGHV status as predictors.  
                      The horizontal bars  on left show the size of fitted coefficients assigned to genetic features, showing those that meet selected cut offs.
                      Matrix above scatter plot indicates patient mutation status for the selected genetic features. 
                      Matrix fields correspond to points in scatter plot (ie patient data is aligned), to indicate how the size 
                      of beta int varies with selected genetic feature."),
                      width = 12))
                
                )
   
      

                
        )
)
        
                    
       
        
  
 ui <- 
   dashboardPage(
     header,
     sidebar,
     body
 )
        

                         

# Define server logic 
server <- function(input, output) { 
  
  #TAB2 OUTPUTS -------------------------------------------------------
  output$drugTarget2 <- renderInfoBox({
    infoBox(
      "Drug Target", 
      if(input$drugP2 == "DMSO"){"Control Treatment"}else{filter(drugs, Name == input$drugP2) %>% select(main_targets)}, 
      icon = icon("bullseye"),
      color = "purple"
    )
  })
  
  output$cytTarget2 <- renderInfoBox({
    infoBox(
      "Stimulus Pathway", if(input$cytP2 == "Not Cytokine"){"No Cytokine"}else{filter(Cytokines, name == input$cytP2) %>% select(pathway)}, 
      icon = icon("dot-circle"),
      color = "yellow"
    )
  })
  
  output$linePlot <- renderPlot({
    
    df.comp1 <- df_patmeta[which(!is.na(df_patmeta[input$omicP2])),]
    
    plotTab2 <- 
      df.comp1 %>%
      dplyr::filter(Drug %in% c("DMSO", input$drugP2)) %>%
      dplyr::filter(Cytokine %in% c("No Cytokine", input$cytP2)) %>%
      dplyr::mutate(treatment_combination =
                      ifelse(Drug=="DMSO" & Cytokine=="No Cytokine","DMSO",
                             ifelse(Drug == input$drugP2 & Cytokine == "No Cytokine",paste0("Drug=", input$drugP2),
                                    ifelse(Drug == "DMSO" & Cytokine == input$cytP2, paste0("Cytokine=", input$cytP2),
                                           paste0("Drug=", input$drugP2, "\nCytokine=", input$cytP2) )))) %>%
      mutate(treatment_combination = factor(treatment_combination, levels=c("DMSO",
                                                                            paste0("Drug=", input$drugP2),
                                                                            paste0("Cytokine=", input$cytP2),
                                                                            paste0("Drug=", input$drugP2, "\nCytokine=", input$cytP2))))
    
    #axis scale
    ymin <- min(plotTab2$Log)-0.1
    
    ymax <- max(plotTab2$Log + 0.1)
    
    #plot
    
    ggplot(plotTab2, aes(treatment_combination, Log, group = PatientID)) +
      
      #plot WT points and then Mutated points
      geom_point(data = plotTab2[which(plotTab2[input$omicP2]==0|plotTab2[input$omicP2]=="M"),], aes_string(colour = input$omicP2), size = 1) +
      geom_line(data = plotTab2[which(plotTab2[input$omicP2]==0|plotTab2[input$omicP2]=="M"),], aes_string(colour = input$omicP2), size = 0.1) +
      geom_point(data = plotTab2[which(plotTab2[input$omicP2]==1|plotTab2[input$omicP2]=="U"),], aes_string(colour = input$omicP2), size = 1) +
      geom_line(data = plotTab2[which(plotTab2[input$omicP2]==1|plotTab2[input$omicP2]=="U"),], aes_string(colour = input$omicP2), size = 0.1) +
      geom_hline(yintercept = 0)+
      
      
      #Sort aesthetics
      ggtitle(paste(input$drugP2, "+", input$cytP2, sep=" ")) + 
      scale_color_manual(values = c("#D0D0CE", palreds[7]), labels =  c("M" = "M", "U" = "U", "1"= input$omicP2, "0"= "WT")) +
      t2 + 
      scale_x_discrete(labels= c( "DMSO", input$drugP2, input$cytP2, paste(input$drugP2, "+", input$cytP2, sep="\n"))) +
      ylab("Log(Viability)") + 
      xlab("")+
      ylim(ymin, ymax) +
      theme(legend.title = element_text(face='bold',hjust = 1, size=11),
            legend.position = c(0.95, 0.9),
            legend.key = element_blank(),
            legend.text = element_text(size=12),
            legend.background = element_rect(color = "black"))
    
    
  })
  
  #plot with mutants only
  output$linePlot_omic1 <- renderPlot({
    df.comp2 <- df_patmeta[which(!is.na(df_patmeta[input$omicP2])),]
    
    plotTab3 <- 
      df.comp2 %>%
      dplyr::filter(Drug %in% c("DMSO", input$drugP2)) %>%
      dplyr::filter(Cytokine %in% c("No Cytokine", input$cytP2)) %>%
      dplyr::mutate(treatment_combination =
                      ifelse(Drug=="DMSO" & Cytokine=="No Cytokine","DMSO",
                             ifelse(Drug == input$drugP2 & Cytokine == "No Cytokine",paste0("Drug=", input$drugP2),
                                    ifelse(Drug == "DMSO" & Cytokine == input$cytP2, paste0("Cytokine=", input$cytP2),
                                           paste0("Drug=", input$drugP2, "\nCytokine=", input$cytP2) )))) %>%
      mutate(treatment_combination = factor(treatment_combination, levels=c("DMSO",
                                                                            paste0("Drug=", input$drugP2),
                                                                            paste0("Cytokine=", input$cytP2),
                                                                            paste0("Drug=", input$drugP2, "\nCytokine=", input$cytP2))))
    
    #axis scale
    ymin <- min(plotTab3$Log)-0.1
    
    ymax <- max(plotTab3$Log + 0.1)
    
    
    plotTab3 <- plotTab3[which(plotTab3[input$omicP2]==1|plotTab3[input$omicP2]=="M"),]
    
   
    #plot
    
    ggplot(plotTab3, aes(treatment_combination, Log, group = PatientID)) +
      
      #plot mutated points only
      geom_point(colour = palreds[7]) +
      geom_line(colour = palreds[7]) +
      geom_hline(yintercept = 0)+
      
      
      #Sort aesthetics
      ggtitle(paste(input$drugP2, "+", input$cytP2, "\nfor", ifelse(input$omicP2 == "IGHV.status", "IGHV Mutated", input$omicP2), "samples", sep=" ")) + 
      t1 + theme( axis.text.x  = element_text(angle = 45, size = 16, face = "bold", hjust = 0.5, vjust = 0.5))+
      scale_x_discrete(labels= c( "DMSO", input$drugP2, input$cytP2, paste(input$drugP2, "+", input$cytP2, sep="\n"))) +
      ylab("Log(Viability)") + 
      xlab("")+
      ylim(ymin, ymax) +
      theme(legend.title = element_text(face='bold', hjust = 1, size=11),
            legend.position = c(0.95, 0.9),
            legend.key = element_blank(),
            legend.text = element_text(size=12),
            legend.background = element_rect(color = "black"))
  })
  
  output$linePlot_omic2 <- renderPlot({
    df.comp3 <- df_patmeta[which(!is.na(df_patmeta[input$omicP2])),]
    
    plotTab4 <- 
      df.comp3 %>%
      dplyr::filter(Drug %in% c("DMSO", input$drugP2)) %>%
      dplyr::filter(Cytokine %in% c("No Cytokine", input$cytP2)) %>%
      dplyr::mutate(treatment_combination =
                      ifelse(Drug=="DMSO" & Cytokine=="No Cytokine","DMSO",
                             ifelse(Drug == input$drugP2 & Cytokine == "No Cytokine",paste0("Drug=", input$drugP2),
                                    ifelse(Drug == "DMSO" & Cytokine == input$cytP2, paste0("Cytokine=", input$cytP2),
                                           paste0("Drug=", input$drugP2, "\nCytokine=", input$cytP2) )))) %>%
      mutate(treatment_combination = factor(treatment_combination, levels=c("DMSO",
                                                                            paste0("Drug=", input$drugP2),
                                                                            paste0("Cytokine=", input$cytP2),
                                                                            paste0("Drug=", input$drugP2, "\nCytokine=", input$cytP2))))
    
    #axis scale
    ymin <- min(plotTab4$Log)-0.1
    
    ymax <- max(plotTab4$Log + 0.1)
    
    
    
    
    plotTab4 <- plotTab4[which(plotTab4[input$omicP2]==0|plotTab4[input$omicP2]=="U"),]
    
                                     
    #plot
    
    ggplot(plotTab4, aes(treatment_combination, Log, group = PatientID)) +
      
      #plot WT points and then Mutated points
      geom_point(colour = "#D0D0CE") +
      geom_line(colour = "#D0D0CE") +
      geom_hline(yintercept = 0)+
      
      
      #Sort aesthetics
      ggtitle(paste(input$drugP2, "+", input$cytP2,"\nfor", ifelse(input$omicP2 == "IGHV.status", "IGHV Unmutated", "WT"),"samples",  sep=" ")) + 
      t1 + theme( axis.text.x  = element_text(angle = 45, size = 16, face = "bold", hjust = 0.5, vjust = 0.5))+
      scale_x_discrete(labels = c( "DMSO", input$drugP2, input$cytP2, paste(input$drugP2, "+", input$cytP2, sep="\n"))) +
      ylab("Log(Viability)") + 
      xlab("") +
      ylim(ymin, ymax) +
      theme(legend.title = element_text(face = 'bold',hjust = 1, size = 11),
            legend.position = c(0.95, 0.9),
            legend.key = element_blank(),
            legend.text = element_text(size = 12),
            legend.background = element_rect(color = "black"))
  })
  
    
#TAB3 OUTPUTS-------------------------------------------------------
    
    output$drugTarget1 <- renderInfoBox({
        infoBox(
            "Drug Target", 
            if(input$drugP1 == "DMSO"){"Control Treatment"}else{filter(drugs, Name == input$drugP1) %>% select(main_targets)}, 
            icon = icon("bullseye"),
            color = "purple"
        )
    })
    
    output$cytTarget1 <- renderInfoBox({
        infoBox(
            "Stimulus Pathway", if(input$cytP1 == "Not Cytokine"){"No Cytokine"}else{filter(Cytokines, name == input$cytP1) %>% select(pathway)}, 
            icon = icon("dot-circle"),
            color = "yellow"
        )
    })
    

    
  
    output$beeswarm <- renderPlot({
      validate(
        need(input$omic2P1 != input$omic1P1, "Please ensure selected genetic features are distinct.")
      )
        if(input$omic2P1 == "NA") {
            
          
          plotTab1 <- df_patmeta %>% 
            dplyr::filter(Cytokine==input$cytP1,
                          Drug == input$drugP1) %>%
            drop_na(input$omic1P1)
          
          
          plotTab1 %>%
            ggplot(aes_string(x=input$omic1P1, y="Log", color=input$omic1P1)) +
            geom_boxplot()+
            geom_beeswarm(cex=1.5) +
            
            stat_compare_means(method = "t.test", 
                               label.x.npc= "centre")+ 
            
            #Aesthetics
            guides(color="none") +
            geom_hline(yintercept = 0) +
            scale_color_manual(values = c(colors[1], colors[6])) +
            ggtitle(paste("Treatment:", input$drugP1, "and", input$cytP1,  sep = " ")) +
            scale_x_discrete(labels=c("0"="WT", "1"= input$omic1P1, "M"= "M", "U"= "U")) +
            xlab("")+
            ylab("Log(Viability)") +
            t2
       
    
    
        } else {
        
          
          plotTab1 <- df_patmeta %>% 
            dplyr::filter(Cytokine==input$cytP1,
                          Drug == input$drugP1) %>%
            drop_na(all_of(input$omic1P1)) %>% 
            drop_na(all_of(input$omic2P1))
          
          #define interaction
          inter <- paste0('interaction(', paste0(c(input$omic1P1, input$omic2P1), collapse = ', ' ),')')
          
      
          
          plotTab1 %>%
            ggplot(aes_string(x=inter, y="Log", color=input$omic1P1)) +
            geom_boxplot()+
            geom_beeswarm(cex=1.5) +
            
            
            #Aesthetics
            guides(color="none") +
            geom_hline(yintercept = 0) +
            scale_color_manual(values = c(colors[1], colors[6])) +
            ggtitle(paste("Treatment:", input$drugP1, "and",input$cytP1,  sep = " ")) +
            scale_x_discrete(labels=c("M.0"="IGHV M\n WT",
                                      "U.0"="IGHV U\n WT",
                                      "M.1"=paste0("IGHV M\n",  input$omic2P1),
                                      "U.1"=paste0("IGHV U\n", input$omic2P1),
                                      "0.M"="WT\nIGHV M",
                                      "0.U"="WT\nIGHV U",
                                      "1.M"=paste0("IGHV M\n",  input$omic1P1),
                                      "1.U"=paste0("IGHV U\n", input$omic1P1), 
                                      "0.0"="WT \n WT",
                                      "0.1"=paste0("WT \n", input$omic2P1),
                                      "1.0"=paste0("WT \n", input$omic1P1),
                                      "1.1"=paste0(input$omic1P1, "\n", input$omic2P1))) +
            xlab("")+
            ylab("Log(Viability)") +
            t2
        }
          
          
          
         
        })
        
 
      #TAB4 OUTPUTS  -------------------------------------------------------
      output$drugTarget3 <- renderInfoBox({
        infoBox(
          "Drug Target", 
          if(input$drugP3 == "DMSO"){"Control Treatment"}else{filter(drugs, Name == input$drugP3) %>% select(main_targets)}, 
          icon = icon("bullseye"),
          color = "purple"
        )
      })
      
      output$cytTarget3 <- renderInfoBox({
        infoBox(
          "Stimulus Pathway", if(input$cytP3 == "Not Cytokine"){"No Cytokine"}else{filter(Cytokines, name == input$cytP3) %>% select(pathway)}, icon = icon("dot-circle"),
          color = "yellow"
        )
      })
      
      output$lassoPlot <- renderPlot({
        test <- paste0(input$drugP3, input$cytP3)
        validate(
          need(
            test != "DMSONo Cytokine", "Treatment cannot be double control. Please chose a different combination.")
        )
        
        
        seaName <- paste(input$drugP3, input$cytP3, sep = ":")
        barValue <- rowMeans(coefficients_lasso[[seaName]])
        freqValue <- rowMeans(abs(sign(coefficients_lasso[[seaName]])))
        barValue <- barValue[abs(barValue) >= input$coefCutP3 & freqValue >= input$freqCutP3]
        barValue <- barValue[order(barValue)]
        
        validate(
          need(length(barValue) != 0, "There are no significant genetic predictors for this drug/stimulus, please choose a different combination.")
        )
        
        
        
        #for the heatmap and scatter plot below the heatmap
        allData <- geneMatrix
        viabValue <- unlist(viabMatrix[seaName,])
        
        tabValue <- allData[, names(barValue),drop=FALSE]
        ord <- order(viabValue)
        viabValue <- viabValue[ord]
        tabValue <- tabValue[ord, ,drop=FALSE]
        sampleIDs <- rownames(tabValue)
        tabValue <- as.tibble(tabValue)
        
        
        tabValue$Sample <- sampleIDs
        
        #Mark different rows for different scaling in heatmap
        #annotate mutations by mutation, methylation or IGHV
        matValue <- gather(tabValue, key = "Var",value = "Value", -Sample)
        matValue$Type <- "mut"
        
        #for Methylation Cluster
        matValue$Type[grep("Methylation",matValue$Var)] <- "meth"
        
        #for IGHV status
        matValue$Type[grep("IGHV",matValue$Var)] <- "ighv"
        
        #change the scale of the value so that IGHV, Methylation and Mutation do not overlap
        matValue[matValue$Type == "mut",]$Value = matValue[matValue$Type == "mut",]$Value + 10
        matValue[matValue$Type == "meth",]$Value = matValue[matValue$Type == "meth",]$Value + 20
        matValue[matValue$Type == "ighv",]$Value = matValue[matValue$Type == "ighv",]$Value + 30
        
        #change continuous to categorical
        matValue$Value <- factor(matValue$Value,levels = sort(unique(matValue$Value)))
        
        
        #change order of heatmap
        matValue$Var <- factor(matValue$Var, levels = names(barValue))
        matValue$Sample <- factor(matValue$Sample, levels = names(viabValue))
        
        #plot the heatmap NB mutated= 11 = black
        p1 <- ggplot(matValue, aes(x=Sample, y=Var)) + 
          geom_tile(aes(fill=Value), color = "white") + #ghost white
          theme_bw()+
          scale_y_discrete(expand=c(0,0)) + 
          theme(axis.title.y = element_text( size=16),
                axis.text.y=element_text(hjust=0, size=14), 
                axis.ticks=element_blank(),
                panel.border=element_rect(colour="gainsboro"),  
                plot.title=element_text(face="bold", size = 18, margin = margin(t = -5, b = 1)), 
                panel.background=element_blank(),
                panel.grid.major=element_blank(), 
                panel.grid.minor=element_blank()) + 
          xlab("Mutation status for each patient") + 
          ylab("") + 
          scale_fill_manual(name="Mutated", 
                            values=c(`10`= offwhite,  #WT
                                     `11`="#373A36", #Mutant
                                     `20`= offwhite, #LP
                                     `20.5`= "#707372", #IP
                                     `21` = "#A8A99E", #HP
                                     `30` = offwhite, #IGHV-U
                                     `31` = "#707372"), #IGHV-M
                            guide=FALSE) + 
          
          ggtitle(ifelse(input$cytP3 =="No Cytokine" , paste(input$drugP3), 
                         ifelse(input$drugP3 == "DMSO", paste(input$cytP3), 
                                paste(input$drugP3, "and", input$cytP3, sep = " "))))
        
        
        
        
        #Plot the bar plot on the left of the heatmap 
        barDF = data.frame(barValue, nm=factor(names(barValue),levels=names(barValue)))
        
        p2 <- ggplot(data=barDF, aes(x=nm, y=barValue)) + 
          geom_bar(stat="identity", 
                   fill=ifelse(barValue<0,
                               palblues[6],palreds[8]), 
                   colour="black", 
                   size=0.3) + 
          scale_x_discrete(expand=c(0,0.5)) + 
          scale_y_continuous(expand=c(0,0)) + 
          coord_flip() + #changed from min(barValue) and max(barValue)
          theme(panel.grid.major=element_blank(), 
                panel.background=element_blank(), 
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.text=element_text(size=11, angle = 45, hjust = 1.5), 
                panel.border=element_blank()) +
          ylab("Size of predictor") + 
          geom_vline(xintercept=c(0.5), 
                     color="black", 
                     size=0.6)
        
        #Plot the scatter plot under the heatmap
        scatterDF = data.frame(X=factor(names(viabValue), 
                                        levels=names(viabValue)), 
                               Y=unlist(viabValue))
        
        p3 <- ggplot(scatterDF, aes(x=X, y=Y)) + 
          geom_point(shape=21, 
                     fill="dimgrey", 
                     colour=colors[1], 
                     size=1.2) + 
          theme_bw() +
          theme(panel.grid.minor=element_blank(), 
                panel.grid.major.x=element_blank(), 
                axis.ticks.x=element_blank(), 
                axis.text.y=element_text(size=12), 
                panel.border=element_rect(colour="dimgrey", size=0.1),
                panel.background=element_rect(fill="white")) +
          xlab("Log(Viability) for each sample, with selected treatment")
        
        
        #Assemble all the plots togehter
        
        # construct the gtable
        wdths = c(0.2, 1.5, 0.2, 0.9*ncol(matValue), 2, 0.2)
        hghts = c(0.2, 0.2, 0.0020*nrow(matValue), 0.16, 1, 0.2, 1.5)
        gt = gtable(widths=unit(wdths, "in"), heights=unit(hghts, "in"))
        
        ## make grobs
        gg1 = ggplotGrob(p1)
        gg2 = ggplotGrob(p2)
        gg3 = ggplotGrob(p3)
        
        
        #make legend grob
        legendFor = c("G", "I", "M")
        
        coldef<-list()
        coldef["I"] <- list(c(offwhite,"#707372")) #U-CLL and M-CLL
        coldef["M"] <- list(c(offwhite, "#707372", "#A8A99E")) #IP, HP, LP
        coldef["G"] <- list(c(offwhite, "#373A36")) #WT and Mutated

      
        
        legends = HGmakelegends(legendFor=c("G","I", "M"),coldef)
        gg4 <- legends[["plot"]]
        
        ## fill in the gtable
        
        #HEATMAP
        #5:1 = "PREDICTORS"
        gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 3, 4) # add heatmap
        gt = gtable_add_grob(gt, gtable_filter(gg1, "panel"), 3, 4) #add legend
        gt = gtable_add_grob(gt, gtable_filter(gg1, "title"), 1, 4) #add title to plot
        gt = gtable_add_grob(gt, gtable_filter(gg1, "axis-l"), 3, 5) # variable names
        gt = gtable_add_grob(gt, gtable_filter(gg1, "xlab-b"), 2, 4) # axis title
        
        #BARPLOT
        gt = gtable_add_grob(gt, gtable_filter(gg2, "panel"), 3, 2) # add barplot
        gt = gtable_add_grob(gt, gtable_filter(gg2, "axis-b"), 4, 2) # y axis for barplot
        gt = gtable_add_grob(gt, gtable_filter(gg2, "xlab-b"), 2, 2) # y lab for barplot
        
        
        #SCATTER PLOT
        gt = gtable_add_grob(gt, gtable_filter(gg3, "panel"), 5, 4) # add scatterplot
        gt = gtable_add_grob(gt, gtable_filter(gg3, "xlab-b"), 6, 4) # x label for scatter plot
        gt = gtable_add_grob(gt, gtable_filter(gg3, "axis-l"), 5, 3) #  axis for scatter plot
        
        #LEGEND
        gt = gtable_add_grob(gt, gg4, 7, 4) #  legend
        
        #plot
        grid.newpage()
        grid.draw(gt)
        
      })
      
      
      
      #TAB5 OUTPUTS  -------------------------------------------------------
      output$drugTarget4 <- renderInfoBox({
        infoBox(
          "Drug Target", 
          if(input$drugP4 == "DMSO"){"DMSO"}else{filter(drugs, Name == input$drugP4) %>% select(main_targets)}, 
          icon = icon("bullseye"),
          color = "purple"
        )
      })
      
      output$cytTarget4 <- renderInfoBox({
        infoBox(
          "Stimulus Pathway", if(input$cytP4 == "Not Cytokine"){"No Cytokine"}else{filter(Cytokines, name == input$cytP4) %>% 
              select(pathway)}, icon = icon("dot-circle"),
          color = "yellow"
        )
      })
      
      output$lassoIntPlot <- renderPlot({
        
        seaName.Int <- paste(input$drugP4, input$cytP4, sep = " + ")
        barValue.Int <- rowMeans(coefficients_int[[seaName.Int]])
        freqValue.Int <- rowMeans(abs(sign(coefficients_int[[seaName.Int]])))
        barValue.Int <- barValue.Int[abs(barValue.Int) >= input$coefCutP4 & freqValue.Int >= input$freqCutP4]
        barValue.Int <- barValue.Int[order(barValue.Int)]
        
        validate(
          need(length(barValue.Int) != 0, "There are no significant genetic predictors for this drug-stimulus interaction, please choose a different combination.")
        )
        
        
        
        
        
        #for the heatmap and scatter plot below the heatmap
        allData.Int <- geneMatrixInt
        viabValue.Int <- unlist(betaMatrix[seaName.Int,])
        
        tabValue.Int <- allData.Int[, names(barValue.Int),drop=FALSE]
        ord.Int <- order(viabValue.Int)
        viabValue.Int <- viabValue.Int[ord.Int]
        tabValue.Int <- tabValue.Int[ord.Int, ,drop=FALSE]
        sampleIDs.Int <- rownames(tabValue.Int)
        tabValue.Int <- as.tibble(tabValue.Int)
        
        
        
        tabValue.Int$Sample <- sampleIDs.Int
        #Mark different rows for different scaling in heatmap
        matValue.Int <- gather(tabValue.Int, key = "Var",value = "Value", -Sample)
        matValue.Int$Type <- "mut"
        
        #for ighv.status
        matValue.Int$Type[grep("IGHV",matValue.Int$Var)] <- "ighv"
        
        #change the scale of the value, let them do not overlap with each other
        matValue.Int[matValue.Int$Type == "mut",]$Value = matValue.Int[matValue.Int$Type == "mut",]$Value + 10
        matValue.Int[matValue.Int$Type == "ighv",]$Value = matValue.Int[matValue.Int$Type == "ighv",]$Value + 30
        
        #change continuous to catagorical
        matValue.Int$Value <- factor(matValue.Int$Value,levels = sort(unique(matValue.Int$Value)))
        
        #change order of heatmap
        matValue.Int$Var <- factor(matValue.Int$Var, levels = names(barValue.Int))
        matValue.Int$Sample <- factor(matValue.Int$Sample, levels = names(viabValue.Int))
        
        #plot the heatmap NB mutated= 11 = black
        p1.Int <- ggplot(matValue.Int, aes(x=Sample, y=Var)) + 
          geom_tile(aes(fill=Value), color = "white") + #ghost white
          theme_bw()+
          scale_y_discrete(expand=c(0,0)) + 
          theme(axis.title.y = element_text( size=16),
                axis.text.y=element_text(hjust=0, size=14), 
                axis.ticks=element_blank(),
                panel.border=element_rect(colour="gainsboro"),  
                plot.title=element_text(face="bold", size = 18, margin = margin(t = -5, b = 1)), 
                panel.background=element_blank(),
                panel.grid.major=element_blank(), 
                panel.grid.minor=element_blank()) + 
          xlab("Mutation status for each patient") + 
          ylab("") + 
          scale_fill_manual(name="Mutated", 
                            values=c(`10`= offwhite,  #WT
                                     `11`= "#373A36", #Mutant
                                     `30` = offwhite, #IGHV-U
                                     `31` = "#707372"), #IGHV-M
                            guide=FALSE) + 
          ggtitle(paste(" Interaction between:\n", input$drugP4, "and", input$cytP4, sep = " "))
        
        
        
        
        #Plot the bar plot on the left of the heatmap 
        barDF.Int = data.frame(barValue.Int, nm = factor(names(barValue.Int),levels = names(barValue.Int)))
        
        p2.Int <- ggplot(data=barDF.Int, aes(x = nm, y = barValue.Int)) + 
          geom_bar(stat="identity", 
                   fill=ifelse(barValue.Int<0,
                               palblues[6],palreds[8]), 
                   colour="black", 
                   size=0.3) + 
          scale_x_discrete(expand = c(0,0.5)) + 
          scale_y_continuous(expand = c(0,0)) + 
          coord_flip() + #changed from min(barValue) and max(barValue)
          theme(panel.grid.major = element_blank(), 
                panel.background = element_blank(), 
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.text = element_text(size = 11, angle = 45, hjust = 1.5), 
                panel.border = element_blank()) +
          ylab("Size of coefficient") + 
          geom_vline(xintercept = c(0.5), 
                     color = "black", 
                     size = 0.6)
        
        #Plot the scatter plot under the heatmap
        scatterDF.Int = data.frame(X=factor(names(viabValue.Int), 
                                            levels=names(viabValue.Int)), 
                                   Y=unlist(viabValue.Int))
        
        p3.Int <- ggplot(scatterDF.Int, aes(x=X, y=Y)) + 
          geom_point(shape=21, 
                     fill="dimgrey", 
                     colour=colors[3], 
                     size=1.2) + 
          theme_bw() +
          theme(panel.grid.minor=element_blank(), 
                panel.grid.major.x=element_blank(), 
                axis.ticks.x=element_blank(), 
                axis.text.y=element_text(size=12), 
                panel.border=element_rect(colour="dimgrey", size=0.1),
                panel.background=element_rect(fill="white")) +
          xlab(expression(paste("Patient-specific ", beta["int"])))
        
        
        #Assemble all the plots together
        
        # construct the gtable
        wdths.Int = c(0.2, 1.5, 0.2, 0.9*ncol(matValue.Int), 2, 0.2)
        hghts.Int = c(0.5, 0.2, 0.0020*nrow(matValue.Int), 0.16, 1, 0.2, 1.5)
        gt.Int = gtable(widths=unit(wdths.Int, "in"), heights=unit(hghts.Int, "in"))
        
        ## make grobs
        gg1.Int = ggplotGrob(p1.Int)
        gg2.Int = ggplotGrob(p2.Int)
        gg3.Int = ggplotGrob(p3.Int)
        
        
        #make legend grob
        
        
        coldef.Int <-list()
        coldef.Int["I"] <- list(c(offwhite, "#707372"))
        coldef.Int["G"] <- list(c(offwhite, "#373A36")) 
      
        
        legends.Int = HGmakelegends(legendFor=c("G","I"),coldef.Int)
        gg4.Int <- legends.Int[["plot"]]
        
        ## fill in the gtable
        
        #HEATMAP
        #5:1 = "PREDICTORS"
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg1.Int, "panel"), 3, 4) # add heatmap
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg1.Int, "panel"), 3, 4) #add legend
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg1.Int, "title"), 1, 4) #add title to plot
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg1.Int, "axis-l"), 3, 5) # variable names
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg1.Int, "xlab-b"), 2, 4) # axis title
        
        #BARPLOT
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg2.Int, "panel"), 3, 2) # add barplot
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg2.Int, "axis-b"), 4, 2) # y axis for barplot
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg2.Int, "xlab-b"), 2, 2) # y lab for barplot
        
        
        #SCATTER PLOT
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg3.Int, "panel"), 5, 4) # add scatterplot
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg3.Int, "xlab-b"), 6, 4) # x label for scatter plot
        gt.Int = gtable_add_grob(gt.Int, gtable_filter(gg3.Int, "axis-l"), 5, 3) #  axis for scatter plot
        
        #LEGEND
        gt.Int = gtable_add_grob(gt.Int, gg4.Int, 7, 4) #  legend
        
        #plot
        grid.newpage()
        grid.draw(gt.Int)
        
        
      })
      
      
    
    
      
    
    }

# Run the application 
shinyApp(ui = ui, server = server)

