### ggplot themes

fontsize=14

## theme for ggplots
t1<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_line(),
  panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_line(size=.4),
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  axis.text.x  = element_text(angle=90, size=16, face="bold", hjust = 1, vjust = 0.4),
  axis.text.y = element_text(size = 18),
  axis.ticks.x = element_line(linetype = "dotted"),
  axis.ticks.length = unit(0.3,"cm"),
  axis.title.x = element_text(face="bold", size=18), 
  axis.title.y = element_text(face="bold", size=18),
  plot.title = element_text(face="bold", size=18, hjust = 0.5),
  strip.text = element_text(size = fontsize)
)

t2<-t1+
  theme( axis.text.x  = element_text(angle=0, size=16, face="bold", hjust = 0.5, vjust = 1))

## theme for legends
t.leg <-  theme(legend.title = element_text(face='bold', 
                                            hjust = 1, size=11),
                legend.key = element_blank(),
                legend.text = element_text(size=12),
                legend.background = element_rect(color = "black"))



### Set colour palettes

#For Categorical: 
colors <- c("#A1BE1F", #green
            "#F4C61F", #yellow
            "#734595", #purple
            "#D41645", #red
            "#3B6FB6", #blue
            "#B65417", #orange
            "#E2E868", #light green
            "#CBA3D8", #light purple
            "#E58F9E", #light purple
            "#8BB8E8", #light blue
            "#F49E17", #light orange
            "#303030", #black
            "#A8A99E", #grey
            "#007B53") #dark green

#For Divergent: 
Divergent <- c("#003DA5", "#2055B0", "#406EBC", "#6086C7", "#809ED2", "#9FB6DD", "#BFCFE9", "#DFE7F4", "white", "white", "white","#F4E0E7", "#E9C2CF", "#DEA3B6", "#D3849E", "#C76586", "#BC476E", "#B12855", "#A6093D")

#for negatives only: 
palblues <- c("#003DA5", "#2055B0", "#406EBC", "#6086C7", "#809ED2", "#9FB6DD", "#BFCFE9", "#DFE7F4")

#for positives only:
palreds <- c("#F4E0E7", "#E9C2CF", "#DEA3B6", "#D3849E", "#C76586", "#BC476E", "#B12855", "#A6093D")

#For mutations: 
Mutant <- c("#b5b5b5","#373A36")
Sex <- c("#707372","#D0D0CE")
IGHV <- c("#373A36","#D0D0CE")
Methylation_cluster <- c("#373A36","#A8A99E","#D0D0CE")

#For drugs: 
drugpal <- c("#734595", "#CBA3D8") #purples

#For cytokines: 
cytpal <- c("#F49E17", "#EFC06E") #yellows


#neutral
offwhite <- "#f8f8ff"
lightergrey <- "#D0D0CE"
darkergrey <- "#707372"


na_color="#f0f0f0"
