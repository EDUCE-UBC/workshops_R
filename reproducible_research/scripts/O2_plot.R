# load packages
library(tidyverse)

# load data
dat <- read.csv("https://raw.githubusercontent.com/EDUCE-UBC/workshops_data_science/master/reproducible_research/data/data.csv")

# create plot of oxygen by depth
O2_plot <- quickplot(data=dat,
                     x=O2_uM, 
                     y=Depth_m, 
                     colour=Season, 
                     main="Saanich Inlet: Seasonal oxygen depth profile")

# save the plot
ggsave("O2_plot.png")