"
lm.function: Calculates p-values from linear models of geochemical data from a single Saanich Inlet cruise

Kim Dill-Mcfarland
kadm@mail.ubc.ca
University of British Columbia

LICENSE
Copyright (C) 2018 Kim Dill-McFarland

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

INPUTS
data: a data.frame or tibble containing variables of interest
cruise: Cruise # for subsetting data, 18-100 available in current data
x: x-variable in lm, for example 'Depth'
y: y-variables in lm, for example c('O2', 'NO3')

OUTPUTS
data frame of p-values named [x]_lm_pvals
"

lm.function <- function(data, cruise, x, y){
  # Load necessary packages
  require(tidyverse)
  
  # Create an empty list to hold results
  pval = list()
  
  # Subset the data to the cruise of interest
  dat.subset <- data %>% filter(Cruise == cruise)
  
  for(y.variable in y){ # Loop through all variables provided in y
    # Fit a linear model 
    model <- lm(dat.subset[[y.variable]] ~ dat.subset[[x]])
    # Summarize the model
    sum <- summary(model)
    # Extract p-values from the summary. Save into the pval list based on the y.variable name
    pval[[y.variable]] <- sum$coefficients[,"Pr(>|t|)"]
    
  }
  # Bind all results into 1 single object
  pval <- as.data.frame(do.call(rbind,pval))
  
  # Create dynamic column names
  col1 <- paste(colnames(pval)[1], "p", sep=".")
  col2 <- paste(x, "p", sep=".")
  table.name <- paste(x, "lm_pvals", sep="_")
  
  # Rename columns
  pval <- pval %>% ###
    rename_at(vars(colnames(pval)), ~c(col1, col2))
  # Rename output table and save to environment
  assign(table.name, pval, envir = .GlobalEnv)
}