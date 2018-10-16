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
"

lm.function <- function(data, cruise, x, y){
  # Load necessary packages
  require(tidyverse)
  # Remove old results file, if exists
  if(file.exists("pval_results.csv")){file.remove("pval_results.csv")}
  
  # Subset the data to the cruise of interest
  dat.subset <- data %>% filter(Cruise == cruise)
  
  for(y.variable in y){ # Loop through all variables provided in y
    # Fit a linear model 
    model <- lm(dat.subset[[y.variable]] ~ dat.subset[[x]])
    # Summarize the model
    sum <- summary(model)
    # Extract p-values from the summary
    # Reformat to a 1x2 data frame
    pval <- as.data.frame(t(sum$coefficients[,"Pr(>|t|)"]))
    # Add y variable name label
    pval$variable <- y.variable
    
    # Print p-values to a table
    write_csv(pval, path="pval_results.csv", append=TRUE)
  }
  # Create dynamic names for columns
  col1 <- paste(colnames(pval)[1], "p", sep=".")
  col2 <- paste(x, "p", sep=".")
  
  # Create dynamic name fo results table
  table.name <- paste(x, "lm_pvals.csv", sep="_")
  
  # Read in p-value results and add column names
  read_csv("pval_results.csv", col_names=FALSE) %>% 
    rename(!!as.name(col1) := X1,
           !!as.name(col2) := X2,
           variable = X3) %>% 
    # Re-write the results, now with column names
    write_csv(path=table.name)
}