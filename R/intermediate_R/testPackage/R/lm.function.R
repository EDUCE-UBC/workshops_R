#' Linear models across multiple variables
#'
#' Estimates p-values from linear models of geochemical data from a single Saanich Inlet cruise
#'
#' @param data data frame
#' @param cruise subsetting variable (numeric)
#' @param x independent variable in linear model (only singular accepted)
#' @param y dependent variable in linear model (list accepted)
#' @return a p-value table
#'
#' @importFrom dplyr filter
#' @export
lm.function <- function(data, cruise, x, y){
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
