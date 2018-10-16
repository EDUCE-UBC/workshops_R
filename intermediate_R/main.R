################################################################################
# load packages
################################################################################
library(tidyverse)
library(lmerTest)
library(devtools)
library(roxygen2)

################################################################################
# loading and cleaning data
# See The R tidyverse workshop for details
################################################################################
raw_data <- read_csv(file="Saanich_Data.csv",
                     col_names=TRUE,
                     na=c("", "NA", "NAN", "ND"))

dat <- raw_data %>%
  select(Cruise, Date, Depth, 
         Temperature, Salinity, Density, 
         WS_O2, WS_NO3, WS_H2S) %>%
  filter(Date >= "2008-02-01") %>%
  rename(O2=WS_O2, NO3=WS_NO3, H2S=WS_H2S) %>%
  mutate(Depth=Depth*1000)

################################################################################
# Classes and attributes of objects
# Vector objects
################################################################################
# Numeric
x <- 1
class(x)
mode(x)

# Integer
x <- 1L
class(x)
mode(x)

# Character
x <- "1L"
class(x)
mode(x)

# Logical
x <- TRUE
class(x)
mode(x)

# Matrix (complex)
x <- c(1,1)
dim(x) <- c(2,1)
class(x)
mode(x)
attributes(x)

### Exercise.
# 1. Assign x the value `"a"`. What are its class and mode?
# 2. Give it dimensions `c(1,1)`. What are its class and mode?

# 1.

# 2.

### End exercise

################################################################################
# Complex objects
################################################################################
# Built in functions
class
class(class)
attributes(class)

# Graphics
## Save a graphic to an object
p1 <- dat %>% 
  ggplot(aes(x=O2, y=Depth)) +
  geom_point() +
  geom_smooth(method="lm") + 
  scale_y_reverse(limits=c(200, 0)) +
  labs(x=expression(O[2]*" "*(mu*g/mL)),
       y="Depth (m)", 
       title="Oxygen decreases with depth and is less variable at lower depths")

p1
class(p1)
attributes(p1)

# Statistical models
m1 <- lm(O2 ~ Depth, data=dat)
m1
class(m1)
attributes(m1)

# Summaries
s1 <- summary(m1)
s1
class(s1)
attributes(s1)

################################################################################
# Data objects
################################################################################
# Data frame
dat.df <- dat[1:5, 1:5]
class(dat.df)

# Matrix
dat.mat <- as.matrix(dat.df)
class(dat.mat)

# Data frame vs. matrix
dat.df
dat.mat

attributes(dat.df)
attributes(dat.mat)

# Listing potential functions for a given object class
methods(class="data.frame")

### Exercise.
# 1. Calculate a summary table of `dat`. What are its class and attributes?
# 2. Read in the raw data table `Saanich_Data.csv` using the base R function `read.table`. What are this object's class and attributes? Are they any different from the object create when we used `read_csv` to read in the same data?

# 1. 

# 2.

### End exercise

################################################################################
# The R list object
################################################################################
# Make a list object
x <- list(data = dat,
          Function = lm,
          model = m1,
          random = list(sample(LETTERS, 5), sample(1:1000, 5)))

# Explore the list object
class(x)
mode(x)
length(x)
attributes(x)
x$random

# Explore an R created list object
m1$coefficients

### Exercise.

# 1. Calculate the `summary()` of `m1` and save it as `m2`.
# 2. What is the class and mode of `m2`?
# 3. Using a single line of code, pull out just the p-values from `m2`.
#     - *Hint*: You will need to use both `$` and `[ ]`.

# 1. 

# 2. 

# 3. 

### End exercise.

################################################################################
# S4 objects
################################################################################
# linear mixed effects model (lmer)
m3 <- lmer(O2 ~ Cruise + (0 + Cruise | Depth), dat)

class(m3)
mode(m3)

# Many of the objects in m3 are very large so we will just look at the names here
names(attributes(m3))

# Accessing objects within S4
m3@beta

# S3 objects within S4
VC <- VarCorr(m3)
class(VC)
attributes(VC)

### Exercise.
# Exercise.

# 1. Compute and store the variance-covariance matrix of `m3` using `vcov()`.
# 2. What class and mode is it?
# 3. What elements does it contain?
# 4. What are the dimensions of `factors` within this object?
  
# 1. 

# 2.

# 3.

# 4.

### End exercise

################################################################################
# Generic functions
################################################################################
# the function class
class(function(){}) 

# A generic function
f <- function(x) 
{ 
  return(x*x)  
} 

# Running a generic function
f(5)

# Storage of a generic function
f

### Exercise.
# 1. Put the following math into a function $$ f(x) = 1 + 2x - 5x^2 + x^3 $$
# 2. Set x to `1:1000/1000*6-1`
# 3. Plot the results with plot(x, f(x) , main="The answer looks like this")

# 1.

# 2.

# 3. 

### End exercise.

################################################################################
# Function arguments
################################################################################
# Example function
f <- function(a,b){
  out = a/b
  return(out)
}

# Arguments in any order with labels
f(a=1, b=2)
f(b=2, a=1)

# Arugments only in correct order if omit names
f(1, 2)
## Does not equal
f(2, 1)

################################################################################
# Triple dot argument
################################################################################
# Example function
f <- function(x, ...)
{ 
  y = cos(x) 
  plot(x, y, ...)
}

# Only using required inputs
f(1:1000/50)
# Adding additional inputs through ...
f(1:1000/50, type='l') 
f(1:1000/50, type='l', col="purple") 

################################################################################
# Scoping
################################################################################
# you can define a variable globally and then use that variable as-is within your function
z <- 1 # Define z globally

f <- function(x) # no local definition of z is given
{ 
  return(x + z) 
}

f(x=2) # Failing to find z in the function or inputs, R looks for it globally and finds z = 1 to result in 2 + 1 = 3

# you can override globally defined variables within your function
z <- 1 # Define z globally

f <- function(x) 
{ 
  z = 2 # local definition of z overrides global definition 
  return(x + z) 
}

f(x=2) # R finds the function's definition of z=2 first and thus gives 2 + 2 = 4

# if your function calls for an input that is defined globally, it will fail.
z <- 1

f <- function(x, z) # function calls for inputs for both x and z 
{ 
  return(x + z) 
}

f(x=2) # Function expects inputs for x and z but not finding z, fails


# if a variable is defined neither in the function nor in the global environment, your function will not run.
y

f = function(x)
{
  y = 2*y # y does not exist anywhere
  return(x + y)
}

f(3) 

### Exercise. For the following, try to determine what the function will return *without running the code* then check yourself by running it.

#1. Remove all instances of x, z, and f( ) from your environment so that you are starting fresh for this exercise.

rm(x)
rm(z)
rm(f)

# 2. What happens when we run `f()`? Why?
f <- function()
{
  return(2*x) 
}

f() 

# 3. What will `f()` return? Why?
x <- 1 

f <- function()
{
  x = 2 
  return(2*x) 
}

f() 


# 4. What does the final `y` call return?
y <- 1

f <- function(x)
{
  y = x+2 
  return(x*y)
}

f(1)

y

# 1. NA

# 2. 

# 3.

# 4. 

### End exercise

################################################################################
# Building a function
# Step 1: Define and test your task
################################################################################
# Subset the data
dat.subset <- dat %>% filter(Cruise == 72)
# Calculate a linear model
model <- lm(dat.subset$O2 ~ dat.subset$Depth)
# Sumamrize the model
sum <- summary(model)
sum
# Extract p-values
pval <- sum$coefficients[,"Pr(>|t|)"]
# View extracted p-values
print(pval)

################################################################################
# Building a function
# Step 2: Turn the task into a function
################################################################################
# Minimum inputs of data and y variable
lm.function <- function(data, y){ ###
  dat.subset <- data %>% filter(Cruise == 72)
  
  model <- lm(dat.subset[[y]] ~ dat.subset$Depth) ###
  
  sum <- summary(model)
  
  pval <- sum$coefficients[,"Pr(>|t|)"]
  
  print(pval)
} ###

# Test the function
lm.function(data=dat, y="O2")

# All possible inputs of data, cruise #, x and Y variables
lm.function <- function(data, cruise, x, y){ ###
  dat.subset <- data %>% filter(Cruise == cruise)
  
  model <- lm(dat.subset[[y]] ~ dat.subset[[x]]) ###
  
  sum <- summary(model)
  
  pval <- sum$coefficients[,"Pr(>|t|)"]
  
  print(pval)
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y="O2")

################################################################################
# Building a function
# Step 3: Add packages
################################################################################
lm.function <- function(data, cruise, x, y){
  require(tidyverse) ###
  
  dat.subset <- data %>% filter(Cruise == cruise)
  
  model <- lm(dat.subset[[y]] ~ dat.subset[[x]])
  
  sum <- summary(model)
  
  pval <- sum$coefficients[,"Pr(>|t|)"]
  
  print(pval)
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y="O2")

################################################################################
# Building a function
# Step 4: Start your documentation
################################################################################
lm.function <- function(data, cruise, x, y){
  # Load necessary packages ###
  require(tidyverse)
  
  # Subset the data to the cruise of interest ###
  dat.subset <- data %>% filter(Cruise == cruise)
  
  # Fit a linear model ###
  model <- lm(dat.subset[[y]] ~ dat.subset[[x]])
  # Summarize the model ###
  sum <- summary(model)
  # Extract p-values from the summary ###
  pval <- sum$coefficients[,"Pr(>|t|)"]
  
  # Print p-values to the console ###
  print(pval)
}


################################################################################
# Building a function
# Step 5: Loop through multiple y-variables
################################################################################
# Example empty loop
for(a in b){
  # Perform some task(s) on each element in b, one after the other
  # A single element in b is represented by a
}

# Example loop
for(year in 2015:2020){
  phrase = paste("The year is", year)
  print(phrase)
}

# Us a different input name in the same example loop
for (date in 2015:2020){
  phrase = paste("The year is", date)
  print(phrase)
}

# Add a loop to the lm.function
lm.function <- function(data, cruise, x, y){
  # Load necessary packages
  require(tidyverse)
  
  # Subset the data to the cruise of interest
  dat.subset <- data %>% filter(Cruise == cruise)
  
  for(y.variable in y){ # Loop through all variables provided in y  ###
    # Fit a linear model 
    model <- lm(dat.subset[[y.variable]] ~ dat.subset[[x]]) ###
    # Summarize the model
    sum <- summary(model)
    # Extract p-values from the summary
    pval <- sum$coefficients[,"Pr(>|t|)"]
    
    # Print p-values to the console
    print(pval)
  }
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2","NO3"))

# Compare results to calculations "by-hand"
dat.subset <- dat %>% filter(Cruise == 72)
# Oxygen
model <- lm(dat.subset$O2 ~ dat.subset$Depth)
sum <- summary(model)
pval <- sum$coefficients[,"Pr(>|t|)"]
print(pval)
# Nitrate
model <- lm(dat.subset$NO3 ~ dat.subset$Depth)
sum <- summary(model)
pval <- sum$coefficients[,"Pr(>|t|)"]
print(pval)

### Exercise
# 1. Apply the current `lm.function` to all the available geochemical variables in the Saanich data set. Which ones appear to be significantly correlated with depth?
# 2. Copy the `lm.function` and alter it to print out the models' adjusted R-squared values instead of p-values. Be sure to run the function with inputs to make sure it works!

# 1.

# 2.

### End exercise

################################################################################
# Building a function
# Step 6.1: Save all iterations of the loop
################################################################################
# Investigate the pval object
class(pval)
attributes(pval)
is.vector(pval)
length(pval)

# Attemp to write p-values to a table on your hard-drive
write_csv(pval, path="pval_results.csv")

# Reformat results so that the above will be successful
write_csv(as.data.frame(pval), path="pval_results.csv")
# Read back in the results table to see how it was saved
read_csv("pval_results.csv")

# Add the save step to the lm.function
# Only saves last loop iteration
lm.function <- function(data, cruise, x, y){
  # Load necessary packages
  require(tidyverse)
  
  # Subset the data to the cruise of interest
  dat.subset <- data %>% filter(Cruise == cruise)
  
  for(y.variable in y){ # Loop through all variables provided in y
    # Fit a linear model 
    model <- lm(dat.subset[[y.variable]] ~ dat.subset[[x]])
    # Summarize the model
    sum <- summary(model)
    # Extract p-values from the summary
    pval <- sum$coefficients[,"Pr(>|t|)"]
    
    # Print p-values to a table
    write_csv(as.data.frame(pval), path="pval_results.csv") ###
  }
}

# Check the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("pval_results.csv")

# Modify save step to save all loop iterations
lm.function <- function(data, cruise, x, y){
  # Load necessary packages
  require(tidyverse)
  
  # Subset the data to the cruise of interest
  dat.subset <- data %>% filter(Cruise == cruise)
  
  for(y.variable in y){ # Loop through all variables provided in y
    # Fit a linear model 
    model <- lm(dat.subset[[y.variable]] ~ dat.subset[[x]])
    # Summarize the model
    sum <- summary(model)
    # Extract p-values from the summary
    pval <- sum$coefficients[,"Pr(>|t|)"]
    
    # Print p-values to a table
    write_csv(as.data.frame(pval), path="pval_results.csv", append=TRUE) ###
  }
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("pval_results.csv")

# Repeated use continues to append data to the table, resulting in duplicate results
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))

read_csv("pval_results.csv")

# Remove old results file before rerunning the function
lm.function <- function(data, cruise, x, y){
  # Load necessary packages
  require(tidyverse)
  # Remove old results file, if exists
  if(file.exists("pval_results.csv")){file.remove("pval_results.csv")}  ###
  
  # Subset the data to the cruise of interest
  dat.subset <- data %>% filter(Cruise == cruise)
  
  for(y.variable in y){ # Loop through all variables provided in y
    # Fit a linear model 
    model <- lm(dat.subset[[y.variable]] ~ dat.subset[[x]])
    # Summarize the model
    sum <- summary(model)
    # Extract p-values from the summary
    pval <- sum$coefficients[,"Pr(>|t|)"]
    
    # Print p-values to a table
    write_csv(as.data.frame(pval), path="pval_results.csv", append=TRUE)
  }
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("pval_results.csv", col_names=FALSE)

################################################################################
# Building a function
# Step 6.2: Beautify the output
# Step 6.2.1: Transpose
################################################################################
# Transpose the pval results table
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
    pval <- t(sum$coefficients[,"Pr(>|t|)"]) ###
    
    # Print p-values to a table
    write_csv(as.data.frame(pval), path="pval_results.csv", append=TRUE)
  }
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("pval_results.csv", col_names=FALSE)


################################################################################
# Building a function
# Step 6.2.2: Loop names
################################################################################
# Add a variable identifying each loop iteration by the variable it is testing
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
    pval <- as.data.frame(t(sum$coefficients[,"Pr(>|t|)"])) ###
    # Add y variable name label
    pval$variable <- y.variable ###
    
    # Print p-values to a table
    write_csv(pval, path="pval_results.csv", append=TRUE)
  }
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("pval_results.csv", col_names=FALSE)

################################################################################
# Building a function
# Step 6.2.3: Column names
################################################################################
# Name the output columns
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
  # Read in p-value results and add column names
  results <- read_csv("pval_results.csv", col_names=FALSE) %>%  ###
    rename(Intercept.p=X1, Depth.p=X2, variable=X3) %>%  ###
    # Re-write the results, now with column names
    write_csv(path="pval_results.csv") ###
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("pval_results.csv")

################################################################################
# Building a function
# Step 6.2.4: Fully dynamic naming
################################################################################
# Make the names dynamic based on the inputs
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
  col1 <- paste(colnames(pval)[1], "p", sep=".") ###
  col2 <- paste(x, "p", sep=".") ###
  
  # Create dynamic name fo results table
  table.name <- paste(x, "lm_pvals.csv", sep="_") ###
  
  # Read in p-value results and add column names
  read_csv("pval_results.csv", col_names=FALSE) %>% 
    rename(!!as.name(col1) := X1, ###
           !!as.name(col2) := X2, ###
           variable = X3) %>%  ###
    # Re-write the results, now with column names
    write_csv(path=table.name) ###
}

# Test the function
lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("Depth_lm_pvals.csv")

################################################################################
# Building a function
# Step 7: Complete documentation (and move to R script and share?)
################################################################################
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

# Load a function from a separate R script
source("lm.function.R")

### Exercise.
# 1. Using our final `lm.function`, determine the linear fits for all geochemical variables for Cruise 12.
# 2. Choose a different x variable and determine if any of the Saanich geochemical variables correlate with it.

# 1.

# 2.

### End exercise.


### CHALLENGE exercise.
#You may not have time in this workshop to complete these challenging exercises. However, we encourage you to complete them on your own to test your knowledge of functions and improve your coding skills!
  
# 1. How would you alter `lm.function` to accept sub-setting to multiple cruise #s? Hint: Think about using `%in%` when filtering the data.
# 2. How would you alter `lm.function` to output FDR corrected p-values?

### End challenge.

################################################################################
# Creating a package
# In a new Rproj
################################################################################
# Explore the NAMESPACE file
# Run in the Terminal, NOT the R console
`head testPackage/NAMESPACE`

# Update the DESCRIPTION and R script files then
# Build your package
devtools::document()

# Explore the updated NAMESPACE file
# Run in the Terminal, NOT the R console
`head testPackage/NAMESPACE`

# Explor the auto generated help page file
# Run in the Terminal, NOT the R console
`head -50 testPackage/man/lm.function.Rd`

################################################################################
# Test your package
################################################################################

# Load your local package into R
devtools::load_all()

# Test that it works on the previous data
library(tidyverse)

dat <- read_csv(file="../Saanich_Data.csv",
                col_names=TRUE,
                na=c("", "NA", "NAN", "ND")) %>% 
  select(Cruise, Date, Depth, 
         Temperature, Salinity, Density, 
         WS_O2, WS_NO3, WS_H2S) %>%
  filter(Date >= "2008-02-01") %>%
  rename(O2=WS_O2, NO3=WS_NO3, H2S=WS_H2S) %>%
  mutate(Depth=Depth*1000)

lm.function(data=dat, cruise=72, x="Depth", y=c("O2", "NO3"))
read_csv("Depth_lm_pvals.csv")

# Run a full check of you package including file paths, dependencies, etc.
devtools::check()