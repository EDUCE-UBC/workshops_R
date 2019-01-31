##################################
# R as a calculator
##################################
# Using R to do basic calculations
1+4
2*28 # Multiplication
5^3 # Exponents

## You would run these in the console, not saved in a script, since you don't need to reuse them

# Assigning variables
## Create a new variable called 'x' and assign the number 5 to it
x <- 5
## This object appears in the Environment (upper right quadrant)
## You can then do calculations on anything saved in the Environment like so
x^2

##################################
# Packages
##################################

## R comes with many packages but there are 1000s more! See https://cran.r-project.org/web/packages/available_packages_by_name.html
## You can install any one of these packages on your own computer with
install.packages("ggplot2")

## This downloads the file(s) on your computer but does not tell R to access them
## To use functions in a package, you much load it *EVERY TIME* you open R
library(ggplot2)

## In general, install.packages() is run in the console once but library() is put in a script (like here) because you will need to run it many times

##################################
# Loading data
##################################
# You can read data into the console with 
read.table("data.csv")
# or assign the data to a variable named 'dat' in the Environment with
dat <- read.table("data.csv")

# We can add additional parameters to format the data correctly.
## sep indicates what character separates the columns
## header indicates if the first row (or header) is column names or not
## stringsAsFactors tells R is strings (a data type) should be read in as factors (another data type) or not. More on data types later

dat <- read.table("data.csv", sep=',', header=TRUE, stringsAsFactors=FALSE)

## Now you can see the full, formatted data
dat

##################################
# Accessing parts of the data
##################################
# Pull out an entire column by name using $
dat$Season
dat$Depth_m

# Or using indexing (numbers)
## Data frames are called as name[rows, columns]
## So calling the second column would be
dat[,2]
## And calling the 4th row, 2nd column would be
dat[4,2]

##################################
# Calculations with data
##################################
# You can do calculations on your data similar to how you would with that first 'x' variable we created
## For example, divide all Depth_m values by 1000 to get Depth in km
dat$Depth_m / 1000

##################################
# Help
##################################
# Access a help page for any function with ?
?read.table
# Or search within all help pages to find the function you want, if you don't know it by name
??"read table"

##################################
# Data types
##################################
# Ask R what the data type is with class()
## Character i.e. words/letters
class(dat$Season)
## Integer i.e. whole numbers
class(dat$Depth_m)
## Numeric i.e. numbers with decimal places
class(dat$O2_uM)
## Logical i.e TRUE/FALSE
class(dat$Add_data)

# Data types matter because R knows it can only run some functions on some data types. For example, you can't add character data since they are words

##################################
# Conditional statements
##################################
# Ask R about your data with conditional statements like 
## Greater than > and less than <
## Equal to ==
## Is missing is.na()

## For example, you can ask R which pieces of the Depth data are greater than 50
dat$Depth_m > 50
## This will return a string of TRUE/FALSE listing whether each Depth is or is not > 50

# If you save this TRUE/FALSE list, you can subset the data. **R will always choose to keep the TRUE values** so be sure to craft your statement accordingly
y <- dat$Depth_m > 50

dat[y,]

## If we compare the dimensions of dat and the y subset, we see that we've reduced the data
dim(dat)
dim(dat[y,])

## Another useful function is to ask R to like all the unique values of a variable like so. This is best used on categorical data like Season
unique(dat$Season)

##################################
# Logical operators
##################################
# These allow you to string together conditional statements like
## & and
## | or

## For example, subset the data to just Depths between 50 and 100
y2 <- dat$Depth_m > 50 & dat$Depth_m < 100

dat[y2,]

dim(dat[y2,])

##################################
# Plotting
##################################
# Use quickplot() from the ggplot2 package to quickly make a plot of oxygen across depths
## Note that you can spread a single function across multiple lines as line as the preceeding line ends with , to tell R to continue looking for inputs on the next line
quickplot(data=dat,
          x=O2_uM,
          y=Depth_m,
          color=Season,
          main="Saanich Inlet")

##################################
# Exercises
##################################
# 2. Using help to identify the necessary arguments for the log function compute the natural logarithm of 4, base 2 logarithm of 4, and base 4 logarithm of 4.
?log

log(4) #or
log(4, base=exp(1))

log2(4) #or
log(4, base=2)

log(4, base=4)

# 3. Using an R function, determine what data type the Depth_m variable is.
class(dat$Depth_m)

# 4. Using indexing and the square bracket operator `[]`:
##  - determine what depth value occurs in the 20th row
dat[20, "Depth_m"] #or
dat[20, 2] #or
dat$Depth_m[20]

##  - return the cell where oxygen equals 91.115
logical.vector <- dat$O2_uM == 91.115
dat[logical.vector, "O2_uM"]

#or all in one
dat[dat$O2_uM == 91.115, "O2_uM"]

# 5. Subset the data to observations where depth equals 100 m. *Hint*: Use a logical vector. 
logical.vector <- dat$Depth_m == 100
dat[logical.vector,]

# 6. Complete the following code to create a stacked scatterplot of oxygen concentrations within the two different seasons, colored by whether or not microbial data are available (example below).

quickplot(data=dat,
          x=Season, 
          y=O2_uM, 
          colour=Add_data, 
          main="Saanich Inlet: Oxygen in Fall vs. Summer")