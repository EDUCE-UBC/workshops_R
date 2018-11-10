# If you don't have it, install ggplot2 package
# install.packages("ggplot2")

# Or run an ifelse statement to install only if not already present on your computer
if(!("ggplot2" %in% installed.packages()[,1])) install.packages("ggplot2")

# Load packages
library(ggplot2)

# Load data
## No formatting
dat <- read.table(file="data.csv")

# Help function
?read.table

# Read in data
## With formatting
dat <- read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

## If there were NAs as many different forms
dat <- read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=FALSE, na.strings=c("NA", "NAN", " ", "."))

## Other sep
# tab = "\t" =/= "  "

# Exploring variables types
is.character(dat$Season)

## Obtain a single column by name
# dat$name_of_column like
dat$Season

# Check data type (class in R)
class(dat$Season)
class(dat$O2_uM)
class(dat$Add_data)

# Maths!
## Variance
var(dat$O2_uM)
# If there were NAs
var(dat$O2_uM, na.rm=TRUE)

# Conditional statements
## Ask R when something is TRUE/FALSE
### Where is oxygen (O2_uM) > 50?
dat$O2_uM > 50

# Other useful statements
## >= or <= greater/less than or equal to
## == equals
## != not equals
## is.na() missing data
## !is.na() not missing

# Other useful functions
unique(dat$Season)
table(dat$Season)
# Summarize counts
table(dat$O2_uM > 50)
table(is.na(dat$O2_uM))
table(is.na(dat))

# Accessing pieces of dat
# data frames are 2 dimensional so they are [rows, columns]
# vectors are 1 dimensional so they are [values]

## Access the 5th row, 3rd column
dat[5,3]
dat[5,"O2_uM"]
dat$O2_uM[5]

# Subsetting data
## Find the value
logical.vector <- dat$O2_uM == 69.828
logical.vector
## Subset to that value
dat[logical.vector,]
dat[logical.vector,"O2_uM"]

# Another example
dat[dat$O2_uM == 69.828,]
dat[dat$O2_uM == 0,]

# Plotting
quickplot(data=dat,
          x=O2_uM, 
          y=Depth_m, 
          colour=Season, 
          main="Saanich Inlet: Seasonal oxygen depth profile") +
  geom_line()

ggsave()

# In base R
plot(x=dat$O2_uM, y=dat$Depth_m)