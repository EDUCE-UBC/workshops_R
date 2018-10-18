# Install packages
## Only need to do once on your machine
# install.packages("ggplot2")

# Load package
library(ggplot2)
# Is the same as
library(package=ggplot2)

# Read in data
## Read into console
read.table(file="data.csv")

# Read the data in AND save it as an object named 'dat' in the R environment
### We want "Season" to remain a character object so we use stringsAsFactors
dat <- read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
# If we wanted R to change the "Season" variable to a factor type object, we would instead set stringsAsFactors to TRUE (or leave it out of the function since this is the default)
dat2 <- read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=TRUE)
# or
dat3 <- read.table(file="data.csv", sep=",", header=TRUE)

# Note: if you had a tab-delimited file, you would use sep="\t"

# Data types
## You can pull a column out of a table by the column name with the $ operator
dat$Season
dat$O2_uM

# Math in R
## You can run any basic math and statistics in R like
### variance
var(dat$O2_uM)
### mean/average
mean(dat$O2_uM)
# median
median(dat$O2_uM)
# You can also do math on entire data frame, but it will only work if all the data within that data frame are numeric (which these data are not)
median(dat)

# Conditions
## You can ask R true/false questions with conditional statements
## Here we ask which values of oxygen are greater than 10
dat$O2_uM > 10

# Counting observations that are TRUE/FALSE
table(dat$O2_uM > 10)

## You can then use this to subset the data where condition(s) are true
# For a data frame, you pull out rows and columns like dat[rows,columns]
## Here we return the rows of dat where oxygen is greater than 10
dat[dat$O2_uM > 10,]

# Some other useful conditions
## Not equals !=
## in %in%
## >= greater or equal to
## exactl equal to ==

# More math
## You can save the output of some maths to another object in R, like depth_km, and then call this new object just as you've been calling the original dat data frame
depth_km = dat$Depth_m / 1000
depth_km

# Find the unique character values of a variable
## Works best for character or factor variables, not numeric
unique(dat$Season)

# Pulling data out of a vector, 1 dimensional
# Similar to a data frame which is dat[rows,columns], you can use [ ] to get data out of a vector, only the vector only has 1 dimension (so there is no comma)
## If we take our save depth in km vector made above,
depth_km
# We can pull out the 5th value
depth_km[5]
# or use conditional statement to pull out all values less than 0.1
depth_km[depth_km < 0.1]

# Summarize all data in table
summary(dat)
summary(dat2)

# Plotting
## Base R
plot(x=dat$O2_uM, y=dat$Depth_m)

## vs. ggplot
quickplot(data=dat,
          x=O2_uM, 
          y=Depth_m)

# vs. pretty ggplot
## automatically created legends are the best!
quickplot(data=dat,
          x=O2_uM, 
          y=Depth_m, 
          colour=Season, 
          main="Saanich Inlet: Seasonal oxygen depth profile",
          xlab="Oxygen", ylab="Depth")

# Just like data, plots can be saved as an object in your R environment
p1 = quickplot(data=dat,
          x=O2_uM, 
          y=Depth_m, 
          colour=Season, 
          main="Saanich Inlet: Seasonal oxygen depth profile",
          xlab="Oxygen", ylab="Depth")

# You can then save the plot to your computer.
## Unless otherwise specified, the dimensions will be the same as those of the current "Plots" window in RStudio
ggsave(plot=p1, "figures/p1.png")
