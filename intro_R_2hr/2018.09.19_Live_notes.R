# Install ggplot2 package
#install.packages("ggplot2")

# Load the package
library(ggplot2)

# Load data
# function(data=, header=TRUE)
read.table(data="data.csv")
# Error =  didn't run
# Warning = it ran but something weird happened
# Nothing = success!

#Simple, tells R nothing about the data
read.table(file="data.csv")

# Comma separated
read.table(file="data.csv", sep=",")

# Column names
read.table(file="data.csv", sep=",", header=TRUE)

# String variables
# Can use lots of spaces
read.table(file = "data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

#Save to environment
dat = read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
#Identical function
dat <- read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

# Ask R what data type it is
class(dat)
class(dat$O2_uM)
class(dat$Add_data)
class(dat$Season) #The is the variable that is incorrectly formatted if stringsasfactors=TRUE

# Why we had to use stringsasfactors=FALSE
dat2 = read.table(file="data.csv", sep=",", header=TRUE)
class(dat2$Add_data)
class(dat2$Season)


# to rename a column
colnames(dat$..XSeason) = "Season"

# Work with data
## calculate variance
var(dat$O2_uM)
var(dat$Season)

# capital O vs zero 0

# Tell me when oxygen is greater than 50
dat$O2_uM > 50

# convert a variable to a new variable
## Depth is in m. What if you want km?
dat$Depth_m / 1000
#  save it for later use
Depth_km = dat$Depth_m / 1000
# save as a new column
dat$"Depth_km" = dat$Depth_m / 1000

# save to computer
write.table(dat, "~/Desktop/Intro_R/dat_new.csv", sep=",", row.names=FALSE)

#Indexing
# to find pieces of data
dat$O2_uM[5]

#from a table/data frame
#data[rows,columns]
dat[5,3]

##subsetting data
#part of tidyverse
filter(Add_data==FALSE)
#with indexing
dat[dat$Add_data==FALSE,4]

#more readable with indexing
logical.vector = dat$Add_data==FALSE
dat[logical.vector,4]

# if data is numeric, 1.00
# x = 1 is TRUE
# x == 1 is FALSE (exactly equals)
# x != 1 (not equals)

#making a plot
#color = colour
quickplot(data=dat, x=O2_uM, y=Depth_m, color=Season, main="Saanich Inlet")

#base R plot
plot(x=dat$O2_uM, y=dat$Depth_m, main="Saanich Inlet")

#making ggplot better
quickplot(data=dat, x=O2_uM, y=Depth_m, color=Season, main="Saanich Inlet") +
  theme_bw()

#even more themes is ggthemes