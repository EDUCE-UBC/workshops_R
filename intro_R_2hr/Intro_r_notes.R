# Load packages
library(ggplot2)

# Load data
read.table(file="/Users/kim/Intro_R/data.csv")
## Is the same as below because "/Users/kim/Intro_R" is my project directory
read.table(file="data.csv")
# But no formatting happens with this code

# Tell R it is in comma separated format with sep
read.table(file="data.csv", sep=",")

# Add column/variable names (header)
read.table(file="data.csv", sep=",", header=TRUE)

# Tells R to not read character data as factors (stringsAsFactors)
read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

# Save as an object in R environment
dat <- read.table(file="data.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)

######################################

# Ask R what type of data these are
class(dat)
class(dat$Season)
class(dat$Depth_m)
class(dat$O2_uM)
class(dat$Add_data)

# Maths!
## Call a specific variable witht the $
var(dat$O2_uM)

dat$O2_uM >50

dat$O2_uM / 1000

unique(dat$Season)

# Selecting using indices
## Vectors and data frames start at index = 1
dat$O2_uM[1]

dat[1,1]

# Subsetting with indices

logical.vector <- dat$O2_uM == 0
logical.vector

# Default pulls out obervations where the logical vector == TRUE
dat[logical.vector,]
dat[logical.vector, 2]

# is the same as
## You don't necessary have to save the logical vector as its own object. You can put the math/question right in the [ ]
dat[dat$O2_uM == 0, 2]

# Useful math
# & and
# | or
# < less than
# > greater than
# <= less than or equal to
# >= greater than or equal to
# != not equal to
################################

# Plotting in ggplot
quickplot(data=dat,
          x=O2_uM, y=Depth_m,
          color=Season,
          main="Saanich Inlet oxygen")

# Packages for Phylogenetics workshop

# tidyverse 
# ape
# seqinr
# vegan
# betapart
# abind
# Matrix
# cowplot

# phyloseq
source('http://bioconductor.org/biocLite.R')
biocLite("phyloseq")

# PhyloFactor
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
biocLite("Biostrings")
install.packages('devtools')
devtools::install_github('reptalex/phylofactor')
