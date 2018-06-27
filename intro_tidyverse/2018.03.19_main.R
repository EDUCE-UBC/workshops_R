################################################################################
# load and install packages
################################################################################
# if the named package is not installed then install it
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("lubridate")) install.packages("lubridate")

# install cowplot but DO NOT load it
if (!("cowplot" %in% installed.packages()[, 1])) install.packages("cowplot")

################################################################################
# loading data with readr
################################################################################
raw_data <- read_csv(file="Saanich_Data.csv",
                     col_names=TRUE,
                     na=c("", "NA", "NAN", "ND"))

# examine the dimensions of the data frame
dim(raw_data)

# look at the first three rows of data
head(raw_data, 3)

################################################################################
# conditional statements and logical operators
################################################################################
# extract the oxygen column from the raw data
oxygen <- raw_data$WS_O2

# subset a few practice values
(oxygen <- oxygen[c(710, 713, 715, 716, 709, 717, 718, 719)])

### conditional statements
# find the indices of the oxygen vector where the value is 204.259
oxygen == 204.259

# find the indices where the value is not 204.259
oxygen != 204.259

# find the indices where the value is greater than 204.259
oxygen > 204.259

# find the indices where the value is less than or equal to 204.259
oxygen <= 204.376

# check if 204.259 is in the oxygen vector
204.259 %in% oxygen

# find the indices where the data is NA
is.na(oxygen)

### logical operators
# find the indices where the data is NOT NA
!is.na(oxygen)

# find the indices where the value is <= 120 AND >= 20
oxygen <= 120 & oxygen >= 20

# find the indices where the value is <= 50 OR >= 150
oxygen <= 50 | oxygen >= 150

### Exercise. Given the depth data below:
# 1. Find the indices where depth is greater or equal to than 0.055
# 2. Check if the value 0.111 is in the depth data
# 3. Find where depth is less than 0.060 OR depth is greater than 0.140
depth <- raw_data$Depth[1:16]

# 1.

# 2.

# 3. 


### End exercise

# remove depth and oxygen from the R Environment
rm(depth, oxygen)

################################################################################
# data wrangling with dplyr
################################################################################
# make a copy of the raw data that we will clean
dat <- raw_data

# select the variables we will need for this workshop
dat <- select(dat, 
              Cruise, Date, Depth, 
              Temperature, Salinity, Density, 
              WS_O2, WS_NO3, WS_H2S)

# or equivalently, select all variables starting with "WS_"
dat <- select(raw_data, 
              Cruise, Date, Depth,
              Temperature, Salinity, Density,
              starts_with("WS_"))

# filter out data prior to February 2008 when a different instrument was used
dat <- filter(dat, Date >= "2008-02-01")

# filter is a very powerful function! For example, filter to retain data 
# collected in June, where Depth is one of 0.1 or 0.2, and nitrate is nonmissing
filter(dat, months(Date) == "June" & Depth %in% c(0.1, 0.2) & !is.na(WS_NO3))

### Exercise. Given the practice data (pdat) below:
# 1. select the Cruise, Date, Depth, PO4, and NO3 variables
# 2. filter the data to retain data on Cruise 72 where Depth is >= to 0.1
# your resulting pdat object should be a 8x5 data frame
pdat <- raw_data

# 1.

# 2. 

### End exercise

# rename the chemical variables so that they are easier to read and type
dat <- rename(dat, O2=WS_O2, NO3=WS_NO3, H2S=WS_H2S)

# arrange can be used to sort all rows by the value of some variable
# for example, we can arrange the data in the order of ascending NO3 data
arrange(dat[1:10, ], NO3)

# mutate Depth to convert its units from kilometres to metres
dat <- mutate(dat, Depth=Depth*1000)

### Exercise. Using the practice data (pdat) below:
# 1. Select the Date, Depth and O2 variables from pdat using select
# 2. Rename the O2 variable to Oxygen using rename
# 3. Keep August observations where Oxygen is nonmissing using filter, months, and !is.na
# 4. Transform Oxygen from micromoles/L to micrograms/L using mutate (multiply Oxygen by 32)
# 5. Run the provided ggplot() code to create a scatterplot of Oxygen vs Depth
pdat <- dat

# 1. 

# 2. 

# 3.

# 4.


ggplot(pdat, aes(x=Oxygen, y=Depth)) +
  geom_point(size=1) +
  geom_smooth(method="loess", se=FALSE) + 
  scale_y_reverse(limits=c(200, 0)) +
  labs(x=expression(O[2]*" "*(mu*g/L)),
       y="Depth (m)", 
       title="Oxygen decreases with depth and is less variable at lower depths")
### End exercise

# remove pdat from the R Environment
rm(pdat)

################################################################################
# cleaning the geochemical data with pipes
################################################################################
dat <- 
  raw_data %>%
  select(Cruise, Date, Depth, 
         Temperature, Salinity, Density, 
         WS_O2, WS_NO3, WS_H2S) %>%
  filter(Date >= "2008-02-01") %>%
  rename(O2=WS_O2, NO3=WS_NO3, H2S=WS_H2S) %>%
  mutate(Depth=Depth*1000)

### Exercise. Uncomment the initial pipe (dat %>%) below and:
# 1. Rewrite your code from the previous exercise using pipes
# 2. Pipe your data into the ggplot function

# dat %>%

### End exercise

################################################################################
# group_by and summarise
################################################################################
# calculate the mean, sd, and sample size of oxygen concentration by depth
dat %>%
  group_by(Depth) %>%
  summarise(Mean_O2=mean(O2, na.rm=TRUE),
            SD_O2=sd(O2, na.rm=TRUE),
            n=n())

### Exercise: Uncomment the initial pipe (dat %>%) below and:
# 1. Calculate median, interquartile range, and sample size of Temperature by depth

# dat %>%

### End exercise

################################################################################
# wide and long data frames -- gather and spread
################################################################################
# Our data is in the wide format
head(dat)

# gather data (except for Cruise, Date, Depth) from wide to long format
dat <- gather(dat, key="Key", value="Value", -Cruise, -Date, -Depth,
              factor_key=TRUE)  # store the variable names as an ordered factor
head(dat)
head(dat$Key)

# we can undo this by spreading the data from long to wide format
dat <- spread(dat, key="Key", value="Value")
head(dat)

################################################################################
# making figures with dplyr and ggplot2
################################################################################
# geom_point: Scatterplot of Oxygen vs Nitrate
dat %>%
  ggplot(aes(x=O2, y=NO3)) +
  geom_point(size=1)

### Exercise. Uncomment the initial pipes (dat %>%) below and:
# 1. Investigate the relationship between O2 and H2S
# 2. Investigate the relationship between NO3 and H2S

# 1.
# dat %>%

# 2.
# dat %>%

### End exercise

# colour aesthetic: Scatterplot of Oxygen vs Nitrate vs H2S
dat %>%
  filter(!is.na(H2S)) %>%
  arrange(H2S) %>%  # arrange the rows by order of ascending H2S value
  
  ggplot(aes(x=O2, y=NO3, colour=H2S)) +
  geom_point(size=2.5) 

# shape aesthetic: Depth vs Oxygen and Hydrogen Sulfide for Cruise 72
dat %>%
  select(Cruise, Depth, O2, H2S) %>%
  filter(Cruise==72) %>%
  gather(key="Chemical", value="Concentration", -Cruise, -Depth) %>%
  
  ggplot(aes(x=Concentration, y=Depth, shape=Chemical)) +
  geom_point() + 
  scale_y_reverse(limits=c(200, 0))

### Exercise: 
# 1. It may be difficult to differentiate between the different shapes in the
#    previous plot so modify the following code to add colours to the shapes as well:

dat %>%
  select(Cruise, Depth, O2, H2S) %>%
  filter(Cruise==72) %>%
  gather(key="Chemical", value="Concentration", -Cruise, -Depth) %>%
  
  ggplot(aes(x=Concentration, y=Depth, shape=Chemical)) +
  geom_point() + 
  scale_y_reverse(limits=c(200, 0))
### End exercise

# geom_line: Time series of O2 at depth 200
dat %>%
  select(Date, Depth, O2) %>%
  filter(Depth == 200 & !is.na(O2)) %>%
  gather(key="Chemical", value="O2 Concentration", -Date, -Depth) %>%
  
  ggplot(aes(x=Date, y=`O2 Concentration`)) +
  geom_point() + 
  geom_line()

# geom_line: Time series of O2 and H2S at depth 200
dat %>%
  select(Date, Depth, O2, H2S) %>%
  filter(Depth == 200 & !is.na(O2) & !is.na(H2S)) %>%
  mutate(H2S=-H2S) %>%
  gather(key="Chemical", value="Concentration", -Date, -Depth) %>%
  
  ggplot(aes(x=Date, y=Concentration, colour=Chemical)) +
  geom_point() + 
  geom_line()

# geom_histogram: Histogram to examine the distribution of Oxygen at depths < 100
dat %>%
  filter(Depth < 100) %>%
  
  ggplot(aes(x=O2)) +
  geom_histogram()

# geom_histogram: Histogram to examine the distribution of Oxygen at depths >= 100
dat %>%
  filter(Depth >= 100) %>%
  
  ggplot(aes(x=O2)) +
  geom_histogram()

### Exercise. Uncomment the initial pipe (dat %>%) below and:
# 1. Investigate the distribution of nitrate across all depths
# 2. Test out different values for the bins argument ("bins=")

# dat %>%

### End exercise

# geom_boxplot: Boxplot to investigate the distribution of 
# chemical concentrations in July at Depths >= 150
dat %>%
  select(-Temperature, -Salinity, -Density) %>%
  filter(months(Date) == "July" & 
           Depth >= 150) %>%
  gather(key="Chemical", value="Concentration", -Cruise, -Date, -Depth, factor_key=TRUE) %>%
  
  ggplot(aes(x=Chemical, y=Concentration)) +
  geom_boxplot()

# facet_wrap: Plot Depth versus Value for all six variables
dat %>% 
  gather(key="Key", value="Value", -Cruise, -Date, -Depth, factor_key=TRUE) %>%
  
  ggplot(aes(x=Value, y=Depth)) +
  geom_point(size=1) +
  scale_y_reverse(limits=c(200, 0)) +
  facet_wrap(~Key, ncol=2, dir="v", scales="free_x")

# labs: changing the labels of your ggplot
dat %>% 
  gather(key="Key", value="Value", -Cruise, -Date, -Depth, factor_key=TRUE) %>%
  
  ggplot(aes(x=Value, y=Depth)) +
  geom_point(size=1) +
  scale_y_reverse(limits=c(200, 0)) +
  facet_wrap(~Key, ncol=2, dir="v", scales="free_x") +
  labs(title="Saanich Inlet Depth Profiles (2008-2014)",
       x="",
       y="Depth (m)")

### Exercise. Uncomment the initial pipe (dat %>%) and:
# 1. Filter to data at depths of 10, 60, 100 or 200
# 2. Plot Oxygen vs Nitrate faceted by Depth without providing arguments for
#    dir or scales

# dat %>% 

### End exercise. 

################################################################################
# fine tuning your ggplots
################################################################################
# create a clean black and white theme
my_theme <- 
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# add the theme to the facetted plot and remove the x-axis label
p1 <- 
  dat %>% 
  gather(key="Key", value="Value", -Cruise, -Date, -Depth, factor_key=TRUE) %>%
  
  ggplot(aes(x=Value, y=Depth)) +
  geom_point(size=1) +
  scale_y_reverse(limits=c(200, 0)) +
  facet_wrap(~Key, ncol=2, dir="v", scales="free_x") + 
  my_theme +
  labs(x="",
       y="Depth (m)")

p1

# add the theme to the scatterplot of O2 vs NO3 vs H2S
p2 <- 
  dat %>%
  filter(!is.na(H2S)) %>%
  arrange(H2S) %>% 
  
  ggplot(aes(x=O2, y=NO3, colour=H2S)) +
  geom_point(size=2) +
  my_theme +
  labs(x="O2 in uM",
       y="NO3 in uM") +
  scale_colour_continuous(name="H2S in uM")

p2

################################################################################
# making a multi-panel figure with cowplot
################################################################################
# arrange p1 and p2 into a multi-panel figure using cowplot's plot_grid
p <- cowplot::plot_grid(p1, p2, labels=c("A", "B"), rel_widths=c(2/5, 3/5))
p

# save the final ggplot to a pdf
ggsave("saanich.pdf", p, width=10, height=6)