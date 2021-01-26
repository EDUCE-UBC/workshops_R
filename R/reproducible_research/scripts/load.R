# load data
dat <- read.csv("https://raw.githubusercontent.com/EDUCE-UBC/educer/main/data-raw/data_intro_ws.csv")

write.csv(dat, "data.csv")
