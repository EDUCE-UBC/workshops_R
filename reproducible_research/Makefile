all: O2_plot.png

data.csv: load.R
	Rscript load.R

O2_plot.png: data.csv plot.R
	Rscript plot.R
