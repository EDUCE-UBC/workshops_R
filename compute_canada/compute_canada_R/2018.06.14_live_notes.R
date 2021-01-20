# Load packages
library(tidyverse) # Manip data
library(phyloseq) # pretty plots

# Load data
## OTU table
OTU <- read.table("data/Saanich.final.opti_mcc.unique_list.shared", sep="\t", header=TRUE, row.names=2)

## Taxonomy table
taxonomy = read.table("data/Saanich.final.opti_mcc.unique_list.0.03.cons.taxonomy", sep="\t", header=TRUE, row.names=1)

## Metadata
metadata = read.table("data/Saanich.metadata.txt", sep="\t", header=TRUE, row.names=1)

## NJ tree
load("data/NJtree.RData")

# Clean up data
## OTU table
OTU.clean = select(OTU, -label, -numOtus)
## Compare this to base R
OTU.clean2 = OTU[,-c(1,2)]

# Keep columns
OTU.test = select(OTU, label, numOtus)

# Taxonomy table
## The long way
taxonomy.clean1 = select(taxonomy, -Size) 

tax.clean2 = separate(taxonomy.clean1,Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

## The short way
taxonomy.clean = select(taxonomy, -Size) %>% 
separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

# phyloseq
OTU.clean.physeq = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)

tax.clean.physeq = tax_table(as.matrix(taxonomy.clean))

metadata.physeq = sample_data(metadata)

# phyloseq object
saanich = phyloseq(OTU.clean.physeq, tax.clean.physeq, metadata.physeq, NJtree)

saanich

# Plot tree
plot_tree(saanich)

plot_tree(saanich, color="Domain")

# Nitrospina
nitrospina = subset_taxa(saanich, Genus == "Nitrospina")

nitrospina

plot_tree(nitrospina)
plot_tree(nitrospina, color="Sample")


# Exercise 1
plot_tree(saanich, color="Phylum")

plot_tree(nitrospina, color="Sample", size="Abundance")


# Bar plots
plot_bar(saanich, fill="Phylum")

# Calculate relative abundance
saanich_percent = transform_sample_counts(saanich, function(x) 100 * x/sum(x))

plot_bar(saanich_percent, fill="Phylum")

plot_bar(saanich_percent, fill="Phylum") +
  geom_bar(aes(fill=Phylum), stat="identity")

# Correlate with metadata
nitrospina_percent = subset_taxa(saanich_percent, Genus == "Nitrospina")

# heatmap
plot_heatmap(nitrospina_percent, method=NULL, sample.label="NO3_uM", sample.order="NO3_uM")
