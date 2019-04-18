################################################################################
# load packages
################################################################################
# WORKING WITH DATA
# Suite of packages for data manipulation and visualization
library(tidyverse)
# Working with raw sequences (.fasta)
library(seqinr)

# PHYLOGENETICS
# Phylogenetic analyses of microbiomes
library(phyloseq)
# Community ecology analyses
library(vegan)
# Phylogenetic factorization
library(phylofactor)
# Working with trees
library(ape)
# Additional ggplot functions for tree visualization
library(ggtree)

# OTHER
# Multi-panel figures
library(cowplot)

# FOR BDTT FUNCTIONS
# Using multidimensional arrays (part of BDTT functions)
library(abind)
# Working with matrices (part of BDTT functions)
library(Matrix)
# Beta-diversity turnover and nestedness
library(betapart)

################################################################################
# Downloading data
################################################################################
download.file(
  "https://raw.githubusercontent.com/EDUCE-UBC/workshop_data/master/Saanich.taxonomy",
  "data/Saanich.taxonomy")

download.file(
  "https://raw.githubusercontent.com/EDUCE-UBC/workshop_data/master/Saanich_OTU.shared",
  "data/Saanich_OTU.shared")

write.csv(
  read.csv("https://raw.githubusercontent.com/EDUCE-UBC/workshop_data/master/Saanich_Data_clean.csv"),
  "data/Saanich_Data_clean.csv", row.names=FALSE)

download.file(
  "https://raw.githubusercontent.com/EDUCE-UBC/workshop_data/master/Saanich_OTU_rep.alignment",
  "data/Saanich.OTU.rep.alignment")

################################################################################
# Loading and cleaning data 
################################################################################
# Taxonomic identity of each sequence
taxonomy <- read_tsv("data/Saanich.taxonomy") %>% 
  # Separate taxa names into columns
  separate(Taxonomy,
           into=c("domain","phylum","class","order","family","genus","species"),
           sep=";") %>% 
  # Remove unused columns
  select(-Size)

# Counts of sequences in each sample
OTU <- read_tsv("data/Saanich_OTU.shared") %>% 
  # Rename sample variable
  mutate(sample=Group) %>% 
  # Remove unused columns
  select(-label, -numOtus, -Group)

# Geochemical data
metadata <- read_csv("data/Saanich_Data_clean.csv") %>% 
  # Filter to only Cruise 72
  filter(Cruise == 72) %>% 
  # Create sample names similar to OTU table
  mutate(sample=ifelse(Depth_m <100, 
                       paste("Saanich_0", Depth_m, sep=""),
                       paste("Saanich_", Depth_m, sep="")))

# Representative sequences
# Reading function from seqinr package
alignment <- read.fasta("data/Saanich.OTU.rep.alignment")
# Change the sequences names to simple OTU numbers
names(alignment)=taxonomy$OTU

################################################################################
# Importing custom R functions
################################################################################
download.file(
"https://raw.githubusercontent.com/EDUCE-UBC/workshops_R/master/phylogenetics/scripts/BDTT_functions.R",
"scripts/BDTT_functions.R"
)

source("scripts/BDTT_functions.R")

################################################################################
# Building a phylogenetic tree (FastTree)
################################################################################

################################################################################
# Tree constraints
################################################################################
# Define contraints
tax.constrain <- taxonomy %>% 
# Remove unclassified domains
filter(domain != "unknown")

align.constrain <- alignment[tax.constrain$OTU]

# set the domain constraints as 1 for Bacteria and 0 for Archaea in the taxonomy table.
tax.constrain <- tax.constrain %>% 
# Create constraint variable
mutate(Bacteria = ifelse(domain == "Bacteria", 1, 0))

# define the phyla constraints, giving each phylum its own column
# List all unique phyla names
phyla.list <- taxonomy %>%
select(phylum) %>%
# Remove unclassified phyla
filter(!grepl("unclassified", phylum)) 
# Keep only unique phyla in a list
phyla.list <- unique(phyla.list$phylum)

# For each phylum, create a column in the constrained taxonomy table with
# a value of 1 for the OTU being that phylum and 0 for any other phyla
for (i in phyla.list){
tax.constrain[[as.character(i)]][tax.constrain$phylum==i]=1
tax.constrain[[as.character(i)]][!tax.constrain$phylum==i]=0
}

# Keep only the 0/1 constraint columns
tax.constrain2 <- tax.constrain %>% 
select(-(OTU:species))

# View result
tax.constrain2

### Save constraints
# Create empty list to hold data
sequences=list()

# For each OTU, list its 0/1 values from the constraints table
for (i in 1:nrow(tax.constrain2)){
sequences[[i]]=tax.constrain2[i,]
}

# save as a fasta.
# Print these 0/1 lists as the "sequence" for each OTU in fasta format
write.fasta(sequences, names=tax.constrain$OTU,
file.out="data/Saanich_tax_constrain.fasta", 
open = "w", nbchar = 60, as.string = FALSE)


# save the modified alignment file where we've trimmed out the unclassified domain.
write.fasta(alignment, names(alignment),
            "data/Saanich_OTU_rep_mod.alignment",
            open = "w", nbchar = 60, as.string = FALSE)

################################################################################
# FastTree (in terminal not R)
################################################################################
# Replace with your project directory
cd /Users/kim/GitHub/workshops_R/phylogenetics/
# Make results folder
mkdir -p results/

# View FastTree help page
/Users/kim/Applications/FastTree/FastTree -help

# run FastTree on our data.
/Users/kim/Applications/FastTree/FastTree \
-gtr -cat 20 \
-constraints data/Saanich_tax_constrain.fasta \
-nt data/Saanich_OTU_rep_mod.alignment \
> results/Saanich_FastTree_constrain

# construct the tree WITHOUT topological constraints
/Users/kim/Applications/FastTree/FastTree \
-gtr -cat 20 \
-nt data/Saanich_OTU_rep_mod.alignment \
> results/Saanich_FastTree

# Download results if needed 
download.file(
"https://raw.githubusercontent.com/EDUCE-UBC/workshops_R/master/phylogenetics/results/Saanich_FastTree",
"results/Saanich_FastTree")

download.file(
"https://raw.githubusercontent.com/EDUCE-UBC/workshops_R/master/phylogenetics/results/Saanich_FastTree_constrain",
"results/Saanich_FastTree_constrain")

################################################################################
# Exercise: FastTree
################################################################################
#hile the above rerunning, compare the two available nucleotide models in FastTree. Discuss the pros and cons of Jukes-Cantor (JT) versus Generalized Time-Reversible (GTR) with your table. You may find this Wikipedia article helpful (https://en.wikipedia.org/wiki/Models_of_DNA_evolution#Most_common_models_of_DNA_evolution).

################################################################################
# Visualizing trees
################################################################################
# read.tree from the ape package
tree <- read.tree('results/Saanich_FastTree')
treeC <- read.tree('results/Saanich_FastTree_constrain')

################################################################################
# phyloseq formatting
################################################################################
# OTU table
OTU.physeq = OTU %>% 
# set sample name as row names
column_to_rownames(var = "sample") %>% 
# Format to matrix
as.matrix() %>% 
# Format to phyloseq OTU table
otu_table(taxa_are_rows=FALSE)

# Taxonomy
tax.physeq = taxonomy %>% 
# set OTU # as row names
column_to_rownames(var = "OTU") %>% 
# Convert to matrix
as.matrix() %>% 
# Convert to phyloseq tax table
tax_table()

# Metadata
metadata.physeq = metadata %>% 
# Copy sample column
mutate(sample2 = sample) %>% 
# set sample names as row names
column_to_rownames(var = "sample2") %>% 
# Convert to phyloseq sample data
sample_data()

# Trees
tree.physeq=phy_tree(tree)
treeC.physeq=phy_tree(treeC)

# assemble them into phyloseq objects
saanich = phyloseq(OTU.physeq, tax.physeq, 
metadata.physeq, tree.physeq) 

saanichC = phyloseq(OTU.physeq, tax.physeq, 
metadata.physeq, treeC.physeq) 

# See objects
saanich
saanichC

################################################################################
# Basic trees
################################################################################
plot_tree(saanich, "treeonly") +
  # Make circular
  coord_polar(theta = "y") +
  # Plot title
  labs(title="Unconstrained tree")

plot_tree(saanichC, "treeonly") +
  # Make circular
  coord_polar(theta = "y") +
  # Plot title
  labs(title="Constrained tree")

################################################################################
# Color by taxa
################################################################################
#Domain
plot_tree(saanich, color="domain") +
# Make circular
coord_polar(theta = "y") +
# Plot title
labs(title="Unconstrained tree") +
# Move legend to bottom
theme(legend.position="bottom")

plot_tree(saanichC, color="domain") +
# Make circular
coord_polar(theta = "y") +
# Plot title
labs(title="Constrained tree") +
# Move legend to bottom
theme(legend.position="bottom")

#Phyla
plot_tree(saanich, color="phylum") +
# Make circular
coord_polar(theta = "y") +
# Plot title
labs(title="Unconstrained tree") +
# Move legend to bottom
theme(legend.position="bottom") +
guides(col = guide_legend(ncol = 2))

plot_tree(saanichC, color="phylum") +
# Make circular
coord_polar(theta = "y") +
# Plot title
labs(title="Constrained tree") +
# Move legend to bottom
theme(legend.position="bottom") +
guides(col = guide_legend(ncol = 2))

################################################################################
# Exercise: Color taxa trees
################################################################################

# 1. Modify the above trees to color by phylum, instead of domain. Your output should be similar to below.
# 2. Compare these two trees. What do you observe? Which do you think better represents these data?
  
################################################################################
# Re-root the tree
################################################################################
# List all OTUs that are Archaea
archaea.list <- taxonomy %>% 
  filter(domain=="Archaea") %>% 
  select(OTU) %>% 
  as.list()

MRCAnode <- getMRCA(phy = treeC, tip = archaea.list$OTU)

#re-root the tree on this branch
treeC.root <- root(phy = treeC, node = MRCAnode, resolve.root = TRUE)

#save the corresponding tree to the disk
write.tree(treeC.root, 'results/Saanich_FastTree_constrain_root')

#in a new phyloseq object.
treeCR.physeq <- phy_tree(treeC.root)

saanichCR <- phyloseq(OTU.physeq, tax.physeq, 
                      metadata.physeq, treeCR.physeq) 

#Plotting this tree and comparing to the unrooted tree
plot_tree(saanichCR, color="phylum") +
  # Make circular
  coord_polar(theta = "y") +
  # Plot title
  labs(title="Constrained, rooted tree") +
  # Move legend to bottom
  theme(legend.position="bottom") +
  guides(col = guide_legend(ncol = 2))

################################################################################
# Calculating beta-diversity
################################################################################
# Bray-Curtis
BC <- vegdist(otu_table(saanichCR), method = "bray")
# Jaccard
Jac  <- vegdist(otu_table(saanichCR), method = "jaccard")
# UniFrac
UF <- UniFrac(saanichCR, weighted=FALSE)
# Weighted UniFrac
wUF <- UniFrac(saanichCR, weighted=TRUE)

#Look at result example
BC

################################################################################
# Principle coordinate analysis (PCoA/PCA)
################################################################################
BC.pcoa = ordinate(saanichCR, "PCoA", distance=BC)
Jac.pcoa = ordinate(saanichCR, "PCoA", distance=Jac)
UF.pcoa = ordinate(saanichCR, "PCoA", distance=UF)
wUF.pcoa = ordinate(saanichCR, "PCoA", distance=wUF)

# Set plot theme for all plots
theme_set(theme_classic())

# Plot each PCoA colored by sample
plot_ordination(saanichCR, BC.pcoa, color="sample") +
  #Make points larger
  geom_point(size = 4) +
  #Remove legend
  theme(legend.position="none") +
  # Add title
  labs(title="Bray-Curtis")

plot_ordination(saanichCR, Jac.pcoa, color="sample") +
  #Make points larger
  geom_point(size = 4) +
  #Remove legend
  theme(legend.position="none") +
  # Add title
  labs(title="Jaccard")

plot_ordination(saanichCR, UF.pcoa, color="sample") +
  #Make points larger
  geom_point(size = 4) +
  #Remove legend
  theme(legend.position="none") +
  # Add title
  labs(title="UniFrac")

plot_ordination(saanichCR, wUF.pcoa, color="sample")+
  #Make points larger
  geom_point(size = 4) +
  #Move legend
  theme(legend.position="bottom") +
  # Add title
  labs(title="weighted UniFrac")

################################################################################
# Non-metric multidimensional scaling (nMDS)
################################################################################
BC.nmds = ordinate(saanichCR, "NMDS", distance=BC)
Jac.nmds = ordinate(saanichCR, "NMDS", distance=Jac)
UF.nmds = ordinate(saanichCR, "NMDS", distance=UF)
wUF.nmds = ordinate(saanichCR, "NMDS", distance=wUF)

# Set plot theme for all plots
theme_set(theme_classic())

# Plot each PCoA colored by sample
plot_ordination(saanichCR, BC.nmds, color="sample") +
  #Make points larger
  geom_point(size = 4) +
  #Remove legend
  theme(legend.position="none") +
  # Add title
  labs(title="Bray-Curtis")

plot_ordination(saanichCR, Jac.nmds, color="sample") +
  #Make points larger
  geom_point(size = 4) +
  #Remove legend
  theme(legend.position="none") +
  # Add title
  labs(title="Jaccard")

plot_ordination(saanichCR, UF.nmds, color="sample") +
  #Make points larger
  geom_point(size = 4) +
  #Remove legend
  theme(legend.position="none") +
  # Add title
  labs(title="UniFrac")

plot_ordination(saanichCR, wUF.nmds, color="sample")+
  #Make points larger
  geom_point(size = 4) +
  #Move legend
  theme(legend.position="bottom") +
  # Add title
  labs(title="weighted UniFrac")

################################################################################
# Exercise: PCoA vs. nMDS
################################################################################
# Consider the two beta-diversity plotting methods above. What are the pros and cons of these methods? Which do you think best represents these data?

################################################################################
# Statistically testing beta-diversity
################################################################################
adonis(BC ~ Depth_m, data=data.frame(sample_data(saanichCR)))
adonis(Jac ~ Depth_m, data=data.frame(sample_data(saanichCR)))
adonis(UF ~ Depth_m, data=data.frame(sample_data(saanichCR)))
adonis(wUF ~ Depth_m, data=data.frame(sample_data(saanichCR)))

adonis(BC ~ O2_uM * NO3_uM * H2S_uM,
data=data.frame(sample_data(saanichCR)))

################################################################################
# Exercise: Randomness and variable order
################################################################################
# 1. Rerun the the last PERMANOVA several times. Do the results change? Why might this be the case and what does this mean for borderline p-value (*e.g.* P near 0.05)?
  
# 2. Next, alter the variable order in this PERMANOVA. Does this change the outcome? Why might this be the case?
  
################################################################################
# Screening by phylogenetic scale
################################################################################
Hnodes <- getHnodes(treeC.root)

hist(Hnodes, n=150)
hist(Hnodes, n=150, xlim=c(0,.5))

#Create incremental list of values from 0 to 0.3, going up by 0.025 each time
slices <- c(seq(from=0, to=0.3, by=0.025)) 

#transpose our OTU table as a matrix.
OTU.mat <- OTU %>% 
# set sample name as row names
column_to_rownames(var = "sample") %>% 
# Transpose
t() %>% 
#Format to matrix
as.matrix()

#run the analysis for Jaccard using our custom functions
multi.Jac <- BDTT(similarity_slices = slices,
tree = treeC.root, sampleOTUs = OTU.mat,
metric = "jac")

#run the same analysis for Bray-Curtis.
multi.BC <- BDTT(similarity_slices = slices,
tree = treeC.root, sampleOTUs = OTU.mat,
metric = "bc")

#Save results
saveRDS(multi.Jac, "results/multi_Jac.RDS")  
saveRDS(multi.BC, "results/multi_BC.RDS")  

#Or you can download and load the results here.
download.file(
"https://github.com/EDUCE-UBC/workshops_R/blob/master/phylogenetics/results/multi_Jac.RDS?raw=true",
"results/multi_Jac.RDS")

multi.Jac <- readRDS("results/multi_Jac.RDS")

###

download.file(
"https://github.com/EDUCE-UBC/workshops_R/blob/master/phylogenetics/results/multi_BC.RDS?raw=true",
"results/multi_BC.RDS")

multi.BC <- readRDS("results/multi_BC.RDS")

################################################################################
# Exploring array results
################################################################################
class(multi.Jac)
dim(multi.Jac)
multi.Jac[1,,]

################################################################################
# Statistical links to metadata
################################################################################
# list variables of interest
predictors <- c("O2_uM", "NO3_uM", "H2S_uM", "Depth_m")
# Create data frame of all slice, variable, metric combinations
StatsRes <- expand.grid(similarity_slices=as.character(slices),
                        predictors=predictors,
                        metric=c("Jac","BC"))
# Add blank columns to hold results
StatsRes[["F.Model"]] = StatsRes[["R2"]] = StatsRes[["Pr(>F)"]]=NA

#View first few rows of data frame
head(StatsRes)

#"fill" it with our results from PERMANOVA of all models constructed in a loop. 
# For each slice
for (i in as.character(slices)){
  # For each variable of interest
  for (j in predictors) {
    # Calculate PERMANOVA for Jaccard
    res <- unlist(adonis(
      formula = multi.Jac[i,,] ~
        data.frame(sample_data(saanichCR))[,j])$aov.tab[1,c(4,5,6)])
    # Add results to table
    StatsRes[(StatsRes$metric=="Jac") & 
               (StatsRes$predictors==j) & 
               (StatsRes$similarity_slices==i), 4:6] = res
    
    # Calculate PERMANOVA for Bray-Curtis
    res <- unlist(adonis(
      formula = multi.BC[i,,] ~
        data.frame(sample_data(saanichCR))[,j])$aov.tab[1,c(4,5,6)])
    # Add results to table
    StatsRes[(StatsRes$metric=="BC") &
               (StatsRes$predictors==j) &
               (StatsRes$similarity_slices==i),4:6] = res
  }
}

# plot the fit profiles using R^2^ along our phylogenetic time scale.
ggplot(StatsRes,
       aes(y=R2, x=similarity_slices, color=predictors,
           group=predictors)) +
  geom_point(size=4) +
  geom_line(size=1) +
  facet_wrap(~metric) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
# Exercise: Interpreting phylogenetic scale
################################################################################
# Consider the above plot. What is the biological relevance to oxygen and depth behaving differently? How does this inform our questions about nutrient gradients in Saanich Inlet? 
  
################################################################################
# Phylogenetic compositional factor analysis (PhyloFactor)
################################################################################
################################################################################
# Data formatting
################################################################################
OTU.mat.sub <- OTU.mat[row.names(OTU.mat) %in% treeC.root$tip.label,]

taxonomy.sub <- taxonomy %>% 
filter(OTU %in% treeC.root$tip.label)

metadata.sub <- metadata %>% 
filter(sample %in% colnames(OTU.mat))

taxonomy.raw <- read_tsv("data/Saanich.taxonomy") %>% 
# subset to only tips in the tree
filter(OTU %in% treeC.root$tip.label) %>% 
#rename OTU column
mutate(OTU_ID = OTU) %>% 
# Add rownames
column_to_rownames(var = "OTU") %>% 
# Remove unused columns and reorder
select(OTU_ID, Taxonomy)

################################################################################
# Detecting microbial variables in the tree
################################################################################
### Simulated data
# Set random seed for reproducibility
set.seed(1)
CellDiam <- rlnorm(length(treeC.root$tip.label))

# Grab a Boolean vector corresponding to the Arcobacter and Actinobacteria tip labels of the Tree object
arco <- treeC.root$tip.label %in%
  taxonomy.sub[taxonomy.sub$genus == 'Arcobacter',]$OTU

actino <- treeC.root$tip.label %in%
  taxonomy.sub[taxonomy.sub$phylum == 'Actinobacteria',]$OTU

# Change those tips to have new values for CellDiam drawn from a different distribution
CellDiam[arco] <- rlnorm(sum(arco))*4
CellDiam[actino] <- rlnorm(sum(actino))/4

# Take the log of CellDiam variable for analysis
logCellDiam <- log(CellDiam)

################################################################################
# PhyloFactor using two-sample tests
################################################################################
pf_twoSample <- twoSampleFactor(logCellDiam, treeC.root, nfactors=5)

pf_twoSample

#Compare to simulated groups
sum(actino)
sum(arco)

#P-values
pf_twoSample$pvals

#plot
cellDiam.tree <- pf.tree(pf_twoSample)

cellDiam.tree$ggplot
#plot colors
cellDiam.tree$legend

################################################################################
# Detecting metadata variables
################################################################################
################################################################################
# PhyloFactor ranked by effect size
################################################################################
pf_depth_F <- PhyloFactor(OTU.mat.sub, treeC.root, 
                          X = metadata.sub, 
                          frmla = Data ~ Depth_m, 
                          nfactors = 5, ncores=2, 
                          stop.early = TRUE, choice='F')

pf_depth_F

#Explore results
pf.taxa(pf_depth_F, taxonomy.raw, factor=1)$group1

ggplot(summary(pf_depth_F, factor=1)$data, 
       aes(x = X.Depth_m, y = Data)) + 
  geom_point(size=4) + 
  geom_smooth(method='lm') +
  labs(title="Factor 1")

pf.taxa(pf_depth_F, taxonomy.raw, factor=2)$group1

ggplot(summary(pf_depth_F, factor=2)$data,
       aes(x = X.Depth_m, y = Data)) + 
  geom_point(size=4) + 
  geom_smooth(method='lm') +
  labs(title="Factor 2")

################################################################################
# PhyloFactor ranked by variance explained
################################################################################
pf_depth_var <- PhyloFactor(OTU.mat.sub, treeC.root,
                            X = metadata.sub,
                            frmla = Data ~ Depth_m,
                            nfactors = 5, ncores=2, 
                            stop.early = TRUE, choice='var')

pf_depth_var

pf.taxa(pf_depth_var, taxonomy.raw, factor=2)$group1

ggplot(summary(pf_depth_var, factor=2)$data, 
       aes(x = X.Depth_m, y = Data)) + 
  geom_point(size=4) + 
  geom_smooth(method='lm') +
  labs(title="Factor 2")

################################################################################
# PhyloFactor across multiple variables
################################################################################
pf_O2_var <- PhyloFactor(OTU.mat.sub, treeC.root, 
X = metadata.sub, 
frmla = Data ~ O2_uM, 
nfactors = 3, ncores=2, 
choice='var')

pf_NO3_var <- PhyloFactor(OTU.mat.sub, treeC.root, 
X = metadata.sub, 
frmla = Data ~ NO3_uM, 
nfactors = 3, ncores=2, 
choice='var')

##For oxygen:
# Draw the PhyloFactors on the phylogeny
O2_tree <- pf.tree(pf_O2_var)
# Save tree plot
O2_tree_plot <- O2_tree$ggplot +
labs(title="Oxygen")

# Plot each of top 3 factors
O2_plot1 <- ggplot(summary(pf_O2_var, factor=1)$data, 
aes(x = X.O2_uM, y = Data)) +
# Add linear best fit
geom_smooth(method='lm', color = O2_tree$legend$colors[1],
se=FALSE) +
# Add data points
geom_point(size=4, color = O2_tree$legend$colors[1]) + 
# Modify text labels
labs(title = "Factor 1", x="", y="ILR abundance")

O2_plot2 <- ggplot(summary(pf_O2_var, factor=2)$data, 
aes(x = X.O2_uM, y = Data)) +
geom_smooth(method='lm', color = O2_tree$legend$colors[2],
se=FALSE) +
geom_point(size=4, color = O2_tree$legend$colors[2]) + 
labs(title = "Factor 2", x="", y="ILR abundance")

O2_plot3 <- ggplot(summary(pf_O2_var, factor=3)$data, 
aes(x = X.O2_uM, y = Data)) +
geom_smooth(method='lm', color = O2_tree$legend$colors[3],
se=FALSE) +
geom_point(size=4, color = O2_tree$legend$colors[3]) + 
labs(title = "Factor 3", x="Oxygen (uM)", y="ILR abundance")

# Combine 3 factors into 1 plot
O2_plots <- plot_grid(O2_plot1, O2_plot2, O2_plot3,
labels = c("B","C","D"), ncol = 1, align = 'v')


##For nitrate:
# Draw the PhyloFactors on the phylogeny
NO3_tree <- pf.tree(pf_NO3_var)
# Save tree plot
NO3_tree_plot <- NO3_tree$ggplot +
labs(title="Nitrate")

# Plot each of top 3 factors
NO3_plot1 <- ggplot(summary(pf_NO3_var, factor=1)$data, 
aes(x = X.NO3_uM, y = Data)) +
geom_smooth(method='lm', color = NO3_tree$legend$colors[1],
se=FALSE) +
geom_point(size=4, color = NO3_tree$legend$colors[1]) + 
labs(title = "Factor 1", x="", y="ILR abundance")

NO3_plot2 <- ggplot(summary(pf_NO3_var, factor=2)$data, 
aes(x = X.NO3_uM, y = Data)) +
geom_smooth(method='lm', color = NO3_tree$legend$colors[2],
se=FALSE) +
geom_point(size=4, color = NO3_tree$legend$colors[2]) + 
labs(title = "Factor 2", x="", y="ILR abundance")

NO3_plot3 <- ggplot(summary(pf_NO3_var, factor=3)$data, 
aes(x = X.NO3_uM, y = Data)) +
geom_smooth(method='lm', color = NO3_tree$legend$colors[3],
se=FALSE) +
geom_point(size=4, color = NO3_tree$legend$colors[3]) + 
labs(title = "Factor 3", x="Nitrate (uM)", y="ILR abundance")

# Combine 3 factors into 1 plot
NO3_plots <- plot_grid(NO3_plot1, NO3_plot2, NO3_plot3,
labels = c("F","G","H"), ncol = 1, align = 'v')


##Combine into 1 plot.
plot_grid(O2_tree_plot, O2_plots,
NO3_tree_plot, NO3_plots,
labels = c("A","","E",""), ncol = 2, align = 'hv',
rel_widths = c(1.5, 1))

################################################################################
# Exercise: Interpreting PhyloFactor
################################################################################
# Consider the results summarized in our last figure. What can you conclude? (*Hint*: Are you confident in the oxygen results given the distribution of data points along the fitted line?)

################################################################################
# Phylogenetic heatmaps
################################################################################
# Actual data
pf.heatmap(pf_O2_var, factors=1:3, 
column.order = order(metadata.sub$O2_uM),
low='purple', high='yellow')

# Predictions
preds <- predict(pf_O2_var)

pf.heatmap(pf_O2_var, factors=1:3, 
Data=preds,
column.order = order(metadata.sub$O2_uM),
low='purple', high='yellow')