##################################################################
######### Analisis de Core Microbiome - Lagunas Costeras #########
##################################################################



############ Elaboró JORGE ALEXANDER ROJAS-VARGAS ######################


# Libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(phyloseq)
library(tidyverse)
library(readr)
library(limma)
library(ComplexHeatmap)
library(ecodist)
library(patchwork)
library(svglite)
library(reshape)
library(microbiome)
 library(microbiomeutilities)
library(gghalves)
#library(fantaxtic)
library(microViz)
library(compositions)
library(RColorBrewer)
library(viridis)
library(edgeR)
library(microbiomeMarker)
library(microbial)
library(UpSetR)
library(lmPerm)
library(ANCOMBC)
library(ampvis2)
library(mixOmics)
library(dendextend)
library(dunn.test)
library(PCAmixdata)
library(RVAideMemoire)
library(permute)
library(jmuOutlier)
library(GGally)
library(CCA)
library(labdsv)
library(FSA)
library(reshape2)
library(tibble)

# Location of work directory
setwd("/Users/jorgerv/Documents/lagunas_costeras")


###########################
##### ABUNDANCE PLOTS #####
###########################

## Run Pavian to visualize the Kraken reports and download count and abundance tables
library(pavian)
pavian::runApp(port=5000)

## Read the count table

#count_class <- read.csv("C_arc_bacteria_count_class.csv", row.names = 1)

# Read the abundance tables
abundance_table_class <- read.csv("C_bacteria_abd_class.csv")
abundance_table_family <- read.csv("C_bacteria_abd_family.csv")
abundance_table_genus <- read.csv("C_bacteria_abd_genus.csv")
abundance_table_species <- read.csv("C_bacteria_abd_species.csv")
abundance_table_annotation <- read.csv("C_bacteria_abd_annotation.csv")

# Filter abundance tables

# Se eliminarán las abundancias menores a 1%. Se puede cambiar este valor si se desea por x > 0.5 o 0
abundance_table_class <- column_to_rownames(abundance_table_class, var = "name")
abundance_table_class_stock <- abundance_table_class[apply(abundance_table_class, 1, function(x) any(x >= 1, na.rm = TRUE)), ]
# Filtrando las taxas cuya sumatoria de abundancia en todas las muestras sea mayor a 0.5%
# abundance_table_class <- abundance_table_class[rowSums(abundance_table_class > 0.5, na.rm = TRUE) > 0, ]
# Filtrando las taxa que tenga al mínimo el 1% en TODAS las muestras
# abundance_table_class <- abundance_table_class[apply(abundance_table_class, 1, function(x) all(x >= 1 | is.na(x))), ]

abundance_table_family <- column_to_rownames(abundance_table_family, var = "name")
abundance_table_family_stock <- abundance_table_family[apply(abundance_table_family, 1, function(x) any(x >= 1, na.rm = TRUE)), ]
# abundance_table_family <- abundance_table_family[rowSums(abundance_table_family > 0.5, na.rm = TRUE) > 0, ]
# abundance_table_family <- abundance_table_family[apply(abundance_table_family, 1, function(x) all(x >= 1 | is.na(x))), ]

abundance_table_genus <- column_to_rownames(abundance_table_genus, var = "name")
abundance_table_genus_stock <- abundance_table_genus[apply(abundance_table_genus, 1, function(x) any(x >= 0.5, na.rm = TRUE)), ]
# abundance_table_genus <- abundance_table_genus[rowSums(abundance_table_genus > 0.5, na.rm = TRUE) > 0, ]
# abundance_table_genus <- abundance_table_genus[apply(abundance_table_genus, 1, function(x) all(x >= 1 | is.na(x))), ]

abundance_table_species <- column_to_rownames(abundance_table_species, var = "name")
abundance_table_species_stock <- abundance_table_species[apply(abundance_table_species, 1, function(x) any(x >= 0.5, na.rm = TRUE)), ]
# abundance_table_species <- abundance_table_species[rowSums(abundance_table_species > 0.5, na.rm = TRUE) > 0, ]
# abundance_table_species <- abundance_table_species[apply(abundance_table_species, 1, function(x) all(x >= 1 | is.na(x))), ]

abundance_table_annotation <- column_to_rownames(abundance_table_annotation, var = "name")


# Read metadata
metadata <- read.csv("metadata_quintin.csv")

require(ggplot2)
require(reshape2)

## PHYLOSEQ files

# Genus
new_column_names <- metadata$ID_1
colnames(abundance_table_genus) <- new_column_names
taxa_table <- row.names(abundance_table_genus)
setdiff(rownames(abundance_table_genus),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(abundance_table_genus, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$ID_1
ps_abund_genus = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_abund_genus <- prune_taxa(taxa_sums(ps_abund_genus)>0, ps_abund_genus)
print(ps_abund_genus)
# Log transform
a_genus_log <- abundance_table_genus_stock-0.001
a_genus_log[a_genus_log < 0] <- 0
a_genus_log <- log((a_genus_log+0.001)/0.001)

# Log transform
a_species_log <- abundance_table_species_stock-0.001
a_species_log[a_species_log < 0] <- 0
a_species_log <- log((a_species_log+0.001)/0.001)

###### CLASS LEVEL ######

# Añado fila de "Other" para recuperar la taxa filtrada
column_sums <- colSums(abundance_table_class_stock)
diff_values <- 100 - column_sums
other_row <- as.data.frame(t(diff_values))
rownames(other_row) <- "Other"
#abundance_table_class_stock <- rbind(abundance_table_class_stock, other_row)
abundance_table_class_stock <- rbind(other_row, abundance_table_class_stock)
new_column_names <- metadata$ID_1
colnames(abundance_table_class_stock) <- new_column_names

prop<-prop.table(data.matrix(abundance_table_class_stock), 2)
dat_m <- melt(prop)
colnames(dat_m)<-c("Class", "Sample", "Abundance")
dat_m$Abundance <- dat_m$Abundance * 100

# Agrego metadatos
date <- metadata$date[match(dat_m$Sample, metadata$ID_1)]
sector <- metadata$Sector[match(dat_m$Sample, metadata$ID_1)]
place <- metadata$place[match(dat_m$Sample, metadata$ID_1)]

dat_m$date <- date
dat_m$sector <- sector
dat_m$place <- place

# Organizar los datos para que queden en orden de menor a mayor según abundancia de X taxa
#X = "Vibrio"
#dat_X <- subset(dat_m, Genus == X)
#sum_abd_X <- aggregate(Abundance ~ Sample, data = dat_X, FUN = sum)
#abundance_X <- sum_abd_X$Abundance
#ordered_Sample <- sum_abd_X$Sample[orde(abundance_X)]
#dat_m$Sample <- factor(dat_m$Sample, levels = ordered_Sample)

# Plot dividido según la fecha
my.colors = color.spectral(19)
n <- 19
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
my.colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#my.colors = createPalette(19,  c("#ff0000", "#00ff00", "#0000ff"))
#swatch(my.colors)

svglite("Abundance_bacteria_class_quintin.svg", width=10, height=8)
ggplot(dat_m, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") + xlab("Sample") + ylab("Abundance") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold")
  ) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  ggtitle("Bacteria class level") +
  scale_fill_manual(
    values = my.colors,
    guide = guide_legend(nrow = 29)
  ) +
  #facet_grid(~ date, scales = "free", space = "free")
  #facet_grid(~ place, scales = "free", space = "free")
  facet_grid(~ sector, scales = "free", space = "free")
dev.off()



# Extraer los datos para la clase "Flaviobacteriia"
class_abundance <- as.numeric(abundance_table_class_stock["Flavobacteriia", ])
class_abundance <- as.numeric(abundance_table_class_stock["Planctomycetia", ])

# Agregar las abundancias a la tabla de metadatos
metadata$Class <- class_abundance

# Filtrar por los hábitats de interés
filtered_data <- metadata[metadata$Habitat %in% c("Seagrass", "Bare"), ]

# Realizar la prueba de Kruskal-Wallis
kruskal_test <- kruskal.test(Class ~ Habitat, data = filtered_data)

# Filtrar por los hábitats de interés
filtered_data <- metadata[metadata$Sector %in% c("A_Inlet", "B_Transition"), ]

# Realizar la prueba de Kruskal-Wallis
kruskal_test <- kruskal.test(Class ~ Sector, data = filtered_data)

# Mostrar los resultados
print(kruskal_test)




###### FAMILY LEVEL ######

# Añado fila de "Other" para recuperar la taxa filtrada
column_sums <- colSums(abundance_table_family_stock)
diff_values <- 100 - column_sums
other_row <- as.data.frame(t(diff_values))
rownames(other_row) <- "Other"
abundance_table_family_stock <- rbind(abundance_table_family_stock, other_row)

prop<-prop.table(data.matrix(abundance_table_family_stock), 2)
dat_m <- melt(prop)
colnames(dat_m)<-c("Family", "Sample", "Abundance")

# Agrego metadatos
date <- metadata$date[match(dat_m$Sample, metadata$sample)]
zone <- metadata$zone[match(dat_m$Sample, metadata$sample)]
place <- metadata$place[match(dat_m$Sample, metadata$sample)]

dat_m$date <- date
dat_m$zone <- zone
dat_m$place <- place

# Organizar los datos para que queden en orden de menor a mayor según abundancia de X taxa
#X = "Vibrio"
#dat_X <- subset(dat_m, Genus == X)
#sum_abd_X <- aggregate(Abundance ~ Sample, data = dat_X, FUN = sum)
#abundance_X <- sum_abd_X$Abundance
#ordered_Sample <- sum_abd_X$Sample[orde(abundance_X)]
#dat_m$Sample <- factor(dat_m$Sample, levels = ordered_Sample)

# Plot dividido según la fecha
svglite("Abundance_bacteria_family_quintin_zone.svg", width=20, height=10)
ggplot(dat_m, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") + xlab("Sample") + ylab("Relative abundance") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold")
  ) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  ggtitle("Bacteria family level") +
  scale_fill_manual(
    values = c("#d7263d","#f46036","#2e294e","#1b998b","#c5d86d","#759aab","#faf2a1","#4d8b31","#ffc800",
               "#1c2541","#3a506b","#5bc0be","#6fffe9","#ff8811","#f4d06f","#fff8f0","#9dd9d2","#eeebd3",
               "#F04C05","#B74211","#7E391D","#442F29","#0B2535","#15640F","#4A5414","#7F4419","#B4331E","#E92323",
               "#EAAD12","#B3A647","#7BA07C","#4499B0","#0C92E5","#5BF30A","#60B835","#667C61","#6B418C","#7005B7",
               "#F44E4E","#B96669","#7D7D83","#42959E","#06ACB8","#EF745C","#D06257","#B15052","#923E4D","#722B47","#531942","#34073D",
               "#074170","#1B4E60","#2F5B51","#436941","#567631","#6A8322","#7E9012","#8ecae6","#219ebc","#023047","#ffb703","#fb8500",
               "#f6d7d7","#df7472","#53032b","#2d0627","#06262c","#72dddf","#0fd780","#2b5403","#a6df72","#c9d70f",
               "#a09ebb","#a8aec1","#b5d2cb","#bfffbc","#a6ffa1","#476a6f","#031a6b","#033860","#004385","#305252",
               "#160f29","#246a73","#368f8b","#f3dfc1","#ddbea8","#ea526f","#e76b74","#d7af70","#c9c19f","#edf7d2",
               "#f94144","#f3722c","#f8961e","#f9844a","#f9c74f","#90be6d","#43aa8b","#4d908e","#577590","#277da1",
               "#f72585","#b5179e","#7209b7","#560bad","#480ca8","#3a0ca3","#3f37c9","#4361ee","#4895ef","#4cc9f0",
               "#d4e09b","#f6f4d2","#cbdfbd","#f19c79","#a44a3f"),
    guide = guide_legend(nrow = 29)
  ) +
  #facet_grid(~ date, scales = "free", space = "free")
  #facet_grid(~ place, scales = "free", space = "free")
  facet_grid(~ zone, scales = "free", space = "free")
dev.off()



###### GENUS LEVEL ######

# Añado fila de "Other" para recuperar la taxa filtrada
column_sums <- colSums(abundance_table_genus_stock)
diff_values <- 100 - column_sums
other_row <- as.data.frame(t(diff_values))
rownames(other_row) <- "Other"
abundance_table_genus_stock <- rbind(abundance_table_genus_stock, other_row)

prop<-prop.table(data.matrix(abundance_table_genus_stock), 2)
dat_m <- melt(prop)
colnames(dat_m)<-c("Genus", "Sample", "Abundance")

# Agrego metadatos
date <- metadata$date[match(dat_m$Sample, metadata$sample)]
zone <- metadata$zone[match(dat_m$Sample, metadata$sample)]
place <- metadata$place[match(dat_m$Sample, metadata$sample)]

dat_m$date <- date
dat_m$zone <- zone
dat_m$place <- place

# Organizar los datos para que queden en orden de menor a mayor según abundancia de X taxa
#X = "Vibrio"
#dat_X <- subset(dat_m, Genus == X)
#sum_abd_X <- aggregate(Abundance ~ Sample, data = dat_X, FUN = sum)
#abundance_X <- sum_abd_X$Abundance
#ordered_Sample <- sum_abd_X$Sample[orde(abundance_X)]
#dat_m$Sample <- factor(dat_m$Sample, levels = ordered_Sample)

# Plot dividido según la fecha
svglite("Abundance_bacteria_genus_quintin_place.svg", width=20, height=10)
ggplot(dat_m, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + xlab("Sample") + ylab("Relative abundance") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold")
  ) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  ggtitle("Bacteria genus level") +
  scale_fill_manual(
    values = c("#d7263d","#f46036","#2e294e","#1b998b","#c5d86d","#759aab","#faf2a1","#4d8b31","#ffc800",
               "#1c2541","#3a506b","#5bc0be","#6fffe9","#ff8811","#f4d06f","#fff8f0","#9dd9d2","#eeebd3",
               "#F04C05","#B74211","#7E391D","#442F29","#0B2535","#15640F","#4A5414","#7F4419","#B4331E","#E92323",
               "#EAAD12","#B3A647","#7BA07C","#4499B0","#0C92E5","#5BF30A","#60B835","#667C61","#6B418C","#7005B7",
               "#F44E4E","#B96669","#7D7D83","#42959E","#06ACB8","#EF745C","#D06257","#B15052","#923E4D","#722B47","#531942","#34073D",
               "#074170","#1B4E60","#2F5B51","#436941","#567631","#6A8322","#7E9012","#8ecae6","#219ebc","#023047","#ffb703","#fb8500",
               "#f6d7d7","#df7472","#53032b","#2d0627","#06262c","#72dddf","#0fd780","#2b5403","#a6df72","#c9d70f",
               "#a09ebb","#a8aec1","#b5d2cb","#bfffbc","#a6ffa1","#476a6f","#031a6b","#033860","#004385","#305252",
               "#160f29","#246a73","#368f8b","#f3dfc1","#ddbea8","#ea526f","#e76b74","#d7af70","#c9c19f","#edf7d2",
               "#f94144","#f3722c","#f8961e","#f9844a","#f9c74f","#90be6d","#43aa8b","#4d908e","#577590","#277da1",
               "#f72585","#b5179e","#7209b7","#560bad","#480ca8","#3a0ca3","#3f37c9","#4361ee","#4895ef","#4cc9f0",
               "#d4e09b","#f6f4d2","#cbdfbd","#f19c79","#a44a3f"),
    guide = guide_legend(nrow = 29)
  ) +
  #facet_grid(~ date, scales = "free", space = "free")
  facet_grid(~ place, scales = "free", space = "free")
  #facet_grid(~ zone, scales = "free", space = "free")
dev.off()


## Heatmap

abund_matrix <- as.matrix(abundance_table_genus_stock)
svglite("heatmap_ab_genus_2hist.svg", width=10, height=10)
heatmap.2(abund_matrix, main="Percentage of Abundances\nBacterial Genus Level", col = bluered(100), 
          trace = "none", density.info = "none", margins=c(10,10), cexRow = 0.7, cexCol = 0.5)
dev.off()

abund_matrix <- as.matrix(a_genus_log)
svglite("heatmap_ab_genus_log_2hist.svg", width=10, height=10)
heatmap.2(abund_matrix, main="Log Transform of Abundances\nBacterial Genus Level", col = bluered(100), 
          trace = "none", density.info = "none", margins=c(10,10), cexRow = 0.7, cexCol = 0.5)
dev.off()

## UpSet plot

genus_comparison <- read.csv("C_bacteria_upset_genus.csv", header = T, sep = ",")
set <- colnames(genus_comparison)[colnames(genus_comparison) != "name"]
svglite("Upset_genera_place.svg", width = 8, height = 8)
#upset(genus_comparison, sets = set,
upset(genus_comparison, sets = c("Inner", "Transition", "Inlet"),
#upset(genus_comparison, sets = c("Bare", "Seagrass"),
#upset(genus_comparison, sets = c("Relaxing", "Upwelling"),
      sets.bar.color="#56B4E9", order.by="freq",empty.intersections = "on",
      mainbar.y.label = "Number of shared and\nnon-shared genera", sets.x.label = "Number of genera",
      point.size = 3, line.size = 0.5, query.legend = "top",
      text.scale = c(2,1.5,1.5,1,1.5,1))
dev.off()

###### SPECIES LEVEL ######

# Añado fila de "Other" para recuperar la taxa filtrada
column_sums <- colSums(abundance_table_species_stock)
diff_values <- 100 - column_sums
other_row <- as.data.frame(t(diff_values))
rownames(other_row) <- "Other"
abundance_table_species_stock <- rbind(abundance_table_species_stock, other_row)

prop<-prop.table(data.matrix(abundance_table_species_stock), 2)
dat_m <- melt(prop)
colnames(dat_m)<-c("Species", "Sample", "Abundance")

# Agrego metadatos
date <- metadata$date[match(dat_m$Sample, metadata$sample)]
zone <- metadata$zone[match(dat_m$Sample, metadata$sample)]
place <- metadata$place[match(dat_m$Sample, metadata$sample)]

dat_m$date <- date
dat_m$zone <- zone
dat_m$place <- place

# Organizar los datos para que queden en orden de menor a mayor según abundancia de X taxa
#X = "Vibrio"
#dat_X <- subset(dat_m, Genus == X)
#sum_abd_X <- aggregate(Abundance ~ Sample, data = dat_X, FUN = sum)
#abundance_X <- sum_abd_X$Abundance
#ordered_Sample <- sum_abd_X$Sample[orde(abundance_X)]
#dat_m$Sample <- factor(dat_m$Sample, levels = ordered_Sample)

# Plot dividido según la fecha
svglite("Abundance_bacteria_species_quintin_date.svg", width=20, height=10)
ggplot(dat_m, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") + xlab("Sample") + ylab("Relative abundance") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold")
  ) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  ggtitle("Bacteria species level") +
  scale_fill_manual(
    values = c("#d7263d","#f46036","#2e294e","#1b998b","#c5d86d","#759aab","#faf2a1","#4d8b31","#ffc800",
               "#1c2541","#3a506b","#5bc0be","#6fffe9","#ff8811","#f4d06f","#fff8f0","#9dd9d2","#eeebd3",
               "#F04C05","#B74211","#7E391D","#442F29","#0B2535","#15640F","#4A5414","#7F4419","#B4331E","#E92323",
               "#EAAD12","#B3A647","#7BA07C","#4499B0","#0C92E5","#5BF30A","#60B835","#667C61","#6B418C","#7005B7",
               "#F44E4E","#B96669","#7D7D83","#42959E","#06ACB8","#EF745C","#D06257","#B15052","#923E4D","#722B47","#531942","#34073D",
               "#074170","#1B4E60","#2F5B51","#436941","#567631","#6A8322","#7E9012","#8ecae6","#219ebc","#023047","#ffb703","#fb8500",
               "#f6d7d7","#df7472","#53032b","#2d0627","#06262c","#72dddf","#0fd780","#2b5403","#a6df72","#c9d70f",
               "#a09ebb","#a8aec1","#b5d2cb","#bfffbc","#a6ffa1","#476a6f","#031a6b","#033860","#004385","#305252",
               "#160f29","#246a73","#368f8b","#f3dfc1","#ddbea8","#ea526f","#e76b74","#d7af70","#c9c19f","#edf7d2",
               "#f94144","#f3722c","#f8961e","#f9844a","#f9c74f","#90be6d","#43aa8b","#4d908e","#577590","#277da1",
               "#f72585","#b5179e","#7209b7","#560bad","#480ca8","#3a0ca3","#3f37c9","#4361ee","#4895ef","#4cc9f0",
               "#d4e09b","#f6f4d2","#cbdfbd","#f19c79","#a44a3f"),
    guide = guide_legend(nrow = 29)
  ) +
  facet_grid(~ date, scales = "free", space = "free")
  #facet_grid(~ place, scales = "free", space = "free")
  #facet_grid(~ zone, scales = "free", space = "free")
dev.off()


## Heatmap

abund_matrix <- as.matrix(abundance_table_species_stock)
svglite("heatmap_ab_species_2hist.svg", width=10, height=10)
heatmap.2(abund_matrix, main="Percentage of Abundances\nBacterial Species Level", col = bluered(100), 
          trace = "none", density.info = "none", margins=c(10,15), cexRow = 0.7, cexCol = 0.5)
dev.off()

abund_matrix <- as.matrix(a_species_log)
svglite("heatmap_ab_species_log_2hist.svg", width=10, height=10)
heatmap.2(abund_matrix, main="Log Transform of Abundances\nBacterial Species Level", col = bluered(100), 
          trace = "none", density.info = "none", margins=c(10,15), cexRow = 0.7, cexCol = 0.5)
dev.off()

## UpSet plot

species_comparison <- read.csv("C_bacteria_upset_species.csv", header = T, sep = ",")
set <- colnames(species_comparison)[colnames(species_comparison) != "name"]
svglite("Upset_species_sector.svg", width = 10, height = 8)
#upset(species_comparison, sets = set,
#upset(species_comparison, sets = c("Inner", "Transition", "Inlet"),
upset(species_comparison, sets = c("Bare", "Seagrass"),
#upset(species_comparison, sets = c("Relaxing", "Upwelling"),
        sets.bar.color="#56B4E9", order.by="freq",empty.intersections = "on",
        mainbar.y.label = "Number of shared and\nnon-shared species", sets.x.label = "Number of species",
        point.size = 3, line.size = 0.5, query.legend = "top",
        text.scale = c(2,1.5,1.5,1,1.5,1))
dev.off()

# Boxplot with specific genus

# Filtrando las taxa que tenga al mínimo el 1% en TODAS las muestras
#abundance_table_class_core <- abundance_table_class[apply(abundance_table_class, 1, function(x) all(x >= 1 | is.na(x))), ]
abundance_table_family_core <- abundance_table_family[apply(abundance_table_family, 1, function(x) all(x > 1 | is.na(x))), ]
#abundance_table_genus_core <- abundance_table_genus[apply(abundance_table_genus, 1, function(x) all(x > 0.1 | is.na(x))), ]
abundance_table_species_core <- abundance_table_species[apply(abundance_table_species, 1, function(x) all(x > 0.1 | is.na(x))), ]

# Filtrar las taxas que cumplen con el criterio
filt_family <- rownames(abundance_table_family_core)
filt_species <- rownames(abundance_table_species_core)

# Almacenar la lista de especies en el vector species
family <- c(filt_family)
species <- c(filt_species)

# Seleccionar géneros a visualizar
#generos <- c("Desulfopila","Bellilinea","Dehalogenimonas","Pseudomonas","Desulfosarcina","Alcanivorax","Dehalogenimonas","Desulfobulbus",
#          "Thioprofundum","Desulfobacterium","Geoalkalibacter","SEEP-SRB1","Candidatus_Magnetananas","Desulfomonile","Desulfotalea",
#          "Desulfovibrio","Blastopirellula","Pelobacter","Desulfonema","Rhodopirellula")

dat_m_family <- dat_m[dat_m$Family %in% family, ]
dat_m_family_21 <- subset(dat_m_family, dat_m_family$date == "2021_oct")
dat_m_family_22 <- subset(dat_m_family, dat_m_family$date == "2022_jun")

#dat_m_species <- dat_m[dat_m$Species %in% species, ]
#dat_m_species_21 <- subset(dat_m_species, dat_m_species$date == "2021_oct")
#dat_m_species_22 <- subset(dat_m_species, dat_m_species$date == "2022_jun")

#svglite("Boxplot_bacterial_family_21.svg", width=15, height=10)
#ggplot(dat_m_family_21, aes(x = Family, y = Abundance, fill = zone)) +
#ggtitle("Comparison of bacterial families present in all samples in 2021/oct") +

svglite("Boxplot_bacterial_family_22.svg", width=15, height=10)
ggplot(dat_m_family_22, aes(x = Family, y = Abundance, fill = zone)) +
ggtitle("Comparison of bacterial families present in all samples in 2022/jun") +

#svglite("Boxplot_bacterial_species_21.svg", width=15, height=10)
#ggplot(dat_m_species_21, aes(x = Species, y = Abundance, fill = zone)) +
#  ggtitle("Comparison of bacterial species in present all samples in 2021/oct") +

svglite("Boxplot_bacterial_species_22.svg", width=15, height=10)
ggplot(dat_m_species_22, aes(x = Species, y = Abundance, fill = zone)) +
  ggtitle("Comparison of bacterial species present in all samples in 2022/jun") +

    geom_boxplot(size = 0.2) +
  xlab("Species") +
  ylab("Relative Abundance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  scale_fill_manual(
    values = c(Boca="#d45500ff",Cabeza="#d4aa00ff",Intermedia="#88aa00ff")
  ) 
dev.off()


############################
#### HEATMAP ABUNDANCES ####
############################

amp_class <- amp_load(otutable = abundance_table_class, metadata = metadata)
class_hm <- amp_heatmap(amp_class,
            group_by = "Habitat",
            facet_by = "Season",
            tax_aggregate = "OTU",
            tax_show = 20,
            color_vector = c("white","darkred"),
            plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right")

amp_genus <- amp_load(otutable = abundance_table_genus, metadata = metadata)
genus_hm <- amp_heatmap(amp_genus,
            group_by = "Season",
            facet_by = "Sector",
            tax_aggregate = "OTU",
            tax_show = 20,
            color_vector = c("white","darkred"),
            plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right")

amp_species <- amp_load(otutable = abundance_table_species, metadata = metadata)
species_hm <- amp_heatmap(amp_species,
            group_by = "Season",
            facet_by = "Sector",
            tax_aggregate = "OTU",
            tax_show = 20,
            color_vector = c("white","darkred"),
            plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right")

svglite("C_bacteria_heatmap_class_season_habitat.svg", width=6, height=6)
class_hm
dev.off()

svglite("C_bacteria_heatmap_genus_sector_season.svg", width=6, height=8)
genus_hm
dev.off()

svglite("C_bacteria_heatmap_species_sector_season.svg", width=7.5, height=8)
species_hm
dev.off()


amp_annot<- amp_load(otutable = abundance_table_annotation, metadata = metadata)
season_hm <- amp_heatmap(amp_annot,
                        group_by = "Habitat",
                        facet_by = "Season",
                        tax_aggregate = "OTU",
                        tax_show = 35,
                        round = 2,
                        color_vector = c("white","darkred"),
                        plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right")

sector_hm <- amp_heatmap(amp_annot,
                        group_by = "Habitat",
                        facet_by = "Sector",
                        tax_aggregate = "OTU",
                        tax_show = 35,
                        round = 2,
                        color_vector = c("white","darkred"),
                        plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right")

both_hm <- amp_heatmap(amp_annot,
                       group_by = "Season",
                       facet_by = "Sector",
                       tax_aggregate = "OTU",
                       tax_show = 35,
                       round = 2,
                       color_vector = c("white","darkred"),
                       plot_colorscale = "sqrt") +
  theme(axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right")

svglite("C_bacteria_heatmap_class_season_habitat.svg", width=6, height=6)
class_hm
dev.off()

svglite("C_bacteria_heatmap_genus_season_habitat.svg", width=6, height=6)
genus_hm
dev.off()

svglite("C_bacteria_heatmap_species_season_habitat.svg", width=7.5, height=6)
species_hm
dev.off()



### Shared groups

# Cargar la tabla
abund_arc_bac <- read.csv("C_archaea_bacteria_abd_genus.csv", row.names = 1)
abund_arc_bac <- read.csv("C_arc_bac_abd_species.csv", row.names = 1)

# Filtrar filas que tengan al menos un valor >= 0.1
abund_arc_bac_filtered <- abund_arc_bac[apply(abund_arc_bac >= 0.1, 1, any), ]

# Definir los grupos de columnas
grupo_inlet <- c("C_17", "C_18", "C_23", "C_24")
grupo_transition <- c("C_19", "C_20", "C_25", "C_26")
grupo_inner <- c("C_21", "C_22","C_27", "C_28")

# Crear subconjuntos de los datos
subgrupo_1 <- abund_arc_bac_filtered[, grupo_inlet]
subgrupo_2 <- abund_arc_bac_filtered[, grupo_transition]
subgrupo_3 <- abund_arc_bac_filtered[, grupo_inner]

# Definir presencia como abundancia ≥ 0.1
presente_1 <- apply(subgrupo_1 >= 0.1, 1, any)  # Filas con al menos un valor ≥ 0.1 en grupo_1
presente_2 <- apply(subgrupo_2 >= 0.1, 1, any)  # Filas con al menos un valor ≥ 0.1 en grupo_2
presente_3 <- apply(subgrupo_3 >= 0.1, 1, any)  # Filas con al menos un valor ≥ 0.1 en grupo_2


# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_2 & presente_3)

# Mostrar resultado
cat("Número de géneros compartidos entre los dos grupos:", generos_compartidos, "\n")

# Obtener los géneros compartidos entre los dos grupos
generos_compartidos <- rownames(abund_arc_bac_filtered)[presente_1 & presente_2 & presente_3]
# Filtrar solo los géneros compartidos
abund_compartidos <- abund_arc_bac_filtered[generos_compartidos, ]
# Sumar las abundancias de los géneros compartidos en cada muestra (columna)
sum_abundancias_por_muestra <- colSums(abund_compartidos, na.rm = TRUE)
# Mostrar el resultado
print(sum_abundancias_por_muestra)


# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_2 & !presente_3)

# Mostrar resultado
cat("Número de géneros compartidos entre los dos grupos:", generos_compartidos, "\n")

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_2 & presente_3 & !presente_1)

# Mostrar resultado
cat("Número de géneros compartidos entre los dos grupos:", generos_compartidos, "\n")

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_3 & !presente_2)

# Mostrar resultado
cat("Número de géneros compartidos entre los dos grupos:", generos_compartidos, "\n")

# Obtener los géneros únicos de cada grupo
unicos_grupo_1 <- sum(presente_1 & !presente_2 & !presente_3)  # Presentes en grupo_1 pero no en grupo_2
unicos_grupo_2 <- sum(presente_2 & !presente_1 & !presente_3)  # Presentes en grupo_2 pero no en grupo_1
unicos_grupo_3 <- sum(presente_3 & !presente_1 & !presente_2)  # Presentes en grupo_2 pero no en grupo_1

# Mostrar los resultados
cat("Número de géneros únicos en el grupo 1:", unicos_grupo_1, "\n")
cat("Número de géneros únicos en el grupo 2:", unicos_grupo_2, "\n")
cat("Número de géneros únicos en el grupo 3:", unicos_grupo_3, "\n")





# Definir los grupos de columnas
grupo_bare <- c("C_17", "C_19", "C_21", "C_23", "C_25", "C_27")
grupo_seagrass <- c("C_18", "C_20", "C_22", "C_24", "C_26", "C_28")
grupo_intense <- c("C_17", "C_18", "C_19", "C_20", "C_21", "C_22")
grupo_relax <- c("C_23", "C_24", "C_25", "C_26", "C_27", "C_28")

# Crear subconjuntos de los datos
subgrupo_1 <- abund_arc_bac_filtered[, grupo_bare]
subgrupo_2 <- abund_arc_bac_filtered[, grupo_seagrass]
subgrupo_3 <- abund_arc_bac_filtered[, grupo_intense]
subgrupo_4 <- abund_arc_bac_filtered[, grupo_relax]

# Definir presencia como abundancia ≥ 0.1
presente_1 <- apply(subgrupo_1 >= 0.1, 1, any)  # Filas con al menos un valor ≥ 0.1 en grupo_1
presente_2 <- apply(subgrupo_2 >= 0.1, 1, any)  # Filas con al menos un valor ≥ 0.1 en grupo_2
presente_3 <- apply(subgrupo_3 >= 0.1, 1, any)  # Filas con al menos un valor ≥ 0.1 en grupo_2
presente_4 <- apply(subgrupo_4 >= 0.1, 1, any)  # Filas con al menos un valor ≥ 0.1 en grupo_2

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_2 & presente_3 & presente_4)
generos_compartidos

# Obtener los géneros compartidos entre los dos grupos
generos_compartidos <- rownames(abund_arc_bac_filtered)[presente_1 & presente_2 & presente_3 & presente_4]
# Filtrar solo los géneros compartidos
abund_compartidos <- abund_arc_bac_filtered[generos_compartidos, ]
# Sumar las abundancias de los géneros compartidos en cada muestra (columna)
sum_abundancias_por_muestra <- colSums(abund_compartidos, na.rm = TRUE)
# Mostrar el resultado
print(sum_abundancias_por_muestra)



# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_3 & presente_4 & !presente_2)
generos_compartidos

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_2 & presente_3 & presente_4 & !presente_1)
generos_compartidos

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_3 & !presente_4 & !presente_2)
generos_compartidos

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_4 & !presente_3 & !presente_2)
generos_compartidos

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_2 & presente_3 & !presente_4 & !presente_1)
generos_compartidos

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_2 & presente_4 & !presente_1 & !presente_3)
generos_compartidos

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_2 & presente_3 & !presente_4)
generos_compartidos

# Contar los géneros compartidos (presentes en ambos grupos)
generos_compartidos <- sum(presente_1 & presente_2 & presente_4 & !presente_3)
generos_compartidos






###########################################
##### DIFFERENTIAL ABUNDANCE ANALYSIS #####
###########################################

## If analysis with archaea and bacteria 

# Cargar las tablas a nivel de clase
count_class_arc_bac <- read.csv("C_arc_bacteria_count_class.csv", row.names = 1)

# Phyloseq object
taxa_table <- row.names(count_class_arc_bac)
setdiff(rownames(count_class_arc_bac),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_class_arc_bac, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$Sample
ps_count_class = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_class <- prune_taxa(taxa_sums(ps_count_class)>0, ps_count_class)
print(ps_count_class)



# Cargar las tablas a nivel de genero
count_genus_arc_bac <- read.csv("C_archaea_bacteria_count_genus.csv", row.names = 1)
abund_genus_arc_bac <- read.csv("C_archaea_bacteria_abd_genus.csv", row.names = 1)


# Read metadata
metadata <- read.csv("metadata_quintin.csv")

# Filtrar filas que tengan al menos un valor >= 0.1
abund_genus_arc_bac_filt <- abund_genus_arc_bac[apply(abund_genus_arc_bac >= 0.1, 1, any), ]

# Filtrar count_genus_arc_bac para que tenga las mismas filas que abund_genus_arc_bac_filt
#count_genus_arc_bac_filt <- count_genus_arc_bac[rownames(abund_genus_arc_bac_filt), , drop = FALSE]

# Phyloseq object
taxa_table <- row.names(count_genus_arc_bac)
setdiff(rownames(count_genus_arc_bac),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_genus_arc_bac, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$Sample
ps_count_genus = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_genus <- prune_taxa(taxa_sums(ps_count_genus)>0, ps_count_genus)
print(ps_count_genus)

save(ps_count_genus, file = "ps_count_genus.RData")

# Cargar las tablas a nivel de especie
count_species_arc_bac <- read.csv("C_archaea_bacteria_count_species.csv", row.names = 1)
abund_species_arc_bac <- read.csv("C_arc_bac_abd_species.csv", row.names = 1)

# Filtrar filas que tengan al menos un valor >= 0.1
#abund_species_arc_bac_filt <- abund_species_arc_bac[apply(abund_species_arc_bac >= 0.1, 1, any), ]

# Filtrar count_genus_arc_bac para que tenga las mismas filas que abund_species_arc_bac_filt
#count_species_arc_bac_filt <- count_species_arc_bac[rownames(abund_species_arc_bac_filt), , drop = FALSE]

taxa_table <- row.names(count_species_arc_bac)
setdiff(rownames(count_species_arc_bac),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_species_arc_bac, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
ps_count_species = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_species <- prune_taxa(taxa_sums(ps_count_species)>0, ps_count_species)
print(ps_count_species)

save(ps_count_species, file = "ps_count_species.RData")


## If analysis with only bacteria

# Read count tables
count_class <- read.csv("C_bacteria_count_class.csv", header = T)
count_class <- column_to_rownames(count_class, var = "name")
count_class <- as.data.frame(count_class)
count_genus <- read.csv("C_bacteria_count_genus.csv", header = T)
count_genus <- column_to_rownames(count_genus, var = "name")
count_genus <- as.data.frame(count_genus)
count_species <- read.csv("C_bacteria_count_species.csv", header = T)
count_species <- column_to_rownames(count_species, var = "name")
count_species <- as.data.frame(count_species)

# Phyloseq objects
taxa_table <- row.names(count_class)
setdiff(rownames(count_class),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_class, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$Sample
ps_count_class = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_class <- prune_taxa(taxa_sums(ps_count_class)>0, ps_count_class)
print(ps_count_class)

taxa_table <- row.names(count_genus)
setdiff(rownames(count_genus),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_genus, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
ps_count_genus = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_genus <- prune_taxa(taxa_sums(ps_count_genus)>0, ps_count_genus)
print(ps_count_genus)

taxa_table <- row.names(count_species)
setdiff(rownames(count_species),taxa_table)
taxa_table =  as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(count_species, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
ps_count_species = phyloseq(OTU_taxa, TAX_taxa, sampledata)
# Keep only taxa with positive sums
ps_count_species <- prune_taxa(taxa_sums(ps_count_species)>0, ps_count_species)
print(ps_count_species)


#ANCOM

#out1 = ancombc(phyloseq = ps_count_class, formula = "Sector + date + place",
#out1 = ancombc(phyloseq = ps_count_species, formula = "Sector + date + place",
out1 = ancombc(phyloseq = ps_count_genus, formula = "Sector + Season + Habitat",
               p_adj_method = "fdr", lib_cut = 0, prv_cut = 0.1,
               group = "Sector", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5, 
               max_iter = 100, conserve = FALSE, alpha = 0.05, global = TRUE)
res1 = out1$res
res1_global = out1$res_global

#LFC
tab_lfc1 = res1$lfc
col_name = c("Taxon","Intercept","TRAN-INL","INN-INL","INTENSE-RELAX","SEAGRASS-BARE")
colnames(tab_lfc1) = col_name
tab_diff1 = res1$diff_abn
colnames(tab_diff1) = col_name
tab_se1 = res1$se
colnames(tab_se1) = col_name
tab_q1 = res1$q_val
colnames(tab_q1) = col_name

df_lfc1 = data.frame(tab_lfc1[,-1] * tab_diff1[,-1], check.names = FALSE) %>%
  mutate(taxon_id = tab_diff1$Taxon) %>%
  dplyr::select(taxon_id, everything())
df_se1 = data.frame(tab_se1[, -1] * tab_diff1[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = tab_diff1$Taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se1)[-1] = paste0(colnames(df_se1)[-1], "SE")

# Filtrar df_lfc1 para que tenga las mismas filas que abund_genus_arc_bac_filt
df_lfc1 <- df_lfc1[df_lfc1$taxon_id %in% rownames(abund_genus_arc_bac_filt), , drop = FALSE]

#ps_small <- phyloseq::subset_samples(ps_count_class, Sector %in% c("B_Transition", "C_Inner"))
ps_small <- phyloseq::subset_samples(ps_count_genus, Sector %in% c("B_Transition", "C_Inner"))
#ps_small <- phyloseq::subset_samples(ps_count_species, Sector %in% c("B_Transition", "C_Inner"))
out2 = ancombc(phyloseq = ps_small, formula = "Sector",
               p_adj_method = "fdr", lib_cut = 0, prv_cut = 0.1,
               group = "Sector", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5, 
               max_iter = 100, conserve = FALSE, alpha = 0.05, global = TRUE)
res2 = out2$res
res2_global = out2$res_global

#LFC
tab_lfc2 = res2$lfc
col_name = c("Taxon","Intercept","INN-TRAN")
colnames(tab_lfc2) = col_name
tab_diff2 = res2$diff_abn
colnames(tab_diff2) = col_name
tab_se2 = res2$se
colnames(tab_se2) = col_name
tab_q2 = res2$q_val
colnames(tab_q2) = col_name

df_lfc2 = data.frame(tab_lfc2[,-1] * tab_diff2[,-1], check.names = FALSE) %>%
  mutate(taxon_id = tab_diff2$Taxon) %>%
  dplyr::select(taxon_id, everything())
df_se2 = data.frame(tab_se2[, -1] * tab_diff2[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = tab_diff2$Taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se2)[-1] = paste0(colnames(df_se2)[-1], "SE")

# Filtrar df_lfc2 para que tenga las mismas filas que abund_genus_arc_bac_filt
df_lfc2 <- df_lfc2[df_lfc1$taxon_id %in% rownames(abund_genus_arc_bac_filt), , drop = FALSE]


# Visualization LFC heatmap

df_fig_Sector = df_lfc1 %>% 
  dplyr::left_join(df_lfc2, by = "taxon_id") %>%
  filter(`TRAN-INL` != 0 | `INN-INL` != 0 | `INN-TRAN` != 0) %>%
  dplyr::filter(`TRAN-INL` < -1.5 | `TRAN-INL` > 1.5 |
                  `INN-INL` < -1.5 | `INN-INL` > 1.5 |
                  `INN-TRAN` < -1.5 | `INN-TRAN` > 1.5) %>%
  transmute(taxon_id, 
            `TRAN-INL` = round(`TRAN-INL`, 2),
            `INN-INL` = round(`INN-INL`, 2),
            `INN-TRAN` = round(`INN-TRAN`, 2)) %>%
  pivot_longer(cols = `TRAN-INL`:`INN-INL`:`INN-TRAN`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_Sector$value))
up = ceiling(max(df_fig_Sector$value))
mid = (lo + up)/2
p_Sector = df_fig_Sector %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(-3, 3),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 3.5) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared - Sector") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1),
        axis.text.y = element_text(size = 10))

#svglite("Log_fold_changes_class_Sector.svg", width=4, height=3)
#svglite("Log_fold_changes_genus_Sector.svg", width=5, height=12)
svglite("Log_fold_changes_species_Sector_250209.svg", width=7, height=7)
p_Sector
dev.off()

df_fig_Date = df_lfc1 %>% 
  filter(`INTENSE-RELAX` != 0) %>%
  dplyr::filter(`INTENSE-RELAX` < -1 | `INTENSE-RELAX` > 1) %>%
  transmute(taxon_id, 
            `INTENSE-RELAX` = round(`INTENSE-RELAX`, 2)) %>%
  pivot_longer(cols = `INTENSE-RELAX`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_Date$value))
up = ceiling(max(df_fig_Date$value))
mid = (lo + up)/2
p_Date = df_fig_Date %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(-3, 3),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 3.5) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared - Date") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

#svglite("Log_fold_changes_class_Date.svg", width=3, height=3)
#svglite("Log_fold_changes_genus_Date.svg", width=3, height=4)
svglite("Log_fold_changes_species_Date_250209.svg", width=4.0, height=3.9)
p_Date
dev.off()


df_fig_Place = df_lfc1 %>% 
  filter(`SEAGRASS-BARE` != 0) %>%
  dplyr::filter(`SEAGRASS-BARE` < -0.64 | `SEAGRASS-BARE` > 0.64) %>%
  transmute(taxon_id, 
            `SEAGRASS-BARE` = round(`SEAGRASS-BARE`, 2)) %>%
  pivot_longer(cols = `SEAGRASS-BARE`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_Place$value))
up = ceiling(max(df_fig_Place$value))
mid = (lo + up)/2
p_Place = df_fig_Place %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(-1, 1),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 3.3) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared - Place") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

#svglite("Log_fold_changes_class_Place.svg", width=3, height=3)
#svglite("Log_fold_changes_genus_Place.svg", width=3, height=3)
svglite("Log_fold_changes_species_Place_250209.svg", width=3.75, height=3)
p_Place
dev.off()


######################################
#### PCA and dispersion analysis ####
######################################

#CLR transform
#ps_c_clr <- microbiome::transform(ps_count_class, "clr")
ps_c_clr <- microbiome::transform(ps_count_genus, "clr")
#ps_c_clr <- microbiome::transform(ps_count_species, "clr")

#Generate Aitchison distance matrix
clr_c_matrix <- phyloseq::distance(ps_c_clr, method = "euclidean")

#ADONIS test
vegan::adonis2(clr_c_matrix ~ phyloseq::sample_data(ps_c_clr)$Sector)
vegan::adonis2(clr_c_matrix ~ phyloseq::sample_data(ps_c_clr)$Habitat)
vegan::adonis2(clr_c_matrix ~ phyloseq::sample_data(ps_c_clr)$Season)
# Pairwise comparisons of microbial communities structure
library(RVAideMemoire)
set.seed(123)
pairwise.perm.manova(clr_c_matrix, phyloseq::sample_data(ps_c_clr)$Sector, nperm = 10000)
pairwise.perm.manova(clr_c_matrix, phyloseq::sample_data(ps_c_clr)$Season, nperm = 10000)
pairwise.perm.manova(clr_c_matrix, phyloseq::sample_data(ps_c_clr)$Habitat, nperm = 10000)

#Dispersion test
dispr <- vegan::betadisper(clr_c_matrix, phyloseq::sample_data(ps_c_clr)$Sector)
dispr
group_levels <- levels(dispr$group)
svglite("Disp_Sector_class_250220.svg", width=4, height=4)
#svglite("Disp_Sector_genus_250220.svg", width=4, height=4)
#svglite("Disp_Sector_species_250220.svg", width=4, height=4)
plot(dispr, main = "", sub = "")
legend("topright", legend = group_levels, col = 1:length(group_levels), pch = 1)
dev.off()
boxplot(dispr, main = "Boxplot Distance to centroid - Sector", xlab = "")
permutest(dispr)
set.seed(123)
permutest(dispr, pairwise = TRUE)

dispr <- vegan::betadisper(clr_c_matrix, phyloseq::sample_data(ps_c_clr)$Season)
dispr
group_levels <- levels(dispr$group)
svglite("Disp_Season_class_250220.svg", width=4, height=4)
#svglite("Disp_Season_genus_250220.svg", width=4, height=4)
#svglite("Disp_Season_species_250220.svg", width=4, height=4)
plot(dispr, main = "", sub = "")
legend("topright", legend = group_levels, col = 1:length(group_levels), pch = 1)
dev.off()
boxplot(dispr, main = "Boxplot Distance to centroid - Sector", xlab = "")
permutest(dispr)
permutest(dispr, pairwise = TRUE)

dispr <- vegan::betadisper(clr_c_matrix, phyloseq::sample_data(ps_c_clr)$Habitat)
dispr
group_levels <- levels(dispr$group)
svglite("Disp_Place_class_250220.svg", width=4, height=4)
#svglite("Disp_Place_genus_250220.svg", width=4, height=4)
#svglite("Disp_Place_species_250220.svg", width=4, height=4)
plot(dispr, main = "", sub = "")
legend("topright", legend = group_levels, col = 1:length(group_levels), pch = 1)
dev.off()
boxplot(dispr, main = "Boxplot Distance to centroid - Sector", xlab = "")
permutest(dispr)
permutest(dispr, pairwise = TRUE)



#######################################
##### ANALYSIS OF ALPHA DIVERSITY #####
#######################################

plot_richness(ps_count_class)

a_my_comparisons <- list(c("A_Inlet", "B_Transition"), c("A_Inlet", "C_Inner"), c("C_Inner", "B_Transition"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
ad_sector <- plot_richness(ps_count_class, x="Sector", measures = c("Observed","Shannon", "Simpson"), color = "Sector") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#2a7fffff","#808080ff","#c8ab37ff"))
  #+ stat_compare_means(method = "kruskal.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Relaxing", "Upwelling"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
ad_season <- plot_richness(ps_count_class, x="Season", measures = c("Observed","Shannon", "Simpson"), color = "Season") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6)+ scale_color_manual(values = c("#ff0000ff","#000080ff")) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.format")

a_my_comparisons <- list(c("Bare", "Seagrass"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
ad_habitat <- plot_richness(ps_count_class, x="Habitat", measures = c("Observed","Shannon", "Simpson"), color = "Habitat") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#ff6600ff","#338000ff"))
  #+ stat_compare_means(method = "kruskal.test", comparisons = a_my_comparisons, label = "p.signif")

svglite("Alpha_Diversity_class.svg", width=6, height=10)
  ad_sector / ad_habitat / ad_season
dev.off()


load("ps_count_genus.RData")

plot_richness(ps_count_genus)

a_my_comparisons <- list(c("A_Inlet", "B_Transition"), c("A_Inlet", "C_Inner"), c("C_Inner", "B_Transition"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
#ad_sector <- plot_richness(ps_count_genus, x="Sector", measures = c("Observed","Shannon", "Simpson"), color = "Sector") + theme_bw() + 
ad_sector <- plot_richness(ps_count_genus, x="Sector", measures = c("Observed","Shannon"), color = "Sector") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#2a7fffff","#808080ff","#c8ab37ff")) 
#  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Relaxing", "Upwelling"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
#ad_season <- plot_richness(ps_count_genus, x="Season", measures = c("Observed","Shannon", "Simpson"), color = "Season") + theme_bw() + 
ad_season <- plot_richness(ps_count_genus, x="Season", measures = c("Observed","Shannon"), color = "Season") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6)+ scale_color_manual(values = c("#ff0000ff","#000080ff"))
#  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Bare", "Seagrass"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
#ad_habitat <- plot_richness(ps_count_genus, x="Habitat", measures = c("Observed","Shannon", "Simpson"), color = "Habitat") + theme_bw() + 
ad_habitat <- plot_richness(ps_count_genus, x="Habitat", measures = c("Observed","Shannon"), color = "Habitat") + theme_bw() + 
    geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#ff6600ff","#338000ff"))
#+ stat_compare_means(method = "kruskal.test", comparisons = a_my_comparisons, label = "p.signif")

svglite("Alpha_Diversity_genus_250316.svg", width=5, height=7)
  ad_sector / ad_habitat / ad_season
dev.off()


load("ps_count_species.RData")

plot_richness(ps_count_species)

a_my_comparisons <- list(c("A_Inlet", "B_Transition"), c("A_Inlet", "C_Inner"), c("C_Inner", "B_Transition"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
ad_sector <- plot_richness(ps_count_species, x="Sector", measures = c("Observed","Shannon", "Simpson"), color = "Sector") + theme_bw() + 
#ad_sector <- plot_richness(ps_count_species, x="Sector", measures = c("Observed","Shannon"), color = "Sector") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#2a7fffff","#808080ff","#c8ab37ff"))
  #stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Relaxing", "Upwelling"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
#ad_season <- plot_richness(ps_count_species, x="Season", measures = c("Observed","Shannon", "Simpson"), color = "Season") + theme_bw() + 
ad_season <- plot_richness(ps_count_species, x="Season", measures = c("Observed","Shannon"), color = "Season") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6)+ scale_color_manual(values = c("#ff0000ff","#000080ff"))
  #stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Bare", "Seagrass"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
#ad_habitat <- plot_richness(ps_count_species, x="Habitat", measures = c("Observed","Shannon", "Simpson"), color = "Habitat") + theme_bw() + 
ad_habitat <- plot_richness(ps_count_species, x="Habitat", measures = c("Observed","Shannon"), color = "Habitat") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#ff6600ff","#338000ff"))
  #stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

svglite("Alpha_Diversity_species_250209.svg", width=5, height=7)
  ad_sector / ad_habitat / ad_season
dev.off()



# With the abundances tables using the vegan package
data_richness <- estimateR(t(count_genus))  
shannon_index <- diversity(abundance_table_genus, index = "shannon")
simpson_index <- diversity(abundance_table_genus, index = "simpson")
chao1_index <- diversity(abundance_table_genus, index = "chao1")
data_alphadiv <- cbind(data_grp, t(data_richness), data_shannon, data_evenness) 



# Using rarefy samples

# Suponiendo que count_genus es una tabla con géneros en filas y muestras en columnas
otu_table_genus <- otu_table(count_genus, taxa_are_rows = TRUE)

# Aplicar rarefacción
Genus_RR <- rarefy_even_depth(ps_count_genus, sample.size = min(sample_sums(otu_table_genus)),
                              rngseed = 711, replace = TRUE, 
                              trimOTUs = TRUE, verbose = TRUE)

a_my_comparisons <- list(c("A_Inlet", "B_Transition"), c("A_Inlet", "C_Inner"), c("C_Inner", "B_Transition"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
ad_sector <- plot_richness(Genus_RR, x="Sector", measures = c("Observed","Shannon"), color = "Sector") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#2a7fffff","#808080ff","#c8ab37ff")) 
  #stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Relaxing", "Upwelling"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
ad_season <- plot_richness(Genus_RR, x="Season", measures = c("Observed","Shannon", "Simpson"), color = "Season") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6)+ scale_color_manual(values = c("#ff0000ff","#000080ff")) +
  #stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Bare", "Seagrass"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))
ad_habitat <- plot_richness(Genus_RR, x="Habitat", measures = c("Observed","Shannon", "Simpson"), color = "Habitat") + theme_bw() + 
  geom_boxplot(lwd = 0.7, alpha = 0.6) + scale_color_manual(values = c("#ff6600ff","#338000ff")) +
  #stat_compare_means(method = "kruskal.test", comparisons = a_my_comparisons, label = "p.signif")

svglite("Alpha_Diversity_rarefied_genus.svg", width=5, height=7)
ad_sector / ad_habitat / ad_season
dev.off()


######################################
##### ANALYSIS OF BETA DIVERSITY #####
######################################


# Read the abundance tables
abundance_table_phylum <- read_tsv("sorted_perc_phylum_matrix.txt")
abundance_table_genus <- read_tsv("sorted_perc_genus_matrix.txt")

# Filter abundance_table_phylum
# Se eliminarán las abundancias menores a 1%. Se puede cambiar este valor si se desea por 0.005 o 0.001
abundance_table_phylum <- column_to_rownames(abundance_table_phylum, var = "phylum")
abundance_table_phylum <- abundance_table_phylum[rowSums(abundance_table_phylum[, -1] > 0) > 0, ]

# Summarize each genus
abundance_table_genus <- aggregate(. ~ genus, data = abundance_table_genus, FUN = sum)
# Filtrar abundance_table_genus
abundance_table_genus <- column_to_rownames(abundance_table_genus, var = "genus")
abundance_table_genus <- abundance_table_genus[rowSums(abundance_table_genus[, -1] > 0) > 0, ]

# Read metadata
metadata <- read_tsv("1_pastos_marinos_metadatos_PM_modif_clasificado.txt")


#### AT GENUS LEVEL ####



# Calcular la disimilitud de Bray-Curtis para cada status
bray_curtis_dist_genus <- vegan::vegdist(t(abundance_table_genus), method = "bray")

bray_curtis_pcoa_genus <- ecodist::pco(bray_curtis_dist_genus)
bray_curtis_pcoa_df_genus <- data.frame(pcoa1 = bray_curtis_pcoa_genus$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa_genus$vectors[,2])
dim(bray_curtis_pcoa_df_genus)
dim(metadata)

bray_curtis_pcoa_df_genus <- cbind(bray_curtis_pcoa_df_genus, metadata)

p1 <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = Temporada), size = 5) +
  stat_ellipse(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, fill = Temporada), geom = "polygon", alpha = 0.2) +
  labs(x = "PC1", y = "PC2", title = "Bray-Curtis PCoA Genera - Temporada") + 
  theme(title = element_text(size = 14, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text = element_text(size = 12))
p1


p2 <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = Z_m), size = 5) +
  stat_ellipse(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, fill = Z_m), geom = "polygon", alpha = 0.2) +
  labs(x = "PC1", y = "PC2", title = "Bray-Curtis PCoA - Z_m") + 
  theme(title = element_text(size = 14, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text = element_text(size = 12))
p2


p3 <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = T_C), size = 5) +
  stat_ellipse(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, fill = T_C), geom = "polygon", alpha = 0.2) +
  labs(x = "PC1", y = "PC2", title = "Bray-Curtis PCoA - T_C") + 
  theme(title = element_text(size = 14, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text = element_text(size = 10))
p3


p4 <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = N_P), size = 5) +
  stat_ellipse(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, fill = N_P), geom = "polygon", alpha = 0.2) +
  labs(x = "PC1", y = "PC2", title = "Bray-Curtis PCoA - NO3.NO2") + 
  theme(title = element_text(size = 14, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text = element_text(size = 12))
p4

p5 <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = Total_aliph), size = 5) +
  stat_ellipse(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, fill = Total_aliph), geom = "polygon", alpha = 0.2) +
  labs(x = "PC1", y = "PC2", title = "Bray-Curtis PCoA - Total_aliph") + 
  theme(title = element_text(size = 14, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text = element_text(size = 12))
p5

p6 <- ggplot() + theme_bw() +
  geom_point(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, color = Total_HAPs), size = 5) +
  stat_ellipse(data = bray_curtis_pcoa_df_genus, aes(x = pcoa1, y = pcoa2, fill = Total_HAPs), geom = "polygon", alpha = 0.2) +
  labs(x = "PC1", y = "PC2", title = "Bray-Curtis PCoA - Total_HAPs") + 
  theme(title = element_text(size = 14, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text = element_text(size = 12))
p6

(p1 + p2) / (p3 + p4)

svglite("Beta_Diversity_genus.svg", width=20, height=10)
  (p1 + p2 + p3) / (p4 + p5 + p6)
dev.off()







##################################
#### PHYSICOCHEMICAL FEATURES ####
##################################

physchem <- read.csv("physchem_v1.csv")
physchem <- column_to_rownames(physchem, var = "UID")
metadata_physchem <- read.csv("physchem_metadata.csv")

abund_arc_bac <- read.csv("C_archaea_bacteria_abd_genus.csv", row.names = 1)
# Filtrar filas que tengan al menos un valor >= 0.1
abund_arc_bac_filtered <- abund_arc_bac[apply(abund_arc_bac >= 0.1, 1, any), ]
# Verificar el resultado
print(dim(abund_arc_bac_filtered))  # Muestra el número de filas y columnas
head(abund_arc_bac_filtered)  # Muestra las primeras filas del dataframe filtrado


# Transponer la tabla de abundancia
physchem_table_t <- physchem
physchem_df <- as.data.frame(physchem_table_t)
physchem_df$UID <- rownames(physchem_df)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(physchem_df, metadata_physchem, by.x = "UID", by.y = "UID")

# Crear una lista de columnas de interés
columnas_de_interes <- c(
  3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

# Mapeo de colores personalizado
my_colors <- c(
  A_Inlet      = "#56b4e9",
  B_Transition = "#808080",
  C_Inner      = "#c8ab37",  
  Intense      = "#FF0000FF",
  Relax        = "#000080FF",
  Bare         = "orange",
  Seagrass     = "darkgreen"
)

# Crear las comparaciones adecuadas
a_my_comparisons <- list(
  c("A_Inlet", "B_Transition"),
  c("A_Inlet", "C_Inner")
)

# Bucle para generar y almacenar las figuras en variables a1, a2, a3, ...
for (i in 1:length(columnas_de_interes)) {
  columna_index <- columnas_de_interes[i]
  columna_de_interes <- colnames(combined_data)[columna_index]
  # Transformación a formato largo
  df_long <- combined_data %>% 
    dplyr::select(.data[[columna_de_interes]], Sector, Season, Habitat) %>%           
    pivot_longer(cols = c(Sector, Season, Habitat),
                 names_to  = "GroupType",
                 values_to = "Group")
  
  # Fijamos el orden de las siete categorías en el eje x
  level_order <- c("A_Inlet", "B_Transition", "C_Inner",
                   "Intense",  "Relax",
                   "Bare",     "Seagrass")
  df_long$Group <- factor(df_long$Group, levels = level_order)
  
  # Crear el boxplot con el t test usando la nueva columna combinada
  figura <- ggplot(df_long, aes(x = Group, y = .data[[columna_de_interes]], fill = Group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
    geom_jitter(width = 0.0, size = 0.8, alpha = 0.5, colour = "black") +
    scale_fill_manual(values = my_colors) +
    stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.format", 
                       position = position_dodge(0.75)) +
    labs(x = "", y = "", fill = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(columna_de_interes)) #+
  #ylim(0, 7.5)  # Fijar los límites del eje y
  
  # Almacenar la figura en una variable dinámica (a1, a2, a3, ...)
  assign(paste0("a", i), figura)
}


svglite("Fig_2.svg", width=20, height=10)
((a1 | a2 | a3 | a4 | a5) / (a6 | a7 | a8 | a9 | a10))
dev.off()

library(ARTool)

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Ajusta el modelo con los tres factores y todas sus interacciones
art_mod <- art(Sand ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
as.data.frame(anova(art_mod)) %>%
  tibble::rownames_to_column("Term") %>%
  mutate(
    label = paste0(
      Term, ": F(", Df, ",", Df.res, ")=", round(`F value`,2),
      ", p=", signif(`Pr(>F)`,2)
    )
  )

# Luego en ggplot, añadir con annotate() o patchwork::plot_annotation()
plot_annotation(
  title = "Sand abundance: ART( Sector×Season×Habitat )",
  subtitle = paste(anova_tab$label, collapse = "  |  ")
)


art_mod <- art(Silt ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(pH ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(TOC_umol ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(TN_umol ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(NH4_g ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(NO3_ug ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(NO2_ug ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(Fe_III_mg ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)
art_mod <- art(Fe_II_mg ~ Sector * Season * Habitat, data = combined_data)
anova(art_mod)

library(dplyr)
library(ARTool)

# 1) Define el vector de variables a analizar
vars <- c(
  "Sand" ,"Silt", "pH", "TOC_umol", "TN_umol",
  "NH4_g", "NO3_ug", "NO2_ug",
  "Fe_III_mg", "Fe_II_mg"
)

# 2) Función para extraer el ANOVA de un modelo ART
get_art_anova <- function(var, data) {
  # Ajusta el modelo
  mod <- art(
    formula = as.formula(paste(var, "~ Sector*Season*Habitat")),
    data    = data
  )
  # Extrae la tabla ANOVA III
  an <- anova(mod)
  # Convierte a data.frame y añade la respuesta
  df <- as.data.frame(an)
  df$Term     <- rownames(df)
  df$Response <- var
  rownames(df) <- NULL
  # Reordena columnas
  df %>%
    select(Response, Term, Df, Df.res, `F value`, `Pr(>F)`)
}

# 3) Itera sobre todas las variables y combina resultados
all_anovas <- bind_rows(
  lapply(vars, get_art_anova, data = combined_data)
)

# 4) Dale formato legible (opcional)
all_anovas <- all_anovas %>%
  mutate(
    `Pr(>F)` = signif(`Pr(>F)`, 3),
    `F value` = round(`F value`, 2)
  )

# 5) Muestra la tabla
print(all_anovas)

write.csv(all_anovas, "all_anova_art_physicochemical.csv")


library(FSA)
dunn_res <- dunnTest(Sand ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(Silt ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(pH ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(TOC_umol ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(TN_umol ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(NH4_g ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(NO3_ug ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(NO2_ug ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(Fe_III_mg ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res
dunn_res <- dunnTest(Fe_II_mg ~ Sector, data = combined_data, method = "bh")$res
# Filtrar comparaciones de interés
dunn_res



library(ggplot2)
library(dplyr)

# Calcular medias para “Sand”
means_df <- combined_data %>%
  group_by(Sector, Season, Habitat) %>%
  summarise(mu = mean(Sand), .groups="drop")

# Facet por Habitat, x = Season, color = Sector
ggplot(means_df, aes(x = Season, y = mu, color = Sector, group = Sector)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Habitat) +
  scale_color_manual(values = my_colors[c("A_Inlet","B_Transition","C_Inner")]) +
  labs(y = "Mean Sand", title = "Interacción Sector × Season por Habitat") +
  theme_minimal()

ggplot(means_df, aes(x = Sector, y = Season, fill = mu)) +
  geom_tile() +
  facet_wrap(~ Habitat) +
  scale_fill_viridis_c(name = "Mean Sand") +
  labs(title = "Sand: interacciones Sector–Season–Habitat") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggpubr)

ggboxplot(combined_data, x = "Sector", y = "Sand",
          color = "Sector", palette = my_colors[c("A_Inlet","B_Transition","C_Inner")]) +
  facet_grid(Season ~ Habitat) +
  stat_compare_means(aes(group = Sector), method = "wilcox.test",
                     label = "p.signif", paired = FALSE) +
  labs(title = "Sand by Sector × Season × Habitat")




## Estadisticas de cada propiedad segun factor

library(dplyr)

# Selecciona las columnas de interés
physicochem_cols <- names(combined_data)[3:12]

# Función para obtener las estadísticas
get_stats <- function(group_var) {
  combined_data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      across(
        all_of(physicochem_cols),
        list(
          mean   = ~ mean(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE),
          sd     = ~ sd(.x, na.rm = TRUE),   # ← desviación estándar
          min    = ~ min(.x, na.rm = TRUE),
          max    = ~ max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}


# Estadísticas por Sector
stats_by_sector <- get_stats("Sector")
print(t(stats_by_sector))
# Estadísticas por Season
stats_by_season <- get_stats("Season")
print(t(stats_by_season))

# Estadísticas por Habitat
stats_by_habitat <- get_stats("Habitat")
print(t(stats_by_habitat))



################################
#### FUNCTIONAL COMPARISONS ####
################################



# Definir el directorio
base_path <- "/Users/jorgerv/Documents/lagunas_costeras/Results_databases/Second_version/17-28/"

# Leer los archivos usando row.names para que la primera columna sean los nombres de las filas
total_cds_quintin <- read.csv(paste0(base_path, "17-28_count_cds.tsv"), sep = "\t", row.names = 1)
valid_samples <- rownames(total_cds_quintin)

## -------------------------
##
## Análisis con genes
##
## -------------------------

metal_genes <- read.csv(paste0(base_path, "count_BACMET_Metal_genes.csv"), sep = ",", row.names = 1)
drug_genes <- read.csv(paste0(base_path, "count_CARD_Drug_Class_genes.csv"), sep = ",", row.names = 1)
drug_class <- read.csv(paste0(base_path, "count_CARD_Drug_Class.csv"), sep = ",", row.names = 1)
nitrogen_genes <- read.csv(paste0(base_path, "count_NITROGEN_pathway_genes.csv"), sep = ",", row.names = 1)
sulfur_genes <- read.csv(paste0(base_path, "count_SULFUR_pathway_genes.csv"), sep = ",", row.names = 1)
virulence_genes <- read.csv(paste0(base_path, "count_VFDB_type_genes.csv"), sep = ",", row.names = 1)

# Función para filtrar las tablas de genes por las muestras en valid_samples
# y eliminar las columnas numéricas que sumen cero, preservando las filas por nombres
filtrar_y_eliminar_columnas_cero <- function(genes, valid_samples) {
  # Filtrar por muestras en valid_samples
  genes_filtered <- genes[rownames(genes) %in% valid_samples, ]
  # Eliminar las columnas cuyo conteo suma cero
  genes_filtered <- genes_filtered[, colSums(genes_filtered) != 0]
  
  return(genes_filtered)
}

# Aplicar la función a cada tabla de genes
metal_genes_filtered <- filtrar_y_eliminar_columnas_cero(metal_genes, valid_samples)
drug_genes_filtered <- filtrar_y_eliminar_columnas_cero(drug_genes, valid_samples)
drug_class_filtered <- filtrar_y_eliminar_columnas_cero(drug_class, valid_samples)
nitrogen_genes_filtered <- filtrar_y_eliminar_columnas_cero(nitrogen_genes, valid_samples)
sulfur_genes_filtered <- filtrar_y_eliminar_columnas_cero(sulfur_genes, valid_samples)
virulence_genes_filtered <- filtrar_y_eliminar_columnas_cero(virulence_genes, valid_samples)

# Renombrar la columna de 'total_cds_quintin' si es necesario
# colnames(total_cds_quintin)[1] <- "sample"  # No es necesario porque ya se usaron rownames

# Función para calcular abundancias con transformación logarítmica
calcular_abundancias_log <- function(genes_filtered, total_cds) {
  # Unir las tablas por las rownames
  genes_abund <- cbind(genes_filtered, Conteo_CDS = total_cds[rownames(genes_filtered), ])
  
  # Dividir los conteos por el total de CDS y aplicar la transformación logarítmica
  genes_abund <- genes_abund %>%
    dplyr::mutate(across(-Conteo_CDS, ~ ((.x / Conteo_CDS)*100))) %>%
    dplyr::select(-Conteo_CDS)
  
  # Eliminar las columnas cuya suma es cero
  genes_abund <- genes_abund %>%
    dplyr::select(where(~ sum(.) != 0))
  
  return(genes_abund)
}

# Aplicar la función a cada tabla de genes
metal_genes_abund <- calcular_abundancias_log(metal_genes_filtered, total_cds_quintin)
drug_genes_abund <- calcular_abundancias_log(drug_genes_filtered, total_cds_quintin)
nitrogen_genes_abund <- calcular_abundancias_log(nitrogen_genes_filtered, total_cds_quintin)
sulfur_genes_abund <- calcular_abundancias_log(sulfur_genes_filtered, total_cds_quintin)
virulence_genes_abund <- calcular_abundancias_log(virulence_genes_filtered, total_cds_quintin)

# Transpose of abundance tables
metal_genes_abund <- as.data.frame(t(metal_genes_abund))
drug_genes_abund <- as.data.frame(t(drug_genes_abund))
nitrogen_genes_abund <- as.data.frame(nitrogen_genes_abund)
sulfur_genes_abund <- as.data.frame(t(sulfur_genes_abund))
virulence_genes_abund <- as.data.frame(t(virulence_genes_abund))

save(metal_genes_abund, file = "genes_abund_metal.RData")
save(drug_genes_abund, file = "genes_abund_drug.RData")
save(nitrogen_genes_abund, file = "genes_abund_nitrogen.RData")
save(sulfur_genes_abund, file = "genes_abund_sulfur.RData")
save(virulence_genes_abund, file = "genes_abund_virulence.RData")

# Most abundant genes
# Calcular el porcentaje promedio de cada columna
metal_genes_abund$mean <- rowMeans(metal_genes_abund)
drug_genes_abund$mean <- rowMeans(drug_genes_abund)
nitrogen_genes_abund$mean <- rowMeans(nitrogen_genes_abund)
sulfur_genes_abund$mean <- rowMeans(sulfur_genes_abund)
virulence_genes_abund$mean <- rowMeans(virulence_genes_abund)


## --------------------------
##
## Análisis con pathways
##
## --------------------------

metal_pathways <- read.csv(paste0(base_path, "count_BACMET_Metal.csv"), sep = ",", row.names = 1)
drug_class <- read.csv(paste0(base_path, "count_CARD_Drug_Class.csv"), sep = ",", row.names = 1)
nitrogen_pathways <- read.csv(paste0(base_path, "count_NITROGEN_pathway.csv"), sep = ",", row.names = 1)
sulfur_pathways <- read.csv(paste0(base_path, "count_SULFUR_pathway.csv"), sep = ",", row.names = 1)
virulence_pathways <- read.csv(paste0(base_path, "count_VFDB_type.csv"), sep = ",", row.names = 1)


# Función para filtrar las tablas de pathways por las muestras en valid_samples
# y eliminar las columnas numéricas que sumen cero, preservando las filas por nombres
filtrar_y_eliminar_columnas_cero <- function(pathways, valid_samples) {
  # Filtrar por muestras en valid_samples
  pathways_filtered <- pathways[rownames(pathways) %in% valid_samples, ]
  # Eliminar las columnas cuyo conteo suma cero
  pathways_filtered <- pathways_filtered[, colSums(pathways_filtered) != 0]
  
  return(pathways_filtered)
}

# Aplicar la función a cada tabla de pathways
metal_pathways_filtered <- filtrar_y_eliminar_columnas_cero(metal_pathways, valid_samples)
drug_pathways_filtered <- filtrar_y_eliminar_columnas_cero(drug_class, valid_samples)
nitrogen_pathways_filtered <- filtrar_y_eliminar_columnas_cero(nitrogen_pathways, valid_samples)
sulfur_pathways_filtered <- filtrar_y_eliminar_columnas_cero(sulfur_pathways, valid_samples)
virulence_pathways_filtered <- filtrar_y_eliminar_columnas_cero(virulence_pathways, valid_samples)

# Renombrar la columna de 'total_cds_quintin' si es necesario
# colnames(total_cds_quintin)[1] <- "sample"  # No es necesario porque ya se usaron rownames

# Función para calcular abundancias con transformación logarítmica FUNCIONO
calcular_abundancias_log <- function(pathways_filtered, total_cds) {
  total_cds <- as.data.frame(total_cds)  # Asegurar que es un dataframe
  pathways_abund <- cbind(pathways_filtered, Conteo_CDS = total_cds[rownames(pathways_filtered), , drop = FALSE])
  
  # Comprobar si 'Conteo_CDS' está en el dataframe
  if (!"Conteo_CDS" %in% colnames(pathways_abund)) {
    stop("Error: 'Conteo_CDS' no está en pathways_abund. Revisa los rownames.")
  }
  
  print("Columnas antes de select:")  
  print(colnames(pathways_abund))  # Debugging
  
  pathways_abund <- pathways_abund %>%
    mutate(across(-Conteo_CDS, ~ ((.x / Conteo_CDS) * 100))) %>%
    dplyr::select(-Conteo_CDS)  # Asegurar que se usa el select de dplyr
  
  pathways_abund <- pathways_abund %>%
    dplyr::select(where(~ sum(.) != 0))
  
  return(pathways_abund)
}


# Aplicar la función a cada tabla de pathways
metal_pathways_abund <- calcular_abundancias_log(metal_pathways_filtered, total_cds_quintin)
drug_pathways_abund <- calcular_abundancias_log(drug_pathways_filtered, total_cds_quintin)
nitrogen_pathways_abund <- calcular_abundancias_log(nitrogen_pathways_filtered, total_cds_quintin)
sulfur_pathways_abund <- calcular_abundancias_log(sulfur_pathways_filtered, total_cds_quintin)
virulence_pathways_abund <- calcular_abundancias_log(virulence_pathways_filtered, total_cds_quintin)

# Transform to data frames the abundance tables
metal_pathways_abund <- as.data.frame(metal_pathways_abund)
drug_pathways_abund <- as.data.frame(drug_pathways_abund)
nitrogen_pathways_abund <- as.data.frame(nitrogen_pathways_abund)
sulfur_pathways_abund <- as.data.frame(sulfur_pathways_abund)
virulence_pathways_abund <- as.data.frame(virulence_pathways_abund)

summary(metal_pathways_abund)
summary(drug_pathways_abund)
summary(nitrogen_pathways_abund)
summary(sulfur_pathways_abund)
summary(virulence_pathways_abund)

# Save data
save(drug_pathways_abund, file = "pathways_abund_drug.RData")
save(metal_pathways_abund, file = "pathways_abund_metal.RData")
save(nitrogen_pathways_abund, file = "pathways_abund_nitrogen.RData")
save(sulfur_pathways_abund, file = "pathways_abund_sulfur.RData")
save(virulence_pathways_abund, file = "pathways_abund_virulence.RData")


load("pathways_abund_metal.RData")
load("pathways_abund_virulence.RData")

summary(metal_pathways_abund)
mean(as.numeric(as.matrix(metal_pathways_abund)))
max(as.numeric(as.matrix(metal_pathways_abund))) - mean(as.numeric(as.matrix(metal_pathways_abund)))
summary(virulence_pathways_abund)
mean(as.numeric(as.matrix(virulence_pathways_abund)))
max(as.numeric(as.matrix(virulence_pathways_abund))) - mean(as.numeric(as.matrix(virulence_pathways_abund)))

## -----------------------------
## 
## BOX PLOTS ABUNDANCE PATHWAYS
##
## -----------------------------


library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
load("pathways_abund_drug.RData")
load("pathways_abund_metal.RData")
load("pathways_abund_nitrogen.RData")
load("pathways_abund_sulfur.RData")
load("pathways_abund_virulence.RData")

metadata <- read.csv("metadata_quintin.csv")

# Transponer la tabla de abundancia
drug_pathways_abund$Sample <- rownames(drug_pathways_abund)
metal_pathways_abund$Sample <- rownames(metal_pathways_abund)
nitrogen_pathways_abund$Sample <- rownames(nitrogen_pathways_abund)
sulfur_pathways_abund$Sample <- rownames(sulfur_pathways_abund)
virulence_pathways_abund$Sample <- rownames(virulence_pathways_abund)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")

# Mapeo de colores personalizado
my_colors <- c(
  A_Inlet      = "#56b4e9",
  B_Transition = "#808080",
  C_Inner      = "#c8ab37",  
  Intense      = "#FF0000FF",
  Relax        = "#000080FF",
  Bare         = "orange",
  Seagrass     = "darkgreen"
)


# Crear una lista de columnas de interés
#Nitrogen
columnas_de_interes <- c(2, 3, 4, 5, 6, 7, 8)
# Metals
columnas_de_interes <- c(2,3,4,8,9,10,11,13,14,15,17,18,19,20,24)
# Sulfur
columnas_de_interes <- c(2, 3, 4, 5, 6, 7, 8)
# Drug
columnas_de_interes <- c(2:28)
# Virulence
columnas_de_interes <- c(2:17)

# Crear las comparaciones adecuadas
a_my_comparisons <- list(
  c("A_Inlet", "B_Transition"),
  c("A_Inlet", "C_Inner")
)

# Bucle para generar y almacenar las figuras en variables a1, a2, a3, ...
for (i in 1:length(columnas_de_interes)) {
  columna_index <- columnas_de_interes[i]
  columna_de_interes <- colnames(combined_data)[columna_index]
  # Transformación a formato largo
  df_long <- combined_data %>% 
    dplyr::select(.data[[columna_de_interes]], Sector, Season, Habitat) %>%           
    pivot_longer(cols = c(Sector, Season, Habitat),
                 names_to  = "GroupType",
                 values_to = "Group")
  
  # Fijamos el orden de las siete categorías en el eje x
  level_order <- c("A_Inlet", "B_Transition", "C_Inner",
                   "Intense",  "Relax",
                   "Bare",     "Seagrass")
  df_long$Group <- factor(df_long$Group, levels = level_order)
  
  # Crear el boxplot con el t test usando la nueva columna combinada
  figura <- ggplot(df_long, aes(x = Group, y = .data[[columna_de_interes]], fill = Group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
    geom_jitter(width = 0.0, size = 0.8, alpha = 0.5, colour = "black") +
    scale_fill_manual(values = my_colors) +
    stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.format", 
                       position = position_dodge(0.75)) +
    labs(x = "", y = "", fill = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(columna_de_interes)) #+
  #ylim(0, 7.5)  # Fijar los límites del eje y
  
  # Almacenar la figura en una variable dinámica (a1, a2, a3, ...)
  assign(paste0("a", i), figura)
}


svglite("annot_nitrogen_metabolism.svg", width=29, height=4)
((a1 | a2 | a3 | a4 | a5 | a6 | a7))
dev.off()

svglite("annot_metal_metabolism.svg", width=29, height=8)
((a1 | a2 | a3 | a4 | a5 | a6 | a7) / (a8 | a9 | a10| a11 | a12 | a13 | a14))
dev.off()

svglite("annot_sulfur_metabolism.svg", width=29, height=4)
((a1 | a2 | a3 | a4 | a5 | a6 | a7))
dev.off()

svglite("annot_drug_metabolism.svg", width=29, height=16)
((a1 | a2 | a3 | a4 | a5 | a6 | a7) / (a8 | a9 | a10| a11 | a12 | a13 | a14) / (a15 | a16 | a17| a18 | a19 | a20 | a21) / (a22 | a23 | a24 | a25 | a26 | a1 | a2))
dev.off()

svglite("annot_virulence_metabolism.svg", width=24.7, height=12)
((a1 | a2 | a3 | a4 | a5 | a6) / (a7 | a8 | a9 | a10 | a11 | a12) / (a13 | a14 | a15 | a16 | a1 | a12))
dev.off()




### ---------------------------
## Diferencias entre pathways
## ----------------------------


# Cargar paquetes necesarios
library(dplyr)
library(tidyr)
library(FSA)        # para dunnTest()
library(rstatix)    # para ajustar p-values si se prefiere

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# 1) Columnas de abundancia
pathways_cols <- combined_data %>%
  select(-Sample, -Sector, -Season, -Habitat, -ID_1, -ID_2) %>%
  names()

# 2) Calcular KW p-values para cada pathway
kw_results <- lapply(pathways_cols, function(var){
  kw <- kruskal.test(combined_data[[var]] ~ combined_data$Sector)
  data.frame(
    Pathway = var,
    P.value = kw$p.value
  )
}) %>%
  bind_rows() %>%
  arrange(P.value)

# 3) Kruskal–Wallis + Dunn para Sector
kw_dunn_results <- lapply(pathways_cols, function(var){
  # 1) Kruskal–Wallis
  kw <- kruskal.test(combined_data[[var]] ~ combined_data$Sector)
  
  # 2) Si KW sale p ≤ 0.05, lanza Dunn
  if(kw$p.value <= 0.1){
    d <- dunnTest(
      x      = combined_data[[var]],
      g      = combined_data$Sector,
      method = "bh"
    )$res
    
    d %>%
      transmute(
        Pathway    = var,
        Comparison = Comparison,
        P.unadj    = P.unadj,
        P.adj      = P.adj
      )
  } else {
    NULL
  }
}) %>%
  bind_rows()

# 4) Wilcoxon para Season y Habitat
wilcox_results <- lapply(c("Season","Habitat"), function(var){
  lapply(pathways_cols, function(col){
    wt <- wilcox.test(
      combined_data[[col]] ~ combined_data[[var]],
      data = combined_data,
      exact = FALSE
    )
    data.frame(
      Pathway  = col,
      Variable = var,
      P.value  = wt$p.value
    )
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  # corrección BH por grupo de Variable
  group_by(Variable) %>%
  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%
  ungroup()

# 5) Mostrar todos los resultados
kw_results
print(wilcox_results, n = 100)

# 6) Filtrar resultados “tendencia” 0.05 < p ≤ 0.1
kw_dunn_sig   <- kw_dunn_results   %>% filter(P.adj <= 0.1)
wilcox_sig    <- wilcox_results    %>% filter(P.adj <= 0.1)

# 7) Listado final
list(
  Dunn_Sector     = kw_dunn_sig,
  Wilcox_SeasonHabitat = wilcox_sig
)








## ----------------------------
## 
## HEATMAPS GENES
##
## ----------------------------


library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
load("genes_abund_metal.RData")
load("genes_abund_drug.RData")
load("genes_abund_nitrogen.RData")
load("genes_abund_sulfur.RData")
load("genes_abund_virulence.RData")

metadata <- read.csv("metadata_quintin.csv")

# Función para convertir matriz en formato largo, unir con metadata, y calcular medias
procesar_genes <- function(genes_abund, metadata) {
  # Convertir la matriz en data.frame
  genes_abund_df <- as.data.frame(t(genes_abund))
  
  # Convertir al formato largo
  genes_long <- genes_abund_df %>%
    rownames_to_column(var = "Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance")
  
  # Unir con metadata
  genes_metadata <- genes_long %>%
    left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calcular medias por Sector
  sector_means <- genes_metadata %>%
    group_by(Gene, Sector) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Sector, values_from = mean_abundance)
  
  # Calcular medias por Season
  season_means <- genes_metadata %>%
    group_by(Gene, Season) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Season, values_from = mean_abundance)
  
  # Calcular medias por Habitat
  habitat_means <- genes_metadata %>%
    group_by(Gene, Habitat) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Habitat, values_from = mean_abundance)
  
  # Unir todas las medias en un único dataframe
  genes_means <- genes_abund_df %>%
    rownames_to_column(var = "Gene") %>%
    left_join(sector_means, by = "Gene") %>%
    left_join(season_means, by = "Gene") %>%
    left_join(habitat_means, by = "Gene")
  
  return(genes_means)
}

# Aplicar la función a cada una de las matrices
metal_means <- procesar_genes(metal_genes_abund, metadata)
drug_means <- procesar_genes(drug_genes_abund, metadata)
nitrogen_means <- procesar_genes(nitrogen_genes_abund, metadata)
sulfur_means <- procesar_genes(sulfur_genes_abund, metadata)
virulence_means <- procesar_genes(virulence_genes_abund, metadata)

# Calculo los valores máximos
#max_value_drug <- max(drug_means[ , (ncol(drug_means) - 6):ncol(drug_means)], na.rm = TRUE)
max_value_nitrogen <- max(nitrogen_means[ , (ncol(nitrogen_means) - 6):ncol(nitrogen_means)], na.rm = TRUE)
max_value_sulfur <- max(sulfur_means[ , (ncol(sulfur_means) - 6):ncol(sulfur_means)], na.rm = TRUE)

# Función para generar los heatmaps de un archivo *_means
crear_heatmaps <- function(means_data, max_val) {
  
  # Convertir a formato largo para Sector, Season y Habitat
  sector_data <- means_data %>%
    select(Gene, starts_with("A_"), starts_with("B_"), starts_with("C_In")) %>%
    pivot_longer(-Gene, names_to = "Sector", values_to = "Abundance")
  
  season_data <- means_data %>%
    select(Gene, starts_with("Intense"), starts_with("Relax")) %>%
    pivot_longer(-Gene, names_to = "Season", values_to = "Abundance")
  
  habitat_data <- means_data %>%
    select(Gene, starts_with("Bare"), starts_with("Seagrass")) %>%
    pivot_longer(-Gene, names_to = "Habitat", values_to = "Abundance")
  
  # Crear heatmap para Sector
  heatmap_sector <- ggplot(sector_data, aes(x = Sector, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Season
  heatmap_season <- ggplot(season_data, aes(x = Season, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los genes en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Habitat
  heatmap_habitat <- ggplot(habitat_data, aes(x = Habitat, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"), 
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los genes en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Unir los tres heatmaps en un solo gráfico utilizando patchwork
  library(patchwork)
  
  # Combinar los tres heatmaps en una figura
  combined_plot <- (heatmap_sector | heatmap_season | heatmap_habitat) + 
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

# Generar heatmaps para cada archivo *_means
#heatmap_drug <- crear_heatmaps(drug_means, max_value_drug)
heatmap_nitrogen <- crear_heatmaps(nitrogen_means, max_value_nitrogen)
heatmap_sulfur <- crear_heatmaps(sulfur_means, max_value_sulfur)

# Mostrar el plot de cada uno
#print(heatmap_drug)
print(heatmap_nitrogen)
print(heatmap_sulfur)

# Definir el número de genes en cada dataset
num_genes_drug <- nrow(drug_means)
num_genes_nitrogen <- nrow(nitrogen_means)
num_genes_sulfur <- nrow(sulfur_means)

# Ajustar la altura dinámicamente
height_factor <- 0.3  # Ajusta este valor según el espaciado que necesites

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_CARD_genes.svg", width = 8, height = max(10, num_genes_drug * height_factor))
heatmap_drug
dev.off()

svglite("heatmap_NITROGEN_genes.svg", width = 8, height = max(8, num_genes_nitrogen * height_factor))
heatmap_nitrogen
dev.off()

svglite("heatmap_SULFUR_genes.svg", width = 8, height = max(8, num_genes_sulfur * height_factor))
heatmap_sulfur
dev.off()


## ----------------------------
## Diferencias entre genes
## ----------------------------


# Cargar paquetes necesarios
library(FSA)      # Para Dunn test
library(dplyr)    # Para manipulación de datos
library(tidyr)    # Para transformar datos
library(rstatix)  # Para múltiples comparaciones ajustadas

# Transponer la tabla de abundancia
drug_genes_abund$Sample <- rownames(drug_genes_abund)
metal_genes_abund$Sample <- rownames(metal_genes_abund)
nitrogen_genes_abund$Sample <- rownames(nitrogen_genes_abund)
sulfur_genes_abund$Sample <- rownames(sulfur_genes_abund)
virulence_genes_abund$Sample <- rownames(virulence_genes_abund)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_genes_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Seleccionar solo las columnas numéricas (excluyendo "Sample", "Sector", "Season", "Habitat")
genes_cols <- combined_data %>%
  select(-Sample, -Sector, -Season, -Habitat, -ID_1, -ID_2)

# Función para aplicar Dunn test (Sector: 3 niveles)
dunn_results <- lapply(names(genes_cols), function(gene) {
  test <- dunnTest(combined_data[[gene]] ~ combined_data$Sector, method = "bh")
  result <- data.frame(Gene = gene, test$res)
  return(result)
}) %>%
  bind_rows() %>%
  filter(P.adj <= 0.1)  # Filtrar p-values ajustados < 0.05

# Función para aplicar Wilcoxon test (Season y Habitat: 2 niveles)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(genes_cols), function(gene) {
    test <- wilcox.test(combined_data[[gene]] ~ combined_data[[variable]])
    data.frame(Gene = gene, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  filter(p.value <= 0.1)  # Filtrar p-values <= 0.05

# Función para aplicar Wilcoxon test (Season y Habitat: 2 niveles)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(pathways_cols), function(pathway) {
    test <- wilcox.test(combined_data[[pathway]] ~ combined_data[[variable]])
    data.frame(Pathway = pathway, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) 


# Asumiendo que wilcox_results es una lista de dos data.frames
for(i in seq_along(wilcox_results)) {
  wilcox_results[[i]]$p.adj <- p.adjust(
    wilcox_results[[i]]$p.value,
    method = "BH"
  )
}

wilcox_results <- wilcox_results %>%
  bind_rows() %>%
  filter(p.adj <= 0.1)  # Filtrar p-adj <= 0.1

# Mostrar los resultados significativos
list(Dunn_Test_Significativo = dunn_results, Wilcoxon_Significativo = wilcox_results)

## NITROGEN

# $Dunn_Test_Significativo
# Gene        Comparison         Z    P.unadj      P.adj
# 1 narC A_Inlet - C_Inner -2.157277 0.03098405 0.09295215

# $Wilcoxon_Significativo
# Gene Variable     p.value
# 1      amoA_A   Season 0.092125104
# 2        napB   Season 0.064935065
# 3        napC   Season 0.093073593
# 4        narJ   Season 0.064935065
# 5        narZ   Season 0.015151515
# 6        narY   Season 0.002164502
# 7        nasB   Season 0.025974026
# 8        narC   Season 0.093073593
# 9        nirD   Season 0.008658009
# 10       nrfC   Season 0.064935065
# 11 gdh_K00261   Season 0.064935065
# 12 gdh_K00260   Season 0.028440662
# 13     amoA_B  Habitat 0.074009970
# 14        hao  Habitat 0.093073593
# 15       nrfD  Habitat 0.074009970


## SULFUR

# $Dunn_Test_Significativo
# Gene             Comparison         Z     P.unadj      P.adj
# 1   sir      A_Inlet - C_Inner  2.157277 0.030984050 0.09295215
# 2  phsA      A_Inlet - C_Inner -2.255336 0.024112275 0.07233682
# 3  phsC B_Transition - C_Inner  2.549510 0.010787449 0.03236235 *
# 4  comC B_Transition - C_Inner -2.451452 0.014228128 0.04268439 *
# 5  dddD      A_Inlet - C_Inner -2.745626 0.006039559 0.01811868 *
# 6  dddT A_Inlet - B_Transition  2.059219 0.039473224 0.05920984
# 7  dmoA B_Transition - C_Inner -2.451452 0.014228128 0.04268439 *
# 8  hpsO      A_Inlet - C_Inner -2.353394 0.018602930 0.05580879
# 9  mddA A_Inlet - B_Transition -2.059219 0.039473224 0.05920984
# 10 mddA      A_Inlet - C_Inner -2.647568 0.008107310 0.02432193 *
# 11  mdh A_Inlet - B_Transition  2.353394 0.018602930 0.05580879
# 12 msmB      A_Inlet - C_Inner -2.135231 0.032742163 0.09822649

# $Wilcoxon_Significativo
# Gene Variable     p.value
# 1   cysJ   Season 0.002164502 *
# 2   cysN   Season 0.015151515 *
# 3   dsrA   Season 0.041125541 *
# 4   dsrB   Season 0.093073593
# 5   dsrD   Season 0.015151515 *
# 6   dsrE   Season 0.015151515 *
# 7   dsrF   Season 0.008658009 *
# 8   dsrJ   Season 0.015151515 *
# 9   dsrT   Season 0.015151515 *
# 10  rdsr   Season 0.064935065
# 11  asrB   Season 0.041125541 *
# 12  shyC   Season 0.025974026 *
# 13  shyD   Season 0.041125541 *
# 14  phsC   Season 0.064935065
# 15  acuI   Season 0.093073593
# 16  acuN   Season 0.064935065
# 17  comD   Season 0.093073593
# 18  comE   Season 0.041125541 *
# 19  dddP   Season 0.093073593
# 20  dmdA   Season 0.092125104
# 21  dmdD   Season 0.002164502 *
# 22  hpsP   Season 0.041125541 *
# 23  isfD   Season 0.025974026 *
# 24  cuyA   Season 0.041125541 *
# 25  cysK   Season 0.093073593
# 26 hdrA1   Season 0.093073593
# 27  mccB   Season 0.041125541 *
# 28  metZ   Season 0.093073593
# 29  qmoC   Season 0.064935065
# 30  iseJ   Season 0.074009970
# 31  psrC   Season 0.028440662 *
# 32  rdlA   Season 0.074009970
# 33  dsrO  Habitat 0.064935065
# 34  phsB  Habitat 0.064935065



## ----------------------------
## 
## HEATMAPS PATHWAYS
##
## ----------------------------


library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
load("pathways_abund_metal.RData")
load("pathways_abund_drug.RData")
load("pathways_abund_nitrogen.RData")
load("pathways_abund_sulfur.RData")
load("pathways_abund_virulence.RData")

metadata <- read.csv("metadata_quintin.csv")

# Función para convertir matriz en formato largo, unir con metadata, y calcular medias
procesar_pathways <- function(pathways_abund, metadata) {
  # Convertir la matriz en data.frame
  pathways_abund_df <- as.data.frame(t(pathways_abund))
  
  # Convertir al formato largo
  pathways_long <- pathways_abund_df %>%
    rownames_to_column(var = "Pathway") %>%
    pivot_longer(-Pathway, names_to = "Sample", values_to = "Abundance")
  
  # Unir con metadata
  pathways_metadata <- pathways_long %>%
    left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calcular medias por Sector
  sector_means <- pathways_metadata %>%
    group_by(Pathway, Sector) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Sector, values_from = mean_abundance)
  
  # Calcular medias por Season
  season_means <- pathways_metadata %>%
    group_by(Pathway, Season) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Season, values_from = mean_abundance)
  
  # Calcular medias por Habitat
  habitat_means <- pathways_metadata %>%
    group_by(Pathway, Habitat) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Habitat, values_from = mean_abundance)
  
  # Unir todas las medias en un único dataframe
  pathways_means <- pathways_abund_df %>%
    rownames_to_column(var = "Pathway") %>%
    left_join(sector_means, by = "Pathway") %>%
    left_join(season_means, by = "Pathway") %>%
    left_join(habitat_means, by = "Pathway")
  
  return(pathways_means)
}

# Aplicar la función a cada una de las matrices
metal_means <- procesar_pathways(metal_pathways_abund, metadata)
drug_means <- procesar_pathways(drug_pathways_abund, metadata)
#nitrogen_means <- procesar_pathways(nitrogen_pathways_abund, metadata)
#sulfur_means <- procesar_pathways(sulfur_pathways_abund, metadata)
virulence_means <- procesar_pathways(virulence_pathways_abund, metadata)

# Calculo los valores máximos
max_value_metal <- max(metal_means[ , (ncol(metal_means) - 6):ncol(metal_means)], na.rm = TRUE)
max_value_drug <- max(drug_means[ , (ncol(drug_means) - 6):ncol(drug_means)], na.rm = TRUE)
#max_value_nitrogen <- max(nitrogen_means[ , (ncol(nitrogen_means) - 6):ncol(nitrogen_means)], na.rm = TRUE)
#max_value_sulfur <- max(sulfur_means[ , (ncol(sulfur_means) - 6):ncol(sulfur_means)], na.rm = TRUE)
max_value_virulence <- max(virulence_means[ , (ncol(virulence_means) - 6):ncol(virulence_means)], na.rm = TRUE)

# Función para Pathways crear los heatmaps de un archivo *_means
crear_heatmaps <- function(means_data, max_val) {
  
  # Convertir a formato largo para Sector, Season y Habitat
  sector_data <- means_data %>%
    select(Pathway, starts_with("A_"), starts_with("B_"), starts_with("C_In")) %>%
    pivot_longer(-Pathway, names_to = "Sector", values_to = "Abundance")
  
  season_data <- means_data %>%
    select(Pathway, starts_with("Intense"), starts_with("Relax")) %>%
    pivot_longer(-Pathway, names_to = "Season", values_to = "Abundance")
  
  habitat_data <- means_data %>%
    select(Pathway, starts_with("Bare"), starts_with("Seagrass")) %>%
    pivot_longer(-Pathway, names_to = "Habitat", values_to = "Abundance")
  
  # Crear heatmap para Sector
  heatmap_sector <- ggplot(sector_data, aes(x = Sector, y = Pathway, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Season
  heatmap_season <- ggplot(season_data, aes(x = Season, y = Pathway, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los pathways en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Habitat
  heatmap_habitat <- ggplot(habitat_data, aes(x = Habitat, y = Pathway, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los pathways en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Unir los tres heatmaps en un solo gráfico utilizando patchwork
  library(patchwork)
  
  # Combinar los tres heatmaps en una figura
  combined_plot <- (heatmap_sector | heatmap_season | heatmap_habitat) + 
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

# Generar heatmaps para cada archivo *_means
heatmap_drug <- crear_heatmaps(drug_means, max_value_drug)
#heatmap_nitrogen <- crear_heatmaps(nitrogen_means, max_value_nitrogen)
#heatmap_sulfur <- crear_heatmaps(sulfur_means, max_value_sulfur)
heatmap_metal <- crear_heatmaps(metal_means, max_value_metal)
heatmap_virulence <- crear_heatmaps(virulence_means, max_value_virulence)

# Mostrar el plot de cada uno
print(heatmap_drug)
#print(heatmap_nitrogen)
#print(heatmap_sulfur)
print(heatmap_metal)
print(heatmap_virulence)

# Definir el número de genes en cada dataset
num_pathways_drug <- nrow(drug_means)
#num_pathways_nitrogen <- nrow(nitrogen_means)
#num_pathways_sulfur <- nrow(sulfur_means)
num_pathways_metal <- nrow(metal_means)
num_pathways_virulence <- nrow(virulence_means)

# Ajustar la altura dinámicamente
height_factor <- 0.3  # Ajusta este valor según el espaciado que necesites

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_CARD_pathway3.svg", width = 8, height = max(10, num_pathways_drug * height_factor))
heatmap_drug
dev.off()

svglite("heatmap_NITROGEN_pathway3.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_nitrogen
dev.off()

svglite("heatmap_SULFUR_pathway3.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_sulfur
dev.off()

svglite("heatmap_METAL_pathway3.svg", width = 8, height = max(8, num_pathways_metal * height_factor))
heatmap_metal
dev.off()

svglite("heatmap_VFDB_pathway3.svg", width = 8, height = max(8, num_pathways_virulence * height_factor))
heatmap_virulence
dev.off()


## ----------------------------
## Diferencias entre pathways
## ----------------------------


# Cargar paquetes necesarios
library(FSA)      # Para Dunn test
library(dplyr)    # Para manipulación de datos
library(tidyr)    # Para transformar datos
library(rstatix)  # Para múltiples comparaciones ajustadas

# Transponer la tabla de abundancia
drug_pathways_abund$Sample <- rownames(drug_pathways_abund)
metal_pathways_abund$Sample <- rownames(metal_pathways_abund)
#nitrogen_pathways_abund$Sample <- rownames(nitrogen_pathways_abund)
#sulfur_pathways_abund$Sample <- rownames(sulfur_pathways_abund)
virulence_pathways_abund$Sample <- rownames(virulence_pathways_abund)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
#combined_data <- merge(nitrogen_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
#combined_data <- merge(sulfur_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Seleccionar solo las columnas numéricas (excluyendo "Sample", "Sector", "Season", "Habitat")
pathways_cols <- combined_data %>%
  select(-Sample, -Sector, -Season, -Habitat, -ID_1, -ID_2)

# Función para aplicar Dunn test (Sector: 3 niveles)
dunn_results <- lapply(names(pathways_cols), function(pathway) {
  test <- dunnTest(combined_data[[pathway]] ~ combined_data$Sector, method = "bh")
  result <- data.frame(Pathway = pathway, test$res)
  return(result)
}) %>%
  bind_rows() %>%
  filter(P.adj <= 0.1)  # Filtrar p-values ajustados < 0.05

library(dplyr)

# Función para aplicar Wilcoxon test (Season y Habitat: 2 niveles)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(pathways_cols), function(pathway) {
    test <- wilcox.test(combined_data[[pathway]] ~ combined_data[[variable]])
    data.frame(Pathway = pathway, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) 


# Asumiendo que wilcox_results es una lista de dos data.frames
for(i in seq_along(wilcox_results)) {
  wilcox_results[[i]]$p.adj <- p.adjust(
    wilcox_results[[i]]$p.value,
    method = "BH"
  )
}

wilcox_results <- wilcox_results %>%
  bind_rows() %>%
  filter(p.adj <= 0.1)  # Filtrar p-adj <= 0.1

# Mostrar los resultados significativos
list(Dunn_Test_Significativo = dunn_results, Wilcoxon_Significativo = wilcox_results)


## ----------------------------
## 
## COMPLETENESS PATHWAYS
##
## ----------------------------

# Load data
load("genes_abund_nitrogen.RData")
load("genes_abund_sulfur.RData")

metadata <- read.csv("metadata_quintin.csv")

# Transponer la tabla de abundancia
nitrogen_genes_abund$Sample <- rownames(nitrogen_genes_abund)
sulfur_genes_abund$Sample <- rownames(sulfur_genes_abund)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(nitrogen_genes_abund, metadata, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_genes_abund, metadata, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

## Para NITROGEN

# 1) Definir la lista de genes por pathway, incluyendo pasajes alternativos
genes_list <- list(
  Anammox = c("hzo", "hzsA", "hzsB", "hzsC", "hdh"),
  
  Nitrogen_fixation = c("anfG", "nifD", "nifH", "nifK", "nifW"),
  
  Nitrification = list(
    groupA = c("amoA_A","amoB_A","amoC_A"),
    groupB = c("amoA_B","amoB_B","amoC_B"),
    common = c("hao","nxrA","nxrB")
  ),
  
  Denitrification = list(
    step1 = list(
      c("napA","napB","napC"),
      c("narH","narJ","narI","narG"),
      c("narV","narW","narY","narZ")
    ),
    step2 = c("nirK","nirS"),
    step3 = c("norB","norC","norZ"),
    step4 = "nosZ"
  ),
  
  Assimilatory_nitrate_reduction = list(
    step1 = list(
      c("narB","narC"),
      "NR",
      c("nasA","nasB")
    ),
    step2 = "nirA"
  ),
  
  Dissimilatory_nitrate_reduction = list(
    step1 = list(
      c("napA","napB","napC"),
      c("narH","narJ","narI","narG"),
      c("narV","narW","narY","narZ")
    ),
    step2 = list(
      c("nirB","nirD"),
      c("nrfA","nrfB","nrfC","nrfD")
    )
  ),
  
  Organic_degradation_and_synthesis = list(
    route1 = c("glsA","glnA","asnB","ansB",
               "gs_K00264","gs_K00265","gs_K00266","gs_K00284"),
    route2 = c("ureA","ureB","ureC","nao","nmo",
               "gdh_K00260","gdh_K00261","gdh_K00262","gdh_K15371")
  )
)

# 2) Cargar datos y unir abundancias con metadata
load("genes_abund_nitrogen.RData")  # debe dejar 'nitrogen_genes_abund' en el entorno
metadata <- read.csv("metadata_quintin.csv")

# Transponer y unir
nitrogen_genes_abund$Sample <- rownames(nitrogen_genes_abund)
combined_data <- merge(nitrogen_genes_abund, metadata, by = "Sample")

# 3) Crear data.frame de presencia/ausencia (0/1) para **todos** los genes 
#    Si un gen no existe en combined_data, se considera 0
desired_genes <- unlist(
  lapply(genes_list, function(x) if(is.list(x)) unlist(x) else x)
) %>% unique()

present_genes <- intersect(desired_genes, names(combined_data))
missing_genes <- setdiff(desired_genes, names(combined_data))

# Construir presence_df con una columna Sample + todos los genes
presence_df <- combined_data %>%
  select(Sample) %>%
  bind_cols(
    combined_data %>%
      select(all_of(present_genes)) %>%
      mutate(across(everything(), ~ as.integer(. > 0)))
  )

# Añadir los genes faltantes como columnas de ceros
for(g in missing_genes){
  presence_df[[g]] <- 0L
}

# Reordenar columnas: Sample + desired_genes
presence_df <- presence_df %>%
  select(Sample, all_of(desired_genes))

# 4) Función auxiliar para pasos OR en pipelines alternativas
step_present <- function(df_row, groups){
  # df_row: vector numérico de 0/1
  # groups: lista de vectores de nombres
  any(sapply(groups, function(g) any(df_row[g] == 1))) * 1
}

# 5) Calcular completitud de cada pathway (%) por muestra
completeness_df <- presence_df %>%
  rowwise() %>%
  mutate(
    Anammox = sum(c_across(all_of(genes_list$Anammox))) / 5 * 100,
    
    Nitrogen_fixation = if_else(
      c_across("anfG") == 1,
      100,
      sum(c_across(setdiff(genes_list$Nitrogen_fixation, "anfG"))) / 4 * 100
    ),
    
    Nitrification = {
      cntA <- sum(c_across(genes_list$Nitrification$groupA))
      cntB <- sum(c_across(genes_list$Nitrification$groupB))
      cntC <- sum(c_across(genes_list$Nitrification$common))
      (max(cntA, cntB) + cntC) / 6 * 100
    },
    
    Denitrification = {
      s1 <- step_present(cur_data(), genes_list$Denitrification$step1)
      s2 <- any(c_across(genes_list$Denitrification$step2)) * 1
      s3 <- any(c_across(genes_list$Denitrification$step3)) * 1
      s4 <- c_across(genes_list$Denitrification$step4)
      (s1 + s2 + s3 + s4) / 4 * 100
    },
    
    Assimilatory_nitrate_reduction = {
      s1 <- step_present(cur_data(), genes_list$Assimilatory_nitrate_reduction$step1)
      s2 <- c_across(genes_list$Assimilatory_nitrate_reduction$step2)
      (s1 + s2) / 2 * 100
    },
    
    Dissimilatory_nitrate_reduction = {
      s1 <- step_present(cur_data(), genes_list$Dissimilatory_nitrate_reduction$step1)
      s2 <- step_present(cur_data(), genes_list$Dissimilatory_nitrate_reduction$step2)
      (s1 + s2) / 2 * 100
    },
    
    Organic_degradation_and_synthesis = {
      c1 <- sum(c_across(genes_list$Organic_degradation_and_synthesis$route1)) /
        length(genes_list$Organic_degradation_and_synthesis$route1) * 100
      c2 <- sum(c_across(genes_list$Organic_degradation_and_synthesis$route2)) /
        length(genes_list$Organic_degradation_and_synthesis$route2) * 100
      max(c1, c2)
    }
  ) %>%
  ungroup() %>%
  select(
    Sample,
    Anammox,
    Nitrogen_fixation,
    Nitrification,
    Denitrification,
    Assimilatory_nitrate_reduction,
    Dissimilatory_nitrate_reduction,
    Organic_degradation_and_synthesis
  )

# 6) Ver resultados
head(completeness_df)
# (Opcional) unirlo a combined_data:
combined_data <- combined_data %>% left_join(completeness_df, by = "Sample")

## Calculo las medias para el completeness

# 2) Función para convertir tu data.frame de completeness en medias por factor
procesar_completeness <- function(completeness_df, metadata) {
  # pivotar a “largo”
  long <- completeness_df %>%
    pivot_longer(
      cols      = -Sample,
      names_to  = "Pathway",
      values_to = "Completeness"
    ) %>%
    left_join(metadata, by = "Sample")
  
  # medias por Sector
  sector_means <- long %>%
    group_by(Pathway, Sector) %>%
    summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Sector, values_from = mean_completeness)
  
  # medias por Season
  season_means <- long %>%
    group_by(Pathway, Season) %>%
    summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Season, values_from = mean_completeness)
  
  # medias por Habitat
  habitat_means <- long %>%
    group_by(Pathway, Habitat) %>%
    summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Habitat, values_from = mean_completeness)
  
  # unir todo
  full_means <- sector_means %>%
    left_join(season_means,  by = "Pathway") %>%
    left_join(habitat_means, by = "Pathway")
  
  return(full_means)
}

# 3) Calcular medias de nitrogen completeness
nitrogen_means <- procesar_completeness(completeness_df, metadata)

# 4) Función para graficar heatmaps de medias de completeness
crear_heatmaps <- function(means_df) {
  # long formato para cada factor
  sec <- means_df %>%
    pivot_longer(-Pathway, names_to="Sector",    values_to="Completeness") %>%
    filter(Sector %in% c("A_Inlet","B_Transition","C_Inner"))
  
  seas <- means_df %>%
    pivot_longer(-Pathway, names_to="Season",    values_to="Completeness") %>%
    filter(Season %in% c("Intense","Relax"))
  
  hab <- means_df %>%
    pivot_longer(-Pathway, names_to="Habitat",   values_to="Completeness") %>%
    filter(Habitat %in% c("Bare","Seagrass"))
  
  # helper
  mk <- function(df, xvar){
    ggplot(df, aes_string(x = xvar, y = "Pathway", fill = "Completeness")) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradientn(
        colors = c("#ffd166","#619b8a","#3e5c76","#6a4c93","#a4161a"),
        limits = c(0,100),
        na.value = "grey80"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = if(xvar=="Sector") element_text() else element_blank(),
        axis.ticks.y = if(xvar=="Sector") element_line() else element_blank()
      ) +
      labs(x = NULL, y = NULL) +
      coord_fixed(ratio = 0.2)
  }
  
  p1 <- mk(sec,  "Sector")
  p2 <- mk(seas,  "Season")
  p3 <- mk(hab,   "Habitat")
  
  (p1 | p2 | p3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# Generar el heatmap de nitrogen completeness
heatmap_nitrogen <- crear_heatmaps(nitrogen_means)
print(heatmap_nitrogen)

# Ajustar la altura dinámicamente
height_factor <- 0.5  # Ajusta este valor según el espaciado que necesites
num_pathways_nitrogen <- 7

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_NITROGEN_completeness.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_nitrogen
dev.off()



## Para SULFUR

# 0) Cargar librerías
library(dplyr)

# 1) Definir la lista de genes por pathway de azufre
genes_list_sulfur <- list(
  Assimilatory_sulfate_reduction = list(
    step1 = list(c("sat"),
                 c("cysN_cysC","cysN","cysD")),
    step2 = c("cysN_cysC","cysC"),
    step3 = list(c("cysH"), c("nrnA")),
    step4 = list(c("cysJ","cysI"), c("sir"))
  ),
  
  Dissimilatory_sulfur_reduction_and_oxidation = list(
    step1 = c("sat"),
    step2 = list(c("aprA","aprB"),
                 c("qmoA","qmoB","qmoC")),
    step3 = c("dsrA","dsrB","dsrL"),
    step4 = c("dsrC"),
    step5 = c("dsrJ","dsrK","dsrM","dsrN","dsrO","dsrP","dsrT"),
    step6 = c("dsrE","dsrF","dsrH"),
    step7 = c("dsrD","rdsr")
  ),
  
  Sulfur_reduction = list(
    step1 = list(c("asrA","asrB","asrC"),
                 c("fsr"),
                 c("mccA")),
    step2 = list(c("hydB","hydD","hydG"),
                 c("shyA","shyB","shyC","shyD")),
    step3 = c("sudA","sudB","rdlA"),
    step4 = list(c("otr"),
                 c("ttrA","ttrB","ttrC")),
    step5 = list(c("psrA","psrB","psrC"),
                 c("sreA","sreB","sreC"))
  ),
  
  SOX_systems = list(
    step1 = list(c("soxA","soxX"),
                 c("soxY","soxZ")),
    step2 = c("soxB"),
    step3 = c("soxC","soxD")
  ),
  
  Sulfur_oxidation = list(
    step1 = list(c("doxA","doxD"),
                 c("tsda","tsdB")),
    step2 = c("fccA","fccB","sqr"),
    step3 = c("glpE","sseA"),
    step4 = list(c("sorA","sorB"),
                 c("soeA","soeB","soeC"))
  ),
  
  Sulfur_disproportionation = list(
    step1 = c("phsA","phsB","phsC"),
    step2 = c("tetH"),
    step3 = c("sor")
  ),
  
  Organic_sulfur_transformation = list(
    route1 = c("dsyB","dddA","dddC","dddD"),
    route2 = c("dsyB","dddK","dddL","dddP","dddQ","dddT","dddW","dddY","prpE"),
    route3 = c("dsyB","dmdA","dmdB","dmdC","dmdD"),
    route4 = list(c("dsyB"),
                  c("dddD"),
                  c("dddK","dddL","dddP","dddQ","dddT","dddW","dddY"),
                  c("dmoA"),
                  c("mddA")),
    route5 = list(c("dsyB"),
                  c("dddD"),
                  c("dddK","dddL","dddP","dddQ","dddT","dddW","dddY"),
                  c("dmsA","dmsB","dmsC"),
                  c("ddhA","ddhB","ddhC")),
    route6 = c("gah","hpsN","hpsO","hpsP","iseJ","isfD","mdh","mtsA","mtsB",
               "pta","sfnG","slcC","slcD","sqdB","sqdD","sqdX","tauX","tauY",
               "tmm","toa","tpa","yihQ")
  ),
  
  Link_between_inorganic_and_organic_sulfur_transformation = c(
    "cuyA","cysE","cysK","cysM","cysO",
    "hdrA1","hdrA2","hdrB1","hdrB2","hdrC1","hdrC2","hdrD","hdrE",
    "mccB","metA","metB","metC","metX","metY","metZ",
    "msmA","msmB","mtoX","ssuD","ssuE","suyA","suyB","tauD",
    "tbuB","tbuC","tmoC","tmoF","touC","touF","xsc"
  )
)

# 2) Preparar presence_df para estos genes
desired_genes <- unique(unlist(genes_list_sulfur))
present_genes <- intersect(desired_genes, names(sulfur_genes_abund))
missing_genes <- setdiff(desired_genes, names(sulfur_genes_abund))

presence_sulfur <- combined_data %>%
  select(Sample) %>%
  bind_cols(
    combined_data %>%
      select(all_of(present_genes)) %>%
      mutate(across(everything(), ~ as.integer(. > 0)))
  )

for(g in missing_genes) {
  presence_sulfur[[g]] <- 0L
}

presence_sulfur <- presence_sulfur %>%
  select(Sample, all_of(desired_genes))

# 3) Función auxiliar para pasos OR
step_present <- function(df_row, groups) {
  any(sapply(groups, function(g) any(df_row[g] == 1))) * 1
}

# 4) Calcular completeness para cada pathway de azufre
completeness_sulfur <- presence_sulfur %>%
  rowwise() %>%
  mutate(
    Assimilatory_sulfate_reduction = {
      s1 <- step_present(cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step1)
      s2 <- all(c_across(genes_list_sulfur$Assimilatory_sulfate_reduction$step2) == 1) * 1
      s3 <- step_present(cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step3)
      s4 <- step_present(cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step4)
      (s1 + s2 + s3 + s4) / 4 * 100
    },
    Dissimilatory_sulfur_reduction_and_oxidation = {
      s1 <- any(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step1)) * 1
      s2 <- step_present(cur_data(), genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step2)
      s3 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step3) == 1) * 1
      s4 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step4) == 1) * 1
      s5 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step5) == 1) * 1
      s6 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step6) == 1) * 1
      s7 <- all(c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step7) == 1) * 1
      (s1 + s2 + s3 + s4 + s5 + s6 + s7) / 7 * 100
    },
    Sulfur_reduction = {
      sum(sapply(genes_list_sulfur$Sulfur_reduction, function(gr) step_present(cur_data(), gr))) / 5 * 100
    },
    SOX_systems = {
      sum(sapply(genes_list_sulfur$SOX_systems, function(gr) step_present(cur_data(), gr))) / 3 * 100
    },
    Sulfur_oxidation = {
      sum(sapply(genes_list_sulfur$Sulfur_oxidation, function(gr) step_present(cur_data(), gr))) / 4 * 100
    },
    Sulfur_disproportionation = {
      # Paso 1: phsA, phsB, phsC todos presentes
      s1 <- all(c_across(genes_list_sulfur$Sulfur_disproportionation$step1) == 1) * 1
      # Paso 2: tetH presente
      s2 <- c_across(genes_list_sulfur$Sulfur_disproportionation$step2)
      # Paso 3: sor presente
      s3 <- c_across(genes_list_sulfur$Sulfur_disproportionation$step3)
      # Completitud
      (s1 + s2 + s3) / 3 * 100
    },
    Organic_sulfur_transformation = {
      # Calcular cada ruta y tomar la máxima
      r1 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route1)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route1) * 100
      r2 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route2)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route2) * 100
      r3 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route3)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route3) * 100
      r4 <- step_present(cur_data(), genes_list_sulfur$Organic_sulfur_transformation$route4) / 5 * 100
      r5 <- step_present(cur_data(), genes_list_sulfur$Organic_sulfur_transformation$route5) / 5 * 100
      r6 <- sum(c_across(genes_list_sulfur$Organic_sulfur_transformation$route6)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route6) * 100
      max(r1,r2,r3,r4,r5,r6)
    },
    Link_between_inorganic_and_organic_sulfur_transformation = {
      sum(c_across(genes_list_sulfur$Link_between_inorganic_and_organic_sulfur_transformation)) /
        length(genes_list_sulfur$Link_between_inorganic_and_organic_sulfur_transformation) * 100
    }
  ) %>%
  ungroup()

# 5) Ver resultados
pathway_cols <- c("Assimilatory_sulfate_reduction", "Dissimilatory_sulfur_reduction_and_oxidation",
                   "Sulfur_reduction", "SOX_systems", "Sulfur_oxidation", "Sulfur_disproportionation",
                   "Organic_sulfur_transformation", "Organic_sulfur_transformation",
                   "Link_between_inorganic_and_organic_sulfur_transformation")
nonzero_pathways <- completeness_sulfur %>%
  select(all_of(pathway_cols)) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  unlist() %>%
  `>`(0) %>%                # vector lógico
  which() %>%               # índices de TRUE
  names()          

completeness_df <- completeness_sulfur %>%
  select(Sample, all_of(nonzero_pathways))


# (Opcional) Unir a combined_data:
combined_data <- combined_data %>% left_join(completeness_sulfur, by="Sample")



# 3) Calcular medias de sulfur completeness
sulfur_means <- procesar_completeness(completeness_df, metadata)

# 4) Función para graficar heatmaps de medias de completeness
crear_heatmaps <- function(means_df) {
  # long formato para cada factor
  sec <- means_df %>%
    pivot_longer(-Pathway, names_to="Sector",    values_to="Completeness") %>%
    filter(Sector %in% c("A_Inlet","B_Transition","C_Inner"))
  
  seas <- means_df %>%
    pivot_longer(-Pathway, names_to="Season",    values_to="Completeness") %>%
    filter(Season %in% c("Intense","Relax"))
  
  hab <- means_df %>%
    pivot_longer(-Pathway, names_to="Habitat",   values_to="Completeness") %>%
    filter(Habitat %in% c("Bare","Seagrass"))
  
  # helper
  mk <- function(df, xvar){
    ggplot(df, aes_string(x = xvar, y = "Pathway", fill = "Completeness")) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradientn(
        colors = c("#ffd166","#619b8a","#3e5c76","#6a4c93","#a4161a"),
        limits = c(0,100),
        na.value = "grey80"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = if(xvar=="Sector") element_text() else element_blank(),
        axis.ticks.y = if(xvar=="Sector") element_line() else element_blank()
      ) +
      labs(x = NULL, y = NULL) +
      coord_fixed(ratio = 0.2)
  }
  
  p1 <- mk(sec,  "Sector")
  p2 <- mk(seas,  "Season")
  p3 <- mk(hab,   "Habitat")
  
  (p1 | p2 | p3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# Generar el heatmap de nitrogen completeness
heatmap_sulfur <- crear_heatmaps(sulfur_means)
print(heatmap_sulfur)

# Ajustar la altura dinámicamente
height_factor <- 0.5  # Ajusta este valor según el espaciado que necesites
num_pathways_sulfur <- 7

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_SULFUR_completeness.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_sulfur
dev.off()

















































### With SUPERFOCUS analysis
# Transponer la tabla de abundancia
#abundance_table_t <- t(abundance_table_annotation)
#abundance_df <- as.data.frame(abundance_table_t)
#abundance_df$SampleID <- rownames(abundance_df)


### With DATABASES analysis
count_nitrogen <- t(nitrogen_pathways_abund)
count_metals <- t(metal_pathways_abund)
count_sulfur <- t(sulfur_pathways_abund)
count_drug <- t(drug_pathways_abund)
count_virulence <- t(virulence_pathways_abund)

count_nitrogen_df <- as.data.frame(count_nitrogen)
count_nitrogen_df$SampleID <- rownames(count_nitrogen_df)

count_metals_df <- as.data.frame(count_metals)
count_metals_df$SampleID <- rownames(count_metals_df)

count_sulfur_df <- as.data.frame(count_sulfur)
count_sulfur_df$SampleID <- rownames(count_sulfur_df)

count_drug_df <- as.data.frame(count_drug)
count_drug_df$SampleID <- rownames(count_drug_df)

count_virulence_df <- as.data.frame(count_virulence)
count_virulence_df$SampleID <- rownames(count_virulence_df)

# Unir las tablas de abundancia y metadatos
#combined_data <- merge(abundance_df, metadata, by.x = "SampleID", by.y = "sample")
combined_data <- merge(count_nitrogen_df, metadata, by.x = "SampleID", by.y = "Sample")
combined_data <- merge(count_metals_df, metadata, by.x = "SampleID", by.y = "Sample")
combined_data <- merge(count_sulfur_df, metadata, by.x = "SampleID", by.y = "Sample")
combined_data <- merge(count_drug_df, metadata, by.x = "SampleID", by.y = "Sample")
combined_data <- merge(count_virulence_df, metadata, by.x = "SampleID", by.y = "Sample")

# Crear una nueva columna combinada de Sector y Season
combined_data$Sector_Season <- paste(combined_data$Sector, combined_data$Season, sep = "_")

# Crear una lista de columnas de interés
#Nitrogen
columnas_de_interes <- c(2, 3, 4, 5, 6, 7, 8)
# Metals
columnas_de_interes <- c(2,3,4,8,9,10,11,12,13,14,15,17,18,19,20)
# Sulfur
columnas_de_interes <- c(2, 3, 4, 5, 6)
# Drug
columnas_de_interes <- c(2:28)
# Virulence
columnas_de_interes <- c(2:17)


# Crear las comparaciones adecuadas
a_my_comparisons <- list(
  c("A_Inlet_Relax", "A_Inlet_Intense"),
  c("B_Transition_Relax", "B_Transition_Intense"),
  c("C_Inner_Relax", "C_Inner_Intense"),
  c("A_Inlet_Intense","B_Transition_Intense"),
  c("A_Inlet_Intense","C_Inner_Intense"),
  c("B_Transition_Intense","C_Inner_Intense"),
  c("A_Inlet_Relax","B_Transition_Relax"),
  c("A_Inlet_Relax","C_Inner_Relax"),
  c("B_Transition_Relax","C_Inner_Relax")
)



# Bucle para generar y almacenar las figuras en variables a1, a2, a3, ...
for (i in 1:length(columnas_de_interes)) {
  columna_index <- columnas_de_interes[i]
  columna_de_interes <- colnames(combined_data)[columna_index]
  
  # Crear el boxplot con el t test usando la nueva columna combinada
  figura <- ggplot(combined_data, aes(x = Sector_Season, y = .data[[columna_de_interes]], fill = Season)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
    labs(x = "", y = "", fill = "Season") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("Relax" = "#000080ff", "Intense" = "#ff0000ff")) +
    #stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.format", 
    #                   position = position_dodge(0.75)) +
    ggtitle(paste(columna_de_interes)) #+
    #ylim(0, 7.5)  # Fijar los límites del eje y
  
  # Almacenar la figura en una variable dinámica (a1, a2, a3, ...)
  assign(paste0("a", i), figura)
  
  # Almacenar la figura en la lista
  #figuras[[paste0("a", i)]] <- figura
  
  # Guardar cada figura en un archivo
  #ggsave(paste0("figura_", i, ".png"), plot = figura, width = 10, height = 7)
}

svglite("annot_nitrogen_metabolism.svg", width=27, height=4)
((a1 | a2 | a3 | a4 | a5 | a6 | a7))
dev.off()

svglite("annot_metal_metabolism.svg", width=27, height=8)
((a1 | a2 | a3 | a4 | a5 | a6 | a7) / (a8 | a9 | a10| a11 | a12 | a13 | a14))
dev.off()

svglite("annot_sulfur_metabolism.svg", width=19.3, height=4)
((a1 | a2 | a3 | a4 | a5))
dev.off()

svglite("annot_drug_metabolism.svg", width=27, height=16)
((a1 | a2 | a3 | a4 | a5 | a6 | a7) / (a8 | a9 | a10| a11 | a12 | a13 | a14) / (a15 | a16 | a17| a18 | a19 | a20 | a21) / (a22 | a23 | a24 | a25 | a26 | a1 | a2))
dev.off()

svglite("annot_virulence_metabolism.svg", width=23, height=12)
((a1 | a2 | a3 | a4 | a5 | a6) / (a7 | a8 | a9 | a10 | a11 | a12) / (a13 | a14 | a15 | a16 | a1 | a12))
dev.off()

summary(metal_pathways_abund)


# Función para calcular promedio, mínimo y máximo
calcular_estadisticas <- function(df) {
  df_matrix <- as.matrix(df)  # Convertir a matriz
  promedio <- mean(df_matrix, na.rm = TRUE)
  minimo <- min(df_matrix, na.rm = TRUE)
  maximo <- max(df_matrix, na.rm = TRUE)
  
  return(c(Promedio = promedio, Minimo = minimo, Maximo = maximo))
}

# Aplicar la función a cada DataFrame
estadisticas_nitrogen <- calcular_estadisticas(count_nitrogen)
estadisticas_metals <- calcular_estadisticas(count_metals)
estadisticas_sulfur <- calcular_estadisticas(count_sulfur)
estadisticas_drug <- calcular_estadisticas(count_drug)
estadisticas_virulence <- calcular_estadisticas(count_virulence)

# Mostrar resultados
estadisticas_nitrogen
estadisticas_metals
estadisticas_sulfur
estadisticas_drug
estadisticas_virulence











### With other DATABASES






#####

# Seleccionar las columnas de interés
cols_interes <- c("A_Inlet", "B_Transition", "C_Inner", "Relax", "Intense", "Seagrass", "Bare")

# Calcular el porcentaje promedio de cada columna
colMeans(nitrogen_means[cols_interes])
colMeans(sulfur_means[cols_interes])




#####
# DIFERENCIAS SIGNIFICATIVAS PORCENTAJES DE GENES

load("genes_abund_metal.RData")
load("genes_abund_drug.RData")
load("genes_abund_nitrogen.RData")
load("genes_abund_sulfur.RData")
load("genes_abund_virulence.RData")

metadata <- read.csv("metadata_quintin.csv")

library(Maaslin2)
metadata.maaslin <- metadata
rownames(metadata.maaslin) <- metadata.maaslin$Sample
metadata.maaslin$Sample <- NULL
fit_dataMaaslin2fit_data <- Maaslin2(
  input_data = t(sulfur_genes_abund),
  input_metadata = metadata.maaslin,
  output = "maaslind2_sulfur_abund",
  fixed_effects = c("Sector"),
  reference = "Sector,A_Inlet"
)


# Read abundance tables
genes_abund_df <- nitrogen_genes_abund
genes_abund_df <- sulfur_genes_abund
genes_abund_df <- t(metal_genes_abund)
genes_abund_df <- t(drug_genes_abund)

# Transform abundance tables
genes_abund_df <- as.data.frame(genes_abund_df)
genes_abund_df$Sample <- rownames(genes_abund_df)

# Unir las tablas de abundancia y metadatos
combined_data_abund <- merge(genes_abund_df, metadata, by.x = "Sample", by.y = "Sample")

# Crear una nueva columna combinada de Sector y Season
combined_data_abund$Sector_Season <- paste(combined_data_abund$Sector, combined_data_abund$Season, sep = "_")
combined_data_abund$Sector_Season <- factor(combined_data_abund$Sector_Season, 
                                      levels = c("A_Inlet_Relax", "A_Inlet_Intense",
                                                 "B_Transition_Relax", "B_Transition_Intense",
                                                 "C_Inner_Relax", "C_Inner_Intense"))

# Cargar paquetes necesarios
library(FSA)      # Para Dunn test
library(dplyr)    # Para manipulación de datos
library(tidyr)    # Para transformar datos
library(rstatix)  # Para múltiples comparaciones ajustadas

# Seleccionar solo las columnas numéricas (excluyendo "Sample", "Sector", "Season", "Habitat")
genes_cols <- combined_data_abund %>%
  select(-Sample, -Sector, -Season, -Habitat, -ID_1, -ID_2)

# Función para aplicar Dunn test (Sector: 3 niveles)
dunn_results <- lapply(names(genes_cols), function(gene) {
  test <- dunnTest(combined_data_abund[[gene]] ~ combined_data_abund$Sector, method = "bh")
  result <- data.frame(Gene = gene, test$res)
  return(result)
}) %>%
  bind_rows() %>%
  filter(P.adj <= 0.05)  # Filtrar p-values ajustados < 0.05
  #filter(P.adj > 0.05 & P.adj <= 0.1)  # Filtrar p-values <= 0.05

# Función para aplicar Wilcoxon test (Season y Habitat: 2 niveles)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(genes_cols), function(gene) {
    test <- wilcox.test(combined_data_abund[[gene]] ~ combined_data_abund[[variable]])
    data.frame(Gene = gene, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  filter(p.value <= 0.05)  # Filtrar p-values <= 0.05
  #filter(p.value > 0.05 & p.value <= 0.1)  # Filtrar p-values <= 0.05

# Mostrar los resultados significativos
list(Dunn_Test_Significativo = dunn_results, Wilcoxon_Significativo = wilcox_results)



####

# Función para convertir matriz en formato largo, unir con metadata, y calcular medias
procesar_genes <- function(genes_abund, metadata) {
  # Convertir la matriz en data.frame
  genes_abund_df <- as.data.frame(genes_abund)
  
  # Convertir al formato largo
  genes_long <- genes_abund_df %>%
    rownames_to_column(var = "Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance")
  
  # Unir con metadata
  genes_metadata <- genes_long %>%
    left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calcular medias por Sector
  sector_means <- genes_metadata %>%
    dplyr::group_by(Gene, Sector) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Sector, values_from = mean_abundance)
  
  # Calcular medias por Season
  season_means <- genes_metadata %>%
    dplyr::group_by(Gene, Season) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Season, values_from = mean_abundance)
  
  # Calcular medias por Habitat
  habitat_means <- genes_metadata %>%
    dplyr::group_by(Gene, Habitat) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Habitat, values_from = mean_abundance)
  
  # Unir todas las medias en un único dataframe
  genes_means <- genes_abund_df %>%
    rownames_to_column(var = "Gene") %>%
    left_join(sector_means, by = "Gene") %>%
    left_join(season_means, by = "Gene") %>%
    left_join(habitat_means, by = "Gene")
  
  return(genes_means)
}

# Aplicar la función a cada una de las matrices
metal_means <- procesar_genes(metal_genes_abund, metadata)
drug_means <- procesar_genes(drug_genes_abund, metadata)
nitrogen_means <- procesar_genes(t(nitrogen_genes_abund), metadata)
sulfur_means <- procesar_genes(sulfur_genes_abund, metadata)
virulence_means <- procesar_genes(virulence_genes_abund, metadata)

# Calculo los valores máximos
max_value_drug <- max(drug_means[ , (ncol(drug_means) - 6):ncol(drug_means)], na.rm = TRUE)
max_value_metal <- max(metal_means[ , (ncol(metal_means) - 6):ncol(metal_means)], na.rm = TRUE)
max_value_nitrogen <- max(nitrogen_means[ , (ncol(nitrogen_means) - 6):ncol(nitrogen_means)], na.rm = TRUE)
max_value_sulfur <- max(sulfur_means[ , (ncol(sulfur_means) - 6):ncol(sulfur_means)], na.rm = TRUE)

# Función para generar los heatmaps de un archivo *_means
crear_heatmaps <- function(means_data, max_val) {
  
  # Convertir a formato largo para Sector, Season y Habitat
  sector_data <- means_data %>%
    dplyr::select(Gene, starts_with("A_"), starts_with("B_"), starts_with("C_In")) %>%
    pivot_longer(-Gene, names_to = "Sector", values_to = "Abundance")
  
  season_data <- means_data %>%
    select(Gene, starts_with("Intense"), starts_with("Relax")) %>%
    pivot_longer(-Gene, names_to = "Season", values_to = "Abundance")
  
  habitat_data <- means_data %>%
    select(Gene, starts_with("Bare"), starts_with("Seagrass")) %>%
    pivot_longer(-Gene, names_to = "Habitat", values_to = "Abundance")
  
  # Crear heatmap para Sector
  heatmap_sector <- ggplot(sector_data, aes(x = Sector, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Season
  heatmap_season <- ggplot(season_data, aes(x = Season, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los genes en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Habitat
  heatmap_habitat <- ggplot(habitat_data, aes(x = Habitat, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"), 
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los genes en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Unir los tres heatmaps en un solo gráfico utilizando patchwork
  library(patchwork)
  
  # Combinar los tres heatmaps en una figura
  combined_plot <- (heatmap_sector | heatmap_season | heatmap_habitat) + 
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

# Generar heatmaps para cada archivo *_means
heatmap_drug <- crear_heatmaps(drug_means, max_value_drug)
heatmap_nitrogen <- crear_heatmaps(nitrogen_means, max_value_nitrogen)
heatmap_sulfur <- crear_heatmaps(sulfur_means, max_value_sulfur)

# Mostrar el plot de cada uno
print(heatmap_drug)
print(heatmap_nitrogen)
print(heatmap_sulfur)

# Definir el número de genes en cada dataset
num_genes_drug <- nrow(drug_means)
num_genes_nitrogen <- nrow(nitrogen_means)
num_genes_sulfur <- nrow(sulfur_means)

# Ajustar la altura dinámicamente
height_factor <- 0.3  # Ajusta este valor según el espaciado que necesites

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_CARD_genes.svg", width = 8, height = max(10, num_genes_drug * height_factor))
heatmap_drug
dev.off()

svglite("heatmap_NITROGEN_genes.svg", width = 8, height = max(8, num_genes_nitrogen * height_factor))
heatmap_nitrogen
dev.off()

svglite("heatmap_SULFUR_genes.svg", width = 8, height = max(8, num_genes_sulfur * height_factor))
heatmap_sulfur
dev.off()


# UpSet plot count tables


load("upset_METAL_genes.RData")
load("upset_DRUG_genes.RData")
load("upset_NITROGENL_genes.RData")
load("upset_SULFUR_genes.RData")
load("upset_VFDB_genes.RData")

# Read metadata
metadata <- read.csv("metadata_quintin.csv")
metadata$Sector <- sub("^[A-Z]_+", "", metadata$Sector)

# Lista de data.frames Upset ya cargados
upsets <- list(
  METAL     = upset_METAL_genes,
  DRUG      = upset_DRUG_genes,
  NITROGEN  = upset_NITROGEN_genes,
  SULFUR    = upset_SULFUR_genes,
  VIRULENCE = upset_VFDB_genes
)

# 2. Función que resume presencia/ausencia por nivel de un factor
summarise_pa <- function(df, meta, factor_name) {
  # Asegura mismo orden de columnas que en meta$Sample
  df <- df[, meta$Sample, drop = FALSE]
  # Niveles del factor
  levels <- unique(meta[[factor_name]])
  # Para cada nivel, presencia = 1 si al menos un sample tiene >0
  pa <- sapply(levels, function(lvl) {
    samp  <- meta$Sample[meta[[factor_name]] == lvl]
    subdf <- df[, samp, drop = FALSE]
    as.integer(rowSums(subdf) > 0)
  })
  # Ajusta filas y columnas
  pa <- t(pa)
  rownames(pa) <- levels
  return(pa)
}

# 3. Genera las tablas de PA para cada grupo y cada factor
pa_tables <- lapply(upsets, function(df) {
  list(
    Sector  = summarise_pa(df, metadata, "Sector"),
    Season  = summarise_pa(df, metadata, "Season"),
    Habitat = summarise_pa(df, metadata, "Habitat")
  )
})

# Ejemplo de acceso:
# pa_tables$METAL$Sector    → matriz PA genes × sectores
# pa_tables$DRUG$Habitat    → matriz PA genes × hábitats
pa_tables$NITROGEN




# Build the presence/absence matrix for METAL
df_pa <- as.data.frame(t(pa_tables$METAL$Sector))
df_pa <- as.data.frame(t(pa_tables$METAL$Season))
df_pa <- as.data.frame(t(pa_tables$METAL$Habitat))

# Plot
svglite("Upset_METAL_Sector.svg", width = 4, height = 4)
#svglite("Upset_METAL_Season.svg", width = 4, height = 4)
#svglite("Upset_METAL_Habitat.svg", width = 4, height = 4)
upset(
  df_pa,
  sets               = c("Inlet","Transition","Inner"),
  sets.bar.color     = c("#56b4e9","#c8ab37","#808080"),
  #sets               = c("Relax","Intense"),
  #sets.bar.color     = c("#ff0000","#0000ff"),
  #sets               = c("Bare","Seagrass"),
  #sets.bar.color     = c("#008000","#ff6600"),
  order.by           = "freq",
  empty.intersections= "on",
  mainbar.y.label    = "Number of shared and\nunique genes\nMetal-associated",
  sets.x.label       = "# genes",
  point.size         = 3,
  line.size          = 0.5,
  query.legend       = "top",
  text.scale         = c(2,2,1.5,1.5,2,1.5)
)
dev.off()


# Build the presence/absence matrix for DRUG
df_pa <- as.data.frame(t(pa_tables$DRUG$Sector))
df_pa <- as.data.frame(t(pa_tables$DRUG$Season))
df_pa <- as.data.frame(t(pa_tables$DRUG$Habitat))

# Plot
svglite("Upset_DRUG_Sector.svg", width = 4, height = 4)
#svglite("Upset_DRUG_Season.svg", width = 4, height = 4)
#svglite("Upset_DRUG_Habitat.svg", width = 4, height = 4)
upset(
  df_pa,
  sets               = c("Inlet","Transition","Inner"),
  sets.bar.color     = c("#56b4e9","#c8ab37","#808080"),
  #sets               = c("Relax","Intense"),
  #sets.bar.color     = c("#ff0000","#0000ff"),
  #sets               = c("Bare","Seagrass"),
  #sets.bar.color     = c("#008000","#ff6600"),
  order.by           = "freq",
  empty.intersections= "on",
  mainbar.y.label    = "Number of shared and\nunique genes\nAntibiotic resistance",
  sets.x.label       = "# genes",
  point.size         = 3,
  line.size          = 0.5,
  query.legend       = "top",
  text.scale         = c(2,2,1.5,1.5,2,1.5)
)
dev.off()


# Build the presence/absence matrix for NITROGEN
df_pa <- as.data.frame(t(pa_tables$NITROGEN$Sector))
df_pa <- as.data.frame(t(pa_tables$NITROGEN$Season))
df_pa <- as.data.frame(t(pa_tables$NITROGEN$Habitat))

# Plot
#svglite("Upset_NITROGEN_Sector.svg", width = 4, height = 4)
#svglite("Upset_NITROGEN_Season.svg", width = 4, height = 4)
svglite("Upset_NITROGEN_Habitat.svg", width = 4, height = 4)
upset(
  df_pa,
  #sets               = c("Inlet","Transition","Inner"),
  #sets.bar.color     = c("#56b4e9","#c8ab37","#808080"),
  #sets               = c("Relax","Intense"),
  #sets.bar.color     = c("#ff0000","#0000ff"),
  sets               = c("Bare","Seagrass"),
  sets.bar.color     = c("#008000","#ff6600"),
  order.by           = "freq",
  empty.intersections= "on",
  mainbar.y.label    = "Number of shared and\nunique genes\nNitrogen metabolism",
  sets.x.label       = "# genes",
  point.size         = 3,
  line.size          = 0.5,
  query.legend       = "top",
  text.scale         = c(2,2,1.5,1.5,2,1.5)
)
dev.off()



# Build the presence/absence matrix for SULFUR
df_pa <- as.data.frame(t(pa_tables$SULFUR$Sector))
df_pa <- as.data.frame(t(pa_tables$SULFUR$Season))
df_pa <- as.data.frame(t(pa_tables$SULFUR$Habitat))

# Plot
svglite("Upset_SULFUR_Sector.svg", width = 4, height = 4)
#svglite("Upset_SULFUR_Season.svg", width = 4, height = 4)
#svglite("Upset_SULFUR_Habitat.svg", width = 4, height = 4)
upset(
  df_pa,
  sets               = c("Inlet","Transition","Inner"),
  sets.bar.color     = c("#56b4e9","#c8ab37","#808080"),
  #sets               = c("Relax","Intense"),
  #sets.bar.color     = c("#ff0000","#0000ff"),
  #sets               = c("Bare","Seagrass"),
  #sets.bar.color     = c("#008000","#ff6600"),
  order.by           = "freq",
  empty.intersections= "on",
  mainbar.y.label    = "Number of shared and\nunique genes\nSulfur metabolism",
  sets.x.label       = "# genes",
  point.size         = 3,
  line.size          = 0.5,
  query.legend       = "top",
  text.scale         = c(2,2,1.5,1.5,2,1.5)
)
dev.off()



# Build the presence/absence matrix for VIRULENCE
df_pa <- as.data.frame(t(pa_tables$VIRULENCE$Sector))
df_pa <- as.data.frame(t(pa_tables$VIRULENCE$Season))
df_pa <- as.data.frame(t(pa_tables$VIRULENCE$Habitat))

# Plot
#svglite("Upset_VIRULENCE_Sector.svg", width = 4, height = 4)
#svglite("Upset_VIRULENCE_Season.svg", width = 4, height = 4)
svglite("Upset_VIRULENCE_Habitat.svg", width = 4, height = 4)
upset(
  df_pa,
  #sets               = c("Inlet","Transition","Inner"),
  #sets.bar.color     = c("#56b4e9","#c8ab37","#808080"),
  #sets               = c("Relax","Intense"),
  #sets.bar.color     = c("#ff0000","#0000ff"),
  sets               = c("Bare","Seagrass"),
  sets.bar.color     = c("#008000","#ff6600"),
  order.by           = "freq",
  empty.intersections= "on",
  mainbar.y.label    = "Number of shared and\nunique genes\nVirulence factors",
  sets.x.label       = "# genes",
  point.size         = 3,
  line.size          = 0.5,
  query.legend       = "top",
  text.scale         = c(2,2,1.5,1.5,2,1.5)
)
dev.off()

















#####
# DIFERENCIAS SIGNIFICATIVAS PORCENTAJES DE GENES

# Read abundance tables
genes_abund_df <- t(nitrogen_pathways_abund)
genes_abund_df <- t(sulfur_pathways_abund)
genes_abund_df <- (metal_pathways_abund)
genes_abund_df <- t(drug_pathways_abund)
genes_abund_df <- t(virulence_pathways_abund)



# Transform abundance tables
genes_abund_df <- as.data.frame(genes_abund_df)
genes_abund_df$Sample <- rownames(genes_abund_df)

# Unir las tablas de abundancia y metadatos
combined_data_abund <- merge(genes_abund_df, metadata, by.x = "Sample", by.y = "Sample")

# Crear una nueva columna combinada de Sector y Season
combined_data_abund$Sector_Season <- paste(combined_data_abund$Sector, combined_data_abund$Season, sep = "_")
combined_data_abund$Sector_Season <- factor(combined_data_abund$Sector_Season, 
                                            levels = c("A_Inlet_Relax", "A_Inlet_Intense",
                                                       "B_Transition_Relax", "B_Transition_Intense",
                                                       "C_Inner_Relax", "C_Inner_Intense"))

# Cargar paquetes necesarios
library(FSA)      # Para Dunn test
library(dplyr)    # Para manipulación de datos
library(tidyr)    # Para transformar datos
library(rstatix)  # Para múltiples comparaciones ajustadas

# Seleccionar solo las columnas numéricas (excluyendo "Sample", "Sector", "Season", "Habitat")
genes_cols <- combined_data_abund %>%
  select(-Sample, -Sector, -Season, -Habitat, -Sector_Season, -ID_1, -ID_2)

# Función para aplicar Dunn test (Sector: 3 niveles)
dunn_results <- lapply(names(genes_cols), function(gene) {
  test <- dunnTest(combined_data_abund[[gene]] ~ combined_data_abund$Sector, method = "bh")
  result <- data.frame(Gene = gene, test$res)
  return(result)
}) %>%
  bind_rows() %>%
  filter(P.adj <= 0.05)  # Filtrar p-values ajustados < 0.05
  #filter(P.adj > 0.05 & P.adj <= 0.1)  # Filtrar p-values <= 0.05
#Pathway        Comparison        Z    P.unadj      P.adj
# Assimilatory_nitrate_reduction A_Inlet - C_Inner 2.647568 0.00810731 0.02432193

# Función para aplicar Wilcoxon test (Season y Habitat: 2 niveles)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(genes_cols), function(gene) {
    test <- wilcox.test(combined_data_abund[[gene]] ~ combined_data_abund[[variable]])
    data.frame(Gene = gene, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  filter(p.value <= 0.05)  # Filtrar p-values <= 0.05
  #filter(p.value > 0.05 & p.value <= 0.1)  # Filtrar p-values <= 0.05

# Mostrar los resultados significativos
list(Dunn_Test_Significativo = dunn_results, Wilcoxon_Significativo = wilcox_results)



####



# Función para convertir matriz en formato largo, unir con metadata, y calcular medias
procesar_pathways <- function(pathways_abund, metadata) {
  # Convertir la matriz en data.frame
  pathways_abund_df <- as.data.frame(pathways_abund)
  
  # Convertir al formato largo
  pathways_long <- pathways_abund_df %>%
    rownames_to_column(var = "Pathways") %>%
    pivot_longer(-Pathways, names_to = "Sample", values_to = "Abundance")
  
  # Unir con metadata
  pathways_metadata <- pathways_long %>%
    left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calcular medias por Sector
  sector_means <- pathways_metadata %>%
    group_by(Pathways, Sector) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Sector, values_from = mean_abundance)
  
  # Calcular medias por Season
  season_means <- pathways_metadata %>%
    group_by(Pathways, Season) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Season, values_from = mean_abundance)
  
  # Calcular medias por Habitat
  habitat_means <- pathways_metadata %>%
    group_by(Pathways, Habitat) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Habitat, values_from = mean_abundance)
  
  # Unir todas las medias en un único dataframe
  pathways_means <- pathways_abund_df %>%
    rownames_to_column(var = "Pathways") %>%
    left_join(sector_means, by = "Pathways") %>%
    left_join(season_means, by = "Pathways") %>%
    left_join(habitat_means, by = "Pathways")
  
  return(pathways_means)
}

# Aplicar la función a cada una de las matrices
metal_means <- procesar_pathways(metal_pathways_abund, metadata)
drug_means <- procesar_pathways(drug_pathways_abund, metadata)
nitrogen_means <- procesar_pathways(nitrogen_pathways_abund, metadata)
sulfur_means <- procesar_pathways(sulfur_pathways_abund, metadata)
virulence_means <- procesar_pathways(virulence_pathways_abund, metadata)

# Calculo los valores máximos
max_value_metal <- max(metal_means[ , (ncol(metal_means) - 6):ncol(metal_means)], na.rm = TRUE)
max_value_drug <- max(drug_means[ , (ncol(drug_means) - 6):ncol(drug_means)], na.rm = TRUE)
max_value_nitrogen <- max(nitrogen_means[ , (ncol(nitrogen_means) - 6):ncol(nitrogen_means)], na.rm = TRUE)
max_value_sulfur <- max(sulfur_means[ , (ncol(sulfur_means) - 6):ncol(sulfur_means)], na.rm = TRUE)
max_value_virulence <- max(virulence_means[ , (ncol(virulence_means) - 6):ncol(virulence_means)], na.rm = TRUE)

# Función para Pathways crear los heatmaps de un archivo *_means
crear_heatmaps <- function(means_data, max_val) {
  
  # Convertir a formato largo para Sector, Season y Habitat
  sector_data <- means_data %>%
    select(Pathways, starts_with("A_"), starts_with("B_"), starts_with("C_In")) %>%
    pivot_longer(-Pathways, names_to = "Sector", values_to = "Abundance")
  
  season_data <- means_data %>%
    select(Pathways, starts_with("Intense"), starts_with("Relax")) %>%
    pivot_longer(-Pathways, names_to = "Season", values_to = "Abundance")
  
  habitat_data <- means_data %>%
    select(Pathways, starts_with("Bare"), starts_with("Seagrass")) %>%
    pivot_longer(-Pathways, names_to = "Habitat", values_to = "Abundance")
  
  # Crear heatmap para Sector
  heatmap_sector <- ggplot(sector_data, aes(x = Sector, y = Pathways, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Season
  heatmap_season <- ggplot(season_data, aes(x = Season, y = Pathways, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los pathways en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Habitat
  heatmap_habitat <- ggplot(habitat_data, aes(x = Habitat, y = Pathways, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los pathways en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Unir los tres heatmaps en un solo gráfico utilizando patchwork
  library(patchwork)
  
  # Combinar los tres heatmaps en una figura
  combined_plot <- (heatmap_sector | heatmap_season | heatmap_habitat) + 
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

# Generar heatmaps para cada archivo *_means
heatmap_drug <- crear_heatmaps(drug_means, max_value_drug)
heatmap_nitrogen <- crear_heatmaps(nitrogen_means, max_value_nitrogen)
heatmap_sulfur <- crear_heatmaps(sulfur_means, max_value_sulfur)
heatmap_metal <- crear_heatmaps(metal_means, max_value_metal)
heatmap_virulence <- crear_heatmaps(virulence_means, max_value_virulence)

# Mostrar el plot de cada uno
print(heatmap_drug)
print(heatmap_nitrogen)
print(heatmap_sulfur)
print(heatmap_metal)
print(heatmap_virulence)

# Definir el número de genes en cada dataset
num_pathways_drug <- nrow(drug_means)
num_pathways_nitrogen <- nrow(nitrogen_means)
num_pathways_sulfur <- nrow(sulfur_means)
num_pathways_metal <- nrow(metal_means)
num_pathways_virulence <- nrow(virulence_means)

# Ajustar la altura dinámicamente
height_factor <- 0.3  # Ajusta este valor según el espaciado que necesites

# Guardar las figuras con dimensiones ajustadas
svglite("heatmap_CARD_pathway2.svg", width = 8, height = max(10, num_pathways_drug * height_factor))
heatmap_drug
dev.off()

svglite("heatmap_NITROGEN_pathway2.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_nitrogen
dev.off()

svglite("heatmap_SULFUR_pathway2.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_sulfur
dev.off()

svglite("heatmap_METAL_pathway2.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_metal
dev.off()

svglite("heatmap_VFDB_pathway2.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_virulence
dev.off()



# UpSet plot count tables

genus_comparison <- read.csv(paste0(base_path, "upset_BACMET_Metal_pathways.csv"), sep = ",", row.names = 1)
svglite("Upset_BACMET_pathways_sector.svg", width = 4, height = 4)
upset(
  genus_comparison, sets = c("Inner", "Transition", "Inlet"),
  #genus_comparison, sets = c("Bare", "Seagrass"),
  #genus_comparison, sets = c("Relaxing", "Upwelling"),
  sets.bar.color="#56B4E9", order.by="freq",empty.intersections = "on",
  mainbar.y.label = "Number of shared and\nnon-shared Metal\nResistance Types (MRTs) ", sets.x.label = "Number of MRTs",
  point.size = 3, line.size = 0.5, query.legend = "top",
  text.scale = c(1.5,1,1,1,1,1))
dev.off()

genus_comparison <- read.csv(paste0(base_path, "upset_CARD_class.csv"), sep = ",", row.names = 1)
svglite("Upset_CARD_pathways_sector.svg", width = 4, height = 4)
upset(
  genus_comparison, sets = c("Inner", "Transition", "Inlet"),
  #genus_comparison, sets = c("Bare", "Seagrass"),
  #genus_comparison, sets = c("Relaxing", "Upwelling"),
  sets.bar.color="#56B4E9", order.by="freq",empty.intersections = "on",
  mainbar.y.label = "Number of shared and\nnon-shared\nAntibiotic class (AC)", sets.x.label = "Number of ACs",
  point.size = 3, line.size = 0.5, query.legend = "top",
  text.scale = c(1.5,1,1,1,1,1))
dev.off()

genus_comparison <- read.csv(paste0(base_path, "upset_NITROGEN_pathways.csv"), sep = ",", row.names = 1)
svglite("Upset_NITROGEN_pathways_season.svg", width = 4, height = 4)
upset(
  #genus_comparison, sets = c("Inner", "Transition", "Inlet"),
  #genus_comparison, sets = c("Bare", "Seagrass"),
  genus_comparison, sets = c("Relaxing", "Upwelling"),
  sets.bar.color="#56B4E9", order.by="freq",empty.intersections = "on",
  mainbar.y.label = "Number of shared and\nnon-shared Nitrogen\nmetabolism pathways", sets.x.label = "Number of pathways",
  point.size = 3, line.size = 0.5, query.legend = "top",
  text.scale = c(1.5,1,1,1,1,1))
dev.off()

genus_comparison <- read.csv(paste0(base_path, "upset_SULFUR_pathways.csv"), sep = ",", row.names = 1)
svglite("Upset_SULFUR_pathways_sector.svg", width = 4, height = 4)
upset(
  genus_comparison, sets = c("Inner", "Transition", "Inlet"),
  #genus_comparison, sets = c("Bare", "Seagrass"),
  #genus_comparison, sets = c("Relaxing", "Upwelling"),
  sets.bar.color="#56B4E9", order.by="freq",empty.intersections = "on",
  mainbar.y.label = "Number of shared and\nnon-shared Sulfur\nmetabolism pathways", sets.x.label = "Number of pathways",
  point.size = 3, line.size = 0.5, query.legend = "top",
  text.scale = c(1.5,1,1,1,1,1))
dev.off()

genus_comparison <- read.csv(paste0(base_path, "upset_VFDB_type.csv"), sep = ",", row.names = 1)
svglite("Upset_VFDB_types_season.svg", width = 4, height = 4)
upset(
  #genus_comparison, sets = c("Inner", "Transition", "Inlet"),
  #genus_comparison, sets = c("Bare", "Seagrass"),
  genus_comparison, sets = c("Relaxing", "Upwelling"),
  sets.bar.color="#56B4E9", order.by="freq",empty.intersections = "on",
  mainbar.y.label = "Number of shared and\nnon-shared Virulence\nfactor types (VFTs)", sets.x.label = "Number of VFTs",
  point.size = 3, line.size = 0.5, query.legend = "top",
  text.scale = c(1.5,1,1,1,1,1))
dev.off()



## Using pathway completeness

nitrogen_completeness <- read.csv(paste0(base_path, "upset_NITROGEN_completeness.csv"), sep = ",", row.names = 1)
sulfur_completeness <- read.csv(paste0(base_path, "upset_SULFUR_completeness.csv"), sep = ",", row.names = 1)

# Función para convertir matriz en formato largo, unir con metadata, y calcular medias
procesar_pathways <- function(pathways_completeness, metadata) {
  # Convertir la matriz en data.frame
  pathways_completeness_df <- as.data.frame(pathways_completeness)
  
  # Convertir al formato largo
  pathways_long <- pathways_completeness_df %>%
    rownames_to_column(var = "Pathways") %>%
    pivot_longer(-Pathways, names_to = "Sample", values_to = "Completeness")
  
  # Unir con metadata
  pathways_metadata <- pathways_long %>%
    left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calcular medias por Sector
  sector_means <- pathways_metadata %>%
    group_by(Pathways, Sector) %>%
    summarize(mean_completeness = mean(Completeness, na.rm = TRUE)) %>%
    pivot_wider(names_from = Sector, values_from = mean_completeness)
  
  # Calcular medias por Season
  season_means <- pathways_metadata %>%
    group_by(Pathways, Season) %>%
    summarize(mean_completeness = mean(Completeness, na.rm = TRUE)) %>%
    pivot_wider(names_from = Season, values_from = mean_completeness)
  
  # Calcular medias por Habitat
  habitat_means <- pathways_metadata %>%
    group_by(Pathways, Habitat) %>%
    summarize(mean_completeness = mean(Completeness, na.rm = TRUE)) %>%
    pivot_wider(names_from = Habitat, values_from = mean_completeness)
  
  # Unir todas las medias en un único dataframe
  pathways_means <- pathways_completeness_df %>%
    rownames_to_column(var = "Pathways") %>%
    left_join(sector_means, by = "Pathways") %>%
    left_join(season_means, by = "Pathways") %>%
    left_join(habitat_means, by = "Pathways")
  
  return(pathways_means)
}

# Aplicar la función a cada una de las matrices
nitrogen_means <- procesar_pathways(nitrogen_completeness, metadata)
sulfur_means <- procesar_pathways(sulfur_completeness, metadata)

# Calculo los valores máximos
max_value_nitrogen <- max(nitrogen_means[ , (ncol(nitrogen_means) - 6):ncol(nitrogen_means)], na.rm = TRUE)
max_value_sulfur <- max(sulfur_means[ , (ncol(sulfur_means) - 6):ncol(sulfur_means)], na.rm = TRUE)

# Función para Pathways crear los heatmaps de un archivo *_means
crear_heatmaps <- function(means_data) {
  
  # Convertir a formato largo para Sector, Season y Habitat
  sector_data <- means_data %>%
    dplyr::select(Pathways, A_Inlet, B_Transition, C_Inner) %>%
    pivot_longer(-Pathways, names_to = "Sector", values_to = "Completeness")
  
  season_data <- means_data %>%
    dplyr::select(Pathways, starts_with("Intense"), starts_with("Relax")) %>%
    pivot_longer(-Pathways, names_to = "Season", values_to = "Completeness")
  
  habitat_data <- means_data %>%
    dplyr::select(Pathways, starts_with("Bare"), starts_with("Seagrass")) %>%
    pivot_longer(-Pathways, names_to = "Habitat", values_to = "Completeness")
  
  # Crear heatmap para Sector
  heatmap_sector <- ggplot(sector_data, aes(x = Sector, y = Pathways, fill = Completeness)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                         na.value = "grey50", limit = c(0, 100),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Season
  heatmap_season <- ggplot(season_data, aes(x = Season, y = Pathways, fill = Completeness)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),   
                         na.value = "grey50", limit = c(0, 100),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los pathways en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Crear heatmap para Habitat
  heatmap_habitat <- ggplot(habitat_data, aes(x = Habitat, y = Pathways, fill = Completeness)) +
    geom_tile(color = "white", linewidth = 0.5) +
    #scale_fill_gradient2(low = "white", high = "darkred", 
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),   
                         na.value = "grey50", limit = c(0, 100),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Ocultar nombres de los pathways en el eje y
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Unir los tres heatmaps en un solo gráfico utilizando patchwork
  library(patchwork)
  
  # Combinar los tres heatmaps en una figura
  combined_plot <- (heatmap_sector | heatmap_season | heatmap_habitat) + 
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

# Generar heatmaps para cada archivo *_means
heatmap_nitrogen <- crear_heatmaps(nitrogen_means)
heatmap_sulfur <- crear_heatmaps(sulfur_means)

# Mostrar el plot de cada uno
print(heatmap_nitrogen)
print(heatmap_sulfur)


# Definir el número de genes en cada dataset

num_pathways_nitrogen <- nrow(nitrogen_means)
num_pathways_sulfur <- nrow(sulfur_means)

# Ajustar la altura dinámicamente

height_factor <- 0.3  # Ajusta este valor según el espaciado que necesites

# Guardar las figuras con dimensiones ajustadas

svglite("heatmap_NITROGEN_completeness2.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_nitrogen
dev.off()

svglite("heatmap_SULFUR_completeness2.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_sulfur
dev.off()












##################################
#### PHYSICOCHEMICAL FEATURES ####
##################################

physchem <- read.csv("physchem_v1.csv")
physchem <- column_to_rownames(physchem, var = "UID")
metadata_physchem <- read.csv("physchem_metadata.csv")

abund_arc_bac <- read.csv("C_archaea_bacteria_abd_genus.csv", row.names = 1)
# Filtrar filas que tengan al menos un valor >= 0.1
abund_arc_bac_filtered <- abund_arc_bac[apply(abund_arc_bac >= 0.1, 1, any), ]
# Verificar el resultado
print(dim(abund_arc_bac_filtered))  # Muestra el número de filas y columnas
head(abund_arc_bac_filtered)  # Muestra las primeras filas del dataframe filtrado


# Transponer la tabla de abundancia
physchem_table_t <- physchem
physchem_df <- as.data.frame(physchem_table_t)
physchem_df$UID <- rownames(physchem_df)

# Unir las tablas de abundancia y metadatos
combined_data <- merge(physchem_df, metadata_physchem, by.x = "UID", by.y = "UID")

# Crear una lista de columnas de interés
columnas_de_interes <- c(
  3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

# Crear nuevas columnas combinadas de Sector-Season y Sector-Habitat
combined_data$Sector_Season <- paste(combined_data$Sector, combined_data$Season, sep = "_")
combined_data$Sector_Habitat <- paste(combined_data$Sector, combined_data$Habitat, sep = "_")

combined_data$Sector_Season <- factor(combined_data$Sector_Season) 
#                                      levels = c("A_Inlet_Relax", "A_Inlet_Intense",
#                                                 "B_Transition_Relax", "B_Transition_Intense",
#                                                 "C_Inner_Relax", "C_Inner_Intense"))
combined_data$Sector_Habitat <- factor(combined_data$Sector_Habitat)
#                                      levels = c("A_Inlet_Bare", "A_Inlet_Seagrass",
#                                                 "B_Transition_Bare", "B_Transition_Seagrass",
#                                                 "C_Inner_Bare", "C_Inner_Seagrass"))

# Crear las comparaciones adecuadas
a_my_comparisons <- list(
  c("A_Inlet_Relax", "A_Inlet_Intense"),
  c("B_Transition_Relax", "B_Transition_Intense"),
  c("C_Inner_Relax", "C_Inner_Intense")
)
  , 
  c("A_Inlet_Intense","B_Transition_Intense"),
  c("A_Inlet_Intense","C_Inner_Intense"),
  c("B_Transition_Intense","C_Inner_Intense"),
  c("A_Inlet_Relax","B_Transition_Relax"),
  c("A_Inlet_Relax","C_Inner_Relax"),
  c("B_Transition_Relax","C_Inner_Relax")
)

b_my_comparisons <- list(
  c("A_Inlet_Bare", "A_Inlet_Seagrass"),
  c("B_Transition_Bare", "B_Transition_Seagrass"),
  c("C_Inner_Bare", "C_Inner_Seagrass")
)
  , 
  c("A_Inlet_Seagrass","B_Transition_Seagrass"),
  c("A_Inlet_Seagrass","C_Inner_Seagrass"),
  c("B_Transition_Seagrass","C_Inner_Seagrass"),
  c("A_Inlet_Bare","B_Transition_Bare"),
  c("A_Inlet_Bare","C_Inner_Bare"),
  c("B_Transition_Bare","C_Inner_Bare")
)


## ------------------------------
## Comparando Sector-Season
## ------------------------------

# Bucle para generar y almacenar las figuras en variables a1, a2, a3, ...
for (i in 1:length(columnas_de_interes)) {
  columna_index <- columnas_de_interes[i]
  columna_de_interes <- colnames(combined_data)[columna_index]
  
  # Crear el boxplot con el t test usando la nueva columna combinada
  figura <- ggplot(combined_data, aes(x = Sector_Season, y = .data[[columna_de_interes]], 
                                      fill = Season)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
    labs(x = "", y = "", fill = "Season") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y= element_text(size=14),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("Relax" = "#000080ff", "Intense" = "#ff0000ff")) +
    stat_compare_means(method = "t.test", comparisons = a_my_comparisons, label = "p.format", 
                       position = position_dodge(0.75)) +
    ggtitle(paste(columna_de_interes))
  
  # Almacenar la figura en una variable dinámica (a1, a2, a3, ...)
  assign(paste0("a", i), figura)
  
  # Almacenar la figura en la lista
  #figuras[[paste0("a", i)]] <- figura
  
  # Guardar cada figura en un archivo
  #ggsave(paste0("figura_", i, ".png"), plot = figura, width = 10, height = 7)
  
  
}

svglite("annot_physchem_v5.svg", width=20, height=10)
((a1 | a2 | a3 | a4 | a5) / ( a6 | a7 | a8 | a9 | a10))
dev.off()


# Con puntos en lugar de boxplot
for (i in seq_along(columnas_de_interes)) {
  idx <- columnas_de_interes[i]
  var <- colnames(combined_data)[idx]
  
  p <- ggplot(combined_data, 
              aes(x = Sector_Season, y = .data[[var]], 
                  colour = Season, shape  = Habitat)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width  = 0.75),
               size = 2, alpha = 0.8) +
    labs(x = NULL, y = NULL, colour = "Season", shape = "Habitat") +
    theme_bw() +
    theme(axis.text.x        = element_text(angle = 45, hjust = 1),
          axis.text.y        = element_text(size = 14),
          plot.title         = element_text(hjust = 0.5),
          legend.title       = element_text(size = 12),
          legend.text        = element_text(size = 10)) +
    scale_color_manual(values = c("Relax"   = "#000080FF", 
                                  "Intense" = "#FF0000FF")) +
    stat_compare_means(method    = "t.test", 
                       comparisons = a_my_comparisons, 
                       label      = "p.format",
                       position   = position_dodge(width = 0.75)) +
    ggtitle(var)
  
  assign(paste0("a", i), p)
}

svglite("annot_physchem_v7.svg", width=20, height=10)
((a1 | a2 | a3 | a4 | a5) / ( a6 | a7 | a8 | a9 | a10))
dev.off()


## Statistical analysis

library(FSA)           # Para la prueba de Dunn
library(multcompView)  # Para generar letras de grupos

# Realizar el test de Kruskal-Wallis
kruskal_result <- kruskal.test(Fe_II_mg ~ Sector_Season, data = combined_data)

# Mostrar resultados
print(kruskal_result)

# Prueba de Dunn con corrección de Benjamini-Hochberg
dunn_result <- dunnTest(Sand ~ Sector_Season, data = combined_data, method = "bh")
dunn_result 

# Acceder a las comparaciones y p-values ajustados
comparisons <- dunn_result$res$Comparison
p_values <- dunn_result$res$P.adj

# Crear una matriz de valores p
grupos <- unique(combined_data$Sector_Season)
p_matrix <- matrix(1, nrow = length(grupos), ncol = length(grupos))
rownames(p_matrix) <- colnames(p_matrix) <- grupos

# Llenar la matriz con los valores p ajustados
for (i in 1:length(comparisons)) {
  pair <- unlist(strsplit(comparisons[i], " - "))
  p_matrix[pair[1], pair[2]] <- p_values[i]
  p_matrix[pair[2], pair[1]] <- p_values[i]
}

# Convertir los valores p en letras de significancia
letters_result <- multcompLetters(p_matrix)$Letters
print(letters_result)






library(rcompanion)
srh <- scheirerRayHare(Fe_II_mg ~ Sector + Season, data = combined_data)
srh <- scheirerRayHare(pH ~ Sector + Season, data = combined_data)
print(srh)


library(ARTool)

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Ajusta el modelo con los tres factores y todas sus interacciones
art_mod <- art(Fe_II_mg ~ Sector * Season * Habitat, data = combined_data)
art_mod <- art(pH ~ Sector * Season * Habitat, data = combined_data)

# Obtén la tabla ANOVA no paramétrica
anova(art_mod)





library(FSA)          # para dunnTest
library(multcompView) # para convertir p‑values a letras

# Dunn post‑hoc entre niveles de Sector y/o Season
dunn_t <- dunnTest(pH ~ Sector,
                        data   = combined_data,
                        method = "bh")

print(dunn_t)    # tabla con Z, p.crudo y P.adj (BH)

# Two‑level post‑hoc for Season
wilcox_t <- wilcox.test(pH ~ Season,
            data = combined_data,
            exact = FALSE)          # exact=FALSE avoids warnings with ties


print(wilcox_t)    # tabla con Z, p.crudo y P.adj (BH)

## ⇒ Obtener letras de significancia ##
comp    <- dunn_sector$res$Comparison
p_adj   <- dunn_sector$res$P.adj
combined_data$Sector <- factor(combined_data$Sector)
sects   <- levels(combined_data$Sector)

# Matriz de p‑values
p_mat <- matrix(1, nrow = length(sects), ncol = length(sects),
                dimnames = list(sects, sects))
for (i in seq_along(comp)) {
  pair <- strsplit(comp[i], " - ")[[1]]
  p_mat[pair[1], pair[2]] <- p_adj[i]
  p_mat[pair[2], pair[1]] <- p_adj[i]
}

letters <- multcompLetters(p_mat)$Letters
print(letters)        # A, B, AB… para anotar en la figura




library(ggplot2)
library(dplyr)
library(tidyr)

# Repetimos los pasos de transformación a formato largo
df_long <- combined_data %>% 
  dplyr::select(Fe_II_mg, Sector, Season, Habitat) %>%           
  pivot_longer(cols = c(Sector, Season, Habitat),
               names_to  = "GroupType",
               values_to = "Group")

# Fijamos el orden de las siete categorías en el eje x
level_order <- c("A_Inlet", "B_Transition", "C_Inner",
                 "Intense",  "Relax",
                 "Bare",     "Seagrass")
df_long$Group <- factor(df_long$Group, levels = level_order)

# Mapeo de colores personalizado
my_colors <- c(
  A_Inlet      = "#56b4e9",
  B_Transition = "#808080",
  C_Inner      = "#c8ab37",  
  Intense      = "#FF0000FF",
  Relax        = "#000080FF",
  Bare         = "orange",
  Seagrass     = "darkgreen"
)

# Plot
ggplot(df_long, aes(x = Group, y = Fe_II_mg, fill = Group)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, colour = "black") +
  scale_fill_manual(values = my_colors) +
  labs(x = NULL, y = "Fe_II_mg", fill = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # quita la leyenda si no la necesitas
  )






## ------------------------------
## Comparando Sector-Habitat
## ------------------------------

# Bucle para generar y almacenar las figuras en variables a1, a2, a3, ...
for (i in 1:length(columnas_de_interes)) {
  columna_index <- columnas_de_interes[i]
  columna_de_interes <- colnames(combined_data)[columna_index]
  
  # Crear el boxplot con el t test usando la nueva columna combinada
  figura <- ggplot(combined_data, aes(x = Sector_Habitat, y = .data[[columna_de_interes]], fill = Habitat)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
    labs(x = "", y = "", fill = "Habitat") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y= element_text(size=14),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("Bare" = "orange", "Seagrass" = "darkgreen")) +
    stat_compare_means(method = "t.test", comparisons = b_my_comparisons, label = "p.format", 
                       position = position_dodge(0.75)) +
    ggtitle(paste(columna_de_interes))
  
  # Almacenar la figura en una variable dinámica (a1, a2, a3, ...)
  assign(paste0("a", i), figura)
  
  # Almacenar la figura en la lista
  #figuras[[paste0("a", i)]] <- figura
  
  # Guardar cada figura en un archivo
  #ggsave(paste0("figura_", i, ".png"), plot = figura, width = 10, height = 7)
  
  
}

svglite("annot_physchem_v6.svg", width=20, height=10)
((a1 | a2 | a3 | a4 | a5) / ( a6 | a7 | a8 | a9 | a10))
dev.off()


## Statistical analysis

library(FSA)           # Para la prueba de Dunn
library(multcompView)  # Para generar letras de grupos

# Realizar el test de Kruskal-Wallis
kruskal_result <- kruskal.test(Sand ~ Sector_Habitat, data = combined_data)

# Mostrar resultados
print(kruskal_result)

# Prueba de Dunn con corrección de Benjamini-Hochberg
dunn_result <- dunnTest(Sand ~ Sector_Habitat, data = combined_data, method = "bh")
dunn_result 

# Acceder a las comparaciones y p-values ajustados
comparisons <- dunn_result$res$Comparison
p_values <- dunn_result$res$P.adj

# Crear una matriz de valores p
grupos <- unique(combined_data$Sector_Habitat)
p_matrix <- matrix(1, nrow = length(grupos), ncol = length(grupos))
rownames(p_matrix) <- colnames(p_matrix) <- grupos

# Llenar la matriz con los valores p ajustados
for (i in 1:length(comparisons)) {
  pair <- unlist(strsplit(comparisons[i], " - "))
  p_matrix[pair[1], pair[2]] <- p_values[i]
  p_matrix[pair[2], pair[1]] <- p_values[i]
}

# Convertir los valores p en letras de significancia
letters_result <- multcompLetters(p_matrix)$Letters
print(letters_result)




## Comparar pH, Habitats

# Realizar el test de Kruskal-Wallis
subset_data <- combined_data %>% filter(Sector %in% c("A_Inlet", "C_Inner"))
kruskal_result <- kruskal.test(pH ~ Sector, data = subset_data)
print(kruskal_result)
kruskal_result <- kruskal.test(pH ~ Habitat, data = subset_data)
print(kruskal_result)
kruskal_result <- kruskal.test(SOC_umol ~ Sector, data = subset_data)
print(kruskal_result)
kruskal_result <- kruskal.test(NT_umol ~ Sector, data = subset_data)
print(kruskal_result)

kruskal_result <- kruskal.test(pH ~ Habitat, data = combined_data)
print(kruskal_result)
kruskal_result <- kruskal.test(SOC_umol ~ Habitat, data = combined_data)
print(kruskal_result)
kruskal_result <- kruskal.test(NH4_g ~ Habitat, data = combined_data)
print(kruskal_result)
kruskal_result <- kruskal.test(NO2_ug ~ Habitat, data = combined_data)
print(kruskal_result)
kruskal_result <- kruskal.test(NO3_ug ~ Habitat, data = combined_data)
print(kruskal_result)
kruskal_result <- kruskal.test(Fe_III_mg ~ Habitat, data = combined_data)
print(kruskal_result)
kruskal_result <- kruskal.test(Fe_II_mg ~ Habitat, data = combined_data)
print(kruskal_result)





#########################
## REDUNDANCY ANALYSIS
#########################


## ------------------------------------------------------
##
## MICROBIOME AND PHYSICOCHEMICAL FEATURES
##
## ------------------------------------------------------


# Load tables
#abund_arc_bac <- read.csv("C_archaea_bacteria_abd_genus.csv", row.names = 1)
#abund_arc_bac <- read.csv("C_arc_bac_abd_species.csv", row.names = 1)
abund_arc_bac <- read.csv("C_archaea_bacteria_count_genus.csv", row.names = 1)
#abund_arc_bac <- read.csv("C_archaea_bacteria_count_species.csv", row.names = 1)

# Read metadata
metadata <- read.csv("metadata_quintin.csv")
#mtd_chem <- read.csv("physchem_metadata.csv", header = T)
data_chem <- read.csv("physchem_v1.csv", header = T)

# summary statistics for all objects (min, mean, max, etc.)
summary(abund_arc_bac)


## ---------------------------------------------------------
##
## Preparing data for analyses
##
## ---------------------------------------------------------

# Filtrar nombres de columnas que no sean "UID", "sample" ni contengan "Rate"
data_chem_num <- data_chem[, !grepl("Rate|UID|Sample", names(data_chem))]
summary(data_chem_num)

vars <- c("pH", "Sand", "Silt", "TOC_umol", "TN_umol",
          "NH4_g", "NO3_ug", "NO2_ug", "Fe_III_mg", "Fe_II_mg")

data_chem_num <- data_chem_num %>% dplyr::select(all_of(vars))


## ---------------------------------------------------------
##
## Pearson correlation to reduce collinearity
##
## ---------------------------------------------------------


# Pairwise Pearson correlations among environmental variables to identify and discard those 
# with a |pairwise correlation| higher than 0.7
cor_matrix_bac <- cor(data_chem_num) # Correlation matrix

# Collinearity in the physicochemical data
# We can visually look for correlations between variables:

#svglite("pearson_phychem_silvia.svg", width = 6, height = 5)
heatmap(abs(cor_matrix_bac), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
#dev.off()

# Filter with Pearson > 0.7
cor_matrix_bac_rm <- cor_matrix_bac
cor_matrix_bac_rm[upper.tri(cor_matrix_bac_rm)] <- 0
diag(cor_matrix_bac_rm) <- 0
data_chem_num.new <- data_chem_num[ , !apply(cor_matrix_bac_rm,    # Remove highly correlated variables
                                             2,
                                             function(x) any(x > 0.7))]

# Less correlated variables (Pearson < 0.7)
less.cor.var <- colnames(data_chem_num.new)
less.cor.var
## Seven parameters: "Sand", "TN_umol", "NH4_g", "NO3_ug", "NO2_ug", "Fe_III_mg", "Fe_II_mg" 

# Como hay triplicados por cada muestra, se calcula el promedio de los parametros
data_chem_tr <- aggregate(data_chem_num.new, by = list(Sample = data_chem$Sample), FUN = mean)
data_chem_tr <- data_chem_tr %>%
  mutate(across(.cols = 2:ncol(data_chem_tr), .fns = ~ as.numeric(as.character(.))))
rownames(data_chem_tr) <- data_chem_tr$Sample
data_chem_tr$Sample <- NULL

# Tambien hago de data_chem_num
data_chem_num_tr <- aggregate(data_chem_num, by = list(Sample = data_chem$Sample), FUN = mean)
data_chem_num_tr <- data_chem_num_tr %>%
  mutate(across(.cols = 2:ncol(data_chem_num_tr), .fns = ~ as.numeric(as.character(.))))
rownames(data_chem_num_tr) <- data_chem_num_tr$Sample
data_chem_num_tr$Sample <- NULL



## ---------------------------------------------------------
##
## Using VIF approach
##
## ---------------------------------------------------------


# Scale and center variables
data_chem_num.new.z <- decostand(data_chem_tr, method = "standardize")

# Apply Hellinger transformation to correct for the double zero problem in abundance table
ab_genus_rda <- t(abund_arc_bac)
ab_genus_hellinger <- as.data.frame(decostand(ab_genus_rda, method = "hellinger"))

# Full RDA model
rda_ab_chem <- rda(ab_genus_hellinger ~ ., data = data_chem_num.new.z)
summary(rda_ab_chem)
##              Inertia Proportion
# Total          0.0885      1.000
# Constrained    0.0727      0.821
# Unconstrained  0.0158      0.179

vif_criteria <- vif.cca(rda_ab_chem)
vif_criteria
# Sand   NT_umol     NH4_g    NO3_ug    NO2_ug Fe_III_mg  Fe_II_mg 
# 8.302     6.810     7.968     2.347     6.523     8.122    10.476 

vars.envfit <- names(which(vif_criteria < 10))

data_chem_tr_2 <- data_chem_tr[, vars.envfit]
  
# Scale and center variables with less correlated variables
data_chem_num.new.z_2 <- decostand(data_chem_tr_2, method = "standardize")

# Full RDA model with less correlated variables
rda_ab_chem_2 <- rda(ab_genus_hellinger ~ ., data = data_chem_num.new.z_2)
summary(rda_ab_chem_2)
#               Inertia Proportion
# Total         0.08855     1.0000
# Constrained   0.06585     0.7436
# Unconstrained 0.02270     0.2564


## ---------------------------------------------------------
##
## Using ordistep approach
##
## ---------------------------------------------------------

# Perform forward and backward selection of explanatory variables
## Based on the workshop: https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
set.seed(123)
fwd.sel <- ordistep(rda(ab_genus_hellinger ~ 1, data = data_chem_num.new.z_2), # lower model limit (simple!)
                    scope = formula(rda_ab_chem_2), # upper model limit (the "full" model)
                    direction = "both",
                    alpha = 0.01,
                    R2scope = FALSE, # can surpass the "full" model's R2
                    pstep = 999,
                    trace = TRUE) # TRUE to see the selection process!

#             Df     AIC      F Pr(>F)   
# NO3_ug     1 -31.060 5.0747  0.010 **
# Fe_III_mg  1 -28.850 2.5393  0.045 * 
# Sand       1 -28.087 1.7671  0.105   
# NO2_ug     1 -27.976 1.6583  0.155   
# TN_umol    1 -26.768 0.5424  0.795   
# NH4_g      1 -26.426 0.2460  0.960  

# Get the formula
fwd.sel$call
# rda(formula = ab_genus_hellinger ~ NO3_ug + Fe_III_mg, data = data_chem_num.new.z)

# Get the selected variables
vars.ordistep <- gsub('^. ', '', rownames(fwd.sel$anova))
vars.ordistep


# Get the statistics adding the collineal variables
core.vars <- less.cor.var
extra.vars  <- setdiff(colnames(data_chem_num_tr), core.vars) 

core.form  <- as.formula(
  paste("ab_genus_hellinger ~", paste(core.vars, collapse = " + "))
)

scope.form <- as.formula(
  paste("ab_genus_hellinger ~", paste(c(core.vars, extra.vars), collapse = " + "))
)

## modelo base sobre data_chem_num_tr  (datos = todas las columnas)
rda_core   <- rda(core.form, data = data_chem_num_tr)

## selección paso a paso
set.seed(123)
step.mod <- ordistep(rda_core,
                     scope      = scope.form,
                     direction  = "both",
                     pstep      = 999,
                     trace      = TRUE)
# Df     AIC      F Pr(>F)
# TOC_umol  1 -39.059 2.0546  0.190
# Silt      1 -37.126 1.3029  0.290
# pH        1 -35.428 0.7349  0.495


## ---------------------------------------------------------
##
## Combining vif and ordistep approaches
##
## ---------------------------------------------------------

# Get the common variables from the vif and ordistep aproaches
vars.bac <- intersect(vars.ordistep, vars.envfit)
vars.bac
# "NO3_ug"    "Fe_III_mg"

# Write our new model
matrix.kept <- data_chem_num.new.z[, vars.bac]
rda_ab_chem.signif <- rda(ab_genus_hellinger ~ ., data = matrix.kept)

vif_criteria <- vif.cca(rda_ab_chem.signif)
vif_criteria

summary(rda_ab_chem.signif)
## Constrained   0.04458     Proportion: 0.5035
## Unconstrained 0.04397     Proportion: 0.4965

# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_ab_chem.signif)
## r.squared = 0.5034835
## adj.r.squared = 0.3931465

set.seed(123)
anova(rda_ab_chem.signif)
## Model     2 0.044583 4.5631  0.002 **

set.seed(123)
anova(rda_ab_chem.signif,by="axis")
## RDA1: p-values = 0.002
## RDA2: p-values = 0.009

set.seed(123)
anova(rda_ab_chem.signif,by="term")
## NO3_ug: p-values = 0.001
## Fe_III_mg: p-values = 0.017

## ---------------------------------------------------------
##
## Ploting RDA analysis
##
## ---------------------------------------------------------

# Fusionar la información de los sectores con los datos transformados
data_chem_num.new.z <- cbind(data_chem_num.new.z, Sector = metadata$Sector)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Habitat = metadata$Habitat)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Season = metadata$Season)

# Definir formas, colores, y tipos de línea
sector_shapes <- c("A_Inlet" = 21, "B_Transition" = 22, "C_Inner" = 23)
habitat_colors <- c("Bare" = "orange", "Seagrass" = "darkgreen")
season_linetype <- c("Intense" = "solid", "Relax" = "dashed")

explained_var <- summary(rda_ab_chem)$cont$importance[2, ] * 100
xlab_text <- paste0("RDA1 (", round(explained_var[1], 2), "%)")
ylab_text <- paste0("RDA2 (", round(explained_var[2], 2), "%)")


# Iniciar la salida en formato SVG
#svglite("rda_silvia_final2.svg", width = 8, height = 5)

# Crear el plot de RDA
plot(rda_ab_chem.signif, type = "n", main = "Redundance Analysis - chemical",
     xlim = c(-2, 2), ylim = c(-1, 1), xlab = xlab_text, ylab = ylab_text)

# Añadir puntos con formas, colores, y líneas
with(data_chem_num.new.z, points(rda_ab_chem.signif, "sites", bg = habitat_colors[Habitat], cex = 2, 
                            pch = sector_shapes[Sector], 
                            lwd = 1.5, # Ajustar el grosor de las líneas
                            col = "black", # Línea negra para bordes
                            lty = season_linetype[Season]))

# Añadir texto y etiquetas
text(rda_ab_chem.signif, display = "sites", cex = .5, col = "gray")
text(rda_ab_chem.signif, display = "cn", cex = .9, col = "blue")

# Añadir la leyenda
legend("topright", legend = c("Bare", "Seagrass"), 
       fill = c("orange", "darkgreen"), title = "Habitat", box.lwd = 0)
legend("topright", legend = c("A_Inlet", "B_Transition", "C_Inner"),
       pch = c(21, 22, 23), title = "Sector", inset = c(0.1, 0.2))
legend("topright", legend = c("Intense", "Relax"), 
       lty = c(1, 2), col = "black", title = "Season", inset = c(0.1, 0.4), box.lwd = 0)

#dev.off()


# Obtener coordenadas de las variables explicativas en el espacio de RDA
env_coords <- scores(rda_ab_chem.signif, display = "bp", scaling = 2)

env_coords

# Función para calcular ángulo entre dos vectores
angle_between <- function(a, b) {
  cos_theta <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  theta_rad <- acos(cos_theta)
  theta_deg <- theta_rad * 180 / pi
  return(theta_deg)
}

# Calcular ángulo entre NO3_ug y Fe_III_mg
vec1 <- env_coords["NO3_ug", ]
vec2 <- env_coords["Fe_III_mg", ]

angle_between(vec1, vec2)










## ------------------------------------------------------
##
## GENES AND PHYSICOCHEMICAL FEATURES
##
## ------------------------------------------------------


# Definir el directorio
base_path <- "/Users/jorgerv/Documents/lagunas_costeras/Results_databases/Second_version/17-28/"

# Load tables
count_metal <- as.data.frame(t(read.csv(paste0(base_path, "count_BACMET_Metal_genes.csv"), sep = ",", row.names = 1)))
count_drug <- as.data.frame(t(read.csv(paste0(base_path, "count_CARD_Drug_Class_genes.csv"), sep = ",", row.names = 1)))
count_nitrogen <- as.data.frame(t(read.csv(paste0(base_path, "count_NITROGEN_pathway_genes.csv"), sep = ",", row.names = 1)))
count_sulfur <- as.data.frame(t(read.csv(paste0(base_path, "count_SULFUR_pathway_genes.csv"), sep = ",", row.names = 1)))
count_virulence <- as.data.frame(t(read.csv(paste0(base_path, "count_VFDB_type_genes.csv"), sep = ",", row.names = 1)))

samples <- c("C_17", "C_18", "C_19", "C_20", "C_21", "C_22", "C_23", "C_24", "C_25", "C_26", "C_27", "C_28")

count_metal <- count_metal[, samples]
count_drug <- count_drug[, samples]
count_nitrogen <- count_nitrogen[, samples]
count_sulfur <- count_sulfur[, samples]
count_virulence <- count_virulence[, samples]

## Tables for UpSet plots
upset_METAL_genes     <- as.data.frame(+ (count_metal     > 0))
upset_DRUG_genes      <- as.data.frame(+ (count_drug      > 0))
upset_NITROGEN_genes  <- as.data.frame(+ (count_nitrogen  > 0))
upset_SULFUR_genes    <- as.data.frame(+ (count_sulfur    > 0))
upset_VFDB_genes <- as.data.frame(+ (count_virulence > 0))

save(upset_METAL_genes, file = "upset_METAL_genes.RData")
save(upset_DRUG_genes, file = "upset_DRUG_genes.RData")
save(upset_NITROGEN_genes, file = "upset_NITROGENL_genes.RData")
save(upset_SULFUR_genes, file = "upset_SULFUR_genes.RData")
save(upset_VFDB_genes, file = "upset_VFDB_genes.RData")

# Bind all tables
count_all <- rbind(count_metal,
                   count_drug,
                   count_nitrogen,
                   count_sulfur,
                   count_virulence)

abund_arc_bac <- count_all

# Read metadata
metadata <- read.csv("metadata_quintin.csv")
#mtd_chem <- read.csv("physchem_metadata.csv", header = T)
data_chem <- read.csv("physchem_v1.csv", header = T)

# summary statistics for all objects (min, mean, max, etc.)
summary(abund_arc_bac)


## ---------------------------------------------------------
##
## Preparing data for analyses
##
## ---------------------------------------------------------

# Filtrar nombres de columnas que no sean "UID", "sample" ni contengan "Rate"
data_chem_num <- data_chem[, !grepl("Rate|UID|Sample", names(data_chem))]
summary(data_chem_num)

vars <- c("pH", "Sand", "Silt", "TOC_umol", "TN_umol",
          "NH4_g", "NO3_ug", "NO2_ug", "Fe_III_mg", "Fe_II_mg")

data_chem_num <- data_chem_num %>% dplyr::select(all_of(vars))


## ---------------------------------------------------------
##
## Pearson correlation to reduce collinearity
##
## ---------------------------------------------------------


# Pairwise Pearson correlations among environmental variables to identify and discard those 
# with a |pairwise correlation| higher than 0.7
cor_matrix_bac <- cor(data_chem_num) # Correlation matrix

# Collinearity in the physicochemical data
# We can visually look for correlations between variables:

#svglite("pearson_phychem_silvia.svg", width = 6, height = 5)
heatmap(abs(cor_matrix_bac), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
#dev.off()

# Filter with Pearson > 0.7
cor_matrix_bac_rm <- cor_matrix_bac
cor_matrix_bac_rm[upper.tri(cor_matrix_bac_rm)] <- 0
diag(cor_matrix_bac_rm) <- 0
data_chem_num.new <- data_chem_num[ , !apply(cor_matrix_bac_rm,    # Remove highly correlated variables
                                             2,
                                             function(x) any(x > 0.7))]

# Less correlated variables (Pearson < 0.7)
less.cor.var <- colnames(data_chem_num.new)
less.cor.var
## Seven parameters: "Sand"      "TN_umol"   "NH4_g"     "NO3_ug"    "NO2_ug"    "Fe_III_mg" "Fe_II_mg" 

# Como hay triplicados por cada muestra, se calcula el promedio de los parametros
data_chem_tr <- aggregate(data_chem_num.new, by = list(Sample = data_chem$Sample), FUN = mean)
data_chem_tr <- data_chem_tr %>%
  mutate(across(.cols = 2:ncol(data_chem_tr), .fns = ~ as.numeric(as.character(.))))
rownames(data_chem_tr) <- data_chem_tr$Sample
data_chem_tr$Sample <- NULL

# Tambien hago de data_chem_num
data_chem_num_tr <- aggregate(data_chem_num, by = list(Sample = data_chem$Sample), FUN = mean)
data_chem_num_tr <- data_chem_num_tr %>%
  mutate(across(.cols = 2:ncol(data_chem_num_tr), .fns = ~ as.numeric(as.character(.))))
rownames(data_chem_num_tr) <- data_chem_num_tr$Sample
data_chem_num_tr$Sample <- NULL





## ---------------------------------------------------------
##
## Using VIF approach
##
## ---------------------------------------------------------


# Scale and center variables
data_chem_num.new.z <- decostand(data_chem_tr, method = "standardize")

# Apply Hellinger transformation to correct for the double zero problem in abundance table
ab_genus_rda <- t(abund_arc_bac)
ab_genus_hellinger <- as.data.frame(decostand(ab_genus_rda, method = "hellinger"))

# Full RDA model
rda_ab_chem <- rda(ab_genus_hellinger ~ ., data = data_chem_num.new.z)
summary(rda_ab_chem)
## With the seven
##              Inertia Proportion

# Total         0.05219      1.000
# Constrained   0.03977      0.762
# Unconstrained 0.01242      0.238

vif_criteria <- vif.cca(rda_ab_chem)
vif_criteria
# Sand   NT_umol     NH4_g    NO3_ug    NO2_ug Fe_III_mg  Fe_II_mg 
# 8.302483  6.810308  7.968274  2.346754  6.522657  8.121830 10.475641

vars.envfit <- names(which(vif_criteria < 10))
vars.envfit
# "Sand"      "TN_umol"   "NH4_g"     "NO3_ug"    "NO2_ug"    "Fe_III_mg"

data_chem_tr_2 <- data_chem_tr[, vars.envfit]

# Scale and center variables with less correlated variables
data_chem_num.new.z_2 <- decostand(data_chem_tr_2, method = "standardize")

# Full RDA model with less correlated variables
rda_ab_chem_2 <- rda(ab_genus_hellinger ~ ., data = data_chem_num.new.z_2)
summary(rda_ab_chem_2)
#               Inertia Proportion
# Total         0.05219     1.0000
# Constrained   0.03553     0.6807
# Unconstrained 0.01666     0.3193

# final model with the variables that passed r and VIF filters
rda_final <- rda(ab_genus_hellinger ~ NO3_ug + Fe_III_mg + Sand +
                   NO2_ug + TN_umol + NH4_g,
                 data = data_chem_num.new.z_2)

# marginal (type = "terms") permutation test, 999 permutations
perm_test <- anova.cca(rda_final,
                       by           = "term",
                       permutations = 999)
print(perm_test)



## ---------------------------------------------------------
##
## Using ordistep approach
##
## ---------------------------------------------------------

# Perform forward and backward selection of explanatory variables
## Based on the workshop: https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
set.seed(123)
fwd.sel <- ordistep(rda(ab_genus_hellinger ~ 1, data = data_chem_num.new.z_2), # lower model limit (simple!)
                    scope = formula(rda_ab_chem_2), # upper model limit (the "full" model)
                    direction = "both",
                    R2scope = FALSE, # can surpass the "full" model's R2
                    pstep = 999,
                    trace = TRUE) # TRUE to see the selection process!

#             Df     AIC      F Pr(>F)   
# NO3_ug     1 -35.352 2.7064  0.005 **
# Fe_III_mg  1 -34.574 1.9085  0.045 * 
# Sand       1 -34.296 1.6353  0.075 . 
# NO2_ug     1 -33.608 0.9877  0.435   
# TN_umol    1 -33.478 0.8691  0.595   
# NH4_g      1 -33.157 0.5823  0.880    

# Get the formula
fwd.sel$call
# rda(formula = ab_genus_hellinger ~ NO3_ug + Fe_III_mg, data = data_chem_num.new.z)

# Get the selected variables
vars.ordistep <- gsub('^. ', '', rownames(fwd.sel$anova))
vars.ordistep
# "NO3_ug"    "Fe_III_mg"



# Get the statistics adding the collineal variables
core.vars <- less.cor.var
extra.vars  <- setdiff(colnames(data_chem_num_tr), core.vars) 

core.form  <- as.formula(
  paste("ab_genus_hellinger ~", paste(core.vars, collapse = " + "))
)

scope.form <- as.formula(
  paste("ab_genus_hellinger ~", paste(c(core.vars, extra.vars), collapse = " + "))
)

## modelo base sobre data_chem_num_tr  (datos = todas las columnas)
rda_core   <- rda(core.form, data = data_chem_num_tr)

## selección paso a paso
set.seed(123)
step.mod <- ordistep(rda_core,
                     scope      = scope.form,
                     direction  = "both",
                     pstep      = 999,
                     trace      = TRUE)

#           Df     AIC      F Pr(>F)
# TOC_umol  1 -40.178 1.3562  0.345
# Silt      1 -38.963 0.9369  0.420
# pH        1 -38.812 0.8875  0.480



## ---------------------------------------------------------
##
## Combining vif and ordistep approaches
##
## ---------------------------------------------------------

# Get the common variables from the vif and ordistep aproaches
vars.bac <- intersect(vars.ordistep, vars.envfit)
vars.bac
# "NO3_ug"    "Fe_III_mg"

# Write our new model
matrix.kept <- data_chem_num.new.z[, vars.bac]
rda_ab_chem.signif <- rda(ab_genus_hellinger ~ ., data = matrix.kept)

summary(rda_ab_chem.signif)
##              Inertia Proportion
# Total         0.05219     1.0000
# Constrained   0.01864     0.3571
# Unconstrained 0.03355     0.6429

# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_ab_chem.signif)
## r.squared = 0.3571451
## adj.r.squared = 0.2142885

set.seed(123)
anova(rda_ab_chem.signif)
##           Df Variance     F Pr(>F)    
## Model     2 0.018640 2.5  0.002 **
## Residual  9 0.033551

set.seed(123)
anova(rda_ab_chem.signif,by="axis")
## RDA1: p-values = 0.003
## RDA2: p-values = 0.025

set.seed(123)
anova(rda_ab_chem.signif,by="term")
## NO3_ug: p-values = 0.002
## Fe_III_mg: p-values = 0.020


## ---------------------------------------------------------
##
## Ploting RDA analysis
##
## ---------------------------------------------------------

# Fusionar la información de los sectores con los datos transformados
data_chem_num.new.z <- cbind(data_chem_num.new.z, Sector = metadata$Sector)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Habitat = metadata$Habitat)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Season = metadata$Season)

# Definir formas, colores, y tipos de línea
sector_shapes <- c("A_Inlet" = 21, "B_Transition" = 22, "C_Inner" = 23)
habitat_colors <- c("Bare" = "orange", "Seagrass" = "darkgreen")
season_linetype <- c("Intense" = "solid", "Relax" = "dashed")

explained_var <- summary(rda_ab_chem)$cont$importance[2, ] * 100
xlab_text <- paste0("RDA1 (", round(explained_var[1], 2), "%)")
ylab_text <- paste0("RDA2 (", round(explained_var[2], 2), "%)")


# Iniciar la salida en formato SVG
#svglite("rda_silvia_genes2.svg", width = 8, height = 5)

# Crear el plot de RDA
plot(rda_ab_chem.signif, type = "n", main = "Redundance Analysis - chemical",
     xlim = c(-2, 2), ylim = c(-1, 1), xlab = xlab_text, ylab = ylab_text)

# Añadir puntos con formas, colores, y líneas
with(data_chem_num.new.z, points(rda_ab_chem.signif, "sites", bg = habitat_colors[Habitat], cex = 2, 
                                 pch = sector_shapes[Sector], 
                                 lwd = 1.5, # Ajustar el grosor de las líneas
                                 col = "black", # Línea negra para bordes
                                 lty = season_linetype[Season]))

# Añadir texto y etiquetas
text(rda_ab_chem.signif, display = "sites", cex = .5, col = "gray")
text(rda_ab_chem.signif, display = "cn", cex = .9, col = "blue")

# Añadir la leyenda
legend("topright", legend = c("Bare", "Seagrass"), 
       fill = c("orange", "darkgreen"), title = "Habitat", box.lwd = 0)
legend("topright", legend = c("A_Inlet", "B_Transition", "C_Inner"),
       pch = c(21, 22, 23), title = "Sector", inset = c(0.1, 0.2))
legend("topright", legend = c("Intense", "Relax"), 
       lty = c(1, 2), col = "black", title = "Season", inset = c(0.1, 0.4), box.lwd = 0)

#dev.off()


# Obtener coordenadas de las variables explicativas en el espacio de RDA
env_coords <- vegan::scores(rda_ab_chem.signif, display = "bp", scaling = 2)

env_coords

# Función para calcular ángulo entre dos vectores
angle_between <- function(a, b) {
  cos_theta <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  theta_rad <- acos(cos_theta)
  theta_deg <- theta_rad * 180 / pi
  return(theta_deg)
}

# Calcular ángulo entre NO3_ug y Fe_III_mg
vec1 <- env_coords["NO3_ug", ]
vec2 <- env_coords["Fe_III_mg", ]

angle_between(vec1, vec2)
# 112.9961




## ----------------------------------------------------------------
## 0.  Packages  (vegan + tidyverse if you want dplyr syntax)
## ----------------------------------------------------------------
library(vegan)
library(dplyr)

## ----------------------------------------------------------------
## 1.  Abundance matrix   (Hellinger-transformed, rows = samples)
## ----------------------------------------------------------------
ab_genus_rda       <- t(abund_arc_bac)                   # samples × genera
ab_genus_hellinger <- decostand(ab_genus_rda, "hellinger")

## ----------------------------------------------------------------
## 2.  Environmental table with ALL candidate variables ----------
##     (already averaged triplicates and numeric!)
## ----------------------------------------------------------------
#  data_chem_tr  -> contains the per-sample means
#  less.cor.var  -> variables that passed the |r| < 0.7 filter
core.vars   <- less.cor.var                                        # 7 low-collinear vars
extra.vars  <- setdiff(colnames(data_chem_tr), core.vars)          # the “collinear” ones

## 2a.  Build the core table and standardise (= z-score)
env_core <- data_chem_tr %>% dplyr::select(all_of(core.vars))
env_core <- decostand(env_core, "standardize")

## 2b.  (Optional) build a table with ALL vars, but keep order:
env_all  <- cbind(env_core,
                  data_chem_tr %>% dplyr::select(all_of(extra.vars)) %>%
                    decostand("standardize"))

## ----------------------------------------------------------------
## 3.  RDA with core variables only
## ----------------------------------------------------------------
rda_core <- rda(ab_genus_hellinger ~ ., data = env_core)
summary(rda_core)

## ----------------------------------------------------------------
## 4.  Stepwise addition of the “extra” variables
##     – starts from the core model
## ----------------------------------------------------------------
set.seed(123)                      # reproducible permutations!
## 1. fórmulas
core.form  <- as.formula(
  paste("ab_genus_hellinger ~", paste(core.vars, collapse = " + "))
)

scope.form <- as.formula(
  paste("ab_genus_hellinger ~", paste(c(core.vars, extra.vars), collapse = " + "))
)

## 2. modelo base sobre env_all  (datos = todas las columnas)
rda_core   <- rda(core.form, data = env_all)

## 3. selección paso a paso
set.seed(123)
step.mod <- ordistep(rda_core,
                     scope      = scope.form,
                     direction  = "both",
                     pstep      = 999,
                     trace      = TRUE)

rda_full <- rda(ab_genus_hellinger ~ ., data = env_all)
anova(rda_full, by = "term")     # F y p para cada variable en el modelo completo

#           Df  Variance      F Pr(>F)  
#Sand       1 0.0132977 3.7266  0.063 .
#TN_umol    1 0.0078939 2.2122  0.182  
#NH4_g      1 0.0070329 1.9709  0.237  
#NO3_ug     1 0.0167737 4.7007  0.034 *
#NO2_ug     1 0.0039565 1.1088  0.436  
#Fe_III_mg  1 0.0168921 4.7339  0.033 *
#Fe_II_mg   1 0.0068775 1.9274  0.244  
#Silt       1 0.0047916 1.3428  0.351  
#TOC_umol   1 0.0039886 1.1178  0.442  
#pH         1 0.0034761 0.9742  0.528  
#Residual   1 0.0035683                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## ----------------------------------------------------------------
## 5.  Inspect the final model
## ----------------------------------------------------------------
summary(step.mod)
vif.cca(step.mod)      # check VIFs of retained predictors
anova(step.mod)        # global test (999 perms)
anova(step.mod, by = "term")   # marginal tests

















###########################
##### CORE MICROBIOME #####
###########################

# Creo los objetos phyloseq para los niveles taxonómicos

taxa_table <- row.names(abundance_table_species)
taxa_table <- as.matrix(taxa_table)
colnames(taxa_table) <- "Species"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(abundance_table_species, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$sample
ps_species = phyloseq(OTU_taxa, TAX_taxa, sampledata)
ps_species_21 <- phyloseq::subset_samples(ps_species, date %in% c("2021_oct"))
ps_species_22 <- phyloseq::subset_samples(ps_species, date %in% c("2022_jun"))
ps_species <- prune_taxa(taxa_sums(ps_species)>0, ps_species)
ps_species <- microbiome::transform(ps_species, "compositional")

taxa_table <- row.names(abundance_table_genus)
taxa_table <- as.matrix(taxa_table)
colnames(taxa_table) <- "Genus"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(abundance_table_genus, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$sample
ps_genus = phyloseq(OTU_taxa, TAX_taxa, sampledata)
ps_genus <- prune_taxa(taxa_sums(ps_genus)>0, ps_genus)
ps_genus <- microbiome::transform(ps_genus, "compositional")

taxa_table <- row.names(abundance_table_family)
taxa_table <- as.matrix(taxa_table)
colnames(taxa_table) <- "Family"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(abundance_table_family, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$sample
ps_family = phyloseq(OTU_taxa, TAX_taxa, sampledata)
ps_family <- prune_taxa(taxa_sums(ps_family)>0, ps_family)
ps_family <- microbiome::transform(ps_family, "compositional")

taxa_table <- row.names(abundance_table_class)
taxa_table <- as.matrix(taxa_table)
colnames(taxa_table) <- "Class"
rownames(taxa_table) <- taxa_table
OTU_taxa <- otu_table(abundance_table_class, taxa_are_rows = TRUE)
TAX_taxa <- tax_table(taxa_table)
sampledata <- sample_data(metadata)
rownames(sampledata) <- metadata$sample
ps_class = phyloseq(OTU_taxa, TAX_taxa, sampledata)
ps_class <- prune_taxa(taxa_sums(ps_class)>0, ps_class)
ps_class <- microbiome::transform(ps_class, "compositional")

# Determino el core microbiome usando "strictly greater by default": 
# con detection 0.5% y prevalencia del 50%

core_species <- core_members(ps_species, detection = 0.5/100, prevalence = 50/100, include.lowest = TRUE)
print(core_species)
ps.core_species <- core(ps_species, detection = 0.5/100, prevalence = 50/100, include.lowest = TRUE)
core_species_21 <- core_members(ps_species_21, detection = 0.5/100, prevalence = 50/100, include.lowest = TRUE)
print(core_species_21)
ps.core_species_21 <- core(ps_species_21, detection = 0.5/100, prevalence = 50/100, include.lowest = TRUE)
core_species_22 <- core_members(ps_species_22, detection = 0.5/100, prevalence = 50/100, include.lowest = TRUE)
print(core_species_22)
ps.core_species_22 <- core(ps_species_22, detection = 0.5/100, prevalence = 50/100, include.lowest = TRUE)

core_genus <- core_members(ps_genus, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)
print(core_genus)
ps.core_genus <- core(ps_genus, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)

core_family <- core_members(ps_family, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)
print(core_family)
ps.core_family <- core(ps_family, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)

core_class <- core_members(ps_class, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)
print(core_class)
ps.core_class <- core(ps_class, detection = 1/100, prevalence = 50/100, include.lowest = TRUE)


prevalences <- seq(5/100, 100/100, 5/100)
detections <- round(10^seq(log10(0.001), log10(.5), length = 10),3)

#p1 <- plot_core(ps.core_species,
#p1 <- plot_core(ps.core_species_21,
p1 <- plot_core(ps.core_species_22,
                #p1 <- plot_core(ps.core_genus,
                #p1 <- plot_core(ps.core_family,
                #p1 <- plot_core(ps.core_class,
                plot.type = "heatmap",
                colours = rev(brewer.pal(5,"RdGy")),
                prevalences = prevalences,
                detections = detections, 
                min.prevalence = 0.01) +
  xlab("Detection Threshold (Relative Abundance)") +
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))
p1 <- p1 + theme_bw() + ylab("Species")
#p1 <- p1 + theme_bw() + ylab("Genus")
#p1 <- p1 + theme_bw() + ylab("Family")
#p1 <- p1 + theme_bw() + ylab("Class")
p1

library(viridis)
print(p1 + scale_fill_viridis())

#svglite("Core_microbiome_of_Species_0.5.svg", width=10, height=5)
#svglite("Core_microbiome_of_Species_0.5_2021.svg", width=10, height=5)
svglite("Core_microbiome_of_Species_0.5_2022.svg", width=10, height=5)
#svglite("Core_microbiome_of_Genus.svg", width=10, height=5)
#svglite("Core_microbiome_of_Family.svg", width=10, height=5)
#svglite("Core_microbiome_of_Class.svg", width=10, height=5)
(p1 + scale_fill_viridis())
dev.off()


#Indices Chao1, Shannon y Simpson

p.shannon <- boxplot_alpha(ps.core, index = "shannon", x_var = "Temporada")
p.shannon <- p.shannon + theme_minimal() + labs(x="Temporada", y="Shannon diversity - Core microbiome") + 
  theme(
    axis.text = element_text(size = 12), 
    axis.title = element_text(size=16),
    legend.text=element_text(size=12),
    legend.title=element_text(size=16)
  )
p.shannon

p.simpson <- boxplot_alpha(ps.core, index = "simpson", x_var = "Temporada")
p.simpson <- p.simpson + theme_minimal() + labs(x="Temporada", y="Simpson diversity - Core microbiome") + 
  theme(
    axis.text = element_text(size = 12), 
    axis.title = element_text(size=16),
    legend.text=element_text(size=12),
    legend.title=element_text(size=16)
  )
p.simpson

# Guardar con patchwork library
svglite("Alpha_diversity_core_microbiome.svg", width=20, height=10)
(p.shannon | p.simpson)
dev.off()

#otra representación con ggpubr library
svglite("Alpha_diversity_core_microbiome.svg", width=20, height=10)
ggarrange(p.shannon,p.simpson, legend="right")
dev.off()

# Alfa diversidad con p-values

#Temporada
#"Q1" = "#0466c8", "Q2" = "#7209b7", "Q3" = "#f26419", "Q4" = "#38b000", "NR" = "#8b8c89"
mycols <- c("#0466c8", "#7209b7", "#f26419", "#38b000", "#8b8c89")
shannon_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_shannon",
                                            group = "Season",
                                            group.order = c("Lluvias", "Nortes"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_shannon <- shannon_pvalue_core + ylab("Shannon Diversity - Core microbiome - Season") + xlab("")

simpson_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_gini_simpson",
                                            group = "Season",
                                            group.order = c("Lluvias", "Nortes"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_simpson <- simpson_pvalue_core + ylab("Simpson Diversity - Core microbiome - Season") + xlab("")

svglite("Alpha_diversity_core_microbiome_p_value_Temporada.svg", width=20, height=10)
(core_shannon | core_simpson)
dev.off()

#Z_m
mycols <- c("brown3", "steelblue", "grey50", "green", "red")
shannon_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_shannon",
                                            group = "Z_m",
                                            group.order = c("Q1", "Q2","Q3", "Q4"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_shannon <- shannon_pvalue_core + ylab("Shannon Diversity - Core microbiome - Z_m") + xlab("")

simpson_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_gini_simpson",
                                            group = "Z_m",
                                            group.order = c("Q1", "Q2","Q3", "Q4"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_simpson <- simpson_pvalue_core + ylab("Simpson Diversity - Core microbiome - Z_m") + xlab("")

svglite("Alpha_diversity_core_microbiome_p_value_Z_m.svg", width=20, height=10)
(core_shannon | core_simpson)
dev.off()

#Total_aliph
mycols <- c("brown3", "steelblue", "grey50", "green", "red")
shannon_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_shannon",
                                            group = "Total_aliph",
                                            group.order = c("Q1", "Q2","Q3", "Q4","NR"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_shannon <- shannon_pvalue_core + ylab("Shannon Diversity - Core microbiome - Total_aliph") + xlab("")

simpson_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_gini_simpson",
                                            group = "Total_aliph",
                                            group.order = c("Q1", "Q2","Q3", "Q4","NR"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_simpson <- simpson_pvalue_core + ylab("Simpson Diversity - Core microbiome -Total_aliph") + xlab("")

svglite("Alpha_diversity_core_microbiome_p_value_Total_aliph.svg", width=20, height=10)
(core_shannon | core_simpson)
dev.off()

#Total_HAPs
mycols <- c("brown3", "steelblue", "grey50", "green", "red")
shannon_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_shannon",
                                            group = "Total_HAPs",
                                            group.order = c("Q1", "Q2","Q3", "Q4","NR"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_shannon <- shannon_pvalue_core + ylab("Shannon Diversity - Core microbiome - Total_HAPs") + xlab("")

simpson_pvalue_core <- plot_diversity_stats(ps.core,
                                            index = "diversity_gini_simpson",
                                            group = "Total_HAPs",
                                            group.order = c("Q1", "Q2","Q3", "Q4","NR"),
                                            group.colors = mycols,
                                            dot.opacity = 0.25,
                                            violin.opacity = 0.25,
                                            stats = TRUE,
                                            label.format = "p.format")
core_simpson <- simpson_pvalue_core + ylab("Simpson Diversity - Core microbiome -Total_HAPs") + xlab("")

svglite("Alpha_diversity_core_microbiome_p_value_Total_HAPs.svg", width=20, height=10)
(core_shannon | core_simpson)
dev.off()


######### END OF CORE MICROBIOME ANALYSIS ###################3






scale_fill_manual(
  values = c("#d7263d","#f46036","#2e294e","#1b998b","#c5d86d","#759aab","#faf2a1","#4d8b31","#ffc800",
             "#1c2541","#3a506b","#5bc0be","#6fffe9","#ff8811","#f4d06f","#fff8f0","#9dd9d2","#eeebd3",
             "#F04C05","#B74211","#7E391D","#442F29","#0B2535","#15640F","#4A5414","#7F4419","#B4331E","#E92323",
             "#EAAD12","#B3A647","#7BA07C","#4499B0","#0C92E5","#5BF30A","#60B835","#667C61","#6B418C","#7005B7",
             "#F44E4E","#B96669","#7D7D83","#42959E","#06ACB8","#EF745C","#D06257","#B15052","#923E4D","#722B47","#531942","#34073D",
             "#074170","#1B4E60","#2F5B51","#436941","#567631","#6A8322","#7E9012","#8ecae6","#219ebc","#023047","#ffb703","#fb8500",
             "#f6d7d7","#df7472","#53032b","#2d0627","#06262c","#72dddf","#0fd780","#2b5403","#a6df72","#c9d70f",
             "#a09ebb","#a8aec1","#b5d2cb","#bfffbc","#a6ffa1","#476a6f","#031a6b","#033860","#004385","#305252",
             "#160f29","#246a73","#368f8b","#f3dfc1","#ddbea8","#ea526f","#e76b74","#d7af70","#c9c19f","#edf7d2",
             "#f94144","#f3722c","#f8961e","#f9844a","#f9c74f","#90be6d","#43aa8b","#4d908e","#577590","#277da1",
             "#f72585","#b5179e","#7209b7","#560bad","#480ca8","#3a0ca3","#3f37c9","#4361ee","#4895ef","#4cc9f0",
             "#d4e09b","#f6f4d2","#cbdfbd","#f19c79","#a44a3f"),
  