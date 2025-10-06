## San Quintin coastal lagoon - Sediment bioinformatic analysis
## Gene abundance and pathway analysis
## Prepared by: Jorge Rojas-Vargas

# Load libraries
library(ggplot2)    # plots
library(dplyr)      # tables
library(tidyr)      # tables
library(tibble)     # tables
library(FSA)        # to do dunnTest()
library(rstatix)    # to adjust p-values
library(patchwork)  # combine plots into a single figure
library(svglite)    # save svg plots

# Load data
load("pathways_abund_drug.RData")
load("pathways_abund_metal.RData")
load("pathways_abund_nitrogen.RData")
load("pathways_abund_sulfur.RData")
load("pathways_abund_virulence.RData")
load("genes_abund_metal.RData")
load("genes_abund_drug.RData")
load("genes_abund_nitrogen.RData")
load("genes_abund_sulfur.RData")
load("genes_abund_virulence.RData")
load("metadata_BSQ.RData")



## -----------------------------
## 
## BOX PLOTS ABUNDANCE PATHWAYS
##
## -----------------------------

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
columnas_de_interes <- c(2:24)
# Sulfur
columnas_de_interes <- c(2, 3, 4, 5, 6, 7, 8)
# Drug
columnas_de_interes <- c(2:34)
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

svglite("annot_metal_metabolism3.svg", width=29, height=8)
((a15 | a16 | a17 | a18 | a19 | a20 | a21) / (a22 | a23 | a24| a25 | a12 | a13 | a14))
dev.off()

svglite("annot_sulfur_metabolism.svg", width=29, height=4)
((a1 | a2 | a3 | a4 | a5 | a6 | a7))
dev.off()

svglite("annot_drug_metabolism2.svg", width=29, height=16)
((a29 | a30 | a31 | a32 | a33 | a6 | a7) / (a8 | a9 | a10| a11 | a12 | a13 | a14) / (a15 | a16 | a17| a18 | a19 | a20 | a21) / (a22 | a23 | a24 | a25 | a26 | a27 | a28))
dev.off()

svglite("annot_virulence_metabolism.svg", width=24.7, height=12)
((a1 | a2 | a3 | a4 | a5 | a6) / (a7 | a8 | a9 | a10 | a11 | a12) / (a13 | a14 | a15 | a16 | a1 | a12))
dev.off()




### ---------------------------
## Diferencias entre pathways
## ----------------------------


# Unir las tablas de abundancia y metadatos
combined_data <- merge(drug_pathways_abund, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abund, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_pathways_abund, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_pathways_abund, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abund, metadata_BSQ, by.x = "Sample", by.y = "Sample")

# Ajusta el modelo con los tres factores y todas sus interacciones
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Columnas de abundancia
pathways_cols <- combined_data %>%
  dplyr::select(-Sample, -Sector, -Season, -Habitat, -ID_1, -ID_2) %>%
  names()

# Calcular KW p-values para cada pathway
kw_results <- lapply(pathways_cols, function(var){
  kw <- kruskal.test(combined_data[[var]] ~ combined_data$Sector)
  data.frame(
    Pathway = var,
    P.value = kw$p.value
  )
}) %>%
  bind_rows() %>%
  arrange(P.value)

# Kruskal–Wallis + Dunn para Sector
kw_dunn_results <- lapply(pathways_cols, function(var){
  # Kruskal–Wallis
  kw <- kruskal.test(combined_data[[var]] ~ combined_data$Sector)
  
  # Si KW sale p ≤ 0.05, lanza Dunn
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

# Wilcoxon para Season y Habitat
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

# Mostrar todos los resultados
kw_results
print(wilcox_results, n = 100)

# Filtrar resultados “tendencia” 0.05 < p ≤ 0.1
kw_dunn_sig   <- kw_dunn_results   %>% filter(P.adj <= 0.1)
wilcox_sig    <- wilcox_results    %>% filter(P.adj <= 0.1)

# Listado final
list(
  Dunn_Sector     = kw_dunn_sig,
  Wilcox_SeasonHabitat = wilcox_sig
)


## ----------------------------
## 
## HEATMAPS GENES
##
## ----------------------------


# Function to convert matrix into long format, join with metadata, and calculate means
process_genes <- function(genes_abundance, metadata) {
  # Convert the matrix into a data.frame
  genes_abundance_df <- as.data.frame(t(genes_abundance))
  
  # Convert to long format
  genes_long <- genes_abundance_df %>%
    rownames_to_column(var = "Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance")
  
  # Join with metadata
  genes_metadata <- genes_long %>%
    left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calculate means by Sector
  sector_means <- genes_metadata %>%
    group_by(Gene, Sector) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Sector, values_from = mean_abundance)
  
  # Calculate means by Season
  season_means <- genes_metadata %>%
    group_by(Gene, Season) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Season, values_from = mean_abundance)
  
  # Calculate means by Habitat
  habitat_means <- genes_metadata %>%
    group_by(Gene, Habitat) %>%
    summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    pivot_wider(names_from = Habitat, values_from = mean_abundance)
  
  # Join all the means into a single dataframe
  genes_means <- genes_abundance_df %>%
    rownames_to_column(var = "Gene") %>%
    left_join(sector_means, by = "Gene") %>%
    left_join(season_means, by = "Gene") %>%
    left_join(habitat_means, by = "Gene")
  
  return(genes_means)
}

# Apply the function to each matrix
metal_means <- process_genes(metal_genes_abundance, metadata_BSQ)
drug_means <- process_genes(drug_genes_abundance, metadata_BSQ)
nitrogen_means <- process_genes(nitrogen_genes_abundance, metadata_BSQ)
sulfur_means <- process_genes(sulfur_genes_abundance, metadata_BSQ)
virulence_means <- process_genes(virulence_genes_abundance, metadata_BSQ)

# Calculate maximum values
max_value_nitrogen <- max(nitrogen_means[ , (ncol(nitrogen_means) - 6):ncol(nitrogen_means)], na.rm = TRUE)
max_value_sulfur <- max(sulfur_means[ , (ncol(sulfur_means) - 6):ncol(sulfur_means)], na.rm = TRUE)

# Function to generate heatmaps from a *_means data frame
create_heatmaps <- function(means_data, max_val) {
  
  # Convert to long format for Sector, Season, and Habitat
  sector_data <- means_data %>%
    select(Gene, starts_with("A_"), starts_with("B_"), starts_with("C_In")) %>%
    pivot_longer(-Gene, names_to = "Sector", values_to = "Abundance")
  
  season_data <- means_data %>%
    select(Gene, starts_with("Intense"), starts_with("Relax")) %>%
    pivot_longer(-Gene, names_to = "Season", values_to = "Abundance")
  
  habitat_data <- means_data %>%
    select(Gene, starts_with("Bare"), starts_with("Seagrass")) %>%
    pivot_longer(-Gene, names_to = "Habitat", values_to = "Abundance")
  
  # Create heatmap for Sector
  heatmap_sector <- ggplot(sector_data, aes(x = Sector, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Create heatmap for Season
  heatmap_season <- ggplot(season_data, aes(x = Season, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Create heatmap for Habitat
  heatmap_habitat <- ggplot(habitat_data, aes(x = Habitat, y = Gene, fill = Abundance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"), 
                         na.value = "grey50", limit = c(0, max_val),
                         name = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = "", fill = "Abundance") +
    coord_fixed(ratio = 0.2) 
  
  # Combine the three heatmaps into one figure
  combined_plot <- (heatmap_sector | heatmap_season | heatmap_habitat) + 
    plot_layout(guides = "collect") & theme(legend.position = "right")
  
  return(combined_plot)
}

# Generate heatmaps for each *_means data frame
heatmap_nitrogen <- create_heatmaps(nitrogen_means, max_value_nitrogen)
heatmap_sulfur <- create_heatmaps(sulfur_means, max_value_sulfur)

# Show each plot
#print(heatmap_drug)
print(heatmap_nitrogen)
print(heatmap_sulfur)

# Define the number of genes in each dataset
num_genes_nitrogen <- nrow(nitrogen_means)
num_genes_sulfur <- nrow(sulfur_means)

# Adjust height dynamically
height_factor <- 0.3  # Adjust this value according to the spacing you need

# Save figures with adjusted dimensions
#svglite("heatmap_NITROGEN_genes.svg", width = 8, height = max(8, num_genes_nitrogen * height_factor))
heatmap_nitrogen
#dev.off()

#svglite("heatmap_SULFUR_genes.svg", width = 8, height = max(8, num_genes_sulfur * height_factor))
heatmap_sulfur
#dev.off()



## ----------------------------
## Differences between genes
## ----------------------------

# Transpose the abundance table
drug_genes_abundance$Sample <- rownames(drug_genes_abundance)
metal_genes_abundance$Sample <- rownames(metal_genes_abundance)
nitrogen_genes_abundance$Sample <- rownames(nitrogen_genes_abundance)
sulfur_genes_abundance$Sample <- rownames(sulfur_genes_abundance)
virulence_genes_abundance$Sample <- rownames(virulence_genes_abundance)

# Merge abundance tables and metadata
combined_data <- merge(drug_genes_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_genes_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_genes_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_genes_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_genes_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")

# Fit the model with the three factors and all their interactions
combined_data$Sector  <- as.factor(combined_data$Sector)
combined_data$Season  <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Select only numeric columns (excluding "Sample", "Sector", "Season", "Habitat")
genes_cols <- combined_data %>%
  select(-Sample, -Sector, -Season, -Habitat, -ID_1, -ID_2)

# Function to apply Dunn test (Sector: 3 levels)
dunn_results <- lapply(names(genes_cols), function(gene) {
  test <- dunnTest(combined_data[[gene]] ~ combined_data$Sector, method = "bh")
  result <- data.frame(Gene = gene, test$res)
  return(result)
}) %>%
  bind_rows() %>%
  filter(P.adj <= 0.1)  # Filter adjusted p-values < 0.05

# Function to apply Wilcoxon test (Season and Habitat: 2 levels)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(genes_cols), function(gene) {
    test <- wilcox.test(combined_data[[gene]] ~ combined_data[[variable]])
    data.frame(Gene = gene, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  filter(p.value <= 0.1)  # Filter p-values <= 0.05

# Function to apply Wilcoxon test (Season and Habitat: 2 levels)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(pathways_cols), function(pathway) {
    test <- wilcox.test(combined_data[[pathway]] ~ combined_data[[variable]])
    data.frame(Pathway = pathway, Variable = variable, p.value = test$p.value)
  }) %>%
    bind_rows()
}) 

# Assuming wilcox_results is a list of two data.frames
for (i in seq_along(wilcox_results)) {
  wilcox_results[[i]]$p.adj <- p.adjust(
    wilcox_results[[i]]$p.value,
    method = "BH"
  )
}

wilcox_results <- wilcox_results %>%
  bind_rows() %>%
  filter(p.adj <= 0.1)  # Filter p-adj <= 0.1

# Show significant results
list(Significant_Dunn_Test = dunn_results, Significant_Wilcoxon = wilcox_results)




## ----------------------------
## 
## HEATMAPS PATHWAYS
##
## ----------------------------

# Function to convert matrix to long format, Join with metadata, and calculate means
process_pathways <- function(pathways_abundance, metadata) {
  # Convert the matrix to a data.frame
  pathways_abundance_df <- as.data.frame(t(pathways_abundance))
  
  # Convert to long format
  pathways_long <- pathways_abundance_df %>%
    tibble::rownames_to_column(var = "Pathway") %>%
    tidyr::pivot_longer(-Pathway, names_to = "Sample", values_to = "Abundance")
  
  # Join with metadata
  pathways_metadata <- pathways_long %>%
    dplyr::left_join(metadata, by = c("Sample" = "Sample"))
  
  # Calculate means by Sector
  sector_means <- pathways_metadata %>%
    dplyr::group_by(Pathway, Sector) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    tidyr::pivot_wider(names_from = Sector, values_from = mean_abundance)
  
  # Calculate means by Season
  season_means <- pathways_metadata %>%
    dplyr::group_by(Pathway, Season) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    tidyr::pivot_wider(names_from = Season, values_from = mean_abundance)
  
  # Calculate means by Habitat
  habitat_means <- pathways_metadata %>%
    dplyr::group_by(Pathway, Habitat) %>%
    dplyr::summarize(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    tidyr::pivot_wider(names_from = Habitat, values_from = mean_abundance)
  
  # Join all means into a single dataframe
  pathways_means <- pathways_abundance_df %>%
    tibble::rownames_to_column(var = "Pathway") %>%
    dplyr::left_join(sector_means,  by = "Pathway") %>%
    dplyr::left_join(season_means,  by = "Pathway") %>%
    dplyr::left_join(habitat_means, by = "Pathway")
  
  return(pathways_means)
}

# Apply the function to each matrix
metal_means      <- process_pathways(metal_pathways_abundance, metadata_BSQ)
drug_means       <- process_pathways(drug_pathways_abundance, metadata_BSQ)
nitrogen_means   <- process_pathways(nitrogen_pathways_abundance, metadata_BSQ)
sulfur_means     <- process_pathways(sulfur_pathways_abundance, metadata_BSQ)
virulence_means  <- process_pathways(virulence_pathways_abundance, metadata_BSQ)

# Calculate maximum values
max_value_metal     <- max(metal_means[    , (ncol(metal_means)     - 6):ncol(metal_means)],     na.rm = TRUE)
max_value_drug      <- max(drug_means[     , (ncol(drug_means)      - 6):ncol(drug_means)],      na.rm = TRUE)
max_value_nitrogen  <- max(nitrogen_means[ , (ncol(nitrogen_means)  - 6):ncol(nitrogen_means)],  na.rm = TRUE)
max_value_sulfur    <- max(sulfur_means[   , (ncol(sulfur_means)    - 6):ncol(sulfur_means)],    na.rm = TRUE)
max_value_virulence <- max(virulence_means[, (ncol(virulence_means) - 6):ncol(virulence_means)], na.rm = TRUE)

# Function to create heatmaps for a *_means data frame (Pathways)
create_heatmaps <- function(means_data, max_val) {
  
  # Convert to long format for Sector, Season, and Habitat
  sector_data <- means_data %>%
    dplyr::select(Pathway, dplyr::starts_with("A_"), dplyr::starts_with("B_"), dplyr::starts_with("C_In")) %>%
    tidyr::pivot_longer(-Pathway, names_to = "Sector", values_to = "Abundance")
  
  season_data <- means_data %>%
    dplyr::select(Pathway, dplyr::starts_with("Intense"), dplyr::starts_with("Relax")) %>%
    tidyr::pivot_longer(-Pathway, names_to = "Season", values_to = "Abundance")
  
  habitat_data <- means_data %>%
    dplyr::select(Pathway, dplyr::starts_with("Bare"), dplyr::starts_with("Seagrass")) %>%
    tidyr::pivot_longer(-Pathway, names_to = "Habitat", values_to = "Abundance")
  
  # Create heatmap for Sector
  heatmap_sector <- ggplot2::ggplot(sector_data, ggplot2::aes(x = Sector, y = Pathway, fill = Abundance)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                                  na.value = "grey50", limit = c(0, max_val),
                                  name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "", fill = "Abundance") +
    ggplot2::coord_fixed(ratio = 0.2) 
  
  # Create heatmap for Season
  heatmap_season <- ggplot2::ggplot(season_data, ggplot2::aes(x = Season, y = Pathway, fill = Abundance)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                                  na.value = "grey50", limit = c(0, max_val),
                                  name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(title = "", fill = "Abundance") +
    ggplot2::coord_fixed(ratio = 0.2) 
  
  # Create heatmap for Habitat
  heatmap_habitat <- ggplot2::ggplot(habitat_data, ggplot2::aes(x = Habitat, y = Pathway, fill = Abundance)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradientn(colors = c("#ffd166", "#619b8a", "#3e5c76", "#6a4c93", "#a4161a"),  
                                  na.value = "grey50", limit = c(0, max_val),
                                  name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(title = "", fill = "Abundance") +
    ggplot2::coord_fixed(ratio = 0.2) 
  
  # Combine the three heatmaps into a single figure
  combined_plot <- (heatmap_sector | heatmap_season | heatmap_habitat) + 
    patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "right")
  
  return(combined_plot)
}

# Generate heatmaps for each *_means data frame
heatmap_drug       <- create_heatmaps(drug_means,       max_value_drug)
heatmap_nitrogen   <- create_heatmaps(nitrogen_means,   max_value_nitrogen)
heatmap_sulfur     <- create_heatmaps(sulfur_means,     max_value_sulfur)
heatmap_metal      <- create_heatmaps(metal_means,      max_value_metal)
heatmap_virulence  <- create_heatmaps(virulence_means,  max_value_virulence)

# Show each plot
print(heatmap_drug)
print(heatmap_nitrogen)
print(heatmap_sulfur)
print(heatmap_metal)
print(heatmap_virulence)

# Define the number of pathways in each dataset
num_pathways_drug      <- nrow(drug_means)
num_pathways_nitrogen  <- nrow(nitrogen_means)
num_pathways_sulfur    <- nrow(sulfur_means)
num_pathways_metal     <- nrow(metal_means)
num_pathways_virulence <- nrow(virulence_means)

# Adjust height dynamically
height_factor <- 0.3  # Adjust this value according to the spacing you need

# Save figures with adjusted dimensions
#svglite("heatmap_CARD_pathway3.svg", width = 8, height = max(10, num_pathways_drug * height_factor))
heatmap_drug
#dev.off()

#svglite("heatmap_NITROGEN_pathway.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_nitrogen
#dev.off()

#svglite("heatmap_SULFUR_pathway3.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_sulfur
#dev.off()

#svglite("heatmap_METAL_pathway3.svg", width = 8, height = max(8, num_pathways_metal * height_factor))
heatmap_metal
#dev.off()

#svglite("heatmap_VFDB_pathway3.svg", width = 8, height = max(8, num_pathways_virulence * height_factor))
heatmap_virulence
#dev.off()


## ----------------------------
## Differences between pathways
## ----------------------------

# Transpose the abundance table
drug_pathways_abundance$Sample      <- rownames(drug_pathways_abundance)
metal_pathways_abundance$Sample     <- rownames(metal_pathways_abundance)
nitrogen_pathways_abundance$Sample  <- rownames(nitrogen_pathways_abundance)
sulfur_pathways_abundance$Sample    <- rownames(sulfur_pathways_abundance)
virulence_pathways_abundance$Sample <- rownames(virulence_pathways_abundance)

# Merge abundance tables and metadata
combined_data <- merge(drug_pathways_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(metal_pathways_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(nitrogen_pathways_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_pathways_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(virulence_pathways_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")

# Fit the model with the three factors and all their interactions
combined_data$Sector  <- as.factor(combined_data$Sector)
combined_data$Season  <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Select only numeric columns (excluding "Sample", "Sector", "Season", "Habitat")
pathways_cols <- combined_data %>%
  dplyr::select(-Sample, -Sector, -Season, -Habitat, -ID_1, -ID_2)

# Function to apply Dunn test (Sector: 3 levels)
dunn_results <- lapply(names(pathways_cols), function(pathway) {
  test <- FSA::dunnTest(combined_data[[pathway]] ~ combined_data$Sector, method = "bh")
  result <- data.frame(Pathway = pathway, test$res)
  return(result)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(P.adj <= 0.1)

library(dplyr)

# Function to apply Wilcoxon test (Season and Habitat: 2 levels)
wilcox_results <- lapply(c("Season", "Habitat"), function(variable) {
  lapply(names(pathways_cols), function(pathway) {
    test <- stats::wilcox.test(combined_data[[pathway]] ~ combined_data[[variable]])
    data.frame(Pathway = pathway, Variable = variable, p.value = test$p.value)
  }) %>%
    dplyr::bind_rows()
})

# Assuming wilcox_results is a list of two data.frames
for (i in seq_along(wilcox_results)) {
  wilcox_results[[i]]$p.adj <- p.adjust(
    wilcox_results[[i]]$p.value,
    method = "BH"
  )
}

wilcox_results <- wilcox_results %>%
  dplyr::bind_rows() %>%
  dplyr::filter(p.adj <= 0.1)

# Show significant results
list(Significant_Dunn_Test = dunn_results, Significant_Wilcoxon = wilcox_results)




## ----------------------------
## 
## COMPLETENESS PATHWAYS
##
## ----------------------------

## ----------------------------
## Differences between pathways
## ----------------------------

# Transpose the abundance table
nitrogen_genes_abundance$Sample <- rownames(nitrogen_genes_abundance)
sulfur_genes_abundance$Sample   <- rownames(sulfur_genes_abundance)

# Merge abundance tables and metadata
combined_data <- merge(nitrogen_genes_abundance, metadata_BSQ, by.x = "Sample", by.y = "Sample")
combined_data <- merge(sulfur_genes_abundance,   metadata_BSQ, by.x = "Sample", by.y = "Sample")

# Fit the model with the three factors and all their interactions
combined_data$Sector  <- as.factor(combined_data$Sector)
combined_data$Season  <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

## -- For NITROGEN --

# Define the gene list per pathway, including alternative routes
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

# Load data and join abundances with metadata
load("genes_abund_nitrogen.RData")  # should leave 'nitrogen_genes_abundance' in the environment
metadata <- read.csv("metadata_quintin.csv")

# Transpose and merge
nitrogen_genes_abundance$Sample <- rownames(nitrogen_genes_abundance)
combined_data <- merge(nitrogen_genes_abundance, metadata, by = "Sample")

# Create presence/absence (0/1) data.frame for **all** genes
#    If a gene does not exist in combined_data, it is considered 0
desired_genes <- unlist(
  lapply(genes_list, function(x) if (is.list(x)) unlist(x) else x)
) %>% unique()

present_genes <- intersect(desired_genes, names(combined_data))
missing_genes <- setdiff(desired_genes, names(combined_data))

# Build presence_df with one Sample column + all genes
presence_df <- combined_data %>%
  dplyr::select(Sample) %>%
  dplyr::bind_cols(
    combined_data %>%
      dplyr::select(dplyr::all_of(present_genes)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(. > 0)))
  )

# Add missing genes as zero columns
for (g in missing_genes) {
  presence_df[[g]] <- 0L
}

# Reorder columns: Sample + desired_genes
presence_df <- presence_df %>%
  dplyr::select(Sample, dplyr::all_of(desired_genes))

# 4) Helper function for OR steps in alternative pipelines
step_present <- function(df_row, groups){
  # df_row: numeric 0/1 vector
  # groups: list of name vectors
  any(sapply(groups, function(g) any(df_row[g] == 1))) * 1
}

# Compute completeness of each pathway (%) per sample
completeness_df <- presence_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Anammox = sum(dplyr::c_across(dplyr::all_of(genes_list$Anammox))) / 5 * 100,
    
    Nitrogen_fixation = dplyr::if_else(
      dplyr::c_across("anfG") == 1,
      100,
      sum(dplyr::c_across(setdiff(genes_list$Nitrogen_fixation, "anfG"))) / 4 * 100
    ),
    
    Nitrification = {
      cntA <- sum(dplyr::c_across(genes_list$Nitrification$groupA))
      cntB <- sum(dplyr::c_across(genes_list$Nitrification$groupB))
      cntC <- sum(dplyr::c_across(genes_list$Nitrification$common))
      (max(cntA, cntB) + cntC) / 6 * 100
    },
    
    Denitrification = {
      s1 <- step_present(dplyr::cur_data(), genes_list$Denitrification$step1)
      s2 <- any(dplyr::c_across(genes_list$Denitrification$step2)) * 1
      s3 <- any(dplyr::c_across(genes_list$Denitrification$step3)) * 1
      s4 <- dplyr::c_across(genes_list$Denitrification$step4)
      (s1 + s2 + s3 + s4) / 4 * 100
    },
    
    Assimilatory_nitrate_reduction = {
      s1 <- step_present(dplyr::cur_data(), genes_list$Assimilatory_nitrate_reduction$step1)
      s2 <- dplyr::c_across(genes_list$Assimilatory_nitrate_reduction$step2)
      (s1 + s2) / 2 * 100
    },
    
    Dissimilatory_nitrate_reduction = {
      s1 <- step_present(dplyr::cur_data(), genes_list$Dissimilatory_nitrate_reduction$step1)
      s2 <- step_present(dplyr::cur_data(), genes_list$Dissimilatory_nitrate_reduction$step2)
      (s1 + s2) / 2 * 100
    },
    
    Organic_degradation_and_synthesis = {
      c1 <- sum(dplyr::c_across(genes_list$Organic_degradation_and_synthesis$route1)) /
        length(genes_list$Organic_degradation_and_synthesis$route1) * 100
      c2 <- sum(dplyr::c_across(genes_list$Organic_degradation_and_synthesis$route2)) /
        length(genes_list$Organic_degradation_and_synthesis$route2) * 100
      max(c1, c2)
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Sample,
    Anammox,
    Nitrogen_fixation,
    Nitrification,
    Denitrification,
    Assimilatory_nitrate_reduction,
    Dissimilatory_nitrate_reduction,
    Organic_degradation_and_synthesis
  )

# See results
head(completeness_df)
# (Optional) join it to combined_data:
combined_data <- combined_data %>% dplyr::left_join(completeness_df, by = "Sample")

## Compute means for completeness

# Function to convert your completeness data.frame into means by factor
process_completeness <- function(completeness_df, metadata) {
  # pivot to “long”
  long <- completeness_df %>%
    tidyr::pivot_longer(
      cols      = -Sample,
      names_to  = "Pathway",
      values_to = "Completeness"
    ) %>%
    dplyr::left_join(metadata, by = "Sample")
  
  # means by Sector
  sector_means <- long %>%
    dplyr::group_by(Pathway, Sector) %>%
    dplyr::summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sector, values_from = mean_completeness)
  
  # means by Season
  season_means <- long %>%
    dplyr::group_by(Pathway, Season) %>%
    dplyr::summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Season, values_from = mean_completeness)
  
  # means by Habitat
  habitat_means <- long %>%
    dplyr::group_by(Pathway, Habitat) %>%
    dplyr::summarise(mean_completeness = mean(Completeness, na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Habitat, values_from = mean_completeness)
  
  # join everything
  full_means <- sector_means %>%
    dplyr::left_join(season_means,  by = "Pathway") %>%
    dplyr::left_join(habitat_means, by = "Pathway")
  
  return(full_means)
}

# Compute means of nitrogen completeness
nitrogen_means <- process_completeness(completeness_df, metadata)

# Function to plot heatmaps of completeness means
create_heatmaps <- function(means_df) {
  # long format for each factor
  sec <- means_df %>%
    tidyr::pivot_longer(-Pathway, names_to = "Sector",  values_to = "Completeness") %>%
    dplyr::filter(Sector %in% c("A_Inlet","B_Transition","C_Inner"))
  
  seas <- means_df %>%
    tidyr::pivot_longer(-Pathway, names_to = "Season",  values_to = "Completeness") %>%
    dplyr::filter(Season %in% c("Intense","Relax"))
  
  hab <- means_df %>%
    tidyr::pivot_longer(-Pathway, names_to = "Habitat", values_to = "Completeness") %>%
    dplyr::filter(Habitat %in% c("Bare","Seagrass"))
  
  # helper
  mk <- function(df, xvar){
    ggplot2::ggplot(df, ggplot2::aes_string(x = xvar, y = "Pathway", fill = "Completeness")) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradientn(
        colors  = c("#ffd166","#619b8a","#3e5c76","#6a4c93","#a4161a"),
        limits  = c(0,100),
        na.value = "grey80"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y = if (xvar == "Sector") ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.y = if (xvar == "Sector") ggplot2::element_line() else ggplot2::element_blank()
      ) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::coord_fixed(ratio = 0.2)
  }
  
  p1 <- mk(sec,  "Sector")
  p2 <- mk(seas, "Season")
  p3 <- mk(hab,  "Habitat")
  
  (p1 | p2 | p3) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
}

# Generate the heatmap of nitrogen completeness
heatmap_nitrogen <- create_heatmaps(nitrogen_means)
print(heatmap_nitrogen)

# Adjust height dynamically
height_factor <- 0.5  # Adjust this value according to the spacing you need
num_pathways_nitrogen <- 7

# Save figures with adjusted dimensions
#svglite("heatmap_NITROGEN_completeness.svg", width = 8, height = max(8, num_pathways_nitrogen * height_factor))
heatmap_nitrogen
#dev.off()



## -- For SULFUR --

# Define the gene list per sulfur pathway
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

# Prepare presence_df for these genes
desired_genes <- unique(unlist(genes_list_sulfur))
present_genes <- intersect(desired_genes, names(sulfur_genes_abundance))
missing_genes <- setdiff(desired_genes, names(sulfur_genes_abundance))

presence_sulfur <- combined_data %>%
  dplyr::select(Sample) %>%
  dplyr::bind_cols(
    combined_data %>%
      dplyr::select(dplyr::all_of(present_genes)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(. > 0)))
  )

for (g in missing_genes) {
  presence_sulfur[[g]] <- 0L
}

presence_sulfur <- presence_sulfur %>%
  dplyr::select(Sample, dplyr::all_of(desired_genes))

# Helper function for OR steps
step_present <- function(df_row, groups) {
  any(sapply(groups, function(g) any(df_row[g] == 1))) * 1
}

# Compute completeness for each sulfur pathway
completeness_sulfur <- presence_sulfur %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Assimilatory_sulfate_reduction = {
      s1 <- step_present(dplyr::cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step1)
      s2 <- all(dplyr::c_across(genes_list_sulfur$Assimilatory_sulfate_reduction$step2) == 1) * 1
      s3 <- step_present(dplyr::cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step3)
      s4 <- step_present(dplyr::cur_data(), genes_list_sulfur$Assimilatory_sulfate_reduction$step4)
      (s1 + s2 + s3 + s4) / 4 * 100
    },
    Dissimilatory_sulfur_reduction_and_oxidation = {
      s1 <- any(dplyr::c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step1)) * 1
      s2 <- step_present(dplyr::cur_data(), genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step2)
      s3 <- all(dplyr::c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step3) == 1) * 1
      s4 <- all(dplyr::c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step4) == 1) * 1
      s5 <- all(dplyr::c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step5) == 1) * 1
      s6 <- all(dplyr::c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step6) == 1) * 1
      s7 <- all(dplyr::c_across(genes_list_sulfur$Dissimilatory_sulfur_reduction_and_oxidation$step7) == 1) * 1
      (s1 + s2 + s3 + s4 + s5 + s6 + s7) / 7 * 100
    },
    Sulfur_reduction = {
      sum(sapply(genes_list_sulfur$Sulfur_reduction, function(gr) step_present(dplyr::cur_data(), gr))) / 5 * 100
    },
    SOX_systems = {
      sum(sapply(genes_list_sulfur$SOX_systems, function(gr) step_present(dplyr::cur_data(), gr))) / 3 * 100
    },
    Sulfur_oxidation = {
      sum(sapply(genes_list_sulfur$Sulfur_oxidation, function(gr) step_present(dplyr::cur_data(), gr))) / 4 * 100
    },
    Sulfur_disproportionation = {
      # Step 1: phsA, phsB, phsC all present
      s1 <- all(dplyr::c_across(genes_list_sulfur$Sulfur_disproportionation$step1) == 1) * 1
      # Step 2: tetH present
      s2 <- dplyr::c_across(genes_list_sulfur$Sulfur_disproportionation$step2)
      # Step 3: sor present
      s3 <- dplyr::c_across(genes_list_sulfur$Sulfur_disproportionation$step3)
      # Completeness
      (s1 + s2 + s3) / 3 * 100
    },
    Organic_sulfur_transformation = {
      # Compute each route and take the maximum
      r1 <- sum(dplyr::c_across(genes_list_sulfur$Organic_sulfur_transformation$route1)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route1) * 100
      r2 <- sum(dplyr::c_across(genes_list_sulfur$Organic_sulfur_transformation$route2)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route2) * 100
      r3 <- sum(dplyr::c_across(genes_list_sulfur$Organic_sulfur_transformation$route3)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route3) * 100
      r4 <- step_present(dplyr::cur_data(), genes_list_sulfur$Organic_sulfur_transformation$route4) / 5 * 100
      r5 <- step_present(dplyr::cur_data(), genes_list_sulfur$Organic_sulfur_transformation$route5) / 5 * 100
      r6 <- sum(dplyr::c_across(genes_list_sulfur$Organic_sulfur_transformation$route6)) /
        length(genes_list_sulfur$Organic_sulfur_transformation$route6) * 100
      max(r1, r2, r3, r4, r5, r6)
    },
    Link_between_inorganic_and_organic_sulfur_transformation = {
      sum(dplyr::c_across(genes_list_sulfur$Link_between_inorganic_and_organic_sulfur_transformation)) /
        length(genes_list_sulfur$Link_between_inorganic_and_organic_sulfur_transformation) * 100
    }
  ) %>%
  dplyr::ungroup()

# See results
pathway_cols <- c("Assimilatory_sulfate_reduction", "Dissimilatory_sulfur_reduction_and_oxidation",
                  "Sulfur_reduction", "SOX_systems", "Sulfur_oxidation", "Sulfur_disproportionation",
                  "Organic_sulfur_transformation", "Organic_sulfur_transformation",
                  "Link_between_inorganic_and_organic_sulfur_transformation")
nonzero_pathways <- completeness_sulfur %>%
  dplyr::select(dplyr::all_of(pathway_cols)) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
  unlist() %>%
  `>`(0) %>%                # logical vector
  which() %>%               # indices of TRUE
  names()          

completeness_df <- completeness_sulfur %>%
  dplyr::select(Sample, dplyr::all_of(nonzero_pathways))

# (Optional) Join to combined_data:
combined_data <- combined_data %>% dplyr::left_join(completeness_sulfur, by = "Sample")

# Compute means of sulfur completeness
sulfur_means <- process_completeness(completeness_df, metadata)

# Function to plot heatmaps of completeness means
create_heatmaps <- function(means_df) {
  # long format for each factor
  sec <- means_df %>%
    tidyr::pivot_longer(-Pathway, names_to = "Sector",  values_to = "Completeness") %>%
    dplyr::filter(Sector %in% c("A_Inlet","B_Transition","C_Inner"))
  
  seas <- means_df %>%
    tidyr::pivot_longer(-Pathway, names_to = "Season",  values_to = "Completeness") %>%
    dplyr::filter(Season %in% c("Intense","Relax"))
  
  hab <- means_df %>%
    tidyr::pivot_longer(-Pathway, names_to = "Habitat", values_to = "Completeness") %>%
    dplyr::filter(Habitat %in% c("Bare","Seagrass"))
  
  # helper
  mk <- function(df, xvar){
    ggplot2::ggplot(df, ggplot2::aes_string(x = xvar, y = "Pathway", fill = "Completeness")) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradientn(
        colors  = c("#ffd166","#619b8a","#3e5c76","#6a4c93","#a4161a"),
        limits  = c(0,100),
        na.value = "grey80"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y = if (xvar == "Sector") ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.y = if (xvar == "Sector") ggplot2::element_line() else ggplot2::element_blank()
      ) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::coord_fixed(ratio = 0.2)
  }
  
  p1 <- mk(sec,  "Sector")
  p2 <- mk(seas, "Season")
  p3 <- mk(hab,  "Habitat")
  
  (p1 | p2 | p3) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
}

# Generate the heatmap of sulfur completeness
heatmap_sulfur <- create_heatmaps(sulfur_means)
print(heatmap_sulfur)

# Adjust height dynamically
height_factor <- 0.5  # Adjust this value according to the spacing you need
num_pathways_sulfur <- 7

# Save figures with adjusted dimensions
#svglite("heatmap_SULFUR_completeness.svg", width = 8, height = max(8, num_pathways_sulfur * height_factor))
heatmap_sulfur
#dev.off()
