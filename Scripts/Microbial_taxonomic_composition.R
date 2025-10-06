## San Quintin coastal lagoon - Sediment bioinformatic analysis
## Microbial taxonomic composition
## Prepared by: Jorge Rojas-Vargas

# Load libraries
library(ggplot2)   # plots
library(ampvis2)   # heatmaps
library(reshape2)  # matrix reshaping
library(svglite)   # export SVG files
library(ANCOMBC)   # differential abundance
library(dplyr)      # tables
library(tidyr)      # tables

# Load data
load("abund_arc_bac_class.RData")
load("abund_arc_bac_genus.RData")
load("abund_arc_bac_species.RData")
load("ps_count_genus.RData")
load("metadata_BSQ.RData")

## ---------------------------------------------
##
## Stacked bar plot at the Class level
##
## ---------------------------------------------

# Remove taxa with abundances < 1% for the stacked bar plot
abund_arc_bac_class_stack   <- abund_arc_bac_class[apply(abund_arc_bac_class,   1, function(x) any(x >= 1, na.rm = TRUE)), ]
abund_arc_bac_genus_stack   <- abund_arc_bac_genus[apply(abund_arc_bac_genus,   1, function(x) any(x >= 1, na.rm = TRUE)), ]
abund_arc_bac_species_stack <- abund_arc_bac_species[apply(abund_arc_bac_species, 1, function(x) any(x >= 1, na.rm = TRUE)), ]

# Add an "Other" row to recover the filtered-out taxa
column_sums <- colSums(abund_arc_bac_class_stack)
diff_values <- 100 - column_sums
other_row <- as.data.frame(t(diff_values))
rownames(other_row) <- "Other"
abund_arc_bac_class_stack <- rbind(other_row, abund_arc_bac_class_stack)
new_column_names <- metadata_BSQ$ID_1
colnames(abund_arc_bac_class_stack) <- new_column_names

# Order classes by their mean abundance (low to high)
abund_arc_bac_class_stack_filt <- apply(abund_arc_bac_class_stack, 1, mean)
top_class <- names(sort(abund_arc_bac_class_stack_filt, decreasing = FALSE)[1:19])
top_class <- c("Other", setdiff(top_class, "Other")) # Ensure "Other" stays first
top_class <- c("Other", "B_Nitrospira", "B_Verrucomicrobiae", "B_Gemmatimonadetes", 
               "B_Cytophagia", "B_Clostridia", "B_Anaerolineae", "B_Spirochaetia", 
               "B_Bacilli", "B_Bacteroidia", "B_Acidimicrobiia", "B_Planctomycetia", 
               "B_Betaproteobacteria", "B_Epsilonproteobacteria", "B_Actinobacteria", 
               "B_Flavobacteriia", "B_Alphaproteobacteria", "B_Deltaproteobacteria", 
               "B_Gammaproteobacteria")

# Order columns alphabetically
abund_arc_bac_class_stack <- abund_arc_bac_class_stack[, order(colnames(abund_arc_bac_class_stack))]

# Reshape table for plotting
dat_m <- melt(data.matrix(abund_arc_bac_class_stack))
colnames(dat_m) <- c("Class", "Sample", "Abundance")

# Add metadata
Sector <- metadata_BSQ$Sector[match(dat_m$Sample, metadata_BSQ$ID_1)]
Season <- metadata_BSQ$Season[match(dat_m$Sample, metadata_BSQ$ID_1)]
dat_m$Sector <- Sector
dat_m$Season <- Season

# Order taxa according to 'top_class'
dat_m$Class <- factor(dat_m$Class, levels = top_class)

# Define color palette for taxa
my_colors <- c(
  "Other" = "#666666",
  "B_Nitrospira" = "#c8ab37",
  "B_Verrucomicrobiae" = "#225500",
  "B_Gemmatimonadetes" = "#677821",
  "B_Cytophagia" = "#89a02c",
  "B_Clostridia" = "#bcd35f",
  "B_Anaerolineae" = "#dde9af",
  "B_Spirochaetia" = "#803300",
  "B_Bacilli"  = "#dd4000",
  "B_Bacteroidia" = "#ff6700",
  "B_Acidimicrobiia" = "#ff9300",
  "B_Planctomycetia" = "#ffcd00",
  "B_Betaproteobacteria" = "#ffeeaa",
  "B_Epsilonproteobacteria" = "#162d50",
  "B_Actinobacteria" = "#214478",
  "B_Flavobacteriia" = "#2c5aa0",
  "B_Alphaproteobacteria" = "#3771c8",
  "B_Deltaproteobacteria" = "#5f8dd3",
  "B_Gammaproteobacteria" = "#afc6e9"
)

# svglite("Fig_2A.svg", width = 10, height = 8)
ggplot(dat_m, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") + xlab("Samples") + ylab("Relative Abundance (%)") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 16), 
    legend.title = element_text(size = 18, face = "bold"),
    strip.text   = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = my_colors,
                    guide = guide_legend(nrow = 29)) +
  # facet_grid(~ Season, scales = "free", space = "free")
  facet_grid(~ Sector, scales = "free", space = "free")
# dev.off()

## ---------------------------------------------
##
## Heatmap of abundances at the Genus and Species levels
##
## ---------------------------------------------

# Genus level
amp_genus <- amp_load(otutable = abund_arc_bac_genus, metadata = metadata_BSQ)

genus_hm_season <- amp_heatmap(
  amp_genus,
  group_by = "Season",
  facet_by = "Sector",
  tax_aggregate = "OTU",
  tax_show = 20,
  color_vector = c("white","darkred"),
  plot_colorscale = "sqrt"
) + theme(
  axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
  axis.text.y = element_text(size = 10),
  legend.position = "right"
)

genus_hm_habitat <- amp_heatmap(
  amp_genus,
  group_by = "Habitat",
  facet_by = "Sector",
  tax_aggregate = "OTU",
  tax_show = 20,
  color_vector = c("white","darkred"),
  plot_colorscale = "sqrt"
) + theme(
  axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
  axis.text.y = element_text(size = 10),
  legend.position = "right"
)

# Species level
amp_species <- amp_load(otutable = abund_arc_bac_species, metadata = metadata_BSQ)

species_hm_season <- amp_heatmap(
  amp_species,
  group_by = "Season",
  facet_by = "Sector",
  tax_aggregate = "OTU",
  tax_show = 20,
  color_vector = c("white","darkred"),
  plot_colorscale = "sqrt"
) + theme(
  axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
  axis.text.y = element_text(size = 10),
  legend.position = "right"
)

species_hm_habitat <- amp_heatmap(
  amp_species,
  group_by = "Habitat",
  facet_by = "Sector",
  tax_aggregate = "OTU",
  tax_show = 20,
  color_vector = c("white","darkred"),
  plot_colorscale = "sqrt"
) + theme(
  axis.text.x = element_text(angle = 45, size = 10, vjust = 1),
  axis.text.y = element_text(size = 10),
  legend.position = "right"
)

# svglite("Fig_S1A.svg", width = 6, height = 6)
genus_hm_season
# dev.off()

# svglite("Fig_S2A.svg", width = 6, height = 8)
genus_hm_habitat
# dev.off()

# svglite("Fig_S1B.svg", width = 7.5, height = 8)
species_hm_season
# dev.off()

# svglite("Fig_S2B.svg", width = 7.5, height = 8)
species_hm_habitat
# dev.off()

## ---------------------------------------------
##
## Heatmap of abundances at the Genus and Species levels
##
## ---------------------------------------------

# Select data
abund_arc_bac <- abund_arc_bac_genus    # if analysis is done at the Genus level
abund_arc_bac <- abund_arc_bac_species  # if analysis is done at the Species level

# Abundance threshold parameter
thr <- 0.1

# Filter rows with at least one value ≥ thr
abund_arc_bac_filtered <- abund_arc_bac[rowSums(abund_arc_bac >= thr, na.rm = TRUE) > 0, , drop = FALSE]

## ---- Helper functions ----
presence_by_group <- function(dat, cols, thr = 0.1) {
  rowSums(dat[, cols, drop = FALSE] >= thr, na.rm = TRUE) > 0
}

# Counts for all combinations, including UNIQUE PER GROUP (k = 1)
combo_counts <- function(pres_list) {
  nm <- names(pres_list); n <- length(pres_list)
  out <- list()
  # Shared across all groups
  out$all_shared <- sum(Reduce(`&`, pres_list))
  # "Only" combinations (present in chosen groups, absent in the rest), now including k = 1
  for (k in 1:n) {
    for (cn in combn(nm, k, simplify = FALSE)) {
      in_vec  <- Reduce(`&`, pres_list[cn])
      others  <- setdiff(nm, cn)
      out_vec <- if (length(others) == 0) rep(TRUE, length(in_vec)) else
        Reduce(`&`, lapply(pres_list[others], function(x) !x))
      out[[paste0("only_", paste(cn, collapse = "_"))]] <- sum(in_vec & out_vec)
    }
  }
  out
}

shared_genera_and_sums <- function(dat, pres_list) {
  shared_all <- Reduce(`&`, pres_list)
  genera <- rownames(dat)[shared_all]
  abund_shared <- dat[genera, , drop = FALSE]
  sums <- colSums(abund_shared, na.rm = TRUE)
  list(genera = genera, sums = sums)
}

print_counts <- function(tag, lst) {
  cat("\n--", tag, "--\n", sep = "")
  for (nm in names(lst)) cat(sprintf("%-30s : %s\n", nm, lst[[nm]]))
}

## =============================================
## Between Sectors
## =============================================
groups_sectors <- list(
  inlet      = c("C_17", "C_18", "C_23", "C_24"),
  transition = c("C_19", "C_20", "C_25", "C_26"),
  inner      = c("C_21", "C_22", "C_27", "C_28")
)

pres_sectors <- lapply(groups_sectors, function(cols)
  presence_by_group(abund_arc_bac_filtered, cols, thr))

counts_sectors <- combo_counts(pres_sectors)  # includes: all_shared, only_inlet, only_transition, only_inner, ...
print_counts("Between Sectors (counts)", counts_sectors)

shared_sectors <- shared_genera_and_sums(abund_arc_bac_filtered, pres_sectors)
# cat("\nShared genera across all three sector groups:\n")
# print(shared_sectors$genera)
cat("\nSum abundance per sample (shared across all sectors):\n")
print(shared_sectors$sums)

## =============================================
## Between Habitats and Seasons
## =============================================
groups_habitats <- list(
  bare     = c("C_17", "C_19", "C_21", "C_23", "C_25", "C_27"),
  seagrass = c("C_18", "C_20", "C_22", "C_24", "C_26", "C_28"),
  intense  = c("C_17", "C_18", "C_19", "C_20", "C_21", "C_22"),
  relax    = c("C_23", "C_24", "C_25", "C_26", "C_27", "C_28")
)

pres_habitats <- lapply(groups_habitats, function(cols)
  presence_by_group(abund_arc_bac_filtered, cols, thr))

counts_habitats <- combo_counts(pres_habitats)
print_counts("Between Habitats & Seasons (counts)", counts_habitats)

shared_habitats <- shared_genera_and_sums(abund_arc_bac_filtered, pres_habitats)
# cat("\nShared genera across all four habitat/season groups:\n")
# print(shared_habitats$genera)
cat("\nSum abundance per sample (shared across all habitats/seasons):\n")
print(shared_habitats$sums)






## ---------------------------------------------
##
## Differential Abundance at genus level
##
## ---------------------------------------------


#ANCOM
out1 = ancombc(phyloseq = ps_count_genus, formula = "Sector + Season + Habitat",
               p_adj_method = "fdr", lib_cut = 0, prv_cut = 0.1,
               group = "Sector", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5, 
               max_iter = 100, conserve = FALSE, alpha = 0.05, global = TRUE)
res1 = out1$res
res1_global = out1$res_global

#LFC
tab_lfc1 = res1$lfc
col_name = c("Taxon","Intercept","TRAN-INL","INN-INL","RELAX-INTENSE","SEAGRASS-BARE")
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

# Filter df_lfc1
df_lfc1 <- df_lfc1[df_lfc1$taxon_id %in% rownames(abund_arc_bac_filtered), , drop = FALSE]

ps_small <- phyloseq::subset_samples(ps_count_genus, Sector %in% c("B_Transition", "C_Inner"))
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

# Filter df_lfc2
df_lfc2 <- df_lfc2[df_lfc2$taxon_id %in% rownames(abund_arc_bac_filtered), , drop = FALSE]


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


#svglite("Fig_S3A.svg", width=5, height=12)
p_Sector
#dev.off()

df_fig_Date = df_lfc1 %>% 
  filter(`RELAX-INTENSE` != 0) %>%
  dplyr::filter(`RELAX-INTENSE` < -1 | `RELAX-INTENSE` > 1) %>%
  transmute(taxon_id, 
            `RELAX-INTENSE` = round(`RELAX-INTENSE`, 2)) %>%
  pivot_longer(cols = `RELAX-INTENSE`, 
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

#svglite("Fig_S3C.svg", width=3, height=4)
p_Date
#dev.off()


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

#svglite("Fig_3C.svg", width=3, height=3)
p_Place
#dev.off()
