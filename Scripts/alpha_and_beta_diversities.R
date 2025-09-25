## San Quintin coastal lagoon - Sediment bioinformatic analysis
## Alpha and Beta diversities
## Prepared by: Jorge Rojas-Vargas

# Load packages
library(phyloseq)
library(microbiome)   # CLR transformation
library(vegan)        # Alpha/beta diversity, PERMANOVA, betadisper
library(RVAideMemoire)# Pairwise permutational MANOVA
library(ggplot2)      # Plot themes/geoms used by plot_richness
library(patchwork)    # Combine ggplots with "/" operator
library(svglite)    # Uncomment to write SVG files

# Phyloseq objects
load("ps_count_class.RData")
load("ps_count_genus.RData")
load("ps_count_species.RData")

## ---------------------------------------------
##
## Alpha diversity
##
## ---------------------------------------------

## -- Genus level --

a_my_comparisons <- list(c("A_Inlet", "B_Transition"), c("A_Inlet", "C_Inner"), c("C_Inner", "B_Transition"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_sector <- plot_richness(ps_count_genus, x = "Sector", measures = c("Observed","Shannon"), color = "Sector") +
  theme_bw() +
  geom_boxplot(lwd = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("#2a7fffff","#808080ff","#c8ab37ff"))
# + ggpubr::stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Relaxing", "Upwelling"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_season <- plot_richness(ps_count_genus, x = "Season", measures = c("Observed","Shannon"), color = "Season") +
  theme_bw() +
  geom_boxplot(lwd = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("#ff0000ff","#000080ff"))
# + ggpubr::stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Bare", "Seagrass"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_habitat <- plot_richness(ps_count_genus, x = "Habitat", measures = c("Observed","Shannon"), color = "Habitat") +
  theme_bw() +
  geom_boxplot(lwd = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("#ff6600ff","#338000ff"))
# + ggpubr::stat_compare_means(method = "kruskal.test", comparisons = a_my_comparisons, label = "p.signif")

# svglite("Fig_3A.svg", width = 5, height = 7)
ad_sector / ad_habitat / ad_season
# dev.off()

## -- Species level --

a_my_comparisons <- list(c("A_Inlet", "B_Transition"), c("A_Inlet", "C_Inner"), c("C_Inner", "B_Transition"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_sector <- plot_richness(ps_count_species, x = "Sector", measures = c("Observed","Shannon"), color = "Sector") +
  theme_bw() +
  geom_boxplot(lwd = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("#2a7fffff","#808080ff","#c8ab37ff"))
# + ggpubr::stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Relaxing", "Upwelling"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_season <- plot_richness(ps_count_species, x = "Season", measures = c("Observed","Shannon"), color = "Season") +
  theme_bw() +
  geom_boxplot(lwd = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("#ff0000ff","#000080ff"))
# + ggpubr::stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

a_my_comparisons <- list(c("Bare", "Seagrass"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****","***","**","*","ns"))

ad_habitat <- plot_richness(ps_count_species, x = "Habitat", measures = c("Observed","Shannon"), color = "Habitat") +
  theme_bw() +
  geom_boxplot(lwd = 0.7, alpha = 0.6) +
  scale_color_manual(values = c("#ff6600ff","#338000ff"))
# + ggpubr::stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

# svglite("Fig_3B.svg", width = 5, height = 7)
ad_sector / ad_habitat / ad_season
# dev.off()

## ---------------------------------------------
##
## Beta diversity - PCoA and dispersion analysis
##
## ---------------------------------------------

set.seed(123)

# Subset a 'dist' object to samples with non-NA groups
subset_dist <- function(d, keep_index) {
  m <- as.matrix(d)
  stats::as.dist(m[keep_index, keep_index, drop = FALSE])
}

# Core analysis for one phyloseq object (one taxonomic level)
analyze_ps <- function(ps, level_name,
                       factors = c("Sector", "Habitat", "Season"),
                       nperm_pairwise = 10000,
                       nperm_adonis  = 999) {
  
  # CLR transform and Aitchison distance
  ps_clr <- microbiome::transform(ps, "clr")
  ait    <- phyloseq::distance(ps_clr, method = "euclidean")
  
  # Metadata
  meta <- as.data.frame(phyloseq::sample_data(ps_clr))
  for (f in factors) if (!is.factor(meta[[f]])) meta[[f]] <- factor(meta[[f]])
  
  # Run stats per factor
  results <- lapply(factors, function(f) {
    idx  <- !is.na(meta[[f]])
    if (!any(idx)) stop("All values are NA for factor: ", f)
    
    d_sub <- subset_dist(ait, idx)
    grp   <- droplevels(factor(meta[[f]][idx]))
    dfgrp <- data.frame(grp = grp)
    
    # PERMANOVA
    adonis_fit <- vegan::adonis2(d_sub ~ grp, data = dfgrp, permutations = nperm_adonis)
    
    # Pairwise MANOVA (permutation)
    pw_fit <- RVAideMemoire::pairwise.perm.manova(d_sub, grp, nperm = nperm_pairwise, p.method = "BH")
    
    # Dispersion and test (PCoA is computed internally by betadisper)
    bd  <- vegan::betadisper(d_sub, grp)
    bdt <- vegan::permutest(bd, permutations = nperm_adonis)
    
    # --- Plot ---
    # svglite::svglite(paste0("Dispersion_", level_name, "_", f, ".svg"), width = 4, height = 4)
    plot(bd, main = paste0("Dispersion - ", level_name, " (", f, ")"), sub = "")
    legend("topright", legend = levels(grp), pch = 1, col = seq_along(levels(grp)))
    # dev.off()
    
    list(
      adonis = adonis_fit,
      pairwise = pw_fit,
      betadisper = bd,
      betadisper_permutest = bdt
    )
  })
  
  names(results) <- factors
  return(results)
}

# Run for all three taxonomic levels
tax_levels <- list(
  class   = ps_count_class,
  genus   = ps_count_genus,
  species = ps_count_species
)

all_results <- lapply(names(tax_levels), function(nm) {
  message(">>> Running analysis for: ", nm)
  analyze_ps(tax_levels[[nm]], level_name = nm)
})
names(all_results) <- names(tax_levels)

# Print concise summaries to console
for (lvl in names(all_results)) {
  cat("\n==============================\nLevel:", lvl, "\n")
  lvl_res <- all_results[[lvl]]
  for (f in names(lvl_res)) {
    cat("\n--- Factor:", f, "---\n")
    cat("\nPERMANOVA (adonis2):\n"); print(lvl_res[[f]]$adonis)
    cat("\nPairwise MANOVA (BH adj):\n"); print(lvl_res[[f]]$pairwise)
    cat("\nDispersion test (betadisper -> permutest):\n"); print(lvl_res[[f]]$betadisper_permutest)
  }
}

# Notes:
# - If zero counts are present, consider applying a pseudocount / zero-replacement
#   before CLR to avoid numerical issues.
# - To save SVG plots, uncomment the svglite() and dev.off() lines where indicated.


