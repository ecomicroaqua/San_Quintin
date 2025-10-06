## San Quintin coastal lagoon - Sediment bioinformatic analysis
## RDA analysis
## Prepared by: Jorge Rojas-Vargas

# Load libraries
library(vegan)    # RDA, vif.cca, decostand, ordistep
library(dplyr)    # Data manipulation (select, mutate, across)
library(svglite)  # SVG export for plots

# Load data
load("abund_arc_bac_genus.RData")
load("count_metal_genes.RData")
load("count_drug_genes.RData")
load("count_nitrogen_genes.RData")
load("count_sulfur_genes.RData")
load("count_virulence_genes.RData")
load("metadata_BSQ.RData")
load("physchememical_characteristics.RData")


## ------------------------------------------------------
##
## MICROBIOME AND PHYSICOCHEMICAL FEATURES
##
## ------------------------------------------------------


## =============================================
## Preparing data for analyses
## =============================================

# Filter column names that are not "UID", "Sample" or contain "Rate"
data_chem_num <- physchem[, !grepl("Rate|UID|Sample", names(physchem))]
summary(data_chem_num)

vars <- c("pH", "Sand", "Silt", "TOC_umol", "TN_umol",
          "NH4_g", "NO3_ug", "NO2_ug", "Fe_III_mg", "Fe_II_mg")

data_chem_num <- data_chem_num %>% dplyr::select(all_of(vars))


## =============================================
## Pearson correlation to reduce collinearity
## =============================================

# Pairwise Pearson correlations among environmental variables
cor_matrix_bac <- cor(data_chem_num) # Correlation matrix

# Collinearity in the physicochemical data
#svglite("pearson_phychem.svg", width = 6, height = 5)
heatmap(abs(cor_matrix_bac), 
        # Compute Pearson correlation (note they are absolute values)
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
data_chem_num.new <- data_chem_num[ , !apply(cor_matrix_bac_rm,
                                             2,
                                             function(x) any(x > 0.7))]

# Less correlated variables (Pearson < 0.7)
less.cor.var <- colnames(data_chem_num.new)
less.cor.var
## Seven parameters: "Sand", "TN_umol", "NH4_g", "NO3_ug", "NO2_ug", "Fe_III_mg", "Fe_II_mg" 

# Since there are triplicates for each sample, the mean of the parameters is calculated
data_chem_tr <- aggregate(data_chem_num.new, by = list(Sample = physchem$Sample), FUN = mean)
data_chem_tr <- data_chem_tr %>%
  mutate(across(.cols = 2:ncol(data_chem_tr), .fns = ~ as.numeric(as.character(.))))
rownames(data_chem_tr) <- data_chem_tr$Sample
data_chem_tr$Sample <- NULL

# Also do it for data_chem_num
data_chem_num_tr <- aggregate(data_chem_num, by = list(Sample = physchem$Sample), FUN = mean)
data_chem_num_tr <- data_chem_num_tr %>%
  mutate(across(.cols = 2:ncol(data_chem_num_tr), .fns = ~ as.numeric(as.character(.))))
rownames(data_chem_num_tr) <- data_chem_num_tr$Sample
data_chem_num_tr$Sample <- NULL


## =============================================
## Using VIF approach
## =============================================

# Scale and center variables
data_chem_num.new.z <- decostand(data_chem_tr, method = "standardize")

# Apply Hellinger transformation to correct for the double zero problem in abundance table
ab_genus_rda <- t(abund_arc_bac_genus)
ab_genus_hellinger <- as.data.frame(decostand(ab_genus_rda, method = "hellinger"))

# Full RDA model
rda_ab_chem <- rda(ab_genus_hellinger ~ ., data = data_chem_num.new.z)
summary(rda_ab_chem)

vif_criteria <- vif.cca(rda_ab_chem)
vif_criteria

vars.envfit <- names(which(vif_criteria < 10))

data_chem_tr_2 <- data_chem_tr[, vars.envfit]

# Scale and center variables with less correlated variables
data_chem_num.new.z_2 <- decostand(data_chem_tr_2, method = "standardize")

# Full RDA model with less correlated variables
rda_ab_chem_2 <- rda(ab_genus_hellinger ~ ., data = data_chem_num.new.z_2)
summary(rda_ab_chem_2)


## =============================================
## Using ordistep approach
## =============================================

# Perform forward and backward selection of explanatory variables
set.seed(123)
fwd.sel <- ordistep(rda(ab_genus_hellinger ~ 1, data = data_chem_num.new.z_2), # lower model limit (simple!)
                    scope = formula(rda_ab_chem_2), # upper model limit (the "full" model)
                    direction = "both",
                    alpha = 0.01,
                    R2scope = FALSE, # can surpass the "full" model's R2
                    pstep = 999,
                    trace = TRUE) # TRUE to see the selection process!

# Get the formula
fwd.sel$call

# Get the selected variables
vars.ordistep <- gsub('^. ', '', rownames(fwd.sel$anova))
vars.ordistep


## =============================================
## Combining vif and ordistep approaches
## =============================================

# Get the common variables from the vif and ordistep approaches
vars.bac <- intersect(vars.ordistep, vars.envfit)
vars.bac

# Write our new model
matrix.kept <- data_chem_num.new.z[, vars.bac]
rda_ab_chem.signif <- rda(ab_genus_hellinger ~ ., data = matrix.kept)

vif_criteria <- vif.cca(rda_ab_chem.signif)
vif_criteria

summary(rda_ab_chem.signif)

# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_ab_chem.signif)

set.seed(123)
anova(rda_ab_chem.signif)

set.seed(123)
anova(rda_ab_chem.signif,by="axis")

set.seed(123)
anova(rda_ab_chem.signif,by="term")


## =============================================
## Plotting RDA analysis
## =============================================

# Merge the sector information with the transformed data
data_chem_num.new.z <- cbind(data_chem_num.new.z, Sector = metadata_BSQ$Sector)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Habitat = metadata_BSQ$Habitat)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Season = metadata_BSQ$Season)

# Define shapes, colors, and line types
sector_shapes <- c("A_Inlet" = 21, "B_Transition" = 22, "C_Inner" = 23)
habitat_colors <- c("Bare" = "orange", "Seagrass" = "darkgreen")
season_linetype <- c("Intense" = "solid", "Relax" = "dashed")

explained_var <- summary(rda_ab_chem)$cont$importance[2, ] * 100
xlab_text <- paste0("RDA1 (", round(explained_var[1], 2), "%)")
ylab_text <- paste0("RDA2 (", round(explained_var[2], 2), "%)")

#svglite("Fig_4D.svg", width = 8, height = 5)
plot(rda_ab_chem.signif, type = "n", main = "Redundancy Analysis - chemical",
     xlim = c(-2, 2), ylim = c(-1, 1), xlab = xlab_text, ylab = ylab_text)

# Add points with shapes, colors, and lines
with(data_chem_num.new.z, points(rda_ab_chem.signif, "sites", bg = habitat_colors[Habitat], cex = 2, 
                                 pch = sector_shapes[Sector], 
                                 lwd = 1.5,
                                 col = "black",
                                 lty = season_linetype[Season]))

# Add text and labels
text(rda_ab_chem.signif, display = "sites", cex = .5, col = "gray")
text(rda_ab_chem.signif, display = "bp", cex = .9, col = "blue")

# Add legend
legend("topright", legend = c("Bare", "Seagrass"), 
       fill = c("orange", "darkgreen"), title = "Habitat", box.lwd = 0)
legend("topright", legend = c("A_Inlet", "B_Transition", "C_Inner"),
       pch = c(21, 22, 23), title = "Sector", inset = c(0.1, 0.2))
legend("topright", legend = c("Intense", "Relax"), 
       lty = c(1, 2), col = "black", title = "Season", inset = c(0.1, 0.4), box.lwd = 0)
#dev.off()

# Get coordinates of explanatory variables in RDA space
env_coords <- scores(rda_ab_chem.signif, display = "bp", scaling = 2)
env_coords

# Function to calculate angle between two vectors
angle_between <- function(a, b) {
  cos_theta <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  theta_rad <- acos(cos_theta)
  theta_deg <- theta_rad * 180 / pi
  return(theta_deg)
}

# Calculate angle between NO3_ug and Fe_III_mg
vec1 <- env_coords["NO3_ug", ]
vec2 <- env_coords["Fe_III_mg", ]

angle_between(vec1, vec2)



## ------------------------------------------------------
##
## GENES AND PHYSICOCHEMICAL FEATURES
##
## ------------------------------------------------------


# Bind all tables
count_all <- rbind(count_metal,
                   count_drug,
                   count_nitrogen,
                   count_sulfur,
                   count_virulence)

# summary statistics for all objects (min, mean, max, etc.)
summary(count_all)


## =============================================
## Using VIF approach
## =============================================

# Apply Hellinger transformation to correct for the double zero problem in abundance table
count_all_rda <- t(count_all)
count_all_hellinger <- as.data.frame(decostand(count_all_rda, method = "hellinger"))

# Scale and center variables
data_chem_num.new.z <- decostand(data_chem_tr, method = "standardize")

# Full RDA model
rda_ab_chem <- rda(count_all_hellinger ~ ., data = data_chem_num.new.z)
summary(rda_ab_chem)

vif_criteria <- vif.cca(rda_ab_chem)
vif_criteria
#      Sand   TN_umol     NH4_g    NO3_ug    NO2_ug Fe_III_mg  Fe_II_mg 
# 8.302483  6.810308  7.968274  2.346754  6.522657  8.121830 10.475641 

vars.envfit <- names(which(vif_criteria < 10))

data_chem_tr_2 <- data_chem_tr[, vars.envfit]

# Scale and center variables with less correlated variables
data_chem_num.new.z_2 <- decostand(data_chem_tr_2, method = "standardize")

# Full RDA model with less correlated variables
rda_ab_chem_2 <- rda(count_all_hellinger ~ ., data = data_chem_num.new.z_2)
summary(rda_ab_chem_2)

## =============================================
## Using ordistep approach
## =============================================

# Perform forward and backward selection of explanatory variables
set.seed(123)
fwd.sel <- ordistep(rda(count_all_hellinger ~ 1, data = data_chem_num.new.z_2), # lower model limit (simple!)
                    scope = formula(rda_ab_chem_2), # upper model limit (the "full" model)
                    direction = "both",
                    alpha = 0.01,
                    R2scope = FALSE, # can surpass the "full" model's R2
                    pstep = 999,
                    trace = TRUE) # TRUE to see the selection process!

# Get the formula
fwd.sel$call

# Get the selected variables
vars.ordistep <- gsub('^. ', '', rownames(fwd.sel$anova))
vars.ordistep


## =============================================
## Combining vif and ordistep approaches
## =============================================

# Get the common variables from the vif and ordistep approaches
vars.bac <- intersect(vars.ordistep, vars.envfit)
vars.bac

# Write our new model
matrix.kept <- data_chem_num.new.z[, vars.bac]
rda_ab_chem.signif <- rda(count_all_hellinger ~ ., data = matrix.kept)

vif_criteria <- vif.cca(rda_ab_chem.signif)
vif_criteria
# NO3_ug Fe_III_mg 
# 1.180111  1.180111

summary(rda_ab_chem.signif)
#               Inertia Proportion
# Total         0.05219     1.0000
# Constrained   0.01864     0.3571
# Unconstrained 0.03355     0.6429

# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_ab_chem.signif)

set.seed(123)
anova(rda_ab_chem.signif)

set.seed(123)
anova(rda_ab_chem.signif,by="axis")

set.seed(123)
anova(rda_ab_chem.signif,by="term")


## =============================================
## Plotting RDA analysis
## =============================================

# Merge the sector information with the transformed data
data_chem_num.new.z <- cbind(data_chem_num.new.z, Sector = metadata_BSQ$Sector)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Habitat = metadata_BSQ$Habitat)
data_chem_num.new.z <- cbind(data_chem_num.new.z, Season = metadata_BSQ$Season)

# Define shapes, colors, and line types
sector_shapes <- c("A_Inlet" = 21, "B_Transition" = 22, "C_Inner" = 23)
habitat_colors <- c("Bare" = "orange", "Seagrass" = "darkgreen")
season_linetype <- c("Intense" = "solid", "Relax" = "dashed")

explained_var <- summary(rda_ab_chem)$cont$importance[2, ] * 100
xlab_text <- paste0("RDA1 (", round(explained_var[1], 2), "%)")
ylab_text <- paste0("RDA2 (", round(explained_var[2], 2), "%)")

#svglite("Fig_4D.svg", width = 8, height = 5)
plot(rda_ab_chem.signif, type = "n", main = "Redundancy Analysis - chemical",
     xlim = c(-2, 2), ylim = c(-1, 1), xlab = xlab_text, ylab = ylab_text)

# Add points with shapes, colors, and lines
with(data_chem_num.new.z, points(rda_ab_chem.signif, "sites", bg = habitat_colors[Habitat], cex = 2, 
                                 pch = sector_shapes[Sector], 
                                 lwd = 1.5,
                                 col = "black",
                                 lty = season_linetype[Season]))

# Add text and labels
text(rda_ab_chem.signif, display = "sites", cex = .5, col = "gray")
text(rda_ab_chem.signif, display = "bp", cex = .9, col = "blue")

# Add legend
legend("topright", legend = c("Bare", "Seagrass"), 
       fill = c("orange", "darkgreen"), title = "Habitat", box.lwd = 0)
legend("topright", legend = c("A_Inlet", "B_Transition", "C_Inner"),
       pch = c(21, 22, 23), title = "Sector", inset = c(0.1, 0.2))
legend("topright", legend = c("Intense", "Relax"), 
       lty = c(1, 2), col = "black", title = "Season", inset = c(0.1, 0.4), box.lwd = 0)
#dev.off()

# Get coordinates of explanatory variables in RDA space
env_coords <- scores(rda_ab_chem.signif, display = "bp", scaling = 2)
env_coords

# Function to calculate angle between two vectors
angle_between <- function(a, b) {
  cos_theta <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  theta_rad <- acos(cos_theta)
  theta_deg <- theta_rad * 180 / pi
  return(theta_deg)
}

# Calculate angle between NO3_ug and Fe_III_mg
vec1 <- env_coords["NO3_ug", ]
vec2 <- env_coords["Fe_III_mg", ]

angle_between(vec1, vec2)
