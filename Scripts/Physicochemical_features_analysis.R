## San Quintin coastal lagoon - Sediment bioinformatic analysis
## Physicochemical characteristics analysis
## Elaborated by: Jorge Rojas-Vargas


# Load libraries
library(ARTool)      # ART test (Aligned Rank Transform models)
library(dplyr)       # Data manipulation
library(FSA)         # Dunn test for multiple comparisons
library(svglite)     # Export plots to SVG format
library(ggplot2)     # Plotting system
library(tidyr)       # Data reshaping
library(patchwork)   # Combine multiple ggplot figures

# Load data
load("physchememical_characteristics.RData")
load("metadata_physchem.RData")

# Merge abundance tables and metadata
combined_data <- merge(physchem, metadata_physchem, by.x = "UID", by.y = "UID")


## ---------------------------------------
##
## Box plots
##
## ---------------------------------------


# Create a list of columns of interest
columnas_de_interes <- c(
  3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

# Custom color mapping
my_colors <- c(
  A_Inlet      = "#56b4e9",
  B_Transition = "#808080",
  C_Inner      = "#c8ab37",  
  Intense      = "#FF0000FF",
  Relax        = "#000080FF",
  Bare         = "orange",
  Seagrass     = "darkgreen"
)

# Loop to generate and store figures in variables a1, a2, a3, ...
for (i in 1:length(columnas_de_interes)) {
  columna_index <- columnas_de_interes[i]
  columna_de_interes <- colnames(combined_data)[columna_index]
  # Transformation to long format
  df_long <- combined_data %>% 
    dplyr::select(.data[[columna_de_interes]], Sector, Season, Habitat) %>%           
    pivot_longer(cols = c(Sector, Season, Habitat),
                 names_to  = "GroupType",
                 values_to = "Group")
  
  # Set the order of the seven categories on the x-axis
  level_order <- c("A_Inlet", "B_Transition", "C_Inner",
                   "Intense",  "Relax",
                   "Bare",     "Seagrass")
  df_long$Group <- factor(df_long$Group, levels = level_order)
  
  # Create the boxplot with t test using the new combined column
  figura <- ggplot(df_long, aes(x = Group, y = .data[[columna_de_interes]], 
                                fill = Group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
    geom_jitter(width = 0.0, size = 0.8, alpha = 0.5, colour = "black") +
    scale_fill_manual(values = my_colors) +
    labs(x = "", y = "", fill = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(columna_de_interes))
  
  # Store the figure in a dynamic variable (a1, a2, a3, ...)
  assign(paste0("a", i), figura)
}


#svglite("Fig_2.svg", width=20, height=10)
((a1 | a2 | a3 | a4 | a5) / (a6 | a7 | a8 | a9 | a10))
#dev.off()




## ---------------------------------------
##
## Statistical analysis
##
## ---------------------------------------

# Adjust model factors
combined_data$Sector <- as.factor(combined_data$Sector)
combined_data$Season <- as.factor(combined_data$Season)
combined_data$Habitat <- as.factor(combined_data$Habitat)

# Define the vector of variables to analyze
vars <- names(combined_data)[3:12]


## -- ART models --

# Function to extract the ANOVA from an ART model
get_art_anova <- function(var, data) {
  # Fit the model
  mod <- art(
    formula = as.formula(paste(var, "~ Sector*Season*Habitat")),
    data    = data
  )
  # Extract the ANOVA type III table
  an <- anova(mod)
  # Convert to data.frame and add the response
  df <- as.data.frame(an)
  df$Term     <- rownames(df)
  df$Response <- var
  rownames(df) <- NULL
  # Reorder columns
  df %>%
    dplyr::select(Response, Term, Df, Df.res, `F value`, `Pr(>F)`)
}

# Iterate over all variables and combine results
all_anovas <- bind_rows(
  lapply(vars, get_art_anova, data = combined_data)
)

# Readable format
all_anovas <- all_anovas %>%
  mutate(
    `Pr(>F)` = signif(`Pr(>F)`, 3),
    `F value` = round(`F value`, 2)
  )

# Show table
print(all_anovas)

#write.csv(all_anovas, "Table_S2.csv")



## -- Dunn test --

# Function to apply Dunn test to multiple variables
run_dunn_tests <- function(vars, data, group_var = "Sector", method = "bh") {
  results <- lapply(vars, function(v) {
    formula <- as.formula(paste(v, "~", group_var))
    res <- dunnTest(formula, data = data, method = method)$res
    res$Variable <- v
    res
  })
  dplyr::bind_rows(results)
}

# Run the function
dunn_results <- run_dunn_tests(vars, combined_data)

# Show results
print(dunn_results)



## -- Wilcox test --

# Wilcoxon for Season and Habitat
wilcox_results <- lapply(c("Season","Habitat"), function(var){
  lapply(vars, function(col){
    wt <- wilcox.test(
      combined_data[[col]] ~ combined_data[[var]],
      data = combined_data,
      exact = FALSE
    )
    data.frame(
      Property  = col,
      Variable = var,
      P.value  = wt$p.value
    )
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  # BH correction by Variable group
  group_by(Variable) %>%
  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%
  ungroup()

# Show all results
print(wilcox_results, n = 20)



## -- Statistical values of properties for each factor --

# Function to obtain the statistics
get_stats <- function(group_var) {
  combined_data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      across(
        all_of(vars),
        list(
          mean   = ~ mean(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE),
          sd     = ~ sd(.x, na.rm = TRUE),
          min    = ~ min(.x, na.rm = TRUE),
          max    = ~ max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}


# Statistics by Sector
stats_by_sector <- get_stats("Sector")
print(t(stats_by_sector))

# Statistics by Season
stats_by_season <- get_stats("Season")
print(t(stats_by_season))

# Statistics by Habitat
stats_by_habitat <- get_stats("Habitat")
print(t(stats_by_habitat))
