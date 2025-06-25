# Figure 8

library(tidyverse)   
library(patchwork)   
library(ggpubr)     
library(cowplot) 

# read the files
df_fullN   <- read_csv("output_fullN.csv")   %>% mutate(N_level = "N sufficient")
df_Nlimited <- read_csv("output_Nlimited.csv") %>% mutate(N_level = "N limited")
df <- bind_rows(df_fullN, df_Nlimited)

# pre process the data
df <- df %>%
  mutate(
    lateral_enabled   = as.character(lateral_enabled),
    lateral_label     = if_else(lateral_enabled == "TRUE", "On", "Off"),
    intercrop         = as.integer(intercrop),
    m3_uptake_cereal  = m3_uptake_cereal * 0.014,
    m3_uptake_cereal  = if_else(intercrop == 1, m3_uptake_cereal * 2, m3_uptake_cereal), # Doubling intercrop values
    N_level           = factor(N_level, levels = c("N sufficient", "N limited")),
    coverage_label    = if_else(abs(coverage - 1.0) < 0.01, "100% roots", "70% roots"),
    cropping_label    = factor(if_else(intercrop == 1, "Cereal intercrop", "Cereal monocrop"), # Make this a factor for consistent order
                               levels = c("Cereal monocrop", "Cereal intercrop")) # Explicitly set order for plotting
  )

# Define grouping variables for the t-tests: we want to compare 'cropping_label'
# for each unique combination of N_level, coverage_label, and lateral_label.
group_vars_for_ttest <- c("N_level", "coverage_label", "lateral_label")

ttest_results <- df %>%
  group_by(across(all_of(group_vars_for_ttest))) %>%
  group_modify(~{
    # Perform t-test comparing 'm3_uptake_cereal' between 'cropping_label' groups
    tst <- t.test(m3_uptake_cereal ~ cropping_label, data = .x)
    tibble(
      p_value = tst$p.value,
      # Store means for potential use, though y.position is dynamic in plot
      mean_monocrop = mean(.x$m3_uptake_cereal[.x$cropping_label == "Cereal monocrop"]),
      mean_intercrop = mean(.x$m3_uptake_cereal[.x$cropping_label == "Cereal intercrop"])
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj   = p.adjust(p_value, method = "BH"), # Benjamini-Hochberg adjustment
    signif  = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"),
    group1  = "Cereal monocrop",  # Matches the first group on x-axis for comparison
    group2  = "Cereal intercrop", # Matches the second group on x-axis for comparison
    label   = signif
  )

plot_combinations <- expand_grid(n_level = levels(df$N_level),
                                 coverage_label = unique(df$coverage_label))

plot_list <- purrr::map2(plot_combinations$n_level,
                         plot_combinations$coverage_label,
                         ~plot_single_butterfly_panel(df, .x, .y))

plot_sufficient <- wrap_plots(plot_list[plot_combinations$n_level == "N sufficient"],
                              ncol = 2, guides = "collect") +
  plot_annotation(title = "N sufficient",
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 14, face = "bold"),
                                legend.position  = "bottom",
                                # Keep axis.text.x = element_blank() here if you don't want redundant x-labels from inner plots
                                axis.text.x = element_blank())) &
  labs(y = expression("Uptake (mg N/" * m^3 * " soil)"),
       fill = "Lateral diffusion")

plot_limited <- wrap_plots(plot_list[plot_combinations$n_level == "N limited"],
                           ncol = 2, guides = "collect") +
  plot_annotation(title = "N limited",
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 14, face = "bold"),
                                legend.position  = "bottom",
                                axis.text.x = element_blank())) &
  labs(y = expression("Uptake (mg N/" * m^3 * " soil)"),
       fill = "Lateral diffusion")

print(plot_sufficient)
print(plot_limited)


plot_single_butterfly_panel <- function(data, n_level_val, coverage_label_val) {
  
  filtered_data <- filter(data, N_level == n_level_val,
                          coverage_label == coverage_label_val)
  
  stats_data    <- filter(ttest_results, N_level == n_level_val,
                          coverage_label == coverage_label_val)
  
  # Compute dynamic y-position and buffer
  plot_max_y <- max(filtered_data$m3_uptake_cereal, na.rm = TRUE)
  y_offset   <- 0.03 * plot_max_y
  y_buffer   <- 0.1 * plot_max_y
  stats_data <- mutate(stats_data, y.position = plot_max_y + y_offset)
  
  ggplot(filtered_data,
         aes(x = cropping_label,
             y = m3_uptake_cereal,
             fill = cropping_label,
             colour = cropping_label)) +
    
    geom_boxplot(width = 0.4,
                 outlier.shape = NA,
                 alpha = 0.8,
                 linewidth = 0.6,
                 position = position_dodge(width = 0.75)) +
    
    stat_pvalue_manual(
      data        = stats_data,
      label       = "label",
      y.position  = "y.position",
      xmin        = "group1",
      xmax        = "group2",
      tip.length  = 0.01,
      step.increase = 0,
      bracket.nudge.y = 0,
      bracket.shorten = 0,
      size        = 4,
      inherit.aes = FALSE,
      position    = position_dodge(width = 0.75)
    ) +
    
    scale_fill_manual(
      values = c("Cereal monocrop" = "#1b9e77", "Cereal intercrop" = "#c03256"),
      name = "Cropping system"
    ) +
    scale_colour_manual(
      values = c("Cereal monocrop" = "#1b9e77", "Cereal intercrop" = "#c03256"),
      guide = "none"
    ) +
    
    labs(title = coverage_label_val,
         x = NULL,
         y = expression("Uptake (mg N/" * m^3 * " soil)")) +
    
    facet_wrap(~ lateral_label, ncol = 2,
               labeller = as_labeller(c(On = "Lateral On", Off = "Lateral Off"))) +
    
    ylim(NA, plot_max_y + y_buffer) +  # Extend y-axis to prevent overlap
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 11))
}


# Prepare combinations
plot_combinations <- expand_grid(n_level = levels(df$N_level),
                                 coverage_label = unique(df$coverage_label))

# Generate all individual panels
plot_list <- purrr::map2(plot_combinations$n_level,
                         plot_combinations$coverage_label,
                         ~plot_single_butterfly_panel(df, .x, .y))

# Split plots by N level
sufficient_plots <- plot_list[plot_combinations$n_level == "N sufficient"]
limited_plots    <- plot_list[plot_combinations$n_level == "N limited"]

# Create centered row titles using cowplot
label_sufficient <- ggdraw() + 
  draw_label("N sufficient", fontface = "bold", size = 13, hjust = 0.5) +
  theme(plot.margin = margin(0, 0, 0, 0))

label_limited <- ggdraw() + 
  draw_label("N limited", fontface = "bold", size = 13, hjust = 0.5) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine all panels first, collect guides and tag A–D
panels_combined <- wrap_plots(c(sufficient_plots, limited_plots),
                              ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")

# Final layout: row label + panels + row label + panels
combined_plot <- (
  label_sufficient / 
    wrap_plots(sufficient_plots, ncol = 2) /
    label_limited / 
    wrap_plots(limited_plots, ncol = 2)
) +
  plot_layout(heights = c(0.06, 1, 0.06, 1)) &
  theme(legend.position = "bottom")

# Now add the tag levels and collect guides on full stack
final_plot <- wrap_elements(full = combined_plot) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

# Display
print(final_plot)
