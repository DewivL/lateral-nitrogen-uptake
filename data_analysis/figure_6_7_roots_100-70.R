# Figure 6, 7

# ---------------------------------------------------------------
# 1. Load libraries
# ---------------------------------------------------------------
library(tidyverse)   # ggplot2, dplyr, readr …
library(patchwork)   # arranging plots
library(ggpubr)      # stat_pvalue_manual()

# ---------------------------------------------------------------
# 2. Load and combine data
# ---------------------------------------------------------------
df_fullN   <- read_csv("output_fullN.csv")   %>% mutate(N_level = "N sufficient")
df_Nlimited <- read_csv("output_Nlimited.csv") %>% mutate(N_level = "N limited")
df <- bind_rows(df_fullN, df_Nlimited)

# ---------------------------------------------------------------
# 3. Pre-process data
# ---------------------------------------------------------------
df <- df %>%
  mutate(
    lateral_enabled   = as.character(lateral_enabled),
    lateral_label     = if_else(lateral_enabled == "TRUE", "On", "Off"),
    intercrop         = as.integer(intercrop),
    m3_uptake_cereal  = m3_uptake_cereal * 0.014,
    m3_uptake_cereal  = if_else(intercrop == 1, m3_uptake_cereal * 2, m3_uptake_cereal), # NEW LINE FOR DOUBLING
    N_level           = factor(N_level, levels = c("N sufficient", "N limited")),
    coverage_label    = if_else(abs(coverage - 1.0) < 0.01, "100% roots", "70% roots"),
    cropping_label    = if_else(intercrop == 1, "Cereal intercrop", "Cereal monocrop")
  )

# ---------------------------------------------------------------
# 4. Run t-tests & build annotation data frame
# ---------------------------------------------------------------
group_vars <- c("N_level", "coverage_label", "cropping_label")

ttest_results <- df %>%
  group_by(across(all_of(group_vars))) %>%
  group_modify(~{
    tst <- t.test(m3_uptake_cereal ~ lateral_label, data = .x)
    tibble(mean_off  = mean(.x$m3_uptake_cereal[.x$lateral_label == "Off"]),
           mean_on   = mean(.x$m3_uptake_cereal[.x$lateral_label == "On"]),
           p_value   = tst$p.value)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj   = p.adjust(p_value, method = "BH"),
    signif  = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"),
    group1  = "Off",
    group2  = "On",
    label   = signif,
    y.position = pmax(mean_off, mean_on) + 0.10 * max(df$m3_uptake_cereal)
  )

# (optional) hide ‘ns’ results
#ttest_results <- filter(ttest_results, label != "ns")

# ---------------------------------------------------------------
# 5. Plotting function for one (N-level × coverage) panel pair
# ---------------------------------------------------------------
plot_single_butterfly_panel <- function(data, n_level_val, coverage_label_val) {
  
  filtered_data <- filter(data, N_level == n_level_val,
                          coverage_label == coverage_label_val)
  
  stats_data    <- filter(ttest_results, N_level == n_level_val,
                          coverage_label == coverage_label_val)
  
  ggplot(filtered_data,
         aes(x = lateral_label,
             y = m3_uptake_cereal,
             fill = lateral_label,
             colour = lateral_label)) +
    
    geom_boxplot(width = 0.4,
                 outlier.shape = NA,
                 alpha = 0.7,
                 linewidth = 0.6) +
    
    stat_pvalue_manual(
      data        = stats_data,
      label       = "label",
      y.position  = "y.position",
      xmin        = "group1",
      xmax        = "group2",
      tip.length  = 0.01,
      size        = 4,
      inherit.aes = FALSE       
    ) +
    
    scale_fill_manual(values = c(Off = "#E69F00", On = "#56B4E9"),
                      name   = "Lateral diffusion") +
    scale_colour_manual(values = c(Off = "#E69F00", On = "#56B4E9"),
                        guide  = "none") +
    labs(title = coverage_label_val,
          x = NULL,
         y = expression("Uptake (mg N/" * m^3 * " soil)")) +
    facet_wrap(~ cropping_label, ncol = 2) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_blank(),,
          legend.position = "bottom",
          
          strip.text  = element_text(face = "bold", size = 11))
}

# ---------------------------------------------------------------
# 6. Generate and display the four panels
# ---------------------------------------------------------------
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


# for figure 5

# Use max values to set appropriate annotation positions
max_A <- max(df_A$m3_uptake_cereal)
max_B <- max(df_B$m3_uptake_cereal)

ttest_A <- ttest_A %>% mutate(y.position = max_A + 0.03)
ttest_B <- ttest_B %>% mutate(y.position = max_B + 0.03)

# --- Plot A: N sufficient ---
plot_A <- ggplot(df_A,
                 aes(x = lateral_label,
                     y = m3_uptake_cereal,
                     fill = lateral_label,
                     colour = lateral_label)) +
  geom_boxplot(width = 0.4,
               outlier.shape = NA,
               alpha = 0.7,
               linewidth = 0.6) +
  stat_pvalue_manual(
    data        = ttest_A,
    label       = "label",
    y.position  = "y.position",
    xmin        = "group1",
    xmax        = "group2",
    tip.length  = 0.01,
    size        = 4,
    inherit.aes = FALSE       
  ) +
  scale_fill_manual(values = c(Off = "#E69F00", On = "#56B4E9"),
                    name   = "Lateral diffusion") +
  scale_colour_manual(values = c(Off = "#E69F00", On = "#56B4E9"),
                      guide  = "none") +
  labs(title = "A   N sufficient",
       x = NULL,
       y = expression("Total N uptake (mg N/" * m^3 * " soil)")) +
  coord_cartesian(ylim = c(min(df_A$m3_uptake_cereal) - 0.02, max_A + 0.06)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# --- Plot B: N limited ---
plot_B <- ggplot(df_B,
                 aes(x = lateral_label,
                     y = m3_uptake_cereal,
                     fill = lateral_label,
                     colour = lateral_label)) +
  geom_boxplot(width = 0.4,
               outlier.shape = NA,
               alpha = 0.7,
               linewidth = 0.6) +
  stat_pvalue_manual(
    data        = ttest_B,
    label       = "label",
    y.position  = "y.position",
    xmin        = "group1",
    xmax        = "group2",
    tip.length  = 0.01,
    size        = 4,
    inherit.aes = FALSE       
  ) +
  scale_fill_manual(values = c(Off = "#E69F00", On = "#56B4E9"),
                    name   = "Lateral diffusion") +
  scale_colour_manual(values = c(Off = "#E69F00", On = "#56B4E9"),
                      guide  = "none") +
  labs(title = "B   N limited",
       x = NULL,
       y = expression("Total N uptake (mg N/" * m^3 * " soil)")) +
  coord_cartesian(ylim = c(min(df_B$m3_uptake_cereal) - 0.02, max_B + 0.06)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# --- Combine and display ---
combined_plot <- plot_A + plot_B + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)
