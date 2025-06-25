# figure 11

library(tidyverse)

# Function to plot N concentration with shaded bands
plot_n_concentration <- function(file_path, plot_title) {
  # Load the data
  df <- read_csv(file_path)
  
  # Pivot + annotate treatments and N level
  df_long <- df %>%
    pivot_longer(
      cols = c(legume_root, cereal_root, no_root),
      names_to = "zone",
      values_to = "N_concentration"
    ) %>%
    mutate(
      zone = recode(zone,
                    legume_root = "Legume",
                    cereal_root = "Cereal",
                    no_root = "no root"
      ),
      treatment_label = case_when(
        treatment_id == 1 ~ "Control",
        treatment_id == 2 ~ "Control - no lat",
        treatment_id == 3 ~ "Intercrop",
        treatment_id == 4 ~ "Intercrop - no lat",
        treatment_id == 5 ~ "70% Control",
        treatment_id == 6 ~ "70% Control - no lat",
        treatment_id == 7 ~ "70% Intercrop",
        treatment_id == 8 ~ "70% Intercrop - no lat",
        TRUE ~ paste("Treatment", treatment_id)
      ),
      Nlevel = ifelse(grepl("Nlimited", file_path), "70%", "100%") # Dynamically set Nlevel
    )
  
  # Calculate summary statistics for the shaded bands
  df_summary <- df_long %>%
    group_by(treatment_label, zone, timestep, Nlevel) %>%
    summarise(
      mean_N = mean(N_concentration, na.rm = TRUE),
      sd_N = sd(N_concentration, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = mean_N - sd_N,
      ymax = mean_N + sd_N
    )
  
  # Plot with labeled facets, lines, and shaded bands
  ggplot(df_summary, aes(x = timestep, y = mean_N, color = zone)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = zone), alpha = 0.2, colour = NA) + # Added colour = NA
    geom_line(linewidth = 1) +
    facet_wrap(~ treatment_label, ncol = 4) +
    labs(
      title = plot_title,
      x = "Timestep (day)",
      y = "Mean soil N concentration (µmol)",
      color = "Zone",
      fill = "Zone" # Add fill legend title
    ) +
    scale_color_manual(
      values = c(
        "Legume" = "red",
        "Cereal" = "orange",
        "no root" = "purple"
      )
    ) +
    scale_fill_manual( # Use scale_fill_manual for the ribbon colors
      values = c(
        "Legume" = "red",
        "Cereal" = "orange",
        "no root" = "blue"
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
}

# Generate plot for N-limited data
plot_n_concentration("zone_timeseries_Nlimited.csv", "Soil N concentration over time by zone (N-limited)")

# Generate plot for full N data
plot_n_concentration("zone_timeseries_fullN.csv", "Soil N concentration over time by zone (N-sufficient)")