# figure 9

library(tidyverse)

# Load and process both files
process_uptake_file <- function(path, label) {
  read_csv(path) %>%
    mutate(
      no_root = as.numeric(no_root),
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
      Nlevel = label
    ) %>%
    pivot_longer(cols = c(legume_root, cereal_root, no_root),
                 names_to = "zone", values_to = "uptake") %>%
    mutate(
      zone = recode(zone,
                    "legume_root" = "Legume root zone",
                    "cereal_root" = "Cereal root zone",
                    "no_root" = "No root zone")
    ) %>%
    filter(zone == "Cereal root zone") %>%
    group_by(treatment_label, replicate, zone, Nlevel) %>%
    arrange(timestep) %>%
    mutate(cumu_uptake = uptake)%>%
    ungroup()
}

# Combine both datasets
df_all <- bind_rows(
  process_uptake_file("zone_cumuptake_Nlimited.csv", "N limited"),
  process_uptake_file("zone_cumuptake_fullN.csv", "N sufficient")
)

# Add coverage and lateral flags
df_all <- df_all %>%
  mutate(
    treatment_type = case_when(
      str_detect(treatment_label, "Intercrop") ~ "Intercrop",
      TRUE ~ "Control"
    ),
    coverage = case_when(
      str_detect(treatment_label, "70%") ~ "70%",
      TRUE ~ "100%"
    ),
    lateral = if_else(str_detect(treatment_label, "no lat"), "No lateral", "Lateral"),
    lateral = factor(lateral, levels = c("Lateral", "No lateral"))  # force facet order
  )

# Color mapping (this is already defined and good to go)
treatment_colors <- c(
  "Control" = "skyblue",
  "Control - no lat" = "darkblue",
  "Intercrop" = "skyblue",
  "Intercrop - no lat" = "darkblue",
  "70% Control" = "green",
  "70% Control - no lat" = "forestgreen",
  "70% Intercrop" = "green",
  "70% Intercrop - no lat" = "forestgreen"
)

# Define key treatments
key_treatments <- c(
  "Control", "Control - no lat",
  "70% Control", "70% Control - no lat"
)

# Filter from df_all (both N levels) and calculate summary for ribbon
df_focus_summary <- df_all %>%
  filter(
    treatment_label %in% key_treatments,
    zone == "Cereal root zone"
  ) %>%
  mutate(cumu_uptake_mmol = cumu_uptake) %>%
  group_by(treatment_label, Nlevel, timestep) %>%
  summarise(
    mean_cumu_uptake = mean(cumu_uptake_mmol, na.rm = TRUE),
    sd_cumu_uptake = sd(cumu_uptake_mmol, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_cumu_uptake - sd_cumu_uptake,
    ymax = mean_cumu_uptake + sd_cumu_uptake
  )


# === Filter only N sufficient for both plots ===
df_all_sufficient <- df_all %>% filter(Nlevel == "N sufficient") # TOGGLE HERE!!!!!!!!!

# === Create linetype info ===
linetype_info <- df_all_sufficient %>%
  distinct(treatment_label, lateral)

# === FIGURE 1A: 100% vs 70% Control ===
key_treatments_control <- c("Control", "Control - no lat", "70% Control", "70% Control - no lat")

df_control_summary <- df_all_sufficient %>%
  filter(
    treatment_label %in% key_treatments_control,
    zone == "Cereal root zone"
  ) %>%
  mutate(cumu_uptake_mmol = cumu_uptake *14 /1000) %>%
  group_by(treatment_label, timestep) %>%
  summarise(
    mean_cumu_uptake = mean(cumu_uptake_mmol, na.rm = TRUE),
    sd_cumu_uptake = sd(cumu_uptake_mmol, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_cumu_uptake - sd_cumu_uptake,
    ymax = mean_cumu_uptake + sd_cumu_uptake
  ) %>%
  left_join(linetype_info, by = "treatment_label")

plot_control <- ggplot(df_control_summary, aes(
  x = timestep,
  y = mean_cumu_uptake,
  color = treatment_label,
  fill = treatment_label,
  linetype = lateral
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  labs(
    title = 'Monocrop',
    x = "Timestep (day)",
    y = "Cumulative N Uptake (mg)",
    color = "Treatment",
    fill = "Treatment",
    linetype = "Lateral Movement"
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  scale_linetype_manual(values = c("Lateral" = "solid", "No lateral" = "dashed")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = c(0.75, 0.25),  # (x, y) in [0,1] space of the plot area
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key.size = unit(0.5, "lines"),  # smaller legend keys
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# === FIGURE 1B: Intercrop ===
key_treatments_intercrop <- c("Intercrop", "Intercrop - no lat", "70% Intercrop", "70% Intercrop - no lat")

df_intercrop_summary <- df_all_sufficient %>%
  filter(
    treatment_label %in% key_treatments_intercrop,
    zone == "Cereal root zone"
  ) %>%
  mutate(cumu_uptake_mmol = cumu_uptake *14 /1000 ) %>%
  group_by(treatment_label, timestep) %>%
  summarise(
    mean_cumu_uptake = mean(cumu_uptake_mmol, na.rm = TRUE),
    sd_cumu_uptake = sd(cumu_uptake_mmol, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_cumu_uptake - sd_cumu_uptake,
    ymax = mean_cumu_uptake + sd_cumu_uptake
  ) %>%
  left_join(linetype_info, by = "treatment_label")

plot_intercrop <- ggplot(df_intercrop_summary, aes(
  x = timestep,
  y = mean_cumu_uptake,
  color = treatment_label,
  fill = treatment_label,
  linetype = lateral
)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  labs(
    title ='Intercrop',
    x = "Timestep (day)",
    y = "Cumulative N Uptake (mg)",
    color = "Treatment",
    fill = "Treatment",
    linetype = "Lateral Movement"
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  scale_linetype_manual(values = c("Lateral" = "solid", "No lateral" = "dashed")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = c(0.75, 0.25),  # (x, y) in [0,1] space of the plot area
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key.size = unit(0.5, "lines"),  # smaller legend keys
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )



final_line_plot <- plot_control + plot_intercrop +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(legend.position = "bottom")
  )

print(final_line_plot)