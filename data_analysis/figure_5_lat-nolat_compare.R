# figure 5
# === Load libraries ===
library(tidyverse)
library(patchwork)

# === Load and merge data ===
df_fullN <- read_csv("output_fullN.csv") %>%
  mutate(N_level = "N sufficient")
df_Nlimited <- read_csv("output_Nlimited.csv") %>%
  mutate(N_level = "N limited")
df <- bind_rows(df_fullN, df_Nlimited)

# === Process data ===
df <- df %>%
  mutate(
    lateral_label = if_else(as.character(lateral_enabled) == "TRUE", "On", "Off"),
    m3_uptake_cereal = m3_uptake_cereal * 0.014,
    N_level = factor(N_level, levels = c("N sufficient", "N limited")),
    coverage_label = case_when(
      abs(coverage - 1.0) < 0.01 ~ "100% roots",
      TRUE ~ "other"
    ),
    cropping_label = if_else(intercrop == 1, "Intercrop", "Monocrop")
  ) %>%
  filter(coverage_label == "100% roots", cropping_label == "Monocrop")

# === Define panel function ===
make_box <- function(n_level_str, tag) {
  df %>% 
    filter(N_level == n_level_str) %>%
    ggplot(aes(x = lateral_label, y = m3_uptake_cereal, fill = lateral_label)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7,
      width = 0.6
    ) +
    scale_fill_manual(values = c("Off" = "#E69F00", "On" = "#56B4E9")) +
    labs(
      title = n_level_str,
      y = expression("Total N uptake (mg N/" * m^3 * " soil)"),
      x = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.text.x = element_blank(),
      legend.position = "none"
    )
}

# === Combine both panels ===
plot_A <- make_box("N sufficient", "A")
plot_B <- make_box("N limited", "B")

final_plot <- plot_A + plot_B + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(legend.position = "bottom")
  ) &
  labs(fill = "Lateral diffusion")

print(final_plot)
