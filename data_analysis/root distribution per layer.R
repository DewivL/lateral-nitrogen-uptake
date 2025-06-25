library(ggplot2)
library(viridis)
library(readr)
library(dplyr)

inferno_trimmed <- (inferno(256)[30:220])

# Read the root pattern file
roots <- read_csv("root_pattern_70.csv")

# calculate total DW per layer
layer_totals <- roots %>%
  group_by(z) %>%
  summarise(layer_DW = sum(DW, na.rm = TRUE)) %>%
  ungroup() # Ungroup after summarizing for easier joining

# join layer_DW back to the roots dataframe
roots <- roots %>%
  left_join(layer_totals, by = "z")

#  Define color_value:
#    If DW > 0 (there is a root), color_value is the total layer_DW.
#    Otherwise (no root), it's NA_real_ (for white/transparent).
roots <- roots %>%
  mutate(color_value = if_else(DW > 0, layer_DW, NA_real_))

# Get the min and max of layer_DW for the color scale limits
min_layer_dw <- min(layer_totals$layer_DW, na.rm = TRUE)
max_layer_dw <- max(layer_totals$layer_DW, na.rm = TRUE)

# Plot
ggplot(roots, aes(x = x, y = y, fill = color_value)) +
  geom_raster() + 
  facet_wrap(~ z, ncol = 5) + 
  scale_fill_gradientn(
    colours = inferno_trimmed,
    na.value = "white",
    # Set limits based on the actual range of layer_DW
    limits = c(min_layer_dw, max_layer_dw),
    name = "Layer Dry Weight (g)" 
  ) +
  coord_fixed() +
  labs(
    title = "Root distribution per layer"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "bold"),
    strip.text = element_text(color = "black", face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_blank(),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )