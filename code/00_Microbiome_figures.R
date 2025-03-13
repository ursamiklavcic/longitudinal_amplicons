# Microbiome FIGURES 

# Choose width based on publication type

width_in <- 170/25.4
height_in <- 225/25.4

# Create the ggplot
  geom_boxplot(size = 0.5) +  # Ensure line width > 0.25 pt
  theme_bw(base_size = 12) +  # Use legible font size
  theme(
    legend.position = "bottom", 
    legend.key.size = unit(10, "pt") # Adjust legend spacing
  ) +
  guides(fill = guide_legend(nrow = 1)) # Force legend into one row

# Save as high-quality PDF (vector format)
ggsave("figure.pdf", plot = p, width = width_in, height = height_in, device = cairo_pdf)

# Save as high-quality TIFF (raster format with LZW compression)
ggsave("figure.tiff", plot = p, width = width_in, height = height_in, dpi = 300, device = "tiff", compression = "lzw")