rm(list=ls())

library(hexSticker)
library(magick)


hexSticker::sticker("data-raw/GGif3jHY.png",
                    package = "Virusparies",       # Adds the title
                    p_color = "#679267",             # Title color
                    h_size = 1.5,
                    h_color = "#679267",
                    h_fill = "white",
                    p_size = 18,                    # Title size for visibility
                    p_y = 0.5,                     # Moves the title below the image
                    p_family = "Roboto-Bold",      # Bold font (ensure it's installed)
                    s_x = 1,                       # Center the image horizontally
                    s_y = 1.2,                       # Keeps the image centered
                    s_width = 0.8,                 # Shrink the image width
                    s_height = 0.8,                # Shrink the image height
                    url = "https://github.com/SergejRuff/Virusparies", # Add URL
                    u_size = 2,                    # Adjust the URL size
                    u_y = 0.05,                    # Position the URL near the bottom
                    u_color = "black",             # URL color
                    filename = "tools/virusparies.png")


jane <- magick::image_read("tools/virusparies.png")
image_scale(jane, "150") %>%
  image_write(path = "tools/logo.png", format = "png")
