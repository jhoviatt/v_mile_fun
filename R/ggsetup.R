
# confirm ggpubr is loaded
library(ggpubr)

# define parameters for saving figures with ggpubr::ggsave()
# MEPS:
# Type area = 169 x 225 mm, column width = 81 mm, gutter width: 7 mm
# single panel figure (with less text) → 81 mm width
# single panel figure (with more text or insets) → 105 mm width
# two panel figure (horizontally aligned) → 169 mm width
# two panel figure (vertically aligned) → 81 mm width
# multiple panel figure → 169 mm width
my_ggsave <- function(filename, plot, 
                      device = 'tiff',
                      path = './Figs/', 
                      dpi = 300, 
                      units = 'mm', 
                      width = 81, 
                      height = 81, 
                      bg = 'white')
{
    ggsave(filename, plot, device = device, path = path, dpi = dpi, units = units, width = width, height = height, bg = bg)
}

# define text size constants
small_text_size <- 6
ax_text_size <- 8
lab_text_size <- 10

# set figure theme
theme_set(theme_classic2(base_size = lab_text_size) +
              theme(text = element_text(size = lab_text_size),
                    axis.text = element_text(size = ax_text_size))) 

