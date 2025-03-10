
# confirm ggpubr is loaded
library(ggpubr)

# define parameters for saving figures with ggpubr::ggsave()
# set for ecography
#Width: 945 (single column), 1476 (1.5 column) or 1961 (double column) pixels (at 300 dpi).
my_ggsave <- function(filename, plot, 
                      device = 'tiff',
                      path = './Figs/', 
                      dpi = 300, 
                      units = 'px', 
                      width = 945, 
                      height = 945, 
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

