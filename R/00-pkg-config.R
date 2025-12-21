#--------------------------------------------
#
# mordor-spillover/R/00-pkg-config.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Packages config file
#
#------------------------------------


# Worker packages ---------------------------------------------------------

library(tidyverse)
library(here)
library(broom)
library(parallel)
library(units)
library(tictoc)
library(renv)
library(osfr)

# Visuals -----------------------------------------------------------------

library(knitr)
library(kableExtra)
library(quarto)

library(patchwork)
library(ggpubr)
require(ggtext)
library(ggspatial)

library(knitr)
library(kableExtra)

# Spatial packages --------------------------------------------------------

library(sf)
sf_use_s2(FALSE)

library(raster)
library(lwgeom)