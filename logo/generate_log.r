# Ubuntu dependencies:
apt-get install --fix-missing \
  cargo \
  libudunits2-dev \
  libgdal-dev \
  gdal-bin

# R dependencies:
# install.packages("tidyverse")
# install.packages("gganimate")
# install.packages("gifski")
# install.packages("png")
install.packages("transformr")
devtools::install_github("thomasp85/transformr")


library(tidyverse)
library(gganimate)
library(transformer)

mtcars %>%
    as_tibble() %>%
    select(cyl, qsec, gear) %>%
    distinct %>%
    ggplot(aes(cyl, qsec)) + 
        geom_line() +
        transition_states(
            gear
        ) +
        ease_aes("linear")


