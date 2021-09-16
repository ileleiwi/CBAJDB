library(tidyverse)

#read in data
path <- "../data/"

metadata <- read_tsv(paste0(path,"metadata.tsv")) %>%
  mutate(salmonella_relative_abundance = round(salmonella_relative_abundance, 4))

set.seed(13)
#all mice in NMDS (n=30)
nmds_mice_salm_abund <- metadata
#nmds mice
nmds_mice_salm_abund %>%
  ggplot(aes(x = treatment, y = salmonella_relative_abundance)) +
  geom_point(shape = 21, size = 3, position = position_jitter(width = .1, seed = 15), fill = "black", alpha = .5) +
  theme_bw() +
  ylim(-.001,1)

#mice with lipocalin (n=12)
lipo_mice_salm_abund <- metadata %>%
  filter(mouse_id %in% c("3_44", "3_45", "3_46", "3_47", "3_48", "3_49", "3_50",
                         "3_13", "3_15", "3_16", "3_20", "3_27", '3_28')) 

lipo_mice_salm_abund %>%
  ggplot(aes(x = treatment, y = salmonella_relative_abundance)) +
  geom_point(shape = 21, size = 3, position = position_jitter(width = .1), fill = "black", alpha = .5) +
  theme_bw() +
  ylim(-.01,1)

