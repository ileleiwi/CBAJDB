library(tidyverse)

#read in data
path <- "../data/"

metaG_samples <- c("13-(11)-S-F-KE", "20-(11)-S-F-KL", "27-(11)-S-F-KS", 
                   "49-(11)-C-F-LO", "47-(11)-C-F-LM", "50-(11)-C-F-LP") 

metadata <- read_tsv(paste0(path,"metadata.tsv")) %>%
  mutate(metaG_factor = case_when(Sample %in% metaG_samples ~ 6, 
                                  TRUE ~ 1)) 

#plot
set.seed(2)
metadata %>%
  ggplot(aes(x = treatment, y = log(ng_Lc_per_g_feces))) +
  geom_boxplot(aes(fill = treatment), outlier.shape = NA) +
  geom_jitter(aes(size = metaG_factor,
                  fill = treatment),
              shape = 21) +
  theme_classic()
