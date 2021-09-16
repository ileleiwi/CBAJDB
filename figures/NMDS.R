library(tidyverse)
library(vegan)
library(plyr)
library(ggtext)

#read in data
path <- "../data/"

metaG_samples <- c("13-(11)-S-F-KE", "20-(11)-S-F-KL", "27-(11)-S-F-KS", 
                   "49-(11)-C-F-LO", "47-(11)-C-F-LM", "50-(11)-C-F-LP")

ft <- read_tsv(paste0(path,"feature_table.tsv")) %>%
  column_to_rownames(var = "asv_id") %>%
  as.matrix() %>%
  t()


metadata <- read_tsv(paste0(path,"metadata.tsv")) %>%
  mutate(metaG_factor = case_when(Sample %in% metaG_samples ~ 6, 
                                  TRUE ~ 1)) 

ft <- ft[metadata$Sample,]

# #NMDS objects
ord_bray<-metaMDS(ft, distance = "bray", autotransform = FALSE, noshare = 0.1, trace = 1)

# #getting scores for NMDS object and saving as df
ord_bray_scrs<-as.data.frame(scores(ord_bray),display="sites")

#perform Shepards goodness test on NMDS and anosim on control vs salmonella
Shepards_goodness_test_results <- goodness(ord_bray)

anosim_results <- anosim(ft, 
                         metadata$treatment, 
                         permutations = 999, 
                         distance = "bray", 
                         strata = NULL,
                         parallel = getOption("mc.cores"))
stress_plot <- stressplot(ord_bray)

#mrpp
#getting grouping information
grouping_variable <- metadata$treatment

mrpp_out <- mrpp(ft,
                 grouping = grouping_variable,
                 distance = "bray")

#make plot dataframe
nmds_plot_df <- ord_bray_scrs %>%
  rownames_to_column(var = "Sample") %>%
  left_join(metadata)

#create hulls for treatment polygons
find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls <- ddply(nmds_plot_df, "treatment", find_hull)

# Bray Curtis NMDS
bray_nmds <- nmds_plot_df %>%
  ggplot(aes(x = NMDS1, 
             y = NMDS2, 
             fill = treatment, 
             size = metaG_factor)) +
  geom_point(aes(fill = treatment),
             shape = 23) +
  geom_polygon(data = hulls,
               aes(x = NMDS1, 
                   y = NMDS2),
               alpha = 0.3,
               size = .5) +
  geom_textbox(data = data.frame(NMDS1 = .75,
                                 NMDS2 = -.6,
                                 treatment = NA,
                                 metaG_factor = 1),
               aes(label = "Anosim
               \nR: 0.9684
                   \nSignificance: 0.001"),
               box.color = "black",
               box.size = 0.5,
               fill = "white",
               size = 4) +
  theme_classic() +
  scale_alpha_identity() +
  scale_size_identity(guide = "legend") +
  labs(title = "Bray-Curtis (Day 11)",
       fill="Treatment") 

bray_nmds
