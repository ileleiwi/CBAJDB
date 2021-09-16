library(tidyverse)

#read in data
path <- "../data/"

ft <- read_tsv(paste0(path,"feature_table.tsv"))

#filter to metaG mice
metaG_samples <- c("13-(11)-S-F-KE", "20-(11)-S-F-KL", "27-(11)-S-F-KS", 
                   "49-(11)-C-F-LO", "47-(11)-C-F-LM", "50-(11)-C-F-LP")

metaG_ft <- ft %>%
  select(asv_id, all_of(metaG_samples)) 

#healthy feature table
healthy_ft <- metaG_ft %>%
  select(asv_id, 5:7) 

healthy_mean <- rowMeans(healthy_ft[2:4])

healthy_ft_mean <- tibble(asv_id = healthy_ft$asv_id,
                          mean_relabund = healthy_mean) %>%
  left_join(taxonomy, by = c("asv_id" = "Feature ID")) %>%
  select(asv_id, Taxon, mean_relabund)

healthy_ft_mean_sorted <- healthy_ft_mean %>%
  mutate(asv_id = factor(asv_id)) %>%
  arrange(desc(mean_relabund)) %>%
  pull(asv_id)

healthy_ft_mean_plot <- healthy_ft_mean %>%
  mutate(asv_id = factor(asv_id, levels = healthy_ft_mean_sorted)) %>%
  arrange(desc(mean_relabund)) %>%
  mutate(number = seq(1,nrow(healthy_ft_mean)))

#infected feature table
infected_ft <- metaG_ft %>%
  select(asv_id, 2:4) 

infected_mean <- rowMeans(infected_ft[2:4])

infected_ft_mean <- tibble(asv_id = infected_ft$asv_id,
                           mean_relabund = infected_mean) %>%
  left_join(taxonomy, by = c("asv_id" = "Feature ID")) %>%
  select(asv_id, Taxon, mean_relabund)

infected_ft_mean_plot <- infected_ft_mean %>%
  mutate(number = seq(1,nrow(infected_ft_mean)),
         asv_id = factor(asv_id, levels = healthy_ft_mean_sorted)) %>%
  arrange(desc(mean_relabund)) 

taxa_levels <- taxonomy %>%
  mutate(`Feature ID` = factor(`Feature ID`, levels = healthy_ft_mean_sorted)) %>%
  arrange(desc(`Feature ID`))

#plots
healthy_ft_mean_plot  %>%
  slice_head(n = 50) %>%
  ggplot(aes(x = asv_id, y = mean_relabund)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ylim(0,0.6) +
  labs(y = "Relative Abundance (%)",
       x = "ASVs")

infected_ft_mean_plot  %>%
  slice_head(n = 50) %>%
  ggplot(aes(x = asv_id, y = mean_relabund)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(y = "Relative Abundance (%)",
       x = "ASVs") +
  geom_text(aes(label = asv_id, angle = 90),position=position_dodge(width=0.1),vjust=1, hjust=-.1, size =3) 
