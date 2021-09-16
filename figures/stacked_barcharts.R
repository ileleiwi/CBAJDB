library(tidyverse)

#read in data
path <- "../data/"

ft <- read_tsv(paste0(path,"feature_table.tsv"))
metadata <- read_tsv(paste0(path,"metadata.tsv"))
taxonomy <- read_tsv(paste0(path,"taxonomy.tsv"))

#filter feature table to high responders (salm abund >= 25%) and control
#drop rows with unassigned taxonomy, create columns for each taxa level

#control mice
metadata_ctrl <- metadata %>%
  filter(treatment == "control")

ctrl_ft <- ft %>%
  select(asv_id, all_of(metadata_ctrl$Sample)) %>%
  left_join(taxonomy, by = c("asv_id" = "Feature ID")) %>%
  select(Taxon, everything(), -Confidence, -asv_id)

ctrl_ft <- ctrl_ft %>%
  separate(col = Taxon,
           into = c("domain",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species"),
           sep = ";") %>%
  mutate(phylum = case_when(is.na(class) ~ "Other",
                            T ~ class)) %>%
  select(-c(domain, phylum, order, family, genus, species)) %>%
  mutate(class = str_replace(class, "D_2__", "")) %>%
  pivot_longer(cols = -class,
               names_to = "sample",
               values_to = "relabund")  %>%
  group_by(class, sample) %>%
  summarise(mean_relabund = mean(relabund)) 

keep_class_ctrl <- ctrl_ft %>%
  na.omit() %>%
  group_by(class) %>%
  summarise(mean = mean(mean_relabund)) %>%
  arrange(desc(mean)) %>%
  pull(class)

keep_class_ctrl_plot <- c(keep_class_ctrl[1:8], "Gammaproteobacteria")

ctrl_ft <- ctrl_ft %>%
  mutate(class = case_when(class %in% keep_class_ctrl_plot ~ class,
                           T ~ "Other")) %>%
  group_by(class, sample) %>%
  summarise(mean_relabund = mean(mean_relabund)) %>%
  mutate(class = factor(class, levels = rev(c("Bacteroidia", "Clostridia", 
                                              "Mollicutes", "Verrucomicrobiae", 
                                              "Erysipelotrichia", "Bacilli", 
                                              "Halanaerobiia", 
                                              "Thermodesulfovibrionia", "Gammaproteobacteria",
                                              "Other"))))


#high responder mice
metadata_hr <- metadata %>%
  filter(treatment == "salmonella")

hr_ft <- ft %>%
  select(asv_id, all_of(metadata_hr$Sample)) %>%
  left_join(taxonomy, by = c("asv_id" = "Feature ID")) %>%
  select(Taxon, everything(), -Confidence, -asv_id)

hr_ft <- hr_ft %>%
  separate(col = Taxon,
           into = c("domain",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species"),
           sep = ";") %>%
  mutate(phylum = case_when(is.na(class) ~ "Other",
                            T ~ class)) %>%
  select(-c(domain, phylum, order, family, genus, species)) %>%
  mutate(class = str_replace(class, "D_2__", "")) %>%
  pivot_longer(cols = -class,
               names_to = "sample",
               values_to = "relabund")  %>%
  group_by(class, sample) %>%
  summarise(mean_relabund = mean(relabund))


keep_class_hr_plot <- rev(c("Gammaproteobacteria", "Bacilli", "Bacteroidia", 
                            "Clostridia", "Verrucomicrobiae", "Mollicutes", 
                            "Erysipelotrichia", "Oxyphotobacteria", "Halanaerobiia", 
                            "Other"))

hr_ft <- hr_ft %>%
  mutate(class = case_when(class %in% keep_class_hr_plot ~ class,
                           T ~ "Other")) %>%
  group_by(class, sample) %>%
  summarise(mean_relabund = mean(mean_relabund)) %>%
  mutate(class = factor(class, levels = keep_class_hr_plot))


#create palettes for plots
#https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

control_taxa <- rev(c("Bacteroidia", "Clostridia", 
                      "Mollicutes", "Verrucomicrobiae", 
                      "Erysipelotrichia", "Bacilli", 
                      "Halanaerobiia", 
                      "Thermodesulfovibrionia", "Gammaproteobacteria",
                      "Other"))

highresponder_taxa <- rev(c("Gammaproteobacteria", "Bacilli", "Bacteroidia", 
                            "Clostridia", "Verrucomicrobiae", "Mollicutes", 
                            "Erysipelotrichia", "Oxyphotobacteria", "Halanaerobiia", 
                            "Other"))

bar_chart_ctrl_pal <- tibble(ctrl = control_taxa,
                             hex = pal[1:10])

bar_chart_hr_pal <- tibble(hr = highresponder_taxa,
                           hex = pal[1:10])

bar_chart_ctrl_pal_reorder <- bar_chart_ctrl_pal %>%
  mutate(ctrl = case_when(ctrl == "Thermodesulfovibrionia" ~ "Oxyphotobacteria",
                          T ~ ctrl),
         ctrl = factor(ctrl, levels = bar_chart_hr_pal$hr)) %>%
  arrange(ctrl) %>%
  mutate(ctrl = as.character(ctrl),
         ctrl = case_when(ctrl == "Oxyphotobacteria" ~ "Thermodesulfovibrionia",
                          T ~ ctrl)) %>%
  pull(hex)

hr_pal <- bar_chart_ctrl_pal_reorder
ctrl_pal <- bar_chart_ctrl_pal$hex
ctrl_pal[3] <- pal[11]

#plot
ctrl_ft %>%
  ggplot(aes(x = sample, y = mean_relabund, fill = class)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(values = ctrl_pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

hr_ft %>%
  ggplot(aes(x = sample, y = mean_relabund, fill = class)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(values = hr_pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
