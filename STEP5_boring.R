source('Utils.R')
load('RData/cluster_anno.RData')



#===> 可视化一波
library(ggrepel)
library(cowplot)
library(ggsci)

meta.data <- seob@meta.data
ggdata <- dplyr::select(meta.data, sample, cell_type) %>%
  group_by(sample, cell_type) %>%
  summarise(count = n()) %>%
  arrange(sample, desc(cell_type)) %>%
  mutate(cumsum = cumsum(count))

ggplot(ggdata, aes(x = sample, 
                   y = count, 
                   fill = cell_type)) +
  geom_col() +
  geom_text(data = dplyr::filter(ggdata, count >300),
            aes(x = sample,
                y = cumsum - 0.5*count,
                label = count)) +
  labs(x = '') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = 'Set3') +
  theme_cowplot()


