source('Utils.R')
rm(list=ls())
load('RData/reducted.RData')



#===> a little problem
#seob的meta.data的mt.percent有些问题,变成了一个list
#怀疑是标准化的过程导致的，想到重新生成一下
colnames(seob@meta.data)
seob@meta.data <- seob@meta.data[,-5]
head(seob@meta.data)



#===> K-NN聚类
# 找细胞最近的细胞
seob <- FindNeighbors(seob,
                      #k.param = 20, #最近的20个细胞
                      dims = 1:30)
# 聚类
seob <- FindClusters(seob,
                     resolution = 0.3, #分辨率值越大，cluster越多
                     random.seed = 1)
p1 <- DimPlot(seob,
              reduction = 'pca',
              group.by = 'seurat_clusters',
              label = T)
p2 <- DimPlot(seob,
              reduction = 'tsne',
              group.by = 'seurat_clusters',
              label = T)
p3 <- DimPlot(seob,
              reduction = 'umap',
              group.by = 'seurat_clusters',
              label = T)
p1 + p2 + p3 + plot_layout(guides = 'collect') &
  theme(legend.position = 'right')



#===> cell annotation
Cell_Marker <- read_csv("data/CellMarker.csv", 
                       col_types = cols(X3 = col_skip())) %>%
  separate_rows(Cell_Marker, sep = ',') %>%
  mutate(Cell_Marker = str_replace_all(Cell_Marker, ' ', '')) %>%
  distinct()

p1 <- DimPlot(seob,
              reduction = 'umap',
              group.by = 'seurat_clusters',
              label = T)

p2 <- DotPlot(seob,
        features = unique(Cell_Marker$Cell_Marker)) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   size = 8))
p1 + p2



#===> cluster2type
# 这里没有仔细划分直接就用的张老师的
# 可以按照这个一个个找如果一起看不太好看的话
DotPlot(seob,
        features = c('LYZ','CD14')) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   size = 8))
cluster2type <- c(
  "0"="CD14+ monocytes",
  "1"="CD8+ T cells",
  "2"="CD8+ T cells",
  "3"="CD8+ T / NK cells",
  "4"="NK cells",
  "5"="CD8+ T cells",
  "6"="B cells",
  "7"="FCGR3A+ monocytes",
  "8"="B cells",
  "9"="Conventional dendritic cells",
  "10"="Megakaryocytes", 
  "11"="CD8+ T cells",
  "12"="Plasmacytoid dendritic cells",
  "13"="Erythrocytes"
)

FeaturePlot(seob,
            reduction = 'umap',
            features = c('GNLY','NKG7','CD3D','CD8A'),
            #split.by = 'sample',
            label = T)
VlnPlot(object = seob,
        features = c('GNLY','NKG7'))


# cluster2type
seob@meta.data$cell_type <- unname(cluster2type[seob@meta.data$seurat_clusters])
meta.data <- seob@meta.data



# save
save(seob, file = 'RData/cluster_anno.RData')







