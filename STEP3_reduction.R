source('Utils.R')
load('RData/normalized2.RData')



#===> PCA分析 基于HVG
seob <- RunPCA(seob)
#看一眼保留多少PC做后续分析，每一个PC都描述细胞之间的差异
#第一种标准化的方法保留前20个主成分就行
#第二种标准化的方法最好多保留点，SCT官网说比较敏感
#每一个主成分都描述对细胞之前差异
ElbowPlot(seob, ndims = 50)
DimPlot(seob,
        reduction = 'pca',
        group.by = 'sample')
DimPlot(seob,
        reduction = 'pca',
        group.by = 'Phase')



#===> T-SNE分析
#不是基于表达举证而是基于PCA
#cluster之间的距离是没有意义的
seob <- RunTSNE(seob,
                dims = 1:30)
DimPlot(seob,
        reduction = 'tsne',
        group.by = 'sample')



#===> UMAP分析
#速度很快
seob <- RunUMAP(seob, 
                dims = 1:30)
DimPlot(seob,
        reduction = 'umap',
        group.by = 'sample')


#===> 可视化三种降维的方法
#主要看样本之间的差异，存不存在批次效应，细胞周期的差异等等
#因为我们主要想看的是细胞类型的差异
#所以需要排除其他差异的影响结果就是最大的差异就是细胞类型的差异
#聚类分析就是将相同的细胞类型进行聚类
p1 <- DimPlot(seob,
              reduction = 'pca',
              group.by = 'sample')
p2 <- DimPlot(seob,
              reduction = 'tsne',
              group.by = 'sample')
p3 <- DimPlot(seob,
              reduction = 'umap',
              group.by = 'sample')
p1 + p2 + p3 + plot_layout(guides = 'collect') &
  theme(legend.position = 'top')


#排除了其他的批次效应因素后的seob
save(seob, file = 'RData/reducted.RData')


