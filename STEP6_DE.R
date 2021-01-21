source('Utils.R')
load('RData/cluster_anno.RData')



#===> 瞅一眼数据
DimPlot(seob,
        group.by = 'cell_type',
        label = T)



#===> 意义不大的差异表达分析
#不同细胞类型之间有哪些差异表达的基因=>细胞功能
Idents(seob) <- seob@meta.data$seurat_clusters
clu1_vs_clu2 <- FindMarkers(seob,
            ident.1 = 1,
            ident.2 = 2)
#pct.1,pct.2表示某个基因在某个clu中有百分之多少的细胞有表达
Idents(seob)


#一个clu与多个clu之间的差别 =>某个细胞特异表达的基因
clu1_vs_all <- FindMarkers(seob,
                           ident.1 = 1)


#===> 合格的new Marker基因鉴定
#Marker基因是指无论环境如何变化，都在特定细胞类型中表达的基因
#所以，它需要在不同的条件下都表达
#一般使用Seurat包的FindConservedMarkers来鉴定

Idents(seob) <- seob@meta.data$cell_type
nk_cell_markers <- FindConservedMarkers(
  seob,
  ident.1 = 'NK cells',  #指定找那个细胞类型的marker
  grouping.var = 'sample'  #指定要在样本内保守
)

DotPlot(seob, 
        features = c('KLRD1','GZMA','CTSW'))



#===> 合格的差异表达分析
#同一种细胞在不同condition之间差异表达的基因

#第一种方法提取某个细胞的seob对象如nk细胞
nk_cell_seob <- subset(seob,
                       subset = cell_type == 'NK cells')
Idents(nk_cell_seob) <- nk_cell_seob@meta.data$sample
nk_cells_stim_vs_ctrl_markers <- FindMarkers(nk_cell_seob,
            ident.1 = 'stim',
            ident.2 = 'ctrl') # 处理这个condition下nk细胞差异表达的基因

Idents(seob) <- seob@meta.data$cell_type
nk_vs_all <- FindMarkers(
  seob,
  ident.1 = 'NK cells')

rownames(nk_cells_stim_vs_ctrl_markers)[rownames(nk_cells_stim_vs_ctrl_markers) 
                                        %!in% rownames(nk_vs_all)]

FeaturePlot(seob,
            split.by = 'sample',
            features = 'SOS2')

#如果批量的话
Idents(seob) <- seob@meta.data$sample
seob_list <- SplitObject(seob, split.by = 'cell_type')
#接着对seob_list中每个细胞类型按照condition差异表达分析
all_cell_trim_vs_ctrl <- list()
for(i in names(seob_list)){
  tmp <- FindMarkers(seob_list[[i]],
                     ident.1 = 'stim',
                     ident.2 = 'ctrl')
  all_cell_trim_vs_ctrl[[i]] <- tmp 
}



DE_marker <- list(
  clu1_vs_clu2 = clu1_vs_clu2,
  clu1_vs_all = clu1_vs_all,
  nk_cell_markers = nk_cell_markers, 
  all_cell_trim_vs_ctrl_list = all_cell_trim_vs_ctrl
)
glimpse(DE_marker)
names(DE_marker)
