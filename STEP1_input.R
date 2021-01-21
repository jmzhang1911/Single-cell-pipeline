source('Utils.R')



#===> input data，构建seurat对象
if(T){
  filepath <- list.dirs('data/GSE96583/', recursive = F)
  seob_list <- list()
  for(i in filepath){
    filen <- basename(i)
    scrna_data <- Read10X(data.dir = i)
    seob <- CreateSeuratObject(counts = scrna_data,
                               min.cells = 3,
                               min.features = 200)
    seob@meta.data$sample = filen
    #seob[['sample']] = sample
    seob_list[[filen]] = seob
  }
  seob <- merge(x = seob_list[[1]],
                y = seob_list[-1],
                add.cell.ids = names(seob_list)) #行名加前缀
  save(seob, file = 'RData/GSE96583_input_seob.RData')
  table(seob@meta.data$sample)
  }



load('RData/GSE96583_input_seob.RData')
#===> 添加线粒体的信息
rownames(seob) #所有gene id
str_subset(rownames(seob), '^MT-') # mt gene id
seob@meta.data$percent.mt <- PercentageFeatureSet(
  seob,
  pattern = '^MT-'
  #features = mt_genelist
)
class(seob)



#===> 细胞周期打分
seob <- CellCycleScoring(
  seob,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)



#===> 简单可视化
# 小提琴箱线图
VlnPlot(seob,
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),
        group.by = 'sample',
        #log = T,
        pt.size = 0.05)
# 山峦图
RidgePlot(object = seob,
          features = c('nFeature_RNA','nCount_RNA','percent.mt'),
          #log = T,
          ncol = 1,
          group.by = 'sample')
# 散点图
p1 <- FeatureScatter(seob,
                     feature1 = 'nCount_RNA',
                     feature2 = 'nFeature_RNA',
                     group.by = 'sample')
p2 <- FeatureScatter(seob,
                     feature1 = 'nCount_RNA',
                     feature2 = 'percent.mt',
                     group.by = 'sample')
p1 + p2



#===> QC
seob <- subset(seob,
               subset = nFeature_RNA > 200 &
                 nFeature_RNA < 2500 &
                 percent.mt < 10)
save(seob, file = 'RData/GSE96583_input_seob_filtered.RData')



#===> 替代图
ggdata <- as.data.frame(seob@meta.data) %>% 
  rownames_to_column(var = 'cell_id') 
ggdata$percent.mt <- as.numeric(ggdata$percent.mt$nCount_RNA)
ggdata <- as.data.frame(ggdata)  
ggdata2 <- ggdata %>%  
  dplyr::tbl_df() %>%
  dplyr::select(cell_id, sample, nCount_RNA, nFeature_RNA ,percent.mt) %>%
  gather(key = 'condition', value = 'count', 3:5)

ggplot(ggdata2, aes(x = count, y = sample)) +
  geom_density_ridges(aes(fill = sample)) +
  facet_wrap(~condition,
             scales = "free_x",nrow = 3) +
  scale_fill_nejm() +
  labs(x = '', y = '') +
  theme_minimal_hgrid()

ggdata <- as.data.frame(seob@meta.data) %>% 
  rownames_to_column(var = 'cell_id') 
ggdata$percent.mt <- as.numeric(ggdata$percent.mt$nCount_RNA)
ggdata <- as.data.frame(ggdata)  
ggdata2 <- ggdata %>%  
  dplyr::tbl_df() %>%
  dplyr::select(cell_id, sample, nCount_RNA, nFeature_RNA ,percent.mt)

model.lm<-lm(formula = nCount_RNA ~ nFeature_RNA, data = ggdata2)
summary(model.lm)
p1 <- ggplot(ggdata2, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(aes(color = sample)) +
  geom_smooth(aes(group = 2), method = 'lm') +
  annotate("text", x = 12300, y = 500, label = "r = 0.8622;pvalue = 2.2e-16") +
  scale_color_nejm() +
  theme_classic()


model.lm<-lm(formula = nCount_RNA ~ percent.mt, data = ggdata2)
summary(model.lm)
p2 <- ggplot(ggdata2, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(aes(color = sample)) +
  geom_smooth(aes(group = 2), method = 'lm') +
  scale_y_continuous(limits = c(0,15)) +
  annotate("text", x = 12300, y = 13, label = "r = 0.1154;pvalue = 2.2e-16") +
  scale_color_nejm() +
  theme_classic()

p1 + p2



