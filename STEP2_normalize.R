source('Utils.R')
load('RData/GSE96583_input_seob_filtered.RData')



seob_list <- SplitObject(seob, split.by = 'sample')
#===> 标准化 第一种方法
#去除测序深度对基因表达的影响，相当于转录组的CPM
# STEP1:将seob对象拆分成列标，因为需要对每个样本先进行标准化

if(T){
  # 第一次标准化
  for(i in 1:length(seob_list)){
    seob <- seob_list[[i]] #也可单括号，感觉如果按照字符取的话需要双括号
    # 可以认为是组内标准化
    seob <- NormalizeData(seob,
                          normalization.method = 'LogNormalize')
    # 筛选表达高度变化的基因, 1000-2000足够
    seob <- FindVariableFeatures(seob,
                                 selection.method = 'vst',
                                 nfeatures = 1000) 
    seob_list[[i]] <- seob
  } 
  # 每组的数据整合
  anchors <- FindIntegrationAnchors(object.list = seob_list,
                                    normalization.method = 'LogNormalize')
  seob <- IntegrateData(anchorset = anchors)
  DefaultAssay(seob) <- 'integrated'
  # 第二次标准化
  seob <- ScaleData(seob, features = rownames(seob))
  save(seob, file = 'RData/normalized1.RData')
}



#===> 标准化第二种方法 SCTransform
if(T){
  for(i in 1:length(seob_list)){
    seob_list[[i]] <- SCTransform(seob_list[[i]],
                                  variable.features.n = 3000,
                                  #vars.to.regress = c('percent.mt','S.Score','G2M.Score'),
                                  verbose = FALSE)
  }
  # 选择用于整合的基因
  features <- SelectIntegrationFeatures(object.list = seob_list,
                                        nfeatures = 3000)
  # 准备整合
  seob_list <- PrepSCTIntegration(object.list = seob_list,
                                  anchor.features = features)
  # 找anchors
  anchors <- FindIntegrationAnchors(object.list = seob_list,
                                    normalization.method = 'SCT',
                                    anchor.features = features)
  # 整合数据
  seob <- IntegrateData(anchorset = anchors,
                        normalization.method = 'SCT')
  DefaultAssay(seob) <- 'integrated'
  save(seob, file = 'RData/normalized2.RData')
}















