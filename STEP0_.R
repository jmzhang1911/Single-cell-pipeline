source('Utils.R')
library(Seurat)


#===> 一般情况下的cellranger结果
scrna_data1 <- Read10X(data.dir = 'data/GSM3972018/')
#构建 Seurat对象
seob1 <- CreateSeuratObject(
  counts = scrna_data1,
  min.cells = 3, #去除小于3个细胞中表达的基因
  min.features = 200) #小于小于200个基因表达的细胞
class(seob1)
str(seob1)

 

#===> 从h5文件构建
scrna_data2 <- Read10X_h5(
  filename = 'data/GSM3489182/GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5'
)
seob2 <- CreateSeuratObject(counts = scrna_data2,
                            min.cells = 3,
                            min.features = 200)



#===> 从表达举证构建
scrna_data3 <- read.table(
  'data/GSM2829942/GSM2829942_HE6W_LA.TPM.txt',
  row.names = 1,
  header = T
)
seob3 <- CreateSeuratObject(counts = scrna_data3,
                            min.cells = 3,
                            min.features = 200)


