source('Utils.R')
load('data/mm_seob.rdata') #小鼠从受精卵的发育



#===> 提取表达矩阵和细胞信息
#看一眼数据发现木有标准化
#因为seob中count和data的数据一样
mm_seob <- NormalizeData(mm_seob)
gene_counts <- as.data.frame(mm_seob@assays$RNA@data)
gene_expr <- as.data.frame(mm_seob@assays$RNA@counts)
gene_counts[1:5,1:5] #标准化前
gene_expr[1:5,1:5] #标准化后

gene_counts <- mm_seob@assays$RNA@data
gene_expr <- mm_seob@assays$RNA@counts
cell_info <- rownames_to_column(mm_seob@meta.data,
                                var = 'cell_ids')



#===> 构建dyno对象
dynob <- wrap_expression(expression = t(as.matrix(gene_expr)), #稀疏或普通
                         counts = t(as.matrix(gene_counts)), #稀疏或普通
                         cell_info = cell_info)
class(dynob);table(cell_info$cell_type)

#===> 提供一些信息如果有的话
#添加一些已知的可知细胞
dynob <- add_prior_information(dynob,
                               start_id = c('zy','zy.1','zy.2','zy.3')
)
#添加分组
dynob <- add_grouping(dynob,
                      cell_info$cell_type)



#===> 轨迹推断
#选择软件
guidelines_shiny(dataset = dynob)

# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = FALSE, 
  expected_topology = NULL, 
  expect_cycles = FALSE, 
  expect_complex_tree = FALSE, 
  n_cells = 268, 
  n_features = 22431, 
  memory = "10GB", 
  prior_information = NULL,  
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)
guidelines$methods_selected # 看一眼候选的软件

#选择好方法使用docker运行，docker管理员运行
model_tmp <- infer_trajectories(dynob,
                                method = 'slingshot')

plot_dimred(model_tmp,
            label_milestones = F,
            grouping = dynob$grouping)





