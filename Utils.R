options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")



library(tidyverse)
library(patchwork)
library(Seurat)
library(dyno)


`%!in%` = Negate(`%in%`)


#===> install packages
if(F){
  BiocManager::install('hdf5r')
  BiocManager::install('Seurat')
  BiocManager::install('metap') 
}

if(F){
  devtools::install("./software/dynutils-master")
  devtools::install("./software/dynparam-master")
  devtools::install("./software/babelwhale-master")
  devtools::install("./software/dynnormaliser-master")
  devtools::install("./software/dynwrap-master")
  devtools::install("./software/dynplot-master")
  devtools::install("./software/dynmethods-master")
  devtools::install("./software/dynguidelines-master")
  devtools::install("./software/dyno-master")
}

