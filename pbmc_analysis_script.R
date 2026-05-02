# 加载包
library(Seurat)
library(SeuratData)
# 安装PBMC3K数据集
InstallData("pbmc3k")
# 加载pbmc3k数据到R环境中
data("pbmc3k")

# 更新Seurat对象到新版本格式
pbmc_updated <- UpdateSeuratObject(pbmc)

# 查看更新后的对象
print("更新后的对象：")
pbmc_updated

# 检查细胞类型分布
print("细胞类型分布：")
table(pbmc_updated$seurat_annotations)

# 如果上面报错，运行这个完整的分析流程：

# 1. 标准化数据
pbmc_updated <- NormalizeData(pbmc_updated)

# 2. 找到高变基因
pbmc_updated <- FindVariableFeatures(pbmc_updated)

# 3. 数据缩放
pbmc_updated <- ScaleData(pbmc_updated)

# 4. PCA降维
pbmc_updated <- RunPCA(pbmc_updated)

# 5. 细胞聚类
pbmc_updated <- FindNeighbors(pbmc_updated, dims = 1:10)
pbmc_updated <- FindClusters(pbmc_updated, resolution = 0.5)

# 6. UMAP可视化
pbmc_updated <- RunUMAP(pbmc_updated, dims = 1:10)

# 7. 绘制聚类图
DimPlot(pbmc_updated, reduction = "umap", label = TRUE)

# 8. 用已知的细胞类型注释来着色
DimPlot(pbmc_updated, reduction = "umap", group.by = "seurat_annotations", label = TRUE)


# 保存整个分析对象
saveRDS(pbmc_updated, file = "pbmc_analysis_results.rds")

# 以后要重新加载时，只需要：
# pbmc_updated <- readRDS("pbmc_analysis_results.rds")

# 提取细胞类型信息并保存
cell_annotations <- pbmc_updated$seurat_annotations
write.csv(cell_annotations, file = "pbmc_cell_type_annotations.csv")