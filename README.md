# 基于 Seurat 的单细胞 RNA 测序 (scRNA-seq) 数据降维与聚类

本项目使用 R 语言与 Seurat 框架，完成了经典的 PBMC 单细胞测序数据的标准分析流水线。

## 📊 分析流程
1. **数据标准化与特征提取**：完成数据的 `NormalizeData` 与高变基因识别。
2. **线性降维**：应用 PCA (Principal Component Analysis) 提取主成分。
3. **非线性流形学习**：应用 **UMAP** 算法将高维稀疏数据映射至二维空间。
4. **细胞聚类与可视化**：成功聚类出包括 Naive CD4 T, B cells, CD14+ Mono 等 8 种核心免疫细胞亚型。

## 🛠 技术栈
* R Language, Seurat, ggplot2
