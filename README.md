# GAD酶基因组挖掘工具

## 简介
这是一个专门用于在全基因组中挖掘GAD酶（谷氨酸脱羧酶）的工具。考虑到GAD酶在不同物种中序列和结构高度不保守的特点，本工具采用了多种互补的挖掘策略，以提高识别准确率。

## 特点
- **增强保守模式识别**：扩展了23个保守模式，解决模式匹配得分为0的问题
- **多策略挖掘**：结合保守模式搜索、结构特征分析、同源性比较和功能域分析
- **灵活配置**：可调节不同特征的权重和阈值
- **无需外部依赖**：简化版只需要Python标准库
- **实时反馈**：提供详细的分析报告和候选序列列表
- **参考序列支持**：内置多种菌株的GAD酶参考序列

## 文件结构
```
gad_miner/
├── optimized_gad_mine.py       # 优化版GAD酶挖掘工具
├── gad_mine_simplified.py      # 简化版GAD酶挖掘工具
├── gad_config.json             # GAD酶特征配置文件
├── real_gad_sequences_example.fasta  # 真实GAD酶序列示例
├── download_real_gad.py        # 下载真实GAD酶序列的脚本
├── README.md                   # 使用说明
└── examples/                  # 示例数据
```

## 使用方法

### 1. 基本用法
```bash
python optimized_gad_mine.py -g protein_sequences.fasta -o results
```

### 2. 参数说明
```
-g, --genome     蛋白质序列文件（FASTA格式）
-o, --output     输出目录（默认：results）
-t, --threshold  筛选阈值（默认0.35）
```

### 3. 示例
```bash
# 使用示例数据
python optimized_gad_mine.py -g real_gad_sequences_example.fasta -o example_results

# 使用自定义阈值
python optimized_gad_mine.py -g my_proteins.fasta -o my_results -t 0.4
```

### 4. 结果文件
- **gad_candidates.tsv**：候选GAD酶列表（TSV格式）
- **gad_report.txt**：详细的分析报告

## 输出格式

### gad_candidates.tsv
```
Protein_ID    Length    Total_Score    Structural_Score    Motif_Score    Homology_Score    Best_Match    Motif_Count    Hydrophobic    Hydrophilic    Charged    Glycine    Proline    Cysteine
```

### gad_report.txt
包含：
- 分析概况
- Top候选列表
- 统计信息
- 参考匹配分布

## 算法原理

### 1. 保守模式搜索（v2版本增强）
基于GAD酶的已知保守区域：
- PLP结合位点（3个变体）：G[ST]G[KR][DE]，G[ST]G[KR]G，G[ST]G[KR]K
- 活性位点（3个变体）：[KR][ST][DE]X[KR]，K[ST]E[KR]，[KR]S[DE][KR]
- 谷氨酸底物结合位点（3个变体）：[FY]X[KR][DE]，Y[KR][DE]X，F[KR][DE]X
- GAD酶的典型保守序列：KK[DE]，KKD，KKE
- 功能域特征（基于Pfam PF00291）：[KR][DE]X[KR][DE]，[FY][KR]X[DE]，G[ST][KR]X[DE]
- 物种特异性保守模式：Ecoli-GAD-motif，Lactobacillus-GAD-motif，Human-GAD-motif
- 高级保守模式：[KR][DE]X[KR][DE]X[KR]，[FY][KR][DE]X[KR][DE]，G[ST][KR][DE]X[KR][DE]

### 2. 结构特征分析
基于GAD酶的物理特性：
- 典型长度范围：400-650氨基酸
- 氨基酸组成特征（疏水性、亲水性、带电氨基酸比例）
- 甘氨酸和脯氨酸含量分析

### 3. 同源性比较
与多种菌株的GAD酶参考序列比较：
- Escherichia coli GadA/B
- Lactobacillus brevis/plantarum Gad
- Bacillus subtilis Gad
- Listeria monocytogenes Gad
- Salmonella enterica Gad
- Human Gad1/Gad2 (GAD67/GAD65)
- Arabidopsis Gad1
- Mouse Gad67

### 4. 功能域分析
基于Pfam、SMART数据库的功能域模式：
- GAD-domain: [KR]{2}[DE]{2}[KR]{2}
- PLP-binding-domain: G[ST]G[KR][DE]X[KR][DE]
- active-site-domain: [KR][ST][DE]X[KR]X[KR][DE]
- glutamate-binding-domain: [FY]X[KR][DE]X[KR][DE]

### 5. 综合评分系统（v2版本优化）
综合所有方法的加权评分：
```
Total Score = Structural_score * 0.35 + Motif_score * 0.25 + Homology_score * 0.25 + Domain_score * 0.15
```

## 针对GAD酶不保守特征的优化（v2版本）

### 1. 扩展保守模式
从8个模式扩展到23个模式，提高模式匹配能力

### 2. 功能域分析
添加基于Pfam数据库的功能域模式匹配

### 3. 权重优化
调整评分权重：结构评分0.35，模式评分0.25，同源性评分0.25，功能域评分0.15

### 4. 核心模式识别
特别标记核心模式匹配（权重≥3），提高识别准确性

### 5. 额外评分机制
如果匹配到核心模式，额外增加分数，提高候选质量

### 6. 宽松的阈值
设置较低的相似度阈值（25%）以识别更多候选

### 7. 参考序列多样性
包含多种不同菌株的参考序列，提高识别能力

## 注意事项

### 1. 蛋白质序列预处理
输入应为蛋白质序列FASTA文件
如果输入是基因组文件，需要使用基因预测工具（如Prodigal）先提取蛋白质序列

### 2. 阈值设置
- 默认阈值0.35适合大多数情况
- 对于更严格的分析，可提高到0.4-0.45
- 对于更宽松的分析，可降低到0.3-0.35

### 3. 结果验证
所有候选序列都需要实验验证
建议结合酶活性测试、功能域分析等方法验证

## 扩展功能

### 1. 自定义参考序列
在 `optimized_gad_mine.py` 的 `reference_gad_sequences` 中添加更多参考序列

### 2. 模式权重调整
在 `gad_config.json` 中调整保守模式的权重

### 3. 添加机器学习分类器
可使用已知GAD和非GAD序列训练机器学习模型，进一步提高准确率

## 性能建议

### 1. 序列数量
适合处理500-5000个蛋白质序列
对于大规模分析（>10000序列），建议分批处理

### 2. 内存使用
每个序列的分析需要少量内存
大型数据集可能需要分批处理

## 开发计划

### 1. 增强版本
- 添加机器学习分类器
- 集成AlphaFold结构预测
- 添加功能域注释（InterProScan）

### 2. Web界面
- 开发Web应用界面
- 支持在线提交和分析

### 3. API接口
- 提供API接口供其他程序调用
- 支持批量处理

## 联系方式
如有问题或建议，请联系开发者。

## 许可证
MIT License