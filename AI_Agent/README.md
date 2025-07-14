# 🔬 Spatial Feature Interpretation Agent

基于GPT的CANVASX空间特征解释AI Agent，用于NSCLC免疫治疗分析。

## 📋 功能特点

- **智能解释**: 基于GPT-4的空间特征生物学解释
- **多线程处理**: 支持批量特征解释，高效并行处理
- **交互式界面**: 友好的命令行交互模式
- **数据驱动**: 基于262个CANVASX特征和10个生物学栖息地
- **结构化输出**: 5维度标准化解释格式

## 🚀 快速开始

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

### 2. 运行演示

系统已预配置API密钥和服务端点，可直接运行：

```bash
python run.py
```

## 💡 使用方法

### 交互式模式

```bash
python spatial_agent.py --interactive
```

支持命令：
- `help` - 显示帮助信息
- `info` - 显示特征统计
- `list` - 列出所有特征
- `list Composition` - 按类别列出特征
- 直接输入特征名进行解释

### 单个特征解释

```bash
python spatial_agent.py --feature "Ripley_K_mean_Habitat08"
```

### 批量处理

```bash
# 从文件读取特征列表
python spatial_agent.py --features_file features.txt --output results.json

# 支持输出格式: .json, .csv, .txt
```

## 📊 特征类别

系统支持262个空间特征，分为6大类：

1. **Composition (10个)** - 栖息地相对丰度
   - `frequency_Habitat01` ~ `frequency_Habitat10`

2. **Diversity (6个)** - 生态复杂性指标
   - `div_Richness`, `div_Shannon`, `div_Simpson`, 等

3. **Spatial metrics (~90个)** - 栖息地内空间组织
   - `Ripley_K_mean_`, `Ripley_L_mean_`, `G_mean_`, 等

4. **Interaction (100个)** - 栖息地间空间耦合
   - `cci_HabitatX_HabitatY`

5. **Distance (55个)** - 栖息地间空间分离度
   - `dis_HabitatX_HabitatY`

6. **Transition (1个)** - 基于熵的空间混合度
   - `SpatialTransitionEntropy`

## 🏠 生物学栖息地

- **H01**: Tumorogenic Core (肿瘤核心)
- **H02**: Macrophage Enriched (巨噬细胞富集)
- **H03**: B-cell Enriched (B细胞富集)
- **H04**: Fibrotic Activity Hub (纤维化中心)
- **H05**: Plasma cell Enriched (浆细胞富集)
- **H06**: Neutrophil Prominent (中性粒细胞)
- **H07**: Tumor Interface (肿瘤界面)
- **H08**: T-lymphonic Enriched (T淋巴细胞)
- **H09**: Pan-immune Active Zone (泛免疫活跃区)
- **H10**: Vasculature Niche (血管利基)

## 📤 输出格式

每个特征解释包含5个维度：

1. **Category** - 特征类别
2. **Cellular Composition** - 相关栖息地的细胞组成
3. **Spatial Property Description** - 空间特性描述
4. **Topological Coupling Tendency** - 拓扑耦合趋势
5. **Biological and Clinical Implication** - 生物学和临床意义

## 🔧 高级配置

```bash
python spatial_agent.py \
    --model gpt-4o \
    --num_workers 8 \
    --timeout 60 \
    --api_url "your-custom-api-url" \
    --interactive
```

## 📁 项目结构

```
YuchenAgent/
├── spatial_agent.py      # 主要Agent类
├── run.py                # 快速启动脚本
├── requirements.txt      # 依赖包
├── README.md            # 说明文档
└── data/                # 数据文件
    ├── Feature_matrix.xlsx
    ├── Feature_annotation.xlsx
    └── Habitat_annotation.docx
```

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个项目！

## �� 许可证

MIT License 