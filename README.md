# Bioinformatic-Analysis-Programs
this respository is made for automatic programs to run by python,these programs will save our time 
# VP2 Programs README

## 程序功能说明

### 1. check_excel.py
**功能**：读取Excel文件，验证列名顺序并打印数据内容
**使用方法**：`python check_excel.py`
**依赖**：pandas
**说明**：默认读取当前目录下的`test_levels_aim_levels_L1_high_L2_medium.xlsx`文件，打印列名顺序和数据内容

### 2. check_known_vp_genes.py
**功能**：检查已知的副溶血弧菌基因
**使用方法**：`python check_known_vp_genes.py`
**说明**：用于验证副溶血弧菌基因组中已知基因的存在情况

### 3. choose.py
**功能**：FNA文件提取工具（早期版本）
**使用方法**：`python choose.py`
**说明**：从NCBI数据中提取FNA文件的早期版本

### 4. choose2.py
**功能**：FNA文件增量提取工具（NCBI数据专用）
**使用方法**：`python choose2.py`
**依赖**：os, shutil, zipfile, hashlib
**说明**：从原始数据目录提取FNA文件到目标目录，支持增量提取和去重

### 5. choose3.py
**功能**：FNA文件提取工具（后续版本）
**使用方法**：`python choose3.py`
**说明**：choose系列工具的后续版本，用于FNA文件提取

### 6. complete_data_verification.py
**功能**：完整的数据验证
**使用方法**：`python complete_data_verification.py`
**说明**：对数据集进行全面验证，确保数据完整性和准确性

### 7. comprehensive_data_check.py
**功能**：综合数据检查
**使用方法**：`python comprehensive_data_check.py`
**说明**：对数据进行综合检查，包括格式、完整性和一致性

### 8. create_gbk_file.py
**功能**：创建GBK编码的测试文件
**使用方法**：`python create_gbk_file.py`
**说明**：生成包含中文内容的GBK编码测试文件

### 9. create_test_excel.py
**功能**：创建测试Excel文件
**使用方法**：`python create_test_excel.py`
**说明**：生成用于测试的Excel文件

### 10. create_test_levels.py
**功能**：创建测试等级文件
**使用方法**：`python create_test_levels.py`
**说明**：生成用于测试的等级数据文件

### 11. data_fix_tool.py
**功能**：数据修复工具
**使用方法**：`python data_fix_tool.py`
**说明**：修复数据中的错误和不一致

### 12. data_status_report.py
**功能**：数据状态报告
**使用方法**：`python data_status_report.py`
**说明**：生成数据状态报告，包括数据量、完整性等信息

### 13. experiment_starter_kit.py
**功能**：实验启动套件
**使用方法**：`python experiment_starter_kit.py`
**说明**：用于启动实验的工具套件

### 14. experiment_starter_kit_fixed.py
**功能**：实验启动套件（修复版）
**使用方法**：`python experiment_starter_kit_fixed.py`
**说明**：修复了bug的实验启动套件版本

### 15. extract_57strains_sequences.py
**功能**：提取57个菌株的序列
**使用方法**：`python extract_57strains_sequences.py`
**说明**：从57个菌株中提取序列信息

### 16. auto_complete_phenotype.py
**功能**：自动补充表型信息
**使用方法**：`python auto_complete_phenotype.py`
**依赖**：requests, csv
**说明**：从NCBI查询血清型并基于毒力基因划分危害等级，自动补充表型数据

### 17. analyze_with_kaptive.py
**功能**：使用Kaptive进行分析
**使用方法**：`python analyze_with_kaptive.py`
**说明**：利用Kaptive工具进行数据分析

### 18. analyze_virulence_genes.py
**功能**：分析毒力基因
**使用方法**：`python analyze_virulence_genes.py`
**说明**：对毒力基因进行分析和注释

### 19. analyze_real_data.py
**功能**：分析真实数据
**使用方法**：`python analyze_real_data.py`
**说明**：对真实数据集进行分析

### 20. analyze_57strains_targets.py
**功能**：基于57个真实GCA菌株的靶标分析
**使用方法**：`python analyze_57strains_targets.py`
**依赖**：pandas, numpy
**说明**：寻找毒力相关靶标，分析已知副溶血弧菌基因

### 21. 01_transfer.py
**功能**：数据传输脚本
**使用方法**：`python 01_transfer.py`
**说明**：用于数据传输的脚本

### 22. extract_sequences_fixed.py
**功能**：提取序列（修复版）
**使用方法**：`python extract_sequences_fixed.py`
**说明**：修复了bug的序列提取工具

### 23. extract_target_sequences.py
**功能**：提取推荐靶标的序列信息
**使用方法**：`python extract_target_sequences.py`
**依赖**：pandas
**说明**：从Roary数据中提取推荐靶标的详细信息

### 24. extract_target_sequences_fixed.py
**功能**：提取靶标序列（修复版）
**使用方法**：`python extract_target_sequences_fixed.py`
**说明**：修复了bug的靶标序列提取工具

### 25. fast_plot.py
**功能**：快速绘图
**使用方法**：`python fast_plot.py`
**说明**：用于快速生成图表的工具

### 26. filter_lines.py
**功能**：过滤文件中的行（基础版）
**使用方法**：`python filter_lines.py`
**说明**：过滤文件中包含特定字符串的行

### 27. filter_lines_improved.py
**功能**：过滤文件中的行（改进版）
**使用方法**：`python filter_lines_improved.py`
**说明**：改进版的行过滤工具，支持更多功能

### 28. filter_lines_latest.py
**功能**：过滤文件中的行（最新版）
**使用方法**：`python filter_lines_latest.py <input_file> <search_content> [options]`
**依赖**：pandas（用于Excel文件）
**说明**：支持多种文件类型（文本、CSV、Excel）的行过滤工具，支持替换功能

### 29. filter_by_field.py
**功能**：精确筛选指定字段的内容
**使用方法**：`python filter_by_field.py <input_file> <field_name> <field_value> [options]`
**依赖**：pandas, openpyxl
**说明**：支持CSV和Excel文件的字段级精确筛选，支持精确匹配、部分匹配、大小写敏感/不敏感、反向匹配等功能
**选项**：
  - `--output`: 输出文件路径
  - `--case-sensitive`: 大小写敏感匹配
  - `--encoding`: 文件编码
  - `--invert`: 反向匹配（显示不匹配的行）
  - `--contains`: 部分匹配（值包含搜索字符串）

### 30. final_recommended_targets.py
**功能**：最终推荐靶标
**使用方法**：`python final_recommended_targets.py`
**说明**：生成最终的推荐靶标列表

### 31. fix_data_consistency.py
**功能**：修复数据一致性
**使用方法**：`python fix_data_consistency.py`
**说明**：修复数据集中的一致性问题

### 32. intelligence.py
**功能**：智能分析工具
**使用方法**：`python intelligence.py`
**说明**：用于智能数据分析的工具

### 33. join_lines.py
**功能**：合并多行文本为单行
**使用方法**：`python join_lines.py [options] [quick_sep]`
**说明**：将多行文本合并为单行，支持自定义分隔符

### 34. organize_folders.py
**功能**：组织文件夹
**使用方法**：`python organize_folders.py`
**说明**：用于组织和管理文件夹结构

### 35. precise_target_mining.py
**功能**：精确靶标挖掘
**使用方法**：`python precise_target_mining.py`
**说明**：进行精确的靶标挖掘分析

### 36. primer_design_guide.py
**功能**：引物设计指南
**使用方法**：`python primer_design_guide.py`
**说明**：提供引物设计的指导和建议

### 37. prioritize_genes.py
**功能**：基因优先级排序
**使用方法**：`python prioritize_genes.py`
**说明**：对基因进行优先级排序

### 38. project_summary_report.py
**功能**：生成项目总结报告
**使用方法**：`python project_summary_report.py`
**依赖**：pandas, datetime
**说明**：生成副溶血弧菌靶标挖掘项目的总结报告

### 39. reassess_targets_real.py
**功能**：重新评估真实靶标
**使用方法**：`python reassess_targets_real.py`
**说明**：对真实靶标进行重新评估

### 40. same_go.py
**功能**：GO分析相关工具
**使用方法**：`python same_go.py`
**说明**：用于GO（基因本体）分析的工具

### 41. verify_t6ss_cluster.py
**功能**：验证T6SS簇
**使用方法**：`python verify_t6ss_cluster.py`
**说明**：验证T6SS（VI型分泌系统）簇的存在和完整性

### 42. virulence_analysis_interactive.py
**功能**：智能菌株危害等级评分+高毒基因筛选工具
**使用方法**：`python virulence_analysis_interactive.py`
**依赖**：pandas, openpyxl, matplotlib, seaborn
**说明**：交互式工具，用于菌株危害等级评分和高毒基因筛选

### 43. test_levels.py
**功能**：测试等级
**使用方法**：`python test_levels.py`
**说明**：用于测试等级划分的工具

## 功能相同的程序版本

### FNA文件提取工具
- choose.py（早期版本）
- choose2.py（当前版本）
- choose3.py（后续版本）

### 行过滤工具
- filter_lines.py（基础版）
- filter_lines_improved.py（改进版）
- filter_lines_latest.py（最新版）

### 序列提取工具
- extract_sequences_fixed.py（修复版）
- extract_target_sequences.py（原始版）
- extract_target_sequences_fixed.py（修复版）

### 实验启动套件
- experiment_starter_kit.py（原始版）
- experiment_starter_kit_fixed.py（修复版）

## 使用建议

1. 对于FNA文件提取，建议使用最新版本的`choose2.py`或`choose3.py`
2. 对于行过滤，建议使用最新版本的`filter_lines_latest.py`
3. 对于序列提取，建议使用修复版的`extract_target_sequences_fixed.py`
4. 对于实验启动，建议使用修复版的`experiment_starter_kit_fixed.py`
5. 对于智能分析，推荐使用`virulence_analysis_interactive.py`，它提供了交互式的分析体验

## 依赖安装

大多数程序需要安装以下依赖：
```bash
pip install pandas numpy matplotlib seaborn openpyxl requests
```

对于特定功能，可能需要额外安装其他依赖，请根据程序的具体需求进行安装。
