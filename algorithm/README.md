# 函数型时间序列断点检测与 Bootstrap 推断

本目录包含用于函数型时间序列（functional time series）断点检测与 Bootstrap 推断的实现、示例和说明。

本 README 致力于说明问题背景、所采用的方法、代码结构、构建与运行步骤，以及输入/输出格式和常见参数说明。

## 问题描述

目标是检测时间序列中结构性断点（regime change）的位置，并对检测到的候选断点进行显著性检验。数据为函数型观测 X_t(s)，t 为时间，s 为函数域上的网格点。我们需要：

- 在观测序列中定位一个或多个断点位置；
- 对每个候选断点计算检验统计量并通过 Moving Block Bootstrap（移动块自助法）估计 p 值或临界值；
- 将结果保存为易解析的文本文件，以便后续可视化和分析。

## 方法概览

主要步骤：

1. Binary Segmentation（二分分割）用于断点检测
   - 使用 Segmented Sum of Generalized Residuals (SSGR) 作为检验统计量，在每个候选分割点计算分段后的残差平方和减少量（或类似的度量）。
   - 采用逐步二分分割 (binary segmentation) 在全样本及子段中递归寻找显著断点，直到满足停机准则（如 BIC 或最小段长度）。

2. Moving Block Bootstrap（移动块自助法）用于推断
   - 对于每个候选断点，固定位置计算原始检验统计量；
   - 使用长度为 `block_length` 的移动块对时间序列重采样（保留时间依赖结构），重复 `B` 次构建检验统计量的经验分布；
   - 计算 p 值或临界值（例如 α=0.05）来判断显著性。

这些方法结合既能定位断点也能对断点显著性给出合理的推断。

## 代码结构

- `C/`：C 语言实现（高性能实现）
  - `breakpoint_detection.c`：二分分割断点检测主程序，编译生成 `break_detect`（或 `break_detect.exe`）。
  - `bootstrap.c`：移动块 bootstrap 实现，计算每个断点的 p 值，编译生成 `bootstrap`（或 `bootstrap.exe`）。
  - `utils.c` / `utils.h`：共享的数据结构与辅助函数（读写二进制数据、计算段均值、SSGR 等）。

- `R/`：R 脚本用于数据生成、可视化和管道控制
  - `01_data_generation.R`：用于生成带断点的函数型时间序列模拟数据（示例 DGP）。
  - `02_visualization.R`：断点检测结果与 bootstrap 结果的可视化（热图、p 值柱状图等）。
  - `03_main_analysis.R`：高层调用脚本，保存数据、调用 C 程序、读取并汇总结果。

- `run_analysis.R`：顶层脚本，执行 `03_main_analysis.R` 的完整流程并生成图形输出。

- `data/`, `results/`：运行时的输入/输出文件位置（如 `data/temp_data.bin`, `data/temp_meta.txt`, `results/breakpoints.txt`）。

## 构建与运行

构建 C 程序（在 `algorithm/code` 下）：

```powershell
cd algorithm\code
make all      # 或 make break_detect; make bootstrap
```

快速运行完整分析（在 workspace 根目录）：

```powershell
cd <workspace_root>   # 如 c:\Users\wei\Desktop\mwe
Rscript algorithm/run_analysis.R
```

`run_analysis.R` 将执行：

- 生成模拟数据并保存到 `algorithm/code/data/`；
- 运行 `break_detect`，将检测结果保存到 `algorithm/code/results/breakpoints.txt`；
- 运行 `bootstrap`（可配置参数 `B` 和 `block_length`），将 bootstrap 结果保存为 `algorithm/code/results/breakpoints.txt.bootstrap`；
- 读取结果并绘图（若 R 能访问图形设备）。

注意：在 Windows 系统上可执行文件为 `*.exe`，脚本中命令已做兼容处理。

## 输入 / 输出 文件格式

- 输入数据（二进制）：
  - `data/temp_data.bin`: 二进制保存的函数型数据矩阵（T × n_grid），按照项目内 `utils` 中定义的格式读写。
  - `data/temp_meta.txt`: 包含 `T`, `n_grid`, `s_min`, `s_max` 等元信息的文本文件。

- 断点检测输出： `results/breakpoints.txt`
  - 文本格式示例：

```
n_breaks=3

position,test_stat,p_value
147,65.7927240965,-1.0000000000
350,83.9392567141,-1.0000000000
451,7.1807956605,-1.0000000000
```

  - `p_value = -1` 表示该文件来自检测程序（未做 bootstrap 推断）。

- Bootstrap 输出： `results/breakpoints.txt.bootstrap`
  - 文本格式示例：

```
n_breaks=3

position,test_stat,p_value
147,65.7927240965,0.0000000000
350,83.9392567141,0.0000000000
451,7.1807956605,0.1841841842
```

  - 该文件由 `bootstrap` 程序自动生成（在原断点文件名后追加 `.bootstrap`）。

## 参数说明

- `max_breaks`：`break_detect` 中允许的最大断点数（默认示例为 10）。
- `min_segment_length`：检测中允许的最小段长度，用于避免过度细分。
- `B`：bootstrap 重采样次数，建议至少 1000 以保证 p 值稳定。
- `block_length`：移动块长度（20-50 常用，或使用 `ceil(T^(1/3))` 作为经验值）。

## 调试提示与已修复问题

- 若 bootstrap 程序意外崩溃，请检查：
  - 传入的 `breakpoint_file` 格式是否正确；
  - 程序内部是否正确初始化结构（`calloc` 用于 BootstrapResult）；
  - bootstrap 生成样本时不要释放共享内存（`s_grid`）导致 double-free。

- 已知改进：
  - 将 `BootstrapResult` 使用 `calloc` 初始化以避免未定义值；
  - 在释放 bootstrap 样本前避免释放共享的 `s_grid` 指针；
  - R 脚本中已修正命令参数顺序与结果解析，使 pipeline 能顺利运行。

## 示例（快速上手）

1. 在 `algorithm/code` 目录构建程序：

```powershell
cd algorithm\code
make all
```

2. 在工作区根目录运行完整分析：

```powershell
Rscript algorithm/run_analysis.R
```

3. 查看结果：

```powershell
type algorithm\code\results\breakpoints.txt
type algorithm\code\results\breakpoints.txt.bootstrap
```

4. 如果需要只运行 bootstrap（快速调试减少重复次数）：

```powershell
cd algorithm\code
./bootstrap ./data/temp_data.bin ./data/temp_meta.txt ./results/breakpoints.txt 100 8
```

## 贡献与联系方式

如果你对算法实现、参数、可视化或文档有建议，请在项目中提交 issue 或 pull request，或直接联系仓库维护者。

---

该 README 提供了足够的信息来理解实现的目的、主要算法思路、代码组织与运行方式。如需更精确的数学公式、伪代码和论文引用列表，我可以继续补充数学推导、伪代码和论文引用列表。
