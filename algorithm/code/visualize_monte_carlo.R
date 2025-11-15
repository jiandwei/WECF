# ============================================================================
# Monte Carlo结果可视化脚本
# ============================================================================

# 切换到 algorithm/code，确保相对路径一致
setwd("algorithm\\code")

# 加载可视化模块
source("R/05_monte_carlo_visualization.R")

# 生成所有可视化
generate_all_visualizations(
    results_dir = "results/monte_carlo",
    output_dir = "figures"
)
