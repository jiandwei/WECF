# 设置工作目录
setwd("algorithm\\code")

# 加载框架
source("R/04_monte_carlo.R")

# 方案1: 快速测试（100次模拟，约5分钟）
stats_quick <- complete_monte_carlo_analysis(n_sim = 100)

# # 方案2: 标准分析（1000次模拟，约30-60分钟）
# stats_full <- complete_monte_carlo_analysis(n_sim = 1000)

# # 方案3: 发表水平（5000次模拟，约3-5小时）
# stats_paper <- complete_monte_carlo_analysis(n_sim = 5000)
