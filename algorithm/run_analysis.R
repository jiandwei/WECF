# ============================================================================
# 主分析脚本
# ============================================================================

# 设置工作目录
setwd("WECF\\algorithm\\code")

# 编译C程序（如果未编译）
if (!file.exists("break_detect.exe")) {
    cat("首次运行，需要编译C程序...\n")
    system("make all")
}

# 加载模块
source("R/01_data_generation.R")
source("R/02_visualization.R")
source("R/03_main_analysis.R")

# 运行完整分析
results <- complete_analysis(
    T = 500,
    breakpoints = c(0.3, 0.7),
    run_bootstrap = TRUE,
    visualize = TRUE
)

# 评估准确性
true_breaks <- results$data$true_break_indices
detected_breaks <- results$breakpoints$breaks$position

cat("\n===== 准确性评估 =====\n")
cat("真实断点:", true_breaks, "\n")
cat("检测断点:", detected_breaks, "\n")

if (length(detected_breaks) == length(true_breaks)) {
    errors <- abs(detected_breaks - true_breaks)
    cat("估计误差:", errors, "\n")
    cat("平均误差:", mean(errors), "\n")
}
