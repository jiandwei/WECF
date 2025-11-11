source("R/01_data_generation.R")
source("R/02_call_c_program_improved.R")

# 生成测试数据
set.seed(42)
data <- generate_functional_breaks(
    T = 300,
    grid_size = 100,
    break_points = c(0.33, 0.67), # 真实断点
    break_type = "mean",
    snr = 3
)
save_data_for_c(data, "data/test")

# 运行改进算法
result <- run_improved_detection(
    input_dir = "data/test",
    significance_level = 0.05, # 5%显著性水平
    bic_penalty = 1.0 # 标准BIC
)

cat("\nDetected breaks:", result$positions, "\n")
cat("P-values:", result$p_values, "\n")
cat("True breaks:", data$true_breaks, "\n")

# 敏感性分析
sensitivity <- sensitivity_analysis("data/test")
print(sensitivity)
