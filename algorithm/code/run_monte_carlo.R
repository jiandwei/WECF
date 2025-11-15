# ============================================================================
# Monte Carlo driver script
# ============================================================================

# 切换到 algorithm/code，确保相对路径一致
setwd("algorithm\\code")
source("R/01_data_generation.R")
source("R/02_visualization.R")
source("R/03_main_analysis.R")
source("R/04_monte_carlo.R")

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
    pattern <- sprintf("^--%s=", name)
    hit <- grep(pattern, args, value = TRUE)
    if (length(hit) == 0) {
        return(default)
    }
    sub(pattern, "", hit[1])
}

get_bool_arg <- function(name, default = FALSE) {
    if (sprintf("--%s", name) %in% args) {
        return(TRUE)
    }
    val <- get_arg(name, NA)
    if (is.na(val)) {
        return(default)
    }
    tolower(val) %in% c("1", "true", "yes", "on")
}

reps <- as.integer(get_arg("reps", 100))
if (is.na(reps) || reps <= 0) {
    stop("--reps 需要为正整数")
}

run_bootstrap <- get_bool_arg("bootstrap", TRUE)

tolerance <- as.integer(get_arg("tolerance", 25))
if (is.na(tolerance) || tolerance <= 0) {
    stop("--tolerance 需要为正整数")
}

base_seed <- as.integer(get_arg("seed", 2024))
if (is.na(base_seed)) {
    stop("--seed 需要为整数")
}

cat("\nMonte Carlo 参数:\n")
cat(sprintf("  重复次数: %d\n", reps))
cat(sprintf("  断点匹配容忍度: %d\n", tolerance))
cat(sprintf("  是否运行 bootstrap: %s\n", ifelse(run_bootstrap, "是", "否")))
cat(sprintf("  基础随机种子: %d\n", base_seed))

# ============================================================================
# 构建全面的断点场景组合
# ============================================================================
# 涵盖以下维度：
# 1. 样本量: T ∈ {300, 500, 700, 1000}
# 2. 断点数量: M ∈ {0, 1, 2, 3, 4}
# 3. 断点位置模式: 早期/中期/晚期/均匀分布/密集/稀疏
# 4. 时间依赖结构: AR(1), MA(1), independent
# 5. 依赖强度: ρ ∈ {0.3, 0.5, 0.7}
# ============================================================================

scenarios <- list()

# ---------- 场景组1: 无断点基准场景 (Null Hypothesis) ----------
scenarios[[length(scenarios) + 1]] <- list(
    name = "S01_T500_no_break_AR_rho05",
    T = 500,
    breakpoints = numeric(0),
    dependence = "AR1",
    rho = 0.5
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S02_T700_no_break_MA_rho03",
    T = 700,
    breakpoints = numeric(0),
    dependence = "MA1",
    rho = 0.3
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S03_T1000_no_break_independent",
    T = 1000,
    breakpoints = numeric(0),
    dependence = "independent"
)

# ---------- 场景组2: 单断点场景 (不同位置) ----------
# 早期断点 (0.2-0.3)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S04_T500_single_early_020_AR_rho05",
    T = 500,
    breakpoints = 0.2,
    dependence = "AR1",
    rho = 0.5
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S05_T700_single_early_025_MA_rho07",
    T = 700,
    breakpoints = 0.25,
    dependence = "MA1",
    rho = 0.7
)

# 中期断点 (0.45-0.55)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S06_T500_single_mid_050_AR_rho03",
    T = 500,
    breakpoints = 0.5,
    dependence = "AR1",
    rho = 0.3
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S07_T1000_single_mid_055_independent",
    T = 1000,
    breakpoints = 0.55,
    dependence = "independent"
)

# 晚期断点 (0.7-0.8)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S08_T700_single_late_075_MA_rho05",
    T = 700,
    breakpoints = 0.75,
    dependence = "MA1",
    rho = 0.5
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S09_T300_single_late_080_AR_rho07",
    T = 300,
    breakpoints = 0.8,
    dependence = "AR1",
    rho = 0.7
)

# ---------- 场景组3: 两断点场景 (不同间距) ----------
# 均匀分布
scenarios[[length(scenarios) + 1]] <- list(
    name = "S10_T500_two_uniform_033_067_AR_rho05",
    T = 500,
    breakpoints = c(0.33, 0.67),
    dependence = "AR1",
    rho = 0.5
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S11_T700_two_uniform_030_070_MA_rho03",
    T = 700,
    breakpoints = c(0.3, 0.7),
    dependence = "MA1",
    rho = 0.3
)

# 密集分布 (间距 < 0.3)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S12_T1000_two_close_040_055_AR_rho07",
    T = 1000,
    breakpoints = c(0.4, 0.55),
    dependence = "AR1",
    rho = 0.7
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S13_T800_two_close_045_060_independent",
    T = 800,
    breakpoints = c(0.45, 0.6),
    dependence = "independent"
)

# 稀疏分布 (间距 > 0.5)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S14_T700_two_sparse_020_075_MA_rho05",
    T = 700,
    breakpoints = c(0.2, 0.75),
    dependence = "MA1",
    rho = 0.5
)

# ---------- 场景组4: 三断点场景 ----------
# 均匀分布
scenarios[[length(scenarios) + 1]] <- list(
    name = "S15_T700_three_uniform_025_050_075_AR_rho05",
    T = 700,
    breakpoints = c(0.25, 0.5, 0.75),
    dependence = "AR1",
    rho = 0.5
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S16_T1000_three_uniform_020_050_080_MA_rho03",
    T = 1000,
    breakpoints = c(0.2, 0.5, 0.8),
    dependence = "MA1",
    rho = 0.3
)

# 不均匀分布 (前密后疏)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S17_T800_three_uneven_020_035_070_AR_rho07",
    T = 800,
    breakpoints = c(0.2, 0.35, 0.7),
    dependence = "AR1",
    rho = 0.7
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S18_T600_three_uneven_030_045_075_independent",
    T = 600,
    breakpoints = c(0.3, 0.45, 0.75),
    dependence = "independent"
)

# ---------- 场景组5: 四断点场景 (高复杂度) ----------
scenarios[[length(scenarios) + 1]] <- list(
    name = "S19_T1000_four_uniform_020_040_060_080_AR_rho05",
    T = 1000,
    breakpoints = c(0.2, 0.4, 0.6, 0.8),
    dependence = "AR1",
    rho = 0.5
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S20_T800_four_mixed_015_035_055_075_MA_rho07",
    T = 800,
    breakpoints = c(0.15, 0.35, 0.55, 0.75),
    dependence = "MA1",
    rho = 0.7
)

scenarios[[length(scenarios) + 1]] <- list(
    name = "S21_T1200_four_dense_030_045_055_070_AR_rho03",
    T = 1200,
    breakpoints = c(0.3, 0.45, 0.55, 0.7),
    dependence = "AR1",
    rho = 0.3
)

# ---------- 场景组6: 极端场景 (边界测试) ----------
# 非常早的断点
scenarios[[length(scenarios) + 1]] <- list(
    name = "S22_T500_single_veryearly_010_AR_rho05",
    T = 500,
    breakpoints = 0.1,
    dependence = "AR1",
    rho = 0.5
)

# 非常晚的断点
scenarios[[length(scenarios) + 1]] <- list(
    name = "S23_T500_single_verylate_090_MA_rho05",
    T = 500,
    breakpoints = 0.9,
    dependence = "MA1",
    rho = 0.5
)

# 小样本 + 多断点
scenarios[[length(scenarios) + 1]] <- list(
    name = "S24_T300_three_025_050_075_AR_rho07",
    T = 300,
    breakpoints = c(0.25, 0.5, 0.75),
    dependence = "AR1",
    rho = 0.7
)

# 大样本 + 密集断点
scenarios[[length(scenarios) + 1]] <- list(
    name = "S25_T1500_five_020_035_050_065_080_MA_rho05",
    T = 1500,
    breakpoints = c(0.2, 0.35, 0.5, 0.65, 0.8),
    dependence = "MA1",
    rho = 0.5
)

# ---------- 场景组7: 不同依赖强度对比 ----------
# 弱依赖 (ρ = 0.3)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S26_T500_two_030_070_AR_rho03_weak",
    T = 500,
    breakpoints = c(0.3, 0.7),
    dependence = "AR1",
    rho = 0.3
)

# 中等依赖 (ρ = 0.5)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S27_T500_two_030_070_AR_rho05_medium",
    T = 500,
    breakpoints = c(0.3, 0.7),
    dependence = "AR1",
    rho = 0.5
)

# 强依赖 (ρ = 0.7)
scenarios[[length(scenarios) + 1]] <- list(
    name = "S28_T500_two_030_070_AR_rho07_strong",
    T = 500,
    breakpoints = c(0.3, 0.7),
    dependence = "AR1",
    rho = 0.7
)

# ---------- 场景组8: 不同样本量对比 ----------
# 小样本
scenarios[[length(scenarios) + 1]] <- list(
    name = "S29_T300_two_033_067_AR_rho05_small",
    T = 300,
    breakpoints = c(0.33, 0.67),
    dependence = "AR1",
    rho = 0.5
)

# 中等样本
scenarios[[length(scenarios) + 1]] <- list(
    name = "S30_T700_two_033_067_AR_rho05_medium",
    T = 700,
    breakpoints = c(0.33, 0.67),
    dependence = "AR1",
    rho = 0.5
)

# 大样本
scenarios[[length(scenarios) + 1]] <- list(
    name = "S31_T1500_two_033_067_AR_rho05_large",
    T = 1500,
    breakpoints = c(0.33, 0.67),
    dependence = "AR1",
    rho = 0.5
)

cat(sprintf("\n共生成 %d 个断点场景\n", length(scenarios)))

results <- run_monte_carlo(
    scenarios = scenarios,
    n_rep = reps,
    detection_params = list(max_breaks = 6, min_segment_length = 40),
    bootstrap_params = list(B = 199),
    tolerance = tolerance,
    run_bootstrap = run_bootstrap,
    base_seed = base_seed,
    output_dir = "results/monte_carlo"
)

cat("\nMonte Carlo 汇总:\n")
print(results$summary)

cat("\n结果文件:\n")
cat(sprintf("  逐次结果: %s\n", results$files$detailed))
cat(sprintf("  汇总结果: %s\n", results$files$summary))
