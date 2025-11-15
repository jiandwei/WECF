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

scenarios <- list(
    list(
        name = "T500_AR_two_breaks",
        T = 500,
        breakpoints = c(0.3, 0.7),
        dependence = "AR1",
        rho = 0.5
    ),
    list(
        name = "T800_MA_single_break",
        T = 800,
        breakpoints = 0.55,
        dependence = "MA1",
        rho = 0.4
    )
)

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
