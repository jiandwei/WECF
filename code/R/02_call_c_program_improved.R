#' 运行改进的断点检测
#'
#' @param input_dir 输入目录
#' @param output_dir 输出目录
#' @param significance_level 显著性水平
#' @param bic_penalty BIC惩罚因子（越大越保守）
run_improved_detection <- function(
    input_dir = "data/mean_break",
    output_dir = "results",
    significance_level = 0.05,
    bic_penalty = 1.0
) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # 编译（如果需要）
    if (
        !file.exists("C/breakpoint_detector_improved.exe") &&
            !file.exists("C/breakpoint_detector_improved")
    ) {
        cat("Compiling improved C program...\n")
        compile_cmd <- "cd C && gcc -o breakpoint_detector_improved core_algorithm_improved.c -lm -fopenmp -O3"
        system(compile_cmd)
    }

    # 运行
    cat(sprintf(
        "Running with α=%.3f, BIC penalty=%.2f...\n",
        significance_level,
        bic_penalty
    ))

    if (.Platform$OS.type == "windows") {
        run_cmd <- sprintf(
            "C\\breakpoint_detector_improved.exe %s %s %.3f %.2f",
            input_dir,
            output_dir,
            significance_level,
            bic_penalty
        )
    } else {
        run_cmd <- sprintf(
            "./C/breakpoint_detector_improved %s %s %.3f %.2f",
            input_dir,
            output_dir,
            significance_level,
            bic_penalty
        )
    }

    system(run_cmd)

    # 读取结果
    result_file <- file.path(output_dir, "detected_breaks.csv")
    if (file.exists(result_file)) {
        result_df <- read.csv(result_file)
        return(list(
            positions = result_df$position,
            p_values = if ("p_value" %in% colnames(result_df)) {
                result_df$p_value
            } else {
                NA
            },
            fractions = if ("fraction" %in% colnames(result_df)) {
                result_df$fraction
            } else {
                NA
            }
        ))
    } else {
        stop("Detection failed!")
    }
}

#' 参数敏感性分析
sensitivity_analysis <- function(
    data_dir,
    significance_levels = c(0.01, 0.05, 0.10),
    bic_penalties = c(0.5, 1.0, 2.0)
) {
    results_grid <- expand.grid(
        alpha = significance_levels,
        penalty = bic_penalties
    )

    results_grid$n_breaks <- NA
    results_grid$positions <- vector("list", nrow(results_grid))

    for (i in 1:nrow(results_grid)) {
        cat(sprintf(
            "\n[%d/%d] α=%.3f, penalty=%.2f\n",
            i,
            nrow(results_grid),
            results_grid$alpha[i],
            results_grid$penalty[i]
        ))

        result <- run_improved_detection(
            input_dir = data_dir,
            output_dir = "results/sensitivity",
            significance_level = results_grid$alpha[i],
            bic_penalty = results_grid$penalty[i]
        )

        results_grid$n_breaks[i] <- length(result$positions)
        results_grid$positions[[i]] <- result$positions
    }

    return(results_grid)
}
