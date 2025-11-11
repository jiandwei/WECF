# ============================================
# 完整工作流程：从数据生成到结果分析
# ============================================

source("R/01_data_generation.R")
source("R/02_call_c_program.R")

#' 运行完整的验证实验
#'
#' @param n_replications 蒙特卡罗重复次数
#' @param T_values 不同的样本量
#' @param break_types 断点类型
run_validation_experiment <- function(
    n_replications = 100,
    T_values = c(200, 300, 500),
    break_types = c("mean", "variance", "distribution")
) {
    library(parallel)
    library(pbapply)

    # 创建结果存储
    all_results <- list()

    for (break_type in break_types) {
        cat("\n", rep("=", 60), "\n", sep = "")
        cat("Break Type:", break_type, "\n")
        cat(rep("=", 60), "\n\n", sep = "")

        for (T in T_values) {
            cat("Sample Size T =", T, "\n")

            # 并行运行多次重复
            replication_results <- pblapply(
                1:n_replications,
                function(rep) {
                    set.seed(42 + rep)

                    # 生成数据
                    data <- generate_functional_breaks(
                        T = T,
                        grid_size = 100,
                        break_points = c(0.33, 0.67),
                        break_type = break_type,
                        snr = 3
                    )

                    # 保存到临时目录
                    temp_dir <- sprintf(
                        "data/temp_%s_T%d_rep%d",
                        break_type,
                        T,
                        rep
                    )
                    save_data_for_c(data, temp_dir)

                    # 运行检测
                    tryCatch(
                        {
                            detected <- run_breakpoint_detection(
                                input_dir = temp_dir,
                                output_dir = "results/temp"
                            )

                            # 评估
                            eval_result <- evaluate_detection(
                                data$true_breaks,
                                detected,
                                tolerance = max(10, floor(T * 0.05)) # 自适应容忍度
                            )

                            # 清理临时文件
                            unlink(temp_dir, recursive = TRUE)

                            return(list(
                                success = TRUE,
                                detected = detected,
                                metrics = eval_result
                            ))
                        },
                        error = function(e) {
                            return(list(
                                success = FALSE,
                                error = as.character(e)
                            ))
                        }
                    )
                },
                cl = detectCores() - 1
            )

            # 汇总结果
            successful <- sapply(replication_results, function(x) x$success)
            cat(sprintf("  Success rate: %.1f%%\n", 100 * mean(successful)))

            if (sum(successful) > 0) {
                metrics <- lapply(replication_results[successful], function(x) {
                    x$metrics
                })

                summary_stats <- data.frame(
                    T = T,
                    break_type = break_type,
                    success_rate = mean(successful),
                    mean_hausdorff = mean(sapply(metrics, function(x) {
                        x$hausdorff
                    })),
                    sd_hausdorff = sd(sapply(metrics, function(x) x$hausdorff)),
                    mean_precision = mean(sapply(metrics, function(x) {
                        x$precision
                    })),
                    mean_recall = mean(sapply(metrics, function(x) x$recall)),
                    mean_f1 = mean(sapply(metrics, function(x) x$f1_score))
                )

                all_results[[sprintf("%s_T%d", break_type, T)]] <- summary_stats

                cat(sprintf(
                    "  Avg Hausdorff: %.2f (SD: %.2f)\n",
                    summary_stats$mean_hausdorff,
                    summary_stats$sd_hausdorff
                ))
                cat(sprintf("  Avg F1 Score: %.3f\n", summary_stats$mean_f1))
            }
        }
    }

    # 合并所有结果
    final_results <- do.call(rbind, all_results)

    return(final_results)
}

#' 可视化验证实验结果
plot_validation_results <- function(results) {
    library(ggplot2)
    library(tidyr)

    # F1分数对比
    p1 <- ggplot(
        results,
        aes(x = T, y = mean_f1, color = break_type, group = break_type)
    ) +
        geom_line(size = 1.2) +
        geom_point(size = 3) +
        labs(
            x = "Sample Size (T)",
            y = "Mean F1 Score",
            title = "Detection Performance vs Sample Size",
            color = "Break Type"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")

    # Hausdorff距离对比
    p2 <- ggplot(
        results,
        aes(x = T, y = mean_hausdorff, color = break_type, group = break_type)
    ) +
        geom_line(size = 1.2) +
        geom_point(size = 3) +
        geom_errorbar(
            aes(
                ymin = mean_hausdorff - sd_hausdorff,
                ymax = mean_hausdorff + sd_hausdorff
            ),
            width = 10
        ) +
        labs(
            x = "Sample Size (T)",
            y = "Mean Hausdorff Distance",
            title = "Localization Accuracy vs Sample Size",
            color = "Break Type"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")

    # Precision vs Recall
    p3 <- ggplot(
        results,
        aes(
            x = mean_recall,
            y = mean_precision,
            color = break_type,
            shape = factor(T)
        )
    ) +
        geom_point(size = 4) +
        geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed",
            color = "gray50"
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(
            x = "Mean Recall",
            y = "Mean Precision",
            title = "Precision-Recall Trade-off",
            color = "Break Type",
            shape = "Sample Size"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")

    combined <- gridExtra::grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

    return(combined)
}

# ============================================
# 运行完整实验
# ============================================

if (FALSE) {
    # 小规模测试
    cat("Running small-scale validation...\n")
    results_small <- run_validation_experiment(
        n_replications = 20,
        T_values = c(200, 300),
        break_types = c("mean", "variance")
    )

    # 可视化
    p <- plot_validation_results(results_small)
    ggsave("results/validation_results.png", p, width = 12, height = 10)

    # 保存数值结果
    write.csv(
        results_small,
        "results/validation_summary.csv",
        row.names = FALSE
    )

    # 打印Markdown表格
    library(knitr)
    cat("\n## Validation Results\n\n")
    print(kable(results_small, digits = 3, format = "markdown"))
}
