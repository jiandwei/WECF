# ============================================
# R调用C程序并分析结果
# ============================================

library(ggplot2)
library(dplyr)

#' 编译并运行C程序
#'
#' @param input_dir 输入数据目录
#' @param output_dir 输出结果目录
#' @return 检测到的断点位置
run_breakpoint_detection <- function(
    input_dir = "data/mean_break",
    output_dir = "results"
) {
    # 确保输出目录存在
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # 编译C程序（如果需要）
    if (
        !file.exists("C/breakpoint_detector.exe") &&
            !file.exists("C/breakpoint_detector")
    ) {
        cat("Compiling C program...\n")
        compile_cmd <- "cd C && gcc -o breakpoint_detector core_algorithm.c -lm -fopenmp -O3"
        system(compile_cmd)
    }

    # 运行C程序
    cat("Running breakpoint detection...\n")
    if (.Platform$OS.type == "windows") {
        run_cmd <- sprintf(
            "C\\breakpoint_detector.exe %s %s",
            input_dir,
            output_dir
        )
    } else {
        run_cmd <- sprintf(
            "./C/breakpoint_detector %s %s",
            input_dir,
            output_dir
        )
    }

    system(run_cmd)

    # 读取结果
    result_file <- file.path(output_dir, "detected_breaks.csv")
    if (file.exists(result_file)) {
        detected <- read.csv(result_file)$position
        return(detected)
    } else {
        stop("Detection failed! No output file found.")
    }
}

#' 评估检测准确性
#'
#' @param true_breaks 真实断点
#' @param detected_breaks 检测到的断点
#' @param tolerance 容忍误差
#' @return 评估指标列表
evaluate_detection <- function(true_breaks, detected_breaks, tolerance = 10) {
    # Hausdorff距离
    if (length(detected_breaks) == 0) {
        hausdorff <- Inf
        precision <- 0
        recall <- 0
    } else if (length(true_breaks) == 0) {
        hausdorff <- Inf
        precision <- 0
        recall <- 1
    } else {
        dist_true_to_detected <- sapply(true_breaks, function(tb) {
            min(abs(tb - detected_breaks))
        })
        dist_detected_to_true <- sapply(detected_breaks, function(db) {
            min(abs(db - true_breaks))
        })
        hausdorff <- max(max(dist_true_to_detected), max(dist_detected_to_true))

        # 正确检测率（召回率）
        correctly_detected <- sum(dist_true_to_detected <= tolerance)
        recall <- correctly_detected / length(true_breaks)

        # 精确率
        true_positives <- sum(dist_detected_to_true <= tolerance)
        precision <- true_positives / length(detected_breaks)
    }

    # F1分数
    if (precision + recall > 0) {
        f1_score <- 2 * (precision * recall) / (precision + recall)
    } else {
        f1_score <- 0
    }

    return(list(
        hausdorff = hausdorff,
        precision = precision,
        recall = recall,
        f1_score = f1_score,
        n_true = length(true_breaks),
        n_detected = length(detected_breaks)
    ))
}

#' 可视化检测结果
#'
#' @param data_dir 数据目录
#' @param detected_breaks 检测到的断点
#' @param output_file 输出图像文件名
plot_detection_results <- function(
    data_dir,
    detected_breaks,
    output_file = NULL
) {
    # 读取原始数据
    X <- as.matrix(read.csv(
        file.path(data_dir, "functional_data.csv"),
        header = FALSE
    ))
    t_grid <- read.csv(file.path(data_dir, "time_grid.csv"), header = FALSE)$V1
    true_breaks <- read.csv(
        file.path(data_dir, "true_breaks.csv"),
        header = FALSE
    )$V1

    # 创建多面板图
    library(gridExtra)

    # 面板1：样本曲线着色
    n_plot <- min(50, nrow(X))
    sample_idx <- sample(1:nrow(X), n_plot)

    # 确定regime标签
    all_breaks <- sort(c(0, detected_breaks, nrow(X)))
    regime_labels <- findInterval(1:nrow(X), all_breaks)

    df1 <- data.frame(
        curve_id = rep(sample_idx, each = length(t_grid)),
        time = rep(t_grid, n_plot),
        value = as.vector(t(X[sample_idx, ])),
        regime = rep(regime_labels[sample_idx], each = length(t_grid))
    )

    p1 <- ggplot(
        df1,
        aes(x = time, y = value, group = curve_id, color = factor(regime))
    ) +
        geom_line(alpha = 0.4, size = 0.3) +
        geom_vline(
            xintercept = t_grid[1] +
                (true_breaks / nrow(X)) * diff(range(t_grid)),
            linetype = "dashed",
            color = "red",
            size = 1,
            alpha = 0.7
        ) +
        geom_vline(
            xintercept = t_grid[1] +
                (detected_breaks / nrow(X)) * diff(range(t_grid)),
            linetype = "solid",
            color = "blue",
            size = 1,
            alpha = 0.7
        ) +
        labs(
            x = "Time",
            y = "X(t)",
            title = "Detected Regimes",
            subtitle = "Red dashed: true breaks | Blue solid: detected breaks",
            color = "Detected\nRegime"
        ) +
        theme_minimal() +
        theme(legend.position = "right")

    # 面板2：均值函数轨迹
    mean_by_regime <- lapply(1:max(regime_labels), function(r) {
        idx <- which(regime_labels == r)
        colMeans(X[idx, , drop = FALSE])
    })

    df2 <- data.frame(
        time = rep(t_grid, max(regime_labels)),
        mean_value = unlist(mean_by_regime),
        regime = rep(1:max(regime_labels), each = length(t_grid))
    )

    p2 <- ggplot(df2, aes(x = time, y = mean_value, color = factor(regime))) +
        geom_line(size = 1.2) +
        geom_vline(
            xintercept = t_grid[1] +
                (true_breaks / nrow(X)) * diff(range(t_grid)),
            linetype = "dashed",
            color = "red",
            size = 0.8,
            alpha = 0.5
        ) +
        labs(
            x = "Time",
            y = "Mean Function",
            title = "Regime-Specific Mean Functions",
            color = "Regime"
        ) +
        theme_minimal() +
        theme(legend.position = "right")

    # 面板3：方差函数轨迹
    var_by_regime <- lapply(1:max(regime_labels), function(r) {
        idx <- which(regime_labels == r)
        apply(X[idx, , drop = FALSE], 2, var)
    })

    df3 <- data.frame(
        time = rep(t_grid, max(regime_labels)),
        var_value = unlist(var_by_regime),
        regime = rep(1:max(regime_labels), each = length(t_grid))
    )

    p3 <- ggplot(df3, aes(x = time, y = var_value, color = factor(regime))) +
        geom_line(size = 1.2) +
        geom_vline(
            xintercept = t_grid[1] +
                (true_breaks / nrow(X)) * diff(range(t_grid)),
            linetype = "dashed",
            color = "red",
            size = 0.8,
            alpha = 0.5
        ) +
        labs(
            x = "Time",
            y = "Variance Function",
            title = "Regime-Specific Variance Functions",
            color = "Regime"
        ) +
        theme_minimal() +
        theme(legend.position = "right")

    # 面板4：检测误差
    if (length(true_breaks) > 0 && length(detected_breaks) > 0) {
        matching <- sapply(true_breaks, function(tb) {
            idx <- which.min(abs(tb - detected_breaks))
            c(
                true = tb,
                detected = detected_breaks[idx],
                error = detected_breaks[idx] - tb
            )
        })

        df4 <- data.frame(
            true_position = matching[1, ],
            detected_position = matching[2, ],
            error = matching[3, ]
        )

        p4 <- ggplot(df4, aes(x = true_position, y = error)) +
            geom_point(size = 4, color = "darkred") +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
            geom_hline(
                yintercept = c(-10, 10),
                linetype = "dotted",
                color = "gray70"
            ) +
            labs(
                x = "True Break Position",
                y = "Detection Error",
                title = "Break Point Detection Accuracy",
                subtitle = sprintf(
                    "Mean absolute error: %.2f",
                    mean(abs(df4$error))
                )
            ) +
            theme_minimal()
    } else {
        # 如果没有匹配的断点，显示提示信息
        p4 <- ggplot() +
            annotate(
                "text",
                x = 0.5,
                y = 0.5,
                label = "No breaks to compare",
                size = 6,
                color = "gray50"
            ) +
            theme_void()
    }

    # 组合图
    combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

    # 保存
    if (!is.null(output_file)) {
        ggsave(output_file, combined_plot, width = 14, height = 10)
        cat("Visualization saved to", output_file, "\n")
    }

    return(combined_plot)
}

# ============================================
# 主测试脚本
# ============================================

if (FALSE) {
    library(knitr)
    library(kableExtra)

    # 测试所有三种类型的断点
    scenarios <- list(
        list(name = "Mean Break", dir = "data/mean_break"),
        list(name = "Variance Break", dir = "data/variance_break"),
        list(name = "Distribution Break", dir = "data/distribution_break")
    )

    results_summary <- data.frame()

    for (scenario in scenarios) {
        cat("\n", rep("=", 60), "\n", sep = "")
        cat("Testing:", scenario$name, "\n")
        cat(rep("=", 60), "\n\n", sep = "")

        # 运行检测
        detected <- run_breakpoint_detection(
            input_dir = scenario$dir,
            output_dir = "results"
        )

        # 读取真实断点
        true_breaks <- read.csv(
            file.path(scenario$dir, "true_breaks.csv"),
            header = FALSE
        )$V1

        cat("\nTrue breaks:", paste(true_breaks, collapse = ", "), "\n")
        cat("Detected breaks:", paste(detected, collapse = ", "), "\n\n")

        # 评估
        eval_result <- evaluate_detection(true_breaks, detected, tolerance = 10)

        cat("Evaluation Metrics:\n")
        cat(sprintf("  Hausdorff distance: %.2f\n", eval_result$hausdorff))
        cat(sprintf("  Precision: %.3f\n", eval_result$precision))
        cat(sprintf("  Recall: %.3f\n", eval_result$recall))
        cat(sprintf("  F1 Score: %.3f\n", eval_result$f1_score))

        # 可视化
        output_file <- sprintf(
            "results/detection_%s.png",
            gsub(" ", "_", tolower(scenario$name))
        )
        plot_detection_results(scenario$dir, detected, output_file)

        # 汇总结果
        results_summary <- rbind(
            results_summary,
            data.frame(
                Scenario = scenario$name,
                N_True = eval_result$n_true,
                N_Detected = eval_result$n_detected,
                Hausdorff = eval_result$hausdorff,
                Precision = eval_result$precision,
                Recall = eval_result$recall,
                F1_Score = eval_result$f1_score
            )
        )
    }

    # 输出汇总表格
    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Summary of All Scenarios\n")
    cat(rep("=", 60), "\n\n", sep = "")

    print(kable(results_summary, digits = 3, format = "markdown"))

    # 保存结果
    write.csv(results_summary, "results/summary_table.csv", row.names = FALSE)
}
