# ============================================================================
# 文件名: 05_monte_carlo_visualization.R
# 功能: Monte Carlo模拟结果可视化
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

#' 加载Monte Carlo结果
#'
#' @param results_dir 结果目录路径
#' @return 列表包含detailed和summary数据框
load_monte_carlo_results <- function(results_dir = "results/monte_carlo") {
    detailed_path <- file.path(results_dir, "monte_carlo_results.csv")
    summary_path <- file.path(results_dir, "monte_carlo_summary.csv")

    if (!file.exists(detailed_path) || !file.exists(summary_path)) {
        stop("未找到Monte Carlo结果文件，请先运行模拟")
    }

    list(
        detailed = read.csv(detailed_path, stringsAsFactors = FALSE),
        summary = read.csv(summary_path, stringsAsFactors = FALSE)
    )
}

#' 提取场景元信息
#'
#' @param scenario_names 场景名称向量
#' @return 数据框包含场景分组信息
extract_scenario_info <- function(scenario_names) {
    # 场景分组规则
    data.frame(
        scenario = scenario_names,
        group = case_when(
            grepl("no_break", scenario_names) ~ "1. 无断点",
            grepl("single", scenario_names) ~ "2. 单断点",
            grepl("two", scenario_names) ~ "3. 两断点",
            grepl("three", scenario_names) ~ "4. 三断点",
            grepl("four|five", scenario_names) ~ "5. 多断点",
            TRUE ~ "其他"
        ),
        n_breaks = case_when(
            grepl("no_break", scenario_names) ~ 0,
            grepl("single", scenario_names) ~ 1,
            grepl("two", scenario_names) ~ 2,
            grepl("three", scenario_names) ~ 3,
            grepl("four", scenario_names) ~ 4,
            grepl("five", scenario_names) ~ 5,
            TRUE ~ NA_real_
        ),
        sample_size = as.numeric(gsub(".*_T(\\d+)_.*", "\\1", scenario_names)),
        dependence = case_when(
            grepl("_AR_", scenario_names) ~ "AR(1)",
            grepl("_MA_", scenario_names) ~ "MA(1)",
            grepl("independent", scenario_names) ~ "Independent",
            TRUE ~ NA_character_
        ),
        rho = case_when(
            grepl("rho03", scenario_names) ~ 0.3,
            grepl("rho05", scenario_names) ~ 0.5,
            grepl("rho07", scenario_names) ~ 0.7,
            TRUE ~ NA_real_
        ),
        stringsAsFactors = FALSE
    )
}

#' 绘制检测精度对比（按断点数量分组）
#'
#' @param summary_df 汇总数据框
#' @param output_file 输出文件路径
plot_accuracy_by_breaks <- function(summary_df, output_file = NULL) {
    scenario_info <- extract_scenario_info(summary_df$scenario)
    df <- cbind(summary_df, scenario_info[, -1])

    # 过滤掉无断点场景
    df_with_breaks <- df %>% filter(n_breaks > 0)

    p <- ggplot(df_with_breaks, aes(x = factor(n_breaks), y = precision)) +
        geom_boxplot(fill = "steelblue", alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
        labs(
            title = "断点检测精度 vs 真实断点数量",
            subtitle = "Precision = TP / (TP + FP)",
            x = "真实断点数量",
            y = "精度 (Precision)"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5)
        ) +
        ylim(0, 1)

    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 8, height = 6, dpi = 300)
    }

    print(p)
    invisible(p)
}

#' 绘制召回率对比
#'
#' @param summary_df 汇总数据框
#' @param output_file 输出文件路径
plot_recall_by_breaks <- function(summary_df, output_file = NULL) {
    scenario_info <- extract_scenario_info(summary_df$scenario)
    df <- cbind(summary_df, scenario_info[, -1])

    df_with_breaks <- df %>% filter(n_breaks > 0)

    p <- ggplot(df_with_breaks, aes(x = factor(n_breaks), y = recall)) +
        geom_boxplot(fill = "coral", alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
        labs(
            title = "断点检测召回率 vs 真实断点数量",
            subtitle = "Recall = TP / (TP + FN)",
            x = "真实断点数量",
            y = "召回率 (Recall)"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5)
        ) +
        ylim(0, 1)

    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 8, height = 6, dpi = 300)
    }

    print(p)
    invisible(p)
}

#' 绘制依赖强度影响分析
#'
#' @param summary_df 汇总数据框
#' @param output_file 输出文件路径
plot_dependence_effect <- function(summary_df, output_file = NULL) {
    scenario_info <- extract_scenario_info(summary_df$scenario)
    df <- cbind(summary_df, scenario_info[, -1])

    # 选择依赖强度对比场景 (S26-S28)
    df_rho <- df %>%
        filter(grepl("S26|S27|S28", scenario)) %>%
        filter(!is.na(rho))

    if (nrow(df_rho) == 0) {
        warning("未找到依赖强度对比场景")
        return(invisible(NULL))
    }

    p1 <- ggplot(
        df_rho,
        aes(x = factor(rho), y = precision, fill = factor(rho))
    ) +
        geom_bar(stat = "identity", alpha = 0.7) +
        labs(
            title = "依赖强度对精度的影响",
            x = "ρ (依赖参数)",
            y = "精度"
        ) +
        scale_fill_brewer(palette = "Set2", name = "ρ") +
        theme_minimal() +
        theme(legend.position = "none")

    p2 <- ggplot(df_rho, aes(x = factor(rho), y = recall, fill = factor(rho))) +
        geom_bar(stat = "identity", alpha = 0.7) +
        labs(
            title = "依赖强度对召回率的影响",
            x = "ρ (依赖参数)",
            y = "召回率"
        ) +
        scale_fill_brewer(palette = "Set2", name = "ρ") +
        theme_minimal() +
        theme(legend.position = "none")

    p <- grid.arrange(p1, p2, ncol = 2)

    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 12, height = 5, dpi = 300)
    }

    invisible(p)
}

#' 绘制样本量影响分析
#'
#' @param summary_df 汇总数据框
#' @param output_file 输出文件路径
plot_sample_size_effect <- function(summary_df, output_file = NULL) {
    scenario_info <- extract_scenario_info(summary_df$scenario)
    df <- cbind(summary_df, scenario_info[, -1])

    # 选择样本量对比场景 (S29-S31)
    df_size <- df %>%
        filter(grepl("S29|S30|S31", scenario))

    if (nrow(df_size) == 0) {
        warning("未找到样本量对比场景")
        return(invisible(NULL))
    }

    p1 <- ggplot(df_size, aes(x = sample_size, y = precision)) +
        geom_line(color = "steelblue", size = 1.2) +
        geom_point(size = 3, color = "steelblue") +
        labs(
            title = "样本量对精度的影响",
            x = "样本量 (T)",
            y = "精度"
        ) +
        theme_minimal() +
        ylim(0, 1)

    p2 <- ggplot(df_size, aes(x = sample_size, y = mean_abs_error)) +
        geom_line(color = "coral", size = 1.2) +
        geom_point(size = 3, color = "coral") +
        labs(
            title = "样本量对断点位置误差的影响",
            x = "样本量 (T)",
            y = "平均绝对误差"
        ) +
        theme_minimal()

    p <- grid.arrange(p1, p2, ncol = 2)

    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 12, height = 5, dpi = 300)
    }

    invisible(p)
}

#' 绘制计算时间分析
#'
#' @param summary_df 汇总数据框
#' @param output_file 输出文件路径
plot_computation_time <- function(summary_df, output_file = NULL) {
    scenario_info <- extract_scenario_info(summary_df$scenario)
    df <- cbind(summary_df, scenario_info[, -1])

    p <- ggplot(
        df,
        aes(x = sample_size, y = detection_time, color = factor(n_breaks))
    ) +
        geom_point(size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        labs(
            title = "计算时间 vs 样本量和断点数量",
            x = "样本量 (T)",
            y = "计算时间 (秒)"
        ) +
        scale_color_brewer(palette = "Set1", name = "断点数") +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 10, height = 6, dpi = 300)
    }

    print(p)
    invisible(p)
}

#' 生成汇总统计表格
#'
#' @param summary_df 汇总数据框
#' @param output_file 输出文件路径
generate_summary_table <- function(summary_df, output_file = NULL) {
    scenario_info <- extract_scenario_info(summary_df$scenario)
    df <- cbind(summary_df, scenario_info[, -1])

    # 按分组汇总
    table_by_group <- df %>%
        group_by(group) %>%
        summarise(
            n_scenarios = n(),
            avg_precision = mean(precision, na.rm = TRUE),
            avg_recall = mean(recall, na.rm = TRUE),
            avg_detection_time = mean(detection_time, na.rm = TRUE),
            avg_mean_abs_error = mean(mean_abs_error, na.rm = TRUE)
        )

    print(table_by_group)

    if (!is.null(output_file)) {
        write.csv(table_by_group, output_file, row.names = FALSE)
    }

    invisible(table_by_group)
}

#' 主函数：生成所有可视化
#'
#' @param results_dir 结果目录
#' @param output_dir 输出目录
generate_all_visualizations <- function(
    results_dir = "results/monte_carlo",
    output_dir = "figures"
) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    cat("加载Monte Carlo结果...\n")
    results <- load_monte_carlo_results(results_dir)

    cat("\n生成可视化图表...\n")

    cat("  [1/6] 精度 vs 断点数量\n")
    plot_accuracy_by_breaks(
        results$summary,
        file.path(output_dir, "mc_precision_by_breaks.png")
    )

    cat("  [2/6] 召回率 vs 断点数量\n")
    plot_recall_by_breaks(
        results$summary,
        file.path(output_dir, "mc_recall_by_breaks.png")
    )

    cat("  [3/6] 依赖强度影响\n")
    plot_dependence_effect(
        results$summary,
        file.path(output_dir, "mc_dependence_effect.png")
    )

    cat("  [4/6] 样本量影响\n")
    plot_sample_size_effect(
        results$summary,
        file.path(output_dir, "mc_sample_size_effect.png")
    )

    cat("  [5/6] 计算时间分析\n")
    plot_computation_time(
        results$summary,
        file.path(output_dir, "mc_computation_time.png")
    )

    cat("  [6/6] 汇总统计表格\n")
    generate_summary_table(
        results$summary,
        file.path(output_dir, "mc_summary_table.csv")
    )

    cat("\n✓ 所有可视化已完成！\n")
    cat(sprintf("  输出目录: %s\n", output_dir))
}

# ============================================================================
# 示例使用
# ============================================================================

if (FALSE) {
    # 切换到代码目录
    setwd("algorithm/code")

    # 生成所有可视化
    generate_all_visualizations(
        results_dir = "results/monte_carlo",
        output_dir = "figures"
    )
}
