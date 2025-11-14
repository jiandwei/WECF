#' ============================================================================
#' 文件名: 04_monte_carlo.R
#' 功能: 蒙特卡洛模拟验证断点检测方法的统计性质
#' 内容:
#'   1. Size测试（检验水平）
#'   2. Power测试（检验功效）
#'   3. 断点位置估计精度
#'   4. 不同场景下的性能比较
#' ============================================================================

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

source("R/01_data_generation.R")
source("R/03_main_analysis.R")

#' ========== 场景定义 ==========

#' 定义模拟场景
define_scenarios <- function() {
    scenarios <- list(
        # ========== Size测试场景 ==========
        size_null = list(
            name = "Size: H0（无断点）",
            T = 300,
            breakpoints = NULL, # 无断点
            n_breaks_true = 0,
            type = "size"
        ),

        # ========== Power测试场景 ==========

        # 1. 单个断点 - 不同位置
        power_single_early = list(
            name = "Power: 单断点（早期，0.2T）",
            T = 300,
            breakpoints = 0.2,
            n_breaks_true = 1,
            type = "power"
        ),

        power_single_mid = list(
            name = "Power: 单断点（中期，0.5T）",
            T = 300,
            breakpoints = 0.5,
            n_breaks_true = 1,
            type = "power"
        ),

        power_single_late = list(
            name = "Power: 单断点（晚期，0.8T）",
            T = 300,
            breakpoints = 0.8,
            n_breaks_true = 1,
            type = "power"
        ),

        # 2. 多个断点
        power_double = list(
            name = "Power: 双断点（0.33T, 0.67T）",
            T = 300,
            breakpoints = c(0.33, 0.67),
            n_breaks_true = 2,
            type = "power"
        ),

        power_triple = list(
            name = "Power: 三断点（0.25T, 0.5T, 0.75T）",
            T = 300,
            breakpoints = c(0.25, 0.5, 0.75),
            n_breaks_true = 3,
            type = "power"
        ),

        # 3. 不同样本量
        power_small_n = list(
            name = "Power: 小样本（T=150）",
            T = 150,
            breakpoints = 0.5,
            n_breaks_true = 1,
            type = "power"
        ),

        power_large_n = list(
            name = "Power: 大样本（T=600）",
            T = 600,
            breakpoints = 0.5,
            n_breaks_true = 1,
            type = "power"
        ),

        # 4. 不同断点幅度（通过修改生成参数）
        power_weak = list(
            name = "Power: 弱断点（小幅度变化）",
            T = 300,
            breakpoints = 0.5,
            n_breaks_true = 1,
            type = "power"
            # 使用默认参数，变化较小
        ),

        power_strong = list(
            name = "Power: 强断点（大幅度变化）",
            T = 300,
            breakpoints = 0.5,
            n_breaks_true = 1,
            type = "power"
            # 使用默认参数，检测应该更容易
        )
    )

    return(scenarios)
}

#' ========== 单次模拟运行 ==========

#' 运行单次模拟
#'
#' @param scenario 场景配置
#' @param max_breaks 最大检测断点数
#' @param min_segment_length 最小segment长度
#' @param sim_id 模拟ID（用于避免文件冲突）
#' @param verbose 是否输出详细信息
#' @return 模拟结果列表
run_single_simulation <- function(
    scenario,
    max_breaks = 10,
    min_segment_length = 30,
    sim_id = 1,
    verbose = FALSE
) {
    # 生成数据
    data <- generate_functional_data(
        T = scenario$T,
        breakpoints = scenario$breakpoints
    )

    # 运行断点检测（不输出详细信息）
    suppressMessages({
        capture.output(
            {
                bp_results <- run_breakpoint_detection(
                    data,
                    method = "binary_seg",
                    params = list(
                        max_breaks = max_breaks,
                        min_segment_length = min_segment_length
                    ),
                    prefix = paste0("sim_", sim_id)
                )
            },
            type = "output"
        )
    })

    # 提取结果
    n_breaks_detected <- bp_results$n_breaks

    # 计算真实断点位置
    if (!is.null(scenario$breakpoints)) {
        true_positions <- round(scenario$breakpoints * scenario$T)
    } else {
        true_positions <- NULL
    }

    # 计算位置估计误差（如果有断点）
    position_errors <- NA
    if (n_breaks_detected > 0 && !is.null(true_positions)) {
        # 匹配检测到的断点与真实断点
        detected_positions <- bp_results$breaks$position

        # 简单匹配：每个检测断点对应最近的真实断点
        errors <- sapply(detected_positions, function(det_pos) {
            min(abs(det_pos - true_positions))
        })
        position_errors <- mean(errors)
    }

    # 计算检测准确性
    is_correct <- FALSE
    if (n_breaks_detected == scenario$n_breaks_true) {
        if (scenario$n_breaks_true == 0) {
            is_correct <- TRUE
        } else if (!is.null(true_positions)) {
            # 检查每个真实断点是否被正确检测（容忍度：±5%）
            tolerance <- ceiling(0.05 * scenario$T)
            detected_positions <- bp_results$position

            all_detected <- sapply(true_positions, function(true_pos) {
                any(abs(detected_positions - true_pos) <= tolerance)
            })

            is_correct <- all(all_detected)
        }
    }

    # 返回结果
    result <- list(
        scenario_name = scenario$name,
        scenario_type = scenario$type,
        n_breaks_true = scenario$n_breaks_true,
        n_breaks_detected = n_breaks_detected,
        is_correct = is_correct,
        position_error = position_errors,
        test_stats = if (n_breaks_detected > 0) bp_results$test_stat else NA
    )

    return(result)
}

#' 运行蒙特卡洛模拟
run_monte_carlo <- function(
    n_sim = 1000,
    scenarios = NULL,
    max_breaks = 10,
    min_segment_length = 30,
    parallel = TRUE,
    n_cores = 8
) {
    cat("========================================\n")
    cat("蒙特卡洛模拟验证\n")
    cat("========================================\n")
    cat("模拟次数:", n_sim, "\n")
    cat("场景数量:", length(scenarios), "\n")

    if (is.null(scenarios)) {
        scenarios <- define_scenarios()
    }

    # 存储结果
    all_results <- list()

    # 对每个场景进行模拟
    for (i in seq_along(scenarios)) {
        scenario <- scenarios[[i]]

        cat("\n----------------------------------------\n")
        cat("场景", i, "/", length(scenarios), ":", scenario$name, "\n")
        cat("----------------------------------------\n")

        start_time <- Sys.time()

        # 运行n_sim次模拟
        if (parallel && requireNamespace("parallel", quietly = TRUE)) {
            cat("使用", n_cores, "个核心并行运行...\n")

            # 获取当前工作目录
            current_wd <- getwd()

            cl <- parallel::makeCluster(n_cores)

            # 在每个工作进程中设置工作目录并加载源文件
            parallel::clusterEvalQ(cl, {
                library(dplyr)
            })

            # 设置工作目录
            parallel::clusterExport(cl, "current_wd", envir = environment())
            parallel::clusterEvalQ(cl, {
                setwd(current_wd)
            })

            # 加载源文件
            parallel::clusterEvalQ(cl, {
                source("R/01_data_generation.R")
                source("R/03_main_analysis.R")
            })

            # 导出必要的对象
            parallel::clusterExport(
                cl,
                c(
                    "run_single_simulation",
                    "scenario",
                    "max_breaks",
                    "min_segment_length"
                ),
                envir = environment()
            )

            results <- parallel::parLapply(cl, 1:n_sim, function(j) {
                run_single_simulation(
                    scenario,
                    max_breaks,
                    min_segment_length,
                    sim_id = j
                )
            })

            parallel::stopCluster(cl)
        } else {
            # 串行运行
            results <- lapply(1:n_sim, function(j) {
                if (j %% 100 == 0) {
                    cat("\r  进度:", j, "/", n_sim)
                }
                run_single_simulation(
                    scenario,
                    max_breaks,
                    min_segment_length,
                    sim_id = j
                )
            })
            cat("\n")
        }

        end_time <- Sys.time()
        time_elapsed <- difftime(end_time, start_time, units = "secs")

        cat("  √ 完成! 用时:", round(time_elapsed, 2), "秒\n")

        # 整理结果
        results_df <- do.call(
            rbind,
            lapply(results, function(r) {
                data.frame(
                    scenario_name = r$scenario_name,
                    scenario_type = r$scenario_type,
                    n_breaks_true = r$n_breaks_true,
                    n_breaks_detected = r$n_breaks_detected,
                    is_correct = r$is_correct,
                    position_error = r$position_error,
                    stringsAsFactors = FALSE
                )
            })
        )

        all_results[[i]] <- results_df
    }

    # 合并所有结果
    final_results <- do.call(rbind, all_results)

    cat("\n========================================\n")
    cat("模拟完成!\n")
    cat("========================================\n")

    return(final_results)
}


#' ========== 结果分析 ==========

#' 计算统计性质
#'
#' @param results 模拟结果数据框
#' @return 统计性质汇总
compute_statistics <- function(results) {
    stats <- results %>%
        group_by(scenario_name, scenario_type, n_breaks_true) %>%
        summarise(
            # Size/Power
            rejection_rate = mean(n_breaks_detected > 0),

            # 正确检测率
            correct_rate = mean(is_correct),

            # 断点数量统计
            mean_n_breaks = mean(n_breaks_detected),
            sd_n_breaks = sd(n_breaks_detected),

            # 位置估计精度（排除NA）
            mean_position_error = mean(position_error, na.rm = TRUE),
            sd_position_error = sd(position_error, na.rm = TRUE),

            # 样本量
            n_sim = n(),

            .groups = "drop"
        )

    # 添加解释性列
    stats <- stats %>%
        mutate(
            # Size = 在H0下的拒绝率
            size = ifelse(scenario_type == "size", rejection_rate, NA),

            # Power = 在H1下的拒绝率
            power = ifelse(scenario_type == "power", rejection_rate, NA)
        )

    return(stats)
}

#' 打印结果摘要
#'
#' @param stats 统计性质汇总
print_summary <- function(stats) {
    cat("\n========================================\n")
    cat("统计性质汇总\n")
    cat("========================================\n")

    # Size测试结果
    size_results <- stats %>% filter(scenario_type == "size")

    if (nrow(size_results) > 0) {
        cat("\n【Size测试】（名义水平 α = 0.05）\n")
        cat("场景: ", size_results$scenario_name, "\n")
        cat("  经验Size: ", round(size_results$size, 4), "\n")
        cat(
            "  标准误: ",
            round(
                sqrt(
                    size_results$size *
                        (1 - size_results$size) /
                        size_results$n_sim
                ),
                4
            ),
            "\n"
        )

        # 判断Size是否合适（95%置信区间）
        se <- sqrt(0.05 * 0.95 / size_results$n_sim)
        ci_lower <- 0.05 - 1.96 * se
        ci_upper <- 0.05 + 1.96 * se

        if (size_results$size >= ci_lower && size_results$size <= ci_upper) {
            cat("  结论: ✓ Size控制良好（在95%置信区间内）\n")
        } else {
            cat("  结论: ✗ Size可能有偏（超出95%置信区间）\n")
        }
    }

    # Power测试结果
    power_results <- stats %>% filter(scenario_type == "power")

    if (nrow(power_results) > 0) {
        cat("\n【Power测试】\n")
        for (i in 1:nrow(power_results)) {
            row <- power_results[i, ]
            cat("\n场景:", row$scenario_name, "\n")
            cat("  Power: ", round(row$power, 4), "\n")
            cat("  正确检测率: ", round(row$correct_rate, 4), "\n")
            cat(
                "  平均检测断点数: ",
                round(row$mean_n_breaks, 2),
                " ± ",
                round(row$sd_n_breaks, 2),
                "\n"
            )
            cat(
                "  平均位置误差: ",
                round(row$mean_position_error, 2),
                " ± ",
                round(row$sd_position_error, 2),
                "\n"
            )
        }
    }

    cat("\n========================================\n")
}

#' ========== 可视化 ==========

#' 绘制Size和Power图
#'
#' @param stats 统计性质汇总
#' @return ggplot对象
plot_size_power <- function(stats) {
    # 准备数据
    size_data <- stats %>%
        filter(scenario_type == "size") %>%
        mutate(metric = "Size", value = size)

    power_data <- stats %>%
        filter(scenario_type == "power") %>%
        mutate(metric = "Power", value = power)

    plot_data <- rbind(
        data.frame(
            scenario_name = size_data$scenario_name,
            metric = size_data$metric,
            value = size_data$value
        ),
        data.frame(
            scenario_name = power_data$scenario_name,
            metric = power_data$metric,
            value = power_data$value
        )
    )

    # Size图
    p1 <- ggplot(size_data, aes(x = scenario_name, y = size)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
        geom_hline(
            yintercept = 0.05,
            color = "red",
            linetype = "dashed",
            linewidth = 1
        ) +
        geom_hline(
            yintercept = 0.05 - 1.96 * sqrt(0.05 * 0.95 / 1000),
            color = "red",
            linetype = "dotted"
        ) +
        geom_hline(
            yintercept = 0.05 + 1.96 * sqrt(0.05 * 0.95 / 1000),
            color = "red",
            linetype = "dotted"
        ) +
        labs(
            title = "Size测试（名义水平 α = 0.05）",
            x = NULL,
            y = "经验Size",
            caption = "红色实线：名义水平\n红色虚线：95%置信区间"
        ) +
        ylim(0, 0.15) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    # Power图
    p2 <- ggplot(
        power_data,
        aes(x = reorder(scenario_name, power), y = power)
    ) +
        geom_bar(stat = "identity", aes(fill = power), width = 0.7) +
        scale_fill_gradient2(
            low = "red",
            mid = "yellow",
            high = "darkgreen",
            midpoint = 0.5,
            limits = c(0, 1)
        ) +
        labs(
            title = "Power测试（不同场景）",
            x = NULL,
            y = "经验Power"
        ) +
        ylim(0, 1.05) +
        coord_flip() +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none"
        )

    # 组合图
    grid.arrange(p1, p2, ncol = 1, heights = c(1, 2))
}

#' 绘制断点数量分布
#'
#' @param results 模拟结果数据框
#' @return ggplot对象
plot_break_distribution <- function(results) {
    # 只显示部分关键场景
    selected_scenarios <- c(
        "Size: H0（无断点）",
        "Power: 单断点（中期，0.5T）",
        "Power: 双断点（0.33T, 0.67T）",
        "Power: 三断点（0.25T, 0.5T, 0.75T）"
    )

    plot_data <- results %>%
        filter(scenario_name %in% selected_scenarios) %>%
        mutate(
            scenario_name = factor(scenario_name, levels = selected_scenarios)
        )

    ggplot(plot_data, aes(x = n_breaks_detected)) +
        geom_histogram(
            binwidth = 1,
            fill = "steelblue",
            color = "white",
            alpha = 0.8
        ) +
        geom_vline(
            aes(xintercept = n_breaks_true),
            color = "red",
            linetype = "dashed",
            linewidth = 1
        ) +
        facet_wrap(~scenario_name, scales = "free_y", ncol = 2) +
        labs(
            title = "检测断点数量分布",
            x = "检测到的断点数",
            y = "频数",
            caption = "红色虚线：真实断点数"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.text = element_text(face = "bold")
        )
}

#' 绘制位置估计精度
#'
#' @param results 模拟结果数据框
#' @return ggplot对象
plot_position_accuracy <- function(results) {
    # 只显示有断点的场景
    plot_data <- results %>%
        filter(scenario_type == "power", !is.na(position_error))

    ggplot(plot_data, aes(x = scenario_name, y = position_error)) +
        geom_boxplot(fill = "steelblue", alpha = 0.7) +
        labs(
            title = "断点位置估计精度",
            x = NULL,
            y = "位置估计误差（时间点）"
        ) +
        coord_flip() +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
}

#' ========== 完整分析流程 ==========

#' 运行完整的蒙特卡洛验证
#'
#' @param n_sim 模拟次数
#' @param save_results 是否保存结果
#' @param output_dir 输出目录
#' @return 统计性质汇总
complete_monte_carlo_analysis <- function(
    n_sim = 1000,
    save_results = TRUE,
    output_dir = "results/monte_carlo"
) {
    # 创建输出目录
    if (save_results && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # 定义场景
    scenarios <- define_scenarios()

    # 运行模拟
    results <- run_monte_carlo(n_sim = n_sim, scenarios = scenarios)

    # 计算统计性质
    stats <- compute_statistics(results)

    # 打印摘要
    print_summary(stats)

    # 生成图形
    cat("\n生成可视化结果...\n")

    # Size和Power图
    png(file.path(output_dir, "size_power.png"), width = 800, height = 800)
    plot_size_power(stats)
    dev.off()

    # 断点数量分布
    p_dist <- plot_break_distribution(results)
    ggsave(
        file.path(output_dir, "break_distribution.png"),
        p_dist,
        width = 10,
        height = 8
    )

    # 位置估计精度
    p_accuracy <- plot_position_accuracy(results)
    ggsave(
        file.path(output_dir, "position_accuracy.png"),
        p_accuracy,
        width = 10,
        height = 8
    )

    # 保存结果
    if (save_results) {
        write.csv(
            results,
            file.path(output_dir, "simulation_results.csv"),
            row.names = FALSE
        )
        write.csv(
            stats,
            file.path(output_dir, "statistics_summary.csv"),
            row.names = FALSE
        )

        cat("\n结果已保存到:", output_dir, "\n")
    }

    return(stats)
}

#' ========== 使用示例 ==========

# 示例1: 快速测试（少量模拟）
# stats <- complete_monte_carlo_analysis(n_sim = 100)

# 示例2: 完整分析（推荐用于论文）
# stats <- complete_monte_carlo_analysis(n_sim = 1000)

# 示例3: 只运行特定场景
# scenarios <- define_scenarios()
# selected <- scenarios[c("size_null", "power_single_mid", "power_double")]
# results <- run_monte_carlo(n_sim = 500, scenarios = selected)
# stats <- compute_statistics(results)
# print_summary(stats)
