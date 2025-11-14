# ============================================================================
# 文件名: 03_main_analysis.R
# 功能: 完整的分析流程，整合所有模块
# ============================================================================

source("R/01_data_generation.R")
source("R/02_visualization.R")

#' 调用C程序进行断点检测
#'
#' @param data functional_timeseries对象
#' @param method 方法 ("binary_seg", "bootstrap")
#' @param params 参数列表

run_breakpoint_detection <- function(
    data,
    method = "binary_seg",
    params = list(),
    prefix = "temp"
) {
    # 默认参数
    default_params <- list(
        max_breaks = 10,
        min_segment_length = 30,
        B = 999,
        block_length = ceiling(nrow(data$data)^(1 / 3))
    )
    params <- modifyList(default_params, params)

    # 保存数据
    files <- save_for_c(data, output_dir = "data", prefix = prefix)

    # 构建命令
    if (method == "binary_seg") {
        cmd <- sprintf(
            "./break_detect %s %s results/%s_breakpoints.txt %d %d",
            files$data_file,
            files$meta_file,
            prefix,
            params$max_breaks,
            params$min_segment_length
        )
    } else if (method == "bootstrap") {
        cmd <- sprintf(
            "./bootstrap %s %s results/%s_breakpoints.txt %d %d",
            files$data_file,
            files$meta_file,
            prefix,
            params$B,
            params$block_length
        )
    } else {
        stop("不支持的method")
    }

    cat("执行命令:", cmd, "\n")

    # 运行C程序
    system_result <- system(cmd, intern = TRUE)
    cat(paste(system_result, collapse = "\n"), "\n")

    # 读取结果
    if (method == "binary_seg") {
        result <- read_breakpoint_results(sprintf(
            "results/%s_breakpoints.txt",
            prefix
        ))
    } else if (method == "bootstrap") {
        result <- read_bootstrap_results(sprintf(
            "results/%s_breakpoints.txt.bootstrap",
            prefix
        ))
    }

    return(result)
}


#' 读取断点检测结果
read_breakpoint_results <- function(filename) {
    lines <- readLines(filename)

    # 解析断点数量
    n_breaks_line <- grep("^n_breaks=", lines, value = TRUE)
    n_breaks <- as.integer(sub("n_breaks=", "", n_breaks_line))

    # 读取断点详情
    header_line <- which(lines == "position,test_stat,p_value")
    breaks_start <- header_line + 1

    if (n_breaks > 0 && length(header_line) > 0) {
        breaks_df <- read.csv(
            text = paste(
                lines[header_line:(breaks_start + n_breaks - 1)],
                collapse = "\n"
            ),
            header = TRUE
        )
    } else {
        breaks_df <- data.frame(
            position = integer(),
            test_stat = numeric(),
            p_value = numeric()
        )
    }

    return(list(
        n_breaks = n_breaks,
        breaks = breaks_df
    ))
}


#' 读取Bootstrap结果
read_bootstrap_results <- function(filename) {
    # Bootstrap结果现在是断点+p值的格式，与breakpoint_results相同
    result <- read_breakpoint_results(filename)

    # 添加一个critical_value字段以保持兼容性
    result$critical_value_005 <- ifelse(
        nrow(result$breaks) > 0,
        min(result$breaks$test_stat[result$breaks$p_value <= 0.05]),
        NA
    )

    return(result)
}


#' 完整分析流程
#'
#' @param T 样本量
#' @param breakpoints 真实断点位置
#' @param run_bootstrap 是否运行bootstrap（耗时）
#' @param visualize 是否生成可视化

complete_analysis <- function(
    T = 500,
    breakpoints = c(0.3, 0.7),
    run_bootstrap = TRUE,
    visualize = TRUE
) {
    cat("======================================\n")
    cat("函数型数据断点检测 - 完整分析流程\n")
    cat("======================================\n\n")

    # ===== 步骤1：生成数据 =====
    cat("[1/4] 生成模拟数据...\n")

    set.seed(2024)
    fts_data <- generate_functional_data(
        T = T,
        breakpoints = breakpoints,
        dependence = "AR1",
        rho = 0.5
    )

    cat(sprintf("  √ 生成 %d 个观测，%d 个断点\n", T, length(breakpoints)))

    # ===== 步骤2：断点检测 =====
    cat("\n[2/4] 执行断点检测...\n")

    break_result <- run_breakpoint_detection(
        fts_data,
        method = "binary_seg",
        params = list(max_breaks = 10, min_segment_length = 30)
    )

    cat(sprintf("  √ 检测到 %d 个断点\n", break_result$n_breaks))

    if (break_result$n_breaks > 0) {
        cat("  断点位置:\n")
        print(break_result$breaks)
    }

    # ===== 步骤3：Bootstrap推断 =====
    boot_result <- NULL

    if (run_bootstrap) {
        cat("\n[3/4] 执行Bootstrap（这可能需要几分钟）...\n")

        boot_result <- run_breakpoint_detection(
            fts_data,
            method = "bootstrap",
            params = list(B = 999, block_length = ceiling(T^(1 / 3)))
        )

        cat(sprintf(
            "  √ Bootstrap完成，检测到 %d 个显著断点 (p < 0.05)\n",
            sum(boot_result$breaks$p_value < 0.05, na.rm = TRUE)
        ))
        if (nrow(boot_result$breaks) > 0) {
            print(boot_result$breaks)
        }
    } else {
        cat("\n[3/4] 跳过Bootstrap\n")
    }

    # ===== 步骤4：可视化 =====
    if (visualize) {
        cat("\n[4/4] 生成可视化...\n")

        # 图1：原始数据
        p1 <- plot_functional_data(fts_data, type = "heatmap")
        print(p1)

        # 图2：检测结果
        plot_breakpoint_results(
            fts_data,
            break_result$breaks$position,
            boot_result
        )

        cat("  √ 可视化完成\n")
    } else {
        cat("\n[4/4] 跳过可视化\n")
    }

    # ===== 返回结果 =====
    cat("\n======================================\n")
    cat("分析完成！\n")
    cat("======================================\n")

    return(list(
        data = fts_data,
        breakpoints = break_result,
        bootstrap = boot_result
    ))
}


# ============================================================================
# 使用示例
# ============================================================================

if (FALSE) {
    # 运行完整分析
    results <- complete_analysis(
        T = 500,
        breakpoints = c(0.3, 0.7),
        run_bootstrap = TRUE,
        visualize = TRUE
    )

    # 查看结果
    print(results$breakpoints)

    # 保存结果
    saveRDS(results, "results/analysis_results.rds")
}
