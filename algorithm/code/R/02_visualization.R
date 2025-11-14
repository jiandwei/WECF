# ============================================================================
# 文件名: 02_visualization.R
# 功能: 可视化函数型数据及断点检测结果
# ============================================================================

library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

#' 绘制函数型时间序列数据
#'
#' @param data functional_timeseries对象
#' @param highlight_breaks 是否高亮显示真实断点
#' @param n_curves 显示的曲线数量（NULL表示全部）
#' @param type 绘图类型 ("curves", "heatmap", "3d")

plot_functional_data <- function(
    data,
    highlight_breaks = TRUE,
    n_curves = 100,
    type = "curves"
) {
    T <- nrow(data$data)
    n_grid <- ncol(data$data)

    if (type == "curves") {
        # ===== 曲线图 =====

        # 选择要显示的曲线
        if (is.null(n_curves) || n_curves >= T) {
            curve_idx <- 1:T
        } else {
            curve_idx <- round(seq(1, T, length.out = n_curves))
        }

        # 准备数据
        plot_df <- data.frame()
        for (t in curve_idx) {
            plot_df <- rbind(
                plot_df,
                data.frame(
                    time = t,
                    s = data$s_grid,
                    value = data$data[t, ],
                    regime = as.factor(data$regime_labels[t])
                )
            )
        }

        # 绘图
        p <- ggplot(
            plot_df,
            aes(x = s, y = value, group = time, color = regime)
        ) +
            geom_line(alpha = 0.3, linewidth = 0.5) +
            scale_color_brewer(palette = "Set1", name = "Regime") +
            labs(
                title = "函数型时间序列",
                subtitle = sprintf(
                    "T = %d, 断点数 = %d",
                    T,
                    length(data$true_breaks)
                ),
                x = "函数域 s",
                y = "X(s)"
            ) +
            theme_minimal(base_size = 14) +
            theme(legend.position = "right")

        return(p)
    } else if (type == "heatmap") {
        # ===== 热图 =====

        plot_df <- melt(data$data)
        colnames(plot_df) <- c("time", "grid_idx", "value")
        plot_df$s <- data$s_grid[plot_df$grid_idx]

        p <- ggplot(plot_df, aes(x = s, y = time, fill = value)) +
            geom_tile() +
            scale_fill_gradient2(
                low = "blue",
                mid = "white",
                high = "red",
                midpoint = 0,
                name = "X(s)"
            ) +
            labs(title = "函数型时间序列热图", x = "函数域 s", y = "时间 t") +
            theme_minimal(base_size = 14)

        # 添加断点线
        if (highlight_breaks && length(data$true_break_indices) > 0) {
            p <- p +
                geom_hline(
                    yintercept = data$true_break_indices,
                    color = "black",
                    linetype = "dashed",
                    linewidth = 1
                )
        }

        return(p)
    } else if (type == "3d") {
        # ===== 3D图（使用plotly） =====

        if (!requireNamespace("plotly", quietly = TRUE)) {
            stop("需要安装plotly包")
        }

        z_matrix <- data$data

        p <- plotly::plot_ly(
            x = data$s_grid,
            y = 1:T,
            z = z_matrix,
            type = "surface",
            colorscale = "Viridis"
        ) %>%
            plotly::layout(
                title = "函数型时间序列3D视图",
                scene = list(
                    xaxis = list(title = "s"),
                    yaxis = list(title = "t"),
                    zaxis = list(title = "X(s)")
                )
            )

        return(p)
    }
}


#' 绘制断点检测结果
#'
#' @param data functional_timeseries对象
#' @param detected_breaks 检测到的断点位置向量
#' @param bootstrap_result bootstrap结果（可选）

plot_breakpoint_results <- function(
    data,
    detected_breaks,
    bootstrap_result = NULL
) {
    plots <- list()

    # ===== 图1：数据 + 真实断点 + 估计断点 =====

    T <- nrow(data$data)
    n_grid <- ncol(data$data)

    plot_df <- melt(data$data)
    colnames(plot_df) <- c("time", "grid_idx", "value")
    plot_df$s <- data$s_grid[plot_df$grid_idx]

    p1 <- ggplot(plot_df, aes(x = s, y = time, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red",
            midpoint = 0
        ) +
        # 真实断点（红色虚线）
        geom_hline(
            yintercept = data$true_break_indices,
            color = "red",
            linetype = "dashed",
            linewidth = 1,
            alpha = 0.8
        ) +
        # 估计断点（绿色实线）
        geom_hline(
            yintercept = detected_breaks,
            color = "green",
            linetype = "solid",
            linewidth = 1,
            alpha = 0.8
        ) +
        labs(
            title = "断点检测结果",
            subtitle = sprintf("真实断点（红虚线）vs 估计断点（绿实线）"),
            x = "函数域 s",
            y = "时间 t"
        ) +
        theme_minimal(base_size = 12)

    plots[[1]] <- p1

    # ===== 图2：断点位置对比 =====

    comparison_df <- data.frame(
        type = c(
            rep("真实", length(data$true_break_indices)),
            rep("估计", length(detected_breaks))
        ),
        position = c(data$true_break_indices, detected_breaks),
        relative = c(data$true_breaks, detected_breaks / T)
    )

    p2 <- ggplot(comparison_df, aes(x = relative, y = type, color = type)) +
        geom_point(size = 4) +
        geom_vline(
            xintercept = c(data$true_breaks, detected_breaks / T),
            linetype = "dotted",
            alpha = 0.3
        ) +
        scale_color_manual(values = c("真实" = "red", "估计" = "green")) +
        xlim(0, 1) +
        labs(title = "断点位置对比", x = "相对位置", y = "") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none")

    plots[[2]] <- p2

    # ===== 图3：Bootstrap分布（如果提供） =====

    if (!is.null(bootstrap_result)) {
        boot_df <- data.frame(
            statistic = bootstrap_result$distribution
        )

        p3 <- ggplot(boot_df, aes(x = statistic)) +
            geom_histogram(
                bins = 50,
                fill = "lightblue",
                color = "black",
                alpha = 0.7
            ) +
            geom_vline(
                xintercept = bootstrap_result$critical_value_005,
                color = "red",
                linetype = "dashed",
                linewidth = 1
            ) +
            annotate(
                "text",
                x = bootstrap_result$critical_value_005,
                y = Inf,
                vjust = 1.5,
                label = sprintf(
                    "CV (α=0.05) = %.2f",
                    bootstrap_result$critical_value_005
                ),
                color = "red"
            ) +
            labs(title = "Bootstrap分布", x = "检验统计量", y = "频数") +
            theme_minimal(base_size = 12)

        plots[[3]] <- p3
    }

    # ===== 图4：每个regime的平均函数 =====

    n_regimes <- length(data$true_breaks) + 1
    mean_df <- data.frame()

    for (j in 1:n_regimes) {
        regime_idx <- which(data$regime_labels == j)
        mean_curve <- colMeans(data$data[regime_idx, , drop = FALSE])

        mean_df <- rbind(
            mean_df,
            data.frame(
                s = data$s_grid,
                mean = mean_curve,
                regime = as.factor(j)
            )
        )
    }

    p4 <- ggplot(
        mean_df,
        aes(x = s, y = mean, color = regime, group = regime)
    ) +
        geom_line(linewidth = 1.5) +
        scale_color_brewer(palette = "Set1", name = "Regime") +
        labs(title = "各Regime的平均函数", x = "函数域 s", y = "E[X(s)]") +
        theme_minimal(base_size = 12)

    plots[[4]] <- p4

    # 组合所有图形
    do.call(grid.arrange, c(plots, ncol = 2))
}


#' 绘制诊断图
#'
#' @param data functional_timeseries对象
#' @param result 断点检测结果

plot_diagnostics <- function(data, result) {
    # 计算残差
    # TODO: 实现残差计算

    # 绘制ACF
    # TODO

    # 绘制QQ图
    # TODO

    message("诊断图功能待完善")
}
