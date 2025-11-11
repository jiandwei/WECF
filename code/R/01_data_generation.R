# ============================================
# 函数型数据生成模块
# ============================================

library(fda)
library(MASS)

#' 生成带断点的函数型数据
#'
#' @param T 样本量
#' @param grid_size 每条曲线的观测点数
#' @param break_points 断点位置（比例）
#' @param break_type 断点类型："mean", "variance", "distribution"
#' @param snr 信噪比
#' @return 包含数据矩阵和真实断点的列表
generate_functional_breaks <- function(
    T = 300,
    grid_size = 100,
    break_points = c(0.33, 0.67),
    break_type = "mean",
    snr = 3
) {
    # 时间网格
    t_grid <- seq(0, 1, length.out = grid_size)

    # 真实断点位置（样本索引）
    true_breaks <- floor(T * break_points)

    # 初始化数据矩阵
    X <- matrix(0, nrow = T, ncol = grid_size)

    # 定义Fourier basis
    nbasis <- 15
    basis <- create.fourier.basis(c(0, 1), nbasis = nbasis)

    # 生成不同regime的数据
    regimes <- c(1, true_breaks, T)
    n_regimes <- length(regimes) - 1

    for (j in 1:n_regimes) {
        start_idx <- regimes[j]
        end_idx <- regimes[j + 1]
        # number of observations in this regime (inclusive)
        n_obs <- end_idx - start_idx + 1

        if (break_type == "mean") {
            # 均值断点：不同regime有不同均值函数
            mean_shift <- (j - 1) * 2 # regime 1: 0, regime 2: 2, regime 3: 4

            # 生成系数
            coef_matrix <- matrix(
                rnorm(n_obs * nbasis, 0, 1 / sqrt(1:nbasis)),
                nrow = n_obs,
                ncol = nbasis
            )

            # 添加均值漂移（主要在第一个基函数上）
            coef_matrix[, 1] <- coef_matrix[, 1] + mean_shift
        } else if (break_type == "variance") {
            # 方差断点：不同regime有不同方差
            var_scale <- c(1, 2, 0.5)[j]

            coef_matrix <- matrix(
                rnorm(n_obs * nbasis, 0, var_scale / sqrt(1:nbasis)),
                nrow = n_obs,
                ncol = nbasis
            )
        } else if (break_type == "distribution") {
            # 分布断点：从正态切换到t分布
            if (j == 1) {
                coef_matrix <- matrix(
                    rnorm(n_obs * nbasis, 0, 1 / sqrt(1:nbasis)),
                    nrow = n_obs,
                    ncol = nbasis
                )
            } else {
                # 使用t分布（自由度=5）模拟厚尾
                coef_matrix <- matrix(
                    rt(n_obs * nbasis, df = 5) / sqrt(1:nbasis),
                    nrow = n_obs,
                    ncol = nbasis
                )
            }
        }

        # 转换为函数型数据
        fd_obj <- fd(t(coef_matrix), basis)
        X[start_idx:end_idx, ] <- t(eval.fd(t_grid, fd_obj))
    }

    # 添加测量误差
    noise_sd <- sd(X) / snr
    X <- X +
        matrix(rnorm(T * grid_size, 0, noise_sd), nrow = T, ncol = grid_size)

    return(list(
        X = X,
        t_grid = t_grid,
        true_breaks = true_breaks,
        regime_labels = rep(1:n_regimes, diff(c(1, true_breaks, T)))
    ))
}

#' 保存数据为CSV供C程序读取
save_data_for_c <- function(data, output_dir = "data") {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    # 保存数据矩阵
    write.table(
        data$X,
        file.path(output_dir, "functional_data.csv"),
        sep = ",",
        row.names = FALSE,
        col.names = FALSE
    )

    # 保存时间网格
    write.table(
        data$t_grid,
        file.path(output_dir, "time_grid.csv"),
        sep = ",",
        row.names = FALSE,
        col.names = FALSE
    )

    # 保存真实断点（用于后续验证）
    write.table(
        data$true_breaks,
        file.path(output_dir, "true_breaks.csv"),
        sep = ",",
        row.names = FALSE,
        col.names = FALSE
    )

    # 保存元数据
    metadata <- data.frame(
        T = nrow(data$X),
        grid_size = ncol(data$X),
        n_breaks = length(data$true_breaks)
    )
    write.csv(
        metadata,
        file.path(output_dir, "metadata.csv"),
        row.names = FALSE
    )

    cat("Data saved to", output_dir, "\n")
}

#' 可视化生成的数据
plot_functional_data <- function(data, max_curves = 50) {
    library(ggplot2)
    library(reshape2)

    # 随机选择max_curves条曲线展示
    n_curves <- min(max_curves, nrow(data$X))
    selected_curves <- sample(1:nrow(data$X), n_curves)

    # 转换为长格式
    df <- data.frame(
        curve_id = rep(selected_curves, each = ncol(data$X)),
        time = rep(data$t_grid, n_curves),
        value = as.vector(t(data$X[selected_curves, ])),
        regime = rep(data$regime_labels[selected_curves], each = ncol(data$X))
    )

    # 绘图
    p <- ggplot(
        df,
        aes(x = time, y = value, group = curve_id, color = factor(regime))
    ) +
        geom_line(alpha = 0.4, size = 0.5) +
        geom_vline(
            xintercept = data$t_grid[1] +
                (data$true_breaks / nrow(data$X)) *
                    (tail(data$t_grid, 1) - data$t_grid[1]),
            linetype = "dashed",
            color = "red",
            size = 1
        ) +
        labs(
            x = "Time",
            y = "X(t)",
            title = "Functional Data with Structural Breaks",
            color = "Regime"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")

    return(p)
}

# ============================================
# 示例：生成测试数据
# ============================================

if (TRUE) {
    # 测试不同类型的断点
    set.seed(42)

    # 1. 均值断点
    data_mean <- generate_functional_breaks(
        T = 300,
        grid_size = 100,
        break_points = c(0.33, 0.67),
        break_type = "mean",
        snr = 3
    )
    save_data_for_c(data_mean, "data/mean_break")
    p1 <- plot_functional_data(data_mean)
    ggsave("results/data_mean_break.png", p1, width = 10, height = 6)

    # 2. 方差断点
    data_var <- generate_functional_breaks(
        T = 300,
        grid_size = 100,
        break_points = c(0.4, 0.7),
        break_type = "variance",
        snr = 3
    )
    save_data_for_c(data_var, "data/variance_break")
    p2 <- plot_functional_data(data_var)
    ggsave("results/data_variance_break.png", p2, width = 10, height = 6)

    # 3. 分布断点
    data_dist <- generate_functional_breaks(
        T = 300,
        grid_size = 100,
        break_points = c(0.5),
        break_type = "distribution",
        snr = 3
    )
    save_data_for_c(data_dist, "data/distribution_break")
    p3 <- plot_functional_data(data_dist)
    ggsave("results/data_distribution_break.png", p3, width = 10, height = 6)
}
