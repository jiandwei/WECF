## 第一部分：R数据生成模块

# ============================================================================
# 文件名: 01_data_generation.R
# 功能: 生成具有结构断点的函数型时间序列数据
# 作者: [Your Name]
# 日期: 2024
# ============================================================================

#' 生成函数型时间序列数据（带断点）
#'
#' @description
#' 该函数生成满足以下DGP的函数型时间序列：
#'   X_t(s) = μ_j(s) + σ_j(s) * ε_t(s)  for t in regime j
#' 其中断点位置由breakpoints参数指定
#'
#' @param T 样本量
#' @param s_grid 函数定义域的离散网格 (默认: [0,1]上的101个点)
#' @param breakpoints 断点位置向量 (比例形式, 如 c(0.3, 0.7))
#' @param mean_functions 均值函数列表 (长度 = length(breakpoints) + 1)
#' @param cov_functions 协方差函数列表
#' @param dependence 时间依赖性类型 ("AR1", "MA1", "independent")
#' @param rho 依赖性参数 (for AR1/MA1)
#' @param seed 随机种子
#'
#' @return 列表包含:
#'   - data: T x length(s_grid) 矩阵
#'   - s_grid: 函数网格
#'   - true_breaks: 真实断点位置
#'   - regime_labels: 每个观测的regime标签
#'
#' @examples
#' # 生成两个断点的数据
#' data <- generate_functional_data(
#'   T = 500,
#'   breakpoints = c(0.3, 0.7),
#'   mean_functions = list(
#'     function(s) sin(2*pi*s),
#'     function(s) cos(2*pi*s),
#'     function(s) s^2
#'   )
#' )

generate_functional_data <- function(
    T = 500,
    s_grid = seq(0, 1, length.out = 101),
    breakpoints = c(0.3, 0.7),
    mean_functions = NULL,
    cov_functions = NULL,
    dependence = "AR1",
    rho = 0.5,
    seed = 123
) {
    set.seed(seed)

    # ========== 参数检查 ==========
    if (any(breakpoints <= 0 | breakpoints >= 1)) {
        stop("断点必须在(0, 1)区间内")
    }
    breakpoints <- sort(breakpoints)
    M <- length(breakpoints) # 断点数量
    n_regimes <- M + 1 # regime数量

    # ========== 默认均值函数 ==========
    if (is.null(mean_functions)) {
        # 提供足够的默认函数
        default_functions <- list(
            function(s) sin(2 * pi * s), # Regime 1
            function(s) cos(2 * pi * s), # Regime 2
            function(s) 2 * s * (1 - s), # Regime 3
            function(s) sin(4 * pi * s), # Regime 4
            function(s) exp(-10 * (s - 0.5)^2), # Regime 5
            function(s) s^2 - s^3 # Regime 6
        )

        if (n_regimes > length(default_functions)) {
            stop(sprintf(
                "需要 %d 个regime，但只定义了 %d 个默认函数。请提供自定义mean_functions。",
                n_regimes,
                length(default_functions)
            ))
        }

        mean_functions <- default_functions[1:n_regimes]
    }

    # ========== 默认协方差函数 ==========
    if (is.null(cov_functions)) {
        # 使用Gaussian过程协方差: C(s,t) = σ^2 * exp(-|s-t|^2 / l^2)
        cov_functions <- lapply(1:n_regimes, function(j) {
            sigma2 <- 0.5 + 0.2 * j # 不同regime的方差
            length_scale <- 0.3
            function(s, t) {
                sigma2 * exp(-((s - t)^2) / (2 * length_scale^2))
            }
        })
    }

    # ========== 计算regime边界 ==========
    break_indices <- round(T * breakpoints)
    regime_bounds <- c(0, break_indices, T)

    # ========== 初始化数据矩阵 ==========
    n_grid <- length(s_grid)
    X <- matrix(NA, nrow = T, ncol = n_grid)
    regime_labels <- integer(T)

    # ========== 逐regime生成数据 ==========
    for (j in 1:n_regimes) {
        cat(sprintf(
            "生成 Regime %d (观测 %d 到 %d)...\n",
            j,
            regime_bounds[j] + 1,
            regime_bounds[j + 1]
        ))

        # 当前regime的时间索引
        t_start <- regime_bounds[j] + 1
        t_end <- regime_bounds[j + 1]
        n_obs <- t_end - t_start + 1

        # 标记regime
        regime_labels[t_start:t_end] <- j

        # 计算均值向量
        mu_j <- sapply(s_grid, mean_functions[[j]])

        # 构建协方差矩阵
        Sigma_j <- outer(s_grid, s_grid, Vectorize(cov_functions[[j]]))

        # 确保协方差矩阵正定
        Sigma_j <- (Sigma_j + t(Sigma_j)) / 2
        eigen_vals <- eigen(
            Sigma_j,
            symmetric = TRUE,
            only.values = TRUE
        )$values
        if (min(eigen_vals) < 1e-6) {
            Sigma_j <- Sigma_j + diag(1e-6, n_grid)
        }

        # 生成创新 ε_t(s)
        if (dependence == "independent") {
            # 独立同分布
            epsilon <- MASS::mvrnorm(
                n = n_obs,
                mu = rep(0, n_grid),
                Sigma = Sigma_j
            )
        } else if (dependence == "AR1") {
            # AR(1)过程: ε_t = ρ * ε_{t-1} + η_t
            epsilon <- matrix(0, nrow = n_obs, ncol = n_grid)
            epsilon[1, ] <- MASS::mvrnorm(
                n = 1,
                mu = rep(0, n_grid),
                Sigma = Sigma_j / (1 - rho^2)
            )

            for (t in 2:n_obs) {
                innovation <- MASS::mvrnorm(
                    n = 1,
                    mu = rep(0, n_grid),
                    Sigma = Sigma_j
                )
                epsilon[t, ] <- rho *
                    epsilon[t - 1, ] +
                    sqrt(1 - rho^2) * innovation
            }
        } else if (dependence == "MA1") {
            # MA(1)过程: ε_t = η_t + θ * η_{t-1}
            theta <- rho
            innovations <- MASS::mvrnorm(
                n = n_obs + 1,
                mu = rep(0, n_grid),
                Sigma = Sigma_j
            )
            epsilon <- innovations[1:n_obs, ] +
                theta * innovations[2:(n_obs + 1), ]
        } else {
            stop("不支持的dependence类型")
        }

        # 生成观测数据: X_t(s) = μ_j(s) + ε_t(s)
        X[t_start:t_end, ] <- sweep(epsilon, 2, mu_j, "+")
    }

    # ========== 返回结果 ==========
    result <- list(
        data = X,
        s_grid = s_grid,
        true_breaks = breakpoints,
        true_break_indices = break_indices,
        regime_labels = regime_labels,
        mean_functions = mean_functions,
        cov_functions = cov_functions,
        parameters = list(
            T = T,
            M = M,
            dependence = dependence,
            rho = rho
        )
    )

    class(result) <- "functional_timeseries"

    return(result)
}


#' 将函数型数据保存为C程序可读的格式
#'
#' @param data functional_timeseries对象
#' @param output_dir 输出目录
#' @param prefix 文件名前缀
#'
#' @return 写入的文件路径列表

save_for_c <- function(
    data,
    output_dir = "WECF\\algorithm\\code\\data",
    prefix = "fts"
) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # ========== 保存主数据矩阵 (二进制格式，提高读取速度) ==========
    data_file <- file.path(output_dir, paste0(prefix, "_data.bin"))
    writeBin(as.vector(t(data$data)), data_file, size = 8) # double precision

    # ========== 保存元数据 (文本格式) ==========
    meta_file <- file.path(output_dir, paste0(prefix, "_meta.txt"))
    writeLines(
        c(
            sprintf("T=%d", nrow(data$data)),
            sprintf("n_grid=%d", ncol(data$data)),
            sprintf("M=%d", length(data$true_breaks)),
            sprintf(
                "breaks=%s",
                paste(data$true_break_indices, collapse = ",")
            ),
            sprintf("s_min=%.6f", min(data$s_grid)),
            sprintf("s_max=%.6f", max(data$s_grid))
        ),
        meta_file
    )

    # ========== 保存s网格 ==========
    grid_file <- file.path(output_dir, paste0(prefix, "_grid.txt"))
    write.table(data$s_grid, grid_file, row.names = FALSE, col.names = FALSE)

    cat(sprintf(
        "数据已保存:\n  - %s\n  - %s\n  - %s\n",
        data_file,
        meta_file,
        grid_file
    ))

    return(list(
        data_file = data_file,
        meta_file = meta_file,
        grid_file = grid_file
    ))
}


# ============================================================================
# 示例使用
# ============================================================================

if (FALSE) {
    # 设置为TRUE以运行示例

    # 生成数据
    set.seed(2024)
    fts_data <- generate_functional_data(
        T = 500,
        breakpoints = c(0.3, 0.7),
        dependence = "AR1",
        rho = 0.5
    )

    # 保存数据
    files <- save_for_c(fts_data, output_dir = "data", prefix = "example")

    # 快速可视化
    matplot(
        fts_data$s_grid,
        t(fts_data$data[1:50, ]),
        type = "l",
        lty = 1,
        col = rgb(0, 0, 0, 0.3),
        xlab = "s",
        ylab = "X(s)",
        main = "前50条函数型观测"
    )
    abline(v = fts_data$true_breaks, col = "red", lwd = 2, lty = 2)
}
