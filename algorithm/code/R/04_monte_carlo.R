# ============================================================================
# 文件名: 04_monte_carlo.R
# 功能: Monte Carlo 模拟工具函数
# ============================================================================

# 假定本文件被 source 之前，01-03 脚本已加载（提供 generate_functional_data 等）

`%||%` <- function(x, y) {
    if (!is.null(x)) x else y
}

#' 对字符串做文件安全化处理
sanitize_name <- function(name) {
    gsub("[^A-Za-z0-9_-]", "_", name)
}

#' 计算断点检测指标
#'
#' @param true_idx 真实断点索引向量
#' @param detected_idx 检测到的断点索引向量
#' @param tolerance 匹配容忍度（以观测为单位）
#'
#' @return 列表，包含 tp/fp/fn 与平均误差
compute_break_metrics <- function(true_idx, detected_idx, tolerance = 25) {
    true_idx <- sort(as.integer(true_idx))
    detected_idx <- sort(as.integer(detected_idx))

    if (length(true_idx) == 0 && length(detected_idx) == 0) {
        return(list(tp = 0L, fp = 0L, fn = 0L, mean_abs_error = NA_real_))
    }

    used_detected <- rep(FALSE, length(detected_idx))
    abs_errors <- numeric(0)
    tp <- 0L

    for (t_break in true_idx) {
        if (length(detected_idx) == 0) {
            next
        }
        diffs <- abs(detected_idx - t_break)
        within_tol <- which(diffs <= tolerance & !used_detected)
        if (length(within_tol) > 0) {
            best_idx <- within_tol[which.min(diffs[within_tol])]
            used_detected[best_idx] <- TRUE
            abs_errors <- c(abs_errors, diffs[best_idx])
            tp <- tp + 1L
        }
    }

    fp <- sum(!used_detected)
    fn <- length(true_idx) - tp
    mean_abs_error <- if (length(abs_errors) > 0) mean(abs_errors) else NA_real_

    list(tp = tp, fp = fp, fn = fn, mean_abs_error = mean_abs_error)
}

#' 运行单次 Monte Carlo replicate
run_single_rep <- function(
    rep_id,
    scenario,
    detection_params,
    bootstrap_params,
    tolerance,
    run_bootstrap,
    base_seed
) {
    scenario_name <- scenario$name %||% sprintf("scenario_%03d", rep_id)
    scenario_args <- scenario
    scenario_args$name <- NULL
    scenario_args$seed <- base_seed + rep_id

    # 生成数据
    fts_data <- do.call(generate_functional_data, scenario_args)

    prefix <- sprintf("%s_rep%04d", sanitize_name(scenario_name), rep_id)

    det_time <- system.time({
        break_result <- run_breakpoint_detection(
            fts_data,
            method = "binary_seg",
            params = detection_params,
            prefix = prefix
        )
    })["elapsed"]

    metrics <- compute_break_metrics(
        fts_data$true_break_indices,
        break_result$breaks$position,
        tolerance
    )

    boot_time <- NA_real_
    sig_breaks <- NA_integer_
    if (run_bootstrap && break_result$n_breaks > 0) {
        boot_time <- system.time({
            boot_result <- run_breakpoint_detection(
                fts_data,
                method = "bootstrap",
                params = bootstrap_params,
                prefix = prefix
            )
        })["elapsed"]
        sig_breaks <- sum(boot_result$breaks$p_value < 0.05, na.rm = TRUE)
    }

    list(
        scenario = scenario_name,
        replicate = rep_id,
        T = nrow(fts_data$data),
        dependence = scenario_args$dependence %||% "AR1",
        rho = scenario_args$rho %||% NA_real_,
        true_breaks = paste(fts_data$true_break_indices, collapse = ";"),
        n_true = length(fts_data$true_break_indices),
        n_detected = break_result$n_breaks,
        tp = metrics$tp,
        fp = metrics$fp,
        fn = metrics$fn,
        precision = if ((metrics$tp + metrics$fp) > 0) {
            metrics$tp / (metrics$tp + metrics$fp)
        } else {
            NA_real_
        },
        recall = if ((metrics$tp + metrics$fn) > 0) {
            metrics$tp / (metrics$tp + metrics$fn)
        } else {
            NA_real_
        },
        mean_abs_error = metrics$mean_abs_error,
        detection_time = as.numeric(det_time),
        bootstrap_time = as.numeric(boot_time),
        significant_breaks = sig_breaks,
        tolerance = tolerance
    )
}

#' 主 Monte Carlo 驱动
run_monte_carlo <- function(
    scenarios,
    n_rep = 100,
    detection_params = list(max_breaks = 5, min_segment_length = 40),
    bootstrap_params = list(B = 199),
    tolerance = 25,
    run_bootstrap = TRUE,
    base_seed = 2024,
    output_dir = "results/monte_carlo"
) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    results <- list()
    idx <- 1L

    for (scenario in scenarios) {
        cat("\n==============================\n")
        cat(sprintf("Scenario: %s\n", scenario$name))
        cat("==============================\n")

        for (rep_id in seq_len(n_rep)) {
            cat(sprintf(
                "  [Scenario %s] Rep %d/%d\n",
                scenario$name,
                rep_id,
                n_rep
            ))
            res <- run_single_rep(
                rep_id = rep_id,
                scenario = scenario,
                detection_params = detection_params,
                bootstrap_params = bootstrap_params,
                tolerance = tolerance,
                run_bootstrap = run_bootstrap,
                base_seed = base_seed
            )
            results[[idx]] <- res
            idx <- idx + 1L
        }
    }

    df_list <- lapply(results, function(res) {
        if (is.null(res)) {
            return(NULL)
        }
        as.data.frame(res, stringsAsFactors = FALSE)
    })
    df_list <- Filter(Negate(is.null), df_list)

    if (length(df_list) == 0) {
        stop("未生成 Monte Carlo 结果，请检查上游流程。")
    }

    results_df <- do.call(
        rbind.data.frame,
        c(df_list, stringsAsFactors = FALSE, make.row.names = FALSE)
    )

    summary_df <- aggregate(
        cbind(
            tp,
            fp,
            fn,
            n_detected,
            n_true,
            detection_time,
            bootstrap_time,
            mean_abs_error,
            significant_breaks
        ) ~ scenario,
        data = results_df,
        FUN = function(x) mean(x, na.rm = TRUE),
        na.action = na.pass
    )
    summary_df$precision <- with(summary_df, tp / (tp + fp))
    summary_df$recall <- with(summary_df, tp / (tp + fn))

    detailed_path <- file.path(output_dir, "monte_carlo_results.csv")
    summary_path <- file.path(output_dir, "monte_carlo_summary.csv")

    write.csv(results_df, detailed_path, row.names = FALSE)
    write.csv(summary_df, summary_path, row.names = FALSE)

    list(
        detailed = results_df,
        summary = summary_df,
        files = list(detailed = detailed_path, summary = summary_path)
    )
}
