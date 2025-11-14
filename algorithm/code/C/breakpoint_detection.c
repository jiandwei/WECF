/*
 * ============================================================================
 * 文件名: breakpoint_detection.c
 * 功能: Binary Segmentation for Functional Data Breakpoint Detection
 * 理论基础: Segmented Sum of Generalized Residuals (SSGR)
 * 编译: gcc -O3 -fopenmp -o break_detect breakpoint_detection.c utils.o -lm
 * 使用: ./break_detect <data_file> <meta_file> <output_file> <max_breaks> <min_segment_length>
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "utils.h" // 使用共享的工具函数

#define MAX_LINE 1024

/* ========== 断点检测特有的数据结构 ========== */

typedef struct
{
    int position;
    double test_stat;
    double p_value;
} Breakpoint;

typedef struct
{
    int n_breaks;
    Breakpoint *breaks;
} BreakpointResult;

/* ========== 断点检测函数声明 ========== */

double compute_test_statistic(FunctionalData *fd, int start, int end, int *best_position);
double compute_bic(FunctionalData *fd, int *breakpoints, int n_breaks);
BreakpointResult *binary_segmentation(FunctionalData *fd, int max_breaks, int min_segment_length);
void save_results(const char *output_file, BreakpointResult *result);
void free_result(BreakpointResult *result);

/* ========== 主函数 ========== */

int main(int argc, char *argv[])
{

    if (argc != 6)
    {
        fprintf(stderr, "用法: %s <data_file> <meta_file> <output_file> <max_breaks> <min_segment_length>\n", argv[0]);
        fprintf(stderr, "  max_breaks: 最大断点数量\n");
        fprintf(stderr, "  min_segment_length: 最小segment长度（推荐 >= 30）\n");
        return 1;
    }

    const char *data_file = argv[1];
    const char *meta_file = argv[2];
    const char *output_file = argv[3];
    int max_breaks = atoi(argv[4]);
    int min_segment_length = atoi(argv[5]);

    printf("================================================\n");
    printf("Binary Segmentation for Breakpoint Detection\n");
    printf("================================================\n");
    printf("最大断点数: %d\n", max_breaks);
    printf("最小segment长度: %d\n", min_segment_length);
    printf("OpenMP线程数: %d\n", omp_get_max_threads());

    // 加载数据
    printf("\n[1/3] 加载数据...\n");
    FunctionalData *fd = load_data(data_file, meta_file);
    if (fd == NULL)
    {
        return 1;
    }
    printf("  √ T = %d, n_grid = %d\n", fd->T, fd->n_grid);

    // 执行断点检测
    printf("\n[2/3] 执行Binary Segmentation...\n");
    double start_time = omp_get_wtime();

    BreakpointResult *result = binary_segmentation(fd, max_breaks, min_segment_length);

    double end_time = omp_get_wtime();
    printf("  √ 检测完成! 用时: %.2f 秒\n", end_time - start_time);
    printf("  检测到 %d 个断点\n", result->n_breaks);

    // 保存结果
    printf("\n[3/3] 保存结果...\n");
    save_results(output_file, result);

    // 清理
    free_result(result);
    free_data(fd);

    printf("\n================================================\n");
    printf("结果已保存到: %s\n", output_file);
    printf("================================================\n");

    return 0;
}

/* ========== Binary Segmentation实现 ========== */

/**
 * Binary Segmentation主算法
 *
 * 算法步骤:
 * 1. 初始化: 整个样本作为一个segment
 * 2. 循环:
 *    a) 对每个segment计算最优分割位置
 *    b) 选择检验统计量最大的segment进行分割
 *    c) 使用BIC准则决定是否接受新断点
 *    d) 直到达到最大断点数或无法找到显著断点
 */
BreakpointResult *binary_segmentation(FunctionalData *fd, int max_breaks, int min_segment_length)
{

    BreakpointResult *result = (BreakpointResult *)malloc(sizeof(BreakpointResult));
    result->breaks = (Breakpoint *)malloc(max_breaks * sizeof(Breakpoint));
    result->n_breaks = 0;

    // 维护当前的segments
    int *seg_starts = (int *)malloc((max_breaks + 1) * sizeof(int));
    int *seg_ends = (int *)malloc((max_breaks + 1) * sizeof(int));
    int n_segments = 1;

    seg_starts[0] = 0;
    seg_ends[0] = fd->T - 1;

    // 维护已找到的断点（用于BIC计算）
    int *current_breaks = (int *)malloc(max_breaks * sizeof(int));

    // Binary Segmentation循环
    for (int iter = 0; iter < max_breaks; iter++)
    {

        printf("  迭代 %d/%d...\n", iter + 1, max_breaks);

        double max_stat = 0.0;
        int best_seg_idx = -1;
        int best_position = -1;

// 对每个segment寻找最优分割点
#pragma omp parallel for
        for (int seg = 0; seg < n_segments; seg++)
        {

            int start = seg_starts[seg];
            int end = seg_ends[seg];

            // 检查segment长度是否足够分割
            if (end - start + 1 < 2 * min_segment_length)
            {
                continue;
            }

            // 计算该segment的最优分割点和检验统计量
            int position;
            double stat = compute_test_statistic(fd, start, end, &position);

// 更新全局最优
#pragma omp critical
            {
                if (stat > max_stat)
                {
                    max_stat = stat;
                    best_seg_idx = seg;
                    best_position = position;
                }
            }
        }

        // 如果没有找到有效分割，停止
        if (best_seg_idx == -1)
        {
            printf("    未找到有效分割点，停止\n");
            break;
        }

        // 临时添加新断点，计算BIC
        current_breaks[result->n_breaks] = best_position;
        double bic = compute_bic(fd, current_breaks, result->n_breaks + 1);

        // BIC准则：是否接受新断点
        double bic_old = (result->n_breaks > 0) ? compute_bic(fd, current_breaks, result->n_breaks) : compute_bic(fd, NULL, 0);

        if (bic > bic_old)
        {
            printf("    BIC增加，拒绝断点 (BIC: %.2f -> %.2f)\n", bic_old, bic);
            break;
        }

        // 接受新断点
        result->breaks[result->n_breaks].position = best_position;
        result->breaks[result->n_breaks].test_stat = max_stat;
        result->breaks[result->n_breaks].p_value = -1.0; // 待Bootstrap计算
        result->n_breaks++;

        printf("    √ 接受断点: position = %d, stat = %.4f\n", best_position, max_stat);

        // 更新segments
        int new_end = seg_ends[best_seg_idx];
        seg_ends[best_seg_idx] = best_position - 1;

        seg_starts[n_segments] = best_position;
        seg_ends[n_segments] = new_end;
        n_segments++;
    }

    // 对断点排序
    for (int i = 0; i < result->n_breaks - 1; i++)
    {
        for (int j = i + 1; j < result->n_breaks; j++)
        {
            if (result->breaks[i].position > result->breaks[j].position)
            {
                Breakpoint temp = result->breaks[i];
                result->breaks[i] = result->breaks[j];
                result->breaks[j] = temp;
            }
        }
    }

    // 清理
    free(seg_starts);
    free(seg_ends);
    free(current_breaks);

    return result;
}

/**
 * 计算给定segment的最优分割位置和检验统计量
 */
double compute_test_statistic(FunctionalData *fd, int start, int end, int *best_position)
{

    double epsilon = 0.1;
    int min_length = (int)(epsilon * (end - start + 1));
    if (min_length < 10)
        min_length = 10;

    double max_stat = 0.0;
    *best_position = -1;

// 并行搜索所有候选分割点
#pragma omp parallel for reduction(max : max_stat)
    for (int tau = start + min_length; tau <= end - min_length; tau++)
    {

        // 计算分割后的SSGR
        double **mean_left = compute_segment_mean(fd, start, tau - 1);
        double **mean_right = compute_segment_mean(fd, tau, end);

        double ssgr_left = compute_ssgr(fd, start, tau - 1, mean_left);
        double ssgr_right = compute_ssgr(fd, tau, end, mean_right);

        // 计算无分割的SSGR
        double **mean_full = compute_segment_mean(fd, start, end);
        double ssgr_full = compute_ssgr(fd, start, end, mean_full);

        // 检验统计量
        double stat = ssgr_full - (ssgr_left + ssgr_right);

#pragma omp critical
        {
            if (stat > max_stat)
            {
                max_stat = stat;
                *best_position = tau;
            }
        }

        // 清理
        free_segment_mean(mean_left, fd->n_grid);
        free_segment_mean(mean_right, fd->n_grid);
        free_segment_mean(mean_full, fd->n_grid);
    }

    return max_stat;
}

/**
 * 计算BIC准则
 * BIC = log(SSR) + k * log(T) / T
 * 其中 k = n_breaks + 1 (segment数量)
 */
double compute_bic(FunctionalData *fd, int *breakpoints, int n_breaks)
{

    // 创建segment边界
    int n_segments = n_breaks + 1;
    int *starts = (int *)malloc(n_segments * sizeof(int));
    int *ends = (int *)malloc(n_segments * sizeof(int));

    if (n_breaks == 0)
    {
        starts[0] = 0;
        ends[0] = fd->T - 1;
    }
    else
    {
        // 第一个segment
        starts[0] = 0;
        ends[0] = breakpoints[0] - 1;

        // 中间segments
        for (int i = 1; i < n_breaks; i++)
        {
            starts[i] = breakpoints[i - 1];
            ends[i] = breakpoints[i] - 1;
        }

        // 最后一个segment
        starts[n_breaks] = breakpoints[n_breaks - 1];
        ends[n_breaks] = fd->T - 1;
    }

    // 计算总SSR
    double total_ssgr = 0.0;

    for (int seg = 0; seg < n_segments; seg++)
    {
        double **mean = compute_segment_mean(fd, starts[seg], ends[seg]);
        double ssgr = compute_ssgr(fd, starts[seg], ends[seg], mean);
        total_ssgr += ssgr;
        free_segment_mean(mean, fd->n_grid);
    }

    // BIC
    double bic = log(total_ssgr / fd->T) + (double)n_segments * log((double)fd->T) / fd->T;

    free(starts);
    free(ends);

    return bic;
}

/**
 * 保存检测结果
 */
void save_results(const char *output_file, BreakpointResult *result)
{

    FILE *fp = fopen(output_file, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "错误: 无法创建输出文件\n");
        return;
    }

    fprintf(fp, "n_breaks=%d\n", result->n_breaks);
    fprintf(fp, "\nposition,test_stat,p_value\n");

    for (int i = 0; i < result->n_breaks; i++)
    {
        fprintf(fp, "%d,%.10f,%.10f\n",
                result->breaks[i].position,
                result->breaks[i].test_stat,
                result->breaks[i].p_value);
    }

    fclose(fp);
}

/**
 * 释放结果
 */
void free_result(BreakpointResult *result)
{
    if (result == NULL)
        return;
    free(result->breaks);
    free(result);
}
