/*
 * ============================================================================
 * 文件名: bootstrap.c
 * 功能: Moving Block Bootstrap for Critical Value Estimation
 * 理论基础: Politis & Romano (1994) - Limit Theorems for Block Bootstrap
 * 编译: gcc -O3 -fopenmp -o bootstrap bootstrap.c -lm
 * 使用: ./bootstrap <data_file> <meta_file> <output_file> <B> <block_length>
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define MAX_LINE 1024

/* ========== 数据结构 ========== */

typedef struct
{
    double **data;
    int T;
    int n_grid;
    double s_min;
    double s_max;
    double *s_grid;
} FunctionalData;

/* ========== 外部函数声明 ========== */
extern FunctionalData *load_data(const char *data_file, const char *meta_file);
extern void free_data(FunctionalData *fd);
extern double compute_ssgr(FunctionalData *fd, int start, int end, double **segment_mean);
extern double **compute_segment_mean(FunctionalData *fd, int start, int end);
extern void free_segment_mean(double **mean, int n_grid);

/* ========== Bootstrap函数声明 ========== */

FunctionalData *generate_bootstrap_sample(FunctionalData *original, int block_length, unsigned int seed);
double compute_test_statistic(FunctionalData *fd);
void moving_block_bootstrap(FunctionalData *fd, int B, int block_length,
                            double *bootstrap_stats, double *critical_values);
void save_bootstrap_results(const char *output_file, double *bootstrap_stats,
                            double *critical_values, int B);

/* ========== 主函数 ========== */

int main(int argc, char *argv[])
{

    if (argc != 6)
    {
        fprintf(stderr, "用法: %s <data_file> <meta_file> <output_file> <B> <block_length>\n", argv[0]);
        fprintf(stderr, "  B: Bootstrap重复次数 (推荐 >= 999)\n");
        fprintf(stderr, "  block_length: Block长度 (推荐 T^{1/3})\n");
        return 1;
    }

    const char *data_file = argv[1];
    const char *meta_file = argv[2];
    const char *output_file = argv[3];
    int B = atoi(argv[4]);
    int block_length = atoi(argv[5]);

    printf("================================================\n");
    printf("Moving Block Bootstrap\n");
    printf("================================================\n");
    printf("Bootstrap重复次数: %d\n", B);
    printf("Block长度: %d\n", block_length);
    printf("OpenMP线程数: %d\n", omp_get_max_threads());

    // 加载数据
    printf("\n[1/3] 加载数据...\n");
    FunctionalData *fd = load_data(data_file, meta_file);
    if (fd == NULL)
    {
        return 1;
    }
    printf("  √ T = %d\n", fd->T);

    // 验证block长度
    if (block_length < 2 || block_length > fd->T / 2)
    {
        fprintf(stderr, "错误: Block长度必须在 [2, T/2] 范围内\n");
        free_data(fd);
        return 1;
    }

    // 分配内存存储bootstrap统计量
    double *bootstrap_stats = (double *)malloc(B * sizeof(double));
    double *critical_values = (double *)malloc(10 * sizeof(double)); // 存储不同α的临界值

    // 执行bootstrap
    printf("\n[2/3] 执行Moving Block Bootstrap...\n");
    double start_time = omp_get_wtime();

    moving_block_bootstrap(fd, B, block_length, bootstrap_stats, critical_values);

    double end_time = omp_get_wtime();
    printf("  √ Bootstrap完成! 用时: %.2f 秒\n", end_time - start_time);

    // 保存结果
    printf("\n[3/3] 保存结果...\n");
    save_bootstrap_results(output_file, bootstrap_stats, critical_values, B);

    // 清理
    free(bootstrap_stats);
    free(critical_values);
    free_data(fd);

    printf("\n================================================\n");
    printf("结果已保存到: %s\n", output_file);
    printf("================================================\n");

    return 0;
}

/* ========== Bootstrap实现 ========== */

/**
 * Moving Block Bootstrap主函数
 *
 * 算法步骤:
 * 1. 对b = 1, ..., B:
 *    a) 生成bootstrap样本 X*_b
 *    b) 计算检验统计量 T*_b
 * 2. 计算临界值: CV_α = quantile(T*_b, 1-α)
 */
void moving_block_bootstrap(FunctionalData *fd, int B, int block_length,
                            double *bootstrap_stats, double *critical_values)
{

    // 计算原始数据的检验统计量（用于对比）
    double original_stat = compute_test_statistic(fd);
    printf("  原始数据检验统计量: %.6f\n", original_stat);

    // 并行执行bootstrap
    printf("  开始bootstrap循环...\n");

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; b++)
    {

        // 生成随机种子（线程安全）
        unsigned int seed = (unsigned int)time(NULL) + b + omp_get_thread_num() * 1000;

        // 生成bootstrap样本
        FunctionalData *boot_sample = generate_bootstrap_sample(fd, block_length, seed);

        // 计算检验统计量
        bootstrap_stats[b] = compute_test_statistic(boot_sample);

        // 释放bootstrap样本
        free_data(boot_sample);

        // 进度显示（每10%显示一次）
        if ((b + 1) % (B / 10) == 0)
        {
#pragma omp critical
            {
                printf("    进度: %d/%d (%.1f%%)\n", b + 1, B, 100.0 * (b + 1) / B);
            }
        }
    }

    // 对bootstrap统计量排序（用于计算分位数）
    printf("  计算临界值...\n");

    // 简单的冒泡排序（对小规模数据足够）
    for (int i = 0; i < B - 1; i++)
    {
        for (int j = 0; j < B - i - 1; j++)
        {
            if (bootstrap_stats[j] > bootstrap_stats[j + 1])
            {
                double temp = bootstrap_stats[j];
                bootstrap_stats[j] = bootstrap_stats[j + 1];
                bootstrap_stats[j + 1] = temp;
            }
        }
    }

    // 计算不同显著性水平的临界值
    double alpha_levels[10] = {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};

    for (int i = 0; i < 10; i++)
    {
        int idx = (int)((1.0 - alpha_levels[i]) * B);
        if (idx >= B)
            idx = B - 1;
        critical_values[i] = bootstrap_stats[idx];

        printf("    α = %.3f: CV = %.6f\n", alpha_levels[i], critical_values[i]);
    }

    // 计算p值
    int count_exceed = 0;
    for (int b = 0; b < B; b++)
    {
        if (bootstrap_stats[b] >= original_stat)
        {
            count_exceed++;
        }
    }
    double p_value = (double)count_exceed / B;
    printf("  p值: %.6f\n", p_value);
}

/**
 * 生成Moving Block Bootstrap样本
 *
 * 算法:
 * 1. 计算需要的block数: n_blocks = ceil(T / block_length)
 * 2. 对每个block，随机选择起始位置
 * 3. 复制连续的block_length个观测
 */
FunctionalData *generate_bootstrap_sample(FunctionalData *original,
                                          int block_length,
                                          unsigned int seed)
{

    // 初始化bootstrap样本
    FunctionalData *boot = (FunctionalData *)malloc(sizeof(FunctionalData));
    boot->T = original->T;
    boot->n_grid = original->n_grid;
    boot->s_min = original->s_min;
    boot->s_max = original->s_max;
    boot->s_grid = original->s_grid; // 共享网格

    // 分配内存
    boot->data = (double **)malloc(boot->T * sizeof(double *));
    for (int t = 0; t < boot->T; t++)
    {
        boot->data[t] = (double *)malloc(boot->n_grid * sizeof(double));
    }

    // 设置随机数生成器
    srand(seed);

    // 计算需要的block数
    int n_blocks = (original->T + block_length - 1) / block_length;

    // 生成bootstrap样本
    int boot_idx = 0;

    for (int block = 0; block < n_blocks; block++)
    {

        // 随机选择block起始位置
        int max_start = original->T - block_length;
        if (max_start < 0)
            max_start = 0;

        int block_start = rand() % (max_start + 1);

        // 复制block
        for (int i = 0; i < block_length && boot_idx < boot->T; i++)
        {
            int orig_idx = block_start + i;
            if (orig_idx >= original->T)
                break;

            memcpy(boot->data[boot_idx], original->data[orig_idx],
                   boot->n_grid * sizeof(double));
            boot_idx++;
        }
    }

    return boot;
}

/**
 * 计算检验统计量
 *
 * 这里使用supF统计量:
 * T = max_{r ∈ [ε, 1-ε]} F_T(r)
 * 其中 F_T(r) 基于SSGR
 */
double compute_test_statistic(FunctionalData *fd)
{

    double epsilon = 0.1; // Trimming参数
    int start_idx = (int)(epsilon * fd->T);
    int end_idx = (int)((1 - epsilon) * fd->T);

    double max_stat = 0.0;

// 搜索所有候选断点位置
#pragma omp parallel for reduction(max : max_stat)
    for (int tau = start_idx; tau <= end_idx; tau++)
    {

        // 计算分割后的SSGR
        double **mean_left = compute_segment_mean(fd, 0, tau);
        double **mean_right = compute_segment_mean(fd, tau + 1, fd->T - 1);

        double ssgr_left = compute_ssgr(fd, 0, tau, mean_left);
        double ssgr_right = compute_ssgr(fd, tau + 1, fd->T - 1, mean_right);

        // 计算无断点的SSGR
        double **mean_full = compute_segment_mean(fd, 0, fd->T - 1);
        double ssgr_full = compute_ssgr(fd, 0, fd->T - 1, mean_full);

        // F统计量
        double F_stat = (ssgr_full - (ssgr_left + ssgr_right)) * tau * (fd->T - tau) / fd->T;

        if (F_stat > max_stat)
        {
            max_stat = F_stat;
        }

        // 清理
        free_segment_mean(mean_left, fd->n_grid);
        free_segment_mean(mean_right, fd->n_grid);
        free_segment_mean(mean_full, fd->n_grid);
    }

    return max_stat;
}

/**
 * 保存bootstrap结果
 */
void save_bootstrap_results(const char *output_file, double *bootstrap_stats,
                            double *critical_values, int B)
{

    FILE *fp = fopen(output_file, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "错误: 无法创建输出文件\n");
        return;
    }

    // 写入临界值
    fprintf(fp, "# Critical Values\n");
    fprintf(fp, "alpha,critical_value\n");

    double alpha_levels[10] = {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
    for (int i = 0; i < 10; i++)
    {
        fprintf(fp, "%.3f,%.10f\n", alpha_levels[i], critical_values[i]);
    }

    // 写入bootstrap分布
    fprintf(fp, "\n# Bootstrap Distribution\n");
    fprintf(fp, "iteration,statistic\n");
    for (int b = 0; b < B; b++)
    {
        fprintf(fp, "%d,%.10f\n", b + 1, bootstrap_stats[b]);
    }

    fclose(fp);
}
