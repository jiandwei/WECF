/*
 * ============================================================================
 * 文件名: bootstrap.c
 * 功能: Moving Block Bootstrap for p-value Calculation
 * 理论基础: Carlstein (1986) Moving Block Bootstrap
 * 编译: gcc -O3 -fopenmp -o bootstrap bootstrap.c utils.o -lm
 * 使用: ./bootstrap <data_file> <meta_file> <breakpoint_file> <n_bootstrap> <block_length>
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
#define MAX_BREAKS 50

/* ========== Bootstrap特有的数据结构 ========== */

typedef struct
{
    int n_breaks;
    int positions[MAX_BREAKS];
    double test_stats[MAX_BREAKS];
    double p_values[MAX_BREAKS];
} BootstrapResult;

/* ========== Bootstrap函数声明 ========== */

FunctionalData *moving_block_bootstrap(FunctionalData *fd, int block_length);
double compute_test_statistic_at_position(FunctionalData *fd, int position);
BootstrapResult *compute_p_values(FunctionalData *fd, const char *breakpoint_file,
                                  int n_bootstrap, int block_length);
void save_bootstrap_results(const char *output_file, BootstrapResult *result);
void free_bootstrap_result(BootstrapResult *result);

/* ========== 主函数 ========== */

int main(int argc, char *argv[])
{

    if (argc != 6)
    {
        fprintf(stderr, "用法: %s <data_file> <meta_file> <breakpoint_file> <n_bootstrap> <block_length>\n", argv[0]);
        fprintf(stderr, "  breakpoint_file: 断点检测结果文件\n");
        fprintf(stderr, "  n_bootstrap: Bootstrap重复次数（推荐 >= 1000）\n");
        fprintf(stderr, "  block_length: 移动块长度（推荐 20-50）\n");
        return 1;
    }

    const char *data_file = argv[1];
    const char *meta_file = argv[2];
    const char *breakpoint_file = argv[3];
    int n_bootstrap = atoi(argv[4]);
    int block_length = atoi(argv[5]);

    printf("================================================\n");
    printf("Moving Block Bootstrap for p-value Calculation\n");
    printf("================================================\n");
    printf("Bootstrap次数: %d\n", n_bootstrap);
    printf("块长度: %d\n", block_length);
    printf("OpenMP线程数: %d\n", omp_get_max_threads());

    // 设置随机数种子
    srand(time(NULL) + omp_get_thread_num());

    // 加载数据
    printf("\n[1/3] 加载数据...\n");
    FunctionalData *fd = load_data(data_file, meta_file);
    if (fd == NULL)
    {
        return 1;
    }
    printf("  √ T = %d, n_grid = %d\n", fd->T, fd->n_grid);

    // 计算p值
    printf("\n[2/3] 执行Bootstrap计算p值...\n");
    double start_time = omp_get_wtime();

    BootstrapResult *result = compute_p_values(fd, breakpoint_file, n_bootstrap, block_length);

    double end_time = omp_get_wtime();
    printf("  √ Bootstrap完成! 用时: %.2f 秒\n", end_time - start_time);

    // 保存结果
    printf("\n[3/3] 保存结果...\n");
    char output_file[256];
    snprintf(output_file, sizeof(output_file), "%s.bootstrap", breakpoint_file);
    save_bootstrap_results(output_file, result);

    // 清理
    free_bootstrap_result(result);
    free_data(fd);

    printf("\n================================================\n");
    printf("结果已保存到: %s\n", output_file);
    printf("================================================\n");

    return 0;
}

/* ========== Moving Block Bootstrap实现 ========== */

/**
 * Moving Block Bootstrap主函数
 *
 * 算法步骤:
 * 1. 对每个检测到的断点:
 *    a) 计算原始检验统计量
 *    b) 生成B个bootstrap样本
 *    c) 对每个bootstrap样本计算检验统计量
 *    d) p值 = #{bootstrap统计量 >= 原始统计量} / B
 */
BootstrapResult *compute_p_values(FunctionalData *fd, const char *breakpoint_file,
                                  int n_bootstrap, int block_length)
{

    BootstrapResult *result = (BootstrapResult *)calloc(1, sizeof(BootstrapResult));

    // 读取断点位置
    FILE *fp = fopen(breakpoint_file, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "错误: 无法打开断点文件\n");
        return NULL;
    }

    char line[MAX_LINE];
    result->n_breaks = 0;

    while (fgets(line, MAX_LINE, fp))
    {
        if (strncmp(line, "n_breaks=", 9) == 0)
        {
            sscanf(line, "n_breaks=%d", &result->n_breaks);
        }
        else if (result->n_breaks > 0 && strchr(line, ',') != NULL &&
                 strncmp(line, "position", 8) != 0)
        {
            int pos;
            double stat;
            if (sscanf(line, "%d,%lf,", &pos, &stat) == 2)
            {
                int idx = 0;
                while (idx < result->n_breaks && result->positions[idx] != 0)
                    idx++;
                if (idx < result->n_breaks)
                {
                    result->positions[idx] = pos;
                    result->test_stats[idx] = stat;
                }
            }
        }
    }
    fclose(fp);

    printf("  读取到 %d 个断点\n", result->n_breaks);

    // 对每个断点计算p值
    for (int b = 0; b < result->n_breaks; b++)
    {

        int position = result->positions[b];
        double original_stat = result->test_stats[b];

        printf("  断点 %d/%d (position = %d)...\n", b + 1, result->n_breaks, position);

        // Verify the position is valid
        if (position <= 10 || position >= fd->T - 10)
        {
            printf("    [ERROR] 无效位置: %d (T=%d)\n", position, fd->T);
            result->p_values[b] = 1.0;
            continue;
        }

        int count_exceed = 0;

        // Bootstrap循环
        // #pragma omp parallel for reduction(+ : count_exceed)
        for (int boot = 0; boot < n_bootstrap; boot++)
        {
            printf("\r    Bootstrap迭代 %d/%d...", boot + 1, n_bootstrap);
            fflush(stdout);

            // 生成bootstrap样本
            FunctionalData *boot_sample = moving_block_bootstrap(fd, block_length);

            if (boot_sample == NULL)
            {
                printf("\n    [ERROR] 无法生成bootstrap样本\n");
                continue;
            }

            // 计算bootstrap检验统计量
            double boot_stat = compute_test_statistic_at_position(boot_sample, position);

            // 计数
            if (boot_stat >= original_stat)
            {
                count_exceed++;
            }

            // 清理 (不要释放共享的s_grid)
            boot_sample->s_grid = NULL; // 避免double-free
            free_data(boot_sample);
        }

        result->p_values[b] = (double)count_exceed / n_bootstrap;
        printf("\n    √ p值 = %.4f\n", result->p_values[b]);
    }

    return result;
}

/**
 * 生成Moving Block Bootstrap样本
 *
 * 算法:
 * 1. 将时间序列分成长度为block_length的块
 * 2. 随机重采样这些块
 * 3. 拼接成新的时间序列
 */
FunctionalData *moving_block_bootstrap(FunctionalData *fd, int block_length)
{

    // 分配bootstrap样本
    FunctionalData *boot = (FunctionalData *)malloc(sizeof(FunctionalData));
    boot->T = fd->T;
    boot->n_grid = fd->n_grid;
    boot->s_min = fd->s_min;
    boot->s_max = fd->s_max;
    boot->s_grid = fd->s_grid; // 共享网格

    boot->data = (double **)malloc(boot->T * sizeof(double *));
    for (int t = 0; t < boot->T; t++)
    {
        boot->data[t] = (double *)malloc(boot->n_grid * sizeof(double));
    }

    // 计算需要的块数
    int n_blocks = (int)ceil((double)fd->T / block_length);

    // 生成bootstrap样本
    int current_t = 0;

    for (int b = 0; b < n_blocks && current_t < boot->T; b++)
    {

        // 随机选择起始位置
        int start = rand() % (fd->T - block_length + 1);

        // 复制块
        for (int offset = 0; offset < block_length && current_t < boot->T; offset++)
        {
            int source_t = start + offset;
            memcpy(boot->data[current_t], fd->data[source_t],
                   boot->n_grid * sizeof(double));
            current_t++;
        }
    }

    return boot;
}

/**
 * 计算给定位置的检验统计量
 */
double compute_test_statistic_at_position(FunctionalData *fd, int position)
{

    // 确保位置有效
    if (position <= 10 || position >= fd->T - 10)
    {
        return 0.0;
    }

    int start = 0;
    int end = fd->T - 1;

    // 计算分割后的SSGR
    double **mean_left = compute_segment_mean(fd, start, position - 1);
    double **mean_right = compute_segment_mean(fd, position, end);

    double ssgr_left = compute_ssgr(fd, start, position - 1, mean_left);
    double ssgr_right = compute_ssgr(fd, position, end, mean_right);

    // 计算无分割的SSGR
    double **mean_full = compute_segment_mean(fd, start, end);
    double ssgr_full = compute_ssgr(fd, start, end, mean_full);

    // 检验统计量
    double stat = ssgr_full - (ssgr_left + ssgr_right);

    // 清理
    free_segment_mean(mean_left, fd->n_grid);
    free_segment_mean(mean_right, fd->n_grid);
    free_segment_mean(mean_full, fd->n_grid);

    return stat;
}

/**
 * 保存Bootstrap结果
 */
void save_bootstrap_results(const char *output_file, BootstrapResult *result)
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
                result->positions[i],
                result->test_stats[i],
                result->p_values[i]);
    }

    fclose(fp);

    printf("  √ 保存完成\n");
}

/**
 * 释放Bootstrap结果
 */
void free_bootstrap_result(BootstrapResult *result)
{
    if (result == NULL)
        return;
    free(result);
}
