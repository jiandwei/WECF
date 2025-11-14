/*
 * ============================================================================
 * 文件名: breakpoint_detection.c
 * 功能: 实现二元分割算法进行断点检测
 * 理论基础: Binary Segmentation with BIC criterion
 * 编译: gcc -O3 -fopenmp -o break_detect breakpoint_detection.c -lm
 * 使用: ./break_detect <data_file> <meta_file> <output_file> <max_breaks> <min_segment_length>
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#define MAX_LINE 1024
#define PI 3.14159265358979323846

/* ========== 数据结构 (与characteristic_function.c相同) ========== */

typedef struct
{
    double **data;
    int T;
    int n_grid;
    double s_min;
    double s_max;
    double *s_grid;
} FunctionalData;

typedef struct
{
    int position;     // 断点位置
    double test_stat; // 检验统计量
    double p_value;   // p值（bootstrap估计）
} Breakpoint;

typedef struct
{
    Breakpoint *breaks;
    int n_breaks;
    double *segment_means; // 每个segment的均值函数
} BreakpointResult;

/* ========== 函数声明 ========== */

FunctionalData *load_data(const char *data_file, const char *meta_file);
void free_data(FunctionalData *fd);

double compute_ssgr(FunctionalData *fd, int start, int end, double **segment_mean);
double **compute_segment_mean(FunctionalData *fd, int start, int end);
void free_segment_mean(double **mean, int n_grid);

BreakpointResult *binary_segmentation(FunctionalData *fd, int max_breaks, int min_seg_len);
int find_best_split(FunctionalData *fd, int start, int end, int min_seg_len, double *stat);

double compute_bic(double ssgr, int n_obs, int n_params);
void save_breakpoints(const char *output_file, BreakpointResult *result);
void free_breakpoint_result(BreakpointResult *result);

/* ========== 主函数 ========== */

int main(int argc, char *argv[])
{

    if (argc != 6)
    {
        fprintf(stderr, "用法: %s <data_file> <meta_file> <output_file> <max_breaks> <min_segment_length>\n", argv[0]);
        return 1;
    }

    const char *data_file = argv[1];
    const char *meta_file = argv[2];
    const char *output_file = argv[3];
    int max_breaks = atoi(argv[4]);
    int min_seg_len = atoi(argv[5]);

    printf("================================================\n");
    printf("断点检测模块 (Binary Segmentation + BIC)\n");
    printf("================================================\n");
    printf("最大断点数: %d\n", max_breaks);
    printf("最小segment长度: %d\n", min_seg_len);
    printf("OpenMP线程数: %d\n", omp_get_max_threads());

    // 加载数据
    printf("\n[1/2] 加载数据...\n");
    FunctionalData *fd = load_data(data_file, meta_file);
    if (fd == NULL)
    {
        return 1;
    }
    printf("  √ T = %d, n_grid = %d\n", fd->T, fd->n_grid);

    // 执行断点检测
    printf("\n[2/2] 执行二元分割算法...\n");
    double start_time = omp_get_wtime();

    BreakpointResult *result = binary_segmentation(fd, max_breaks, min_seg_len);

    double end_time = omp_get_wtime();
    printf("  √ 检测完成! 用时: %.2f 秒\n", end_time - start_time);
    printf("  √ 检测到 %d 个断点\n", result->n_breaks);

    // 保存结果
    save_breakpoints(output_file, result);

    // 清理
    free_breakpoint_result(result);
    free_data(fd);

    printf("\n================================================\n");
    printf("结果已保存到: %s\n", output_file);
    printf("================================================\n");

    return 0;
}

/* ========== 核心算法实现 ========== */

/**
 * 二元分割算法主流程
 *
 * 算法步骤:
 * 1. 初始化活跃segment集合 S = {[1, T]}
 * 2. 对每个segment in S:
 *    a) 寻找最优分割点 τ* = argmin_{τ} SSGR(segment)
 *    b) 计算BIC(1 break) vs BIC(0 break)
 *    c) 如果BIC改进且统计量显著，则分割
 * 3. 重复直到无法继续分割或达到max_breaks
 */
BreakpointResult *binary_segmentation(FunctionalData *fd, int max_breaks, int min_seg_len)
{

    // 初始化结果结构
    BreakpointResult *result = (BreakpointResult *)malloc(sizeof(BreakpointResult));
    result->breaks = (Breakpoint *)malloc(max_breaks * sizeof(Breakpoint));
    result->n_breaks = 0;

    // 活跃segment队列 (用简单数组实现)
    int *seg_starts = (int *)malloc((max_breaks + 1) * sizeof(int));
    int *seg_ends = (int *)malloc((max_breaks + 1) * sizeof(int));
    int n_segments = 1;

    seg_starts[0] = 0;
    seg_ends[0] = fd->T - 1;

    // 当前BIC
    double **current_mean = compute_segment_mean(fd, 0, fd->T - 1);
    double current_ssgr = compute_ssgr(fd, 0, fd->T - 1, current_mean);
    double current_bic = compute_bic(current_ssgr, fd->T, fd->n_grid);

    free_segment_mean(current_mean, fd->n_grid);

    printf("  初始BIC = %.4f\n", current_bic);

    // 迭代分割
    for (int iter = 0; iter < max_breaks; iter++)
    {

        printf("\n  --- 迭代 %d ---\n", iter + 1);

        // 寻找所有segment中最优的分割点
        int best_seg_idx = -1;
        int best_split_pos = -1;
        double best_improvement = 0.0;
        double best_test_stat = 0.0;

        for (int seg = 0; seg < n_segments; seg++)
        {
            int start = seg_starts[seg];
            int end = seg_ends[seg];
            int seg_length = end - start + 1;

            if (seg_length < 2 * min_seg_len)
            {
                continue; // segment太短，无法分割
            }

            // 寻找该segment的最优分割点
            double test_stat;
            int split_pos = find_best_split(fd, start, end, min_seg_len, &test_stat);

            if (split_pos < 0)
            {
                continue; // 未找到有效分割点
            }

            // 计算分割后的BIC改进
            double **mean_left = compute_segment_mean(fd, start, split_pos);
            double **mean_right = compute_segment_mean(fd, split_pos + 1, end);

            double ssgr_left = compute_ssgr(fd, start, split_pos, mean_left);
            double ssgr_right = compute_ssgr(fd, split_pos + 1, end, mean_right);

            double new_ssgr = ssgr_left + ssgr_right;
            double new_bic = compute_bic(new_ssgr, fd->T, (result->n_breaks + 2) * fd->n_grid);

            double improvement = current_bic - new_bic;

            free_segment_mean(mean_left, fd->n_grid);
            free_segment_mean(mean_right, fd->n_grid);

            if (improvement > best_improvement)
            {
                best_improvement = improvement;
                best_seg_idx = seg;
                best_split_pos = split_pos;
                best_test_stat = test_stat;
            }
        }

        // 判断是否继续分割
        if (best_seg_idx < 0 || best_improvement <= 0)
        {
            printf("  √ 无进一步改进，算法终止\n");
            break;
        }

        // 执行分割
        printf("  √ 在segment [%d, %d] 的位置 %d 处分割\n",
               seg_starts[best_seg_idx], seg_ends[best_seg_idx], best_split_pos);
        printf("     BIC改进: %.4f, 检验统计量: %.4f\n", best_improvement, best_test_stat);

        // 记录断点
        result->breaks[result->n_breaks].position = best_split_pos;
        result->breaks[result->n_breaks].test_stat = best_test_stat;
        result->breaks[result->n_breaks].p_value = -1.0; // 需要bootstrap估计
        result->n_breaks++;

        // 更新segment列表
        int old_start = seg_starts[best_seg_idx];
        int old_end = seg_ends[best_seg_idx];

        seg_ends[best_seg_idx] = best_split_pos;
        seg_starts[n_segments] = best_split_pos + 1;
        seg_ends[n_segments] = old_end;
        n_segments++;

        current_bic -= best_improvement;
    }

    // 清理
    free(seg_starts);
    free(seg_ends);

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

    return result;
}

/**
 * 在给定segment中寻找最优分割点
 * 返回: 分割位置 (若无有效分割则返回-1)
 */
int find_best_split(FunctionalData *fd, int start, int end, int min_seg_len, double *stat)
{

    int seg_length = end - start + 1;
    if (seg_length < 2 * min_seg_len)
    {
        return -1;
    }

    int best_pos = -1;
    double min_ssgr = DBL_MAX;

// 并行搜索所有候选分割点
#pragma omp parallel
    {
        double local_min_ssgr = DBL_MAX;
        int local_best_pos = -1;

#pragma omp for
        for (int tau = start + min_seg_len; tau <= end - min_seg_len; tau++)
        {

            // 计算左右两个segment的SSGR
            double **mean_left = compute_segment_mean(fd, start, tau);
            double **mean_right = compute_segment_mean(fd, tau + 1, end);

            double ssgr_left = compute_ssgr(fd, start, tau, mean_left);
            double ssgr_right = compute_ssgr(fd, tau + 1, end, mean_right);
            double total_ssgr = ssgr_left + ssgr_right;

            free_segment_mean(mean_left, fd->n_grid);
            free_segment_mean(mean_right, fd->n_grid);

            if (total_ssgr < local_min_ssgr)
            {
                local_min_ssgr = total_ssgr;
                local_best_pos = tau;
            }
        }

// 归约
#pragma omp critical
        {
            if (local_min_ssgr < min_ssgr)
            {
                min_ssgr = local_min_ssgr;
                best_pos = local_best_pos;
            }
        }
    }

    // 计算检验统计量 (简化版: CUSUM-type)
    if (best_pos > 0)
    {
        double **mean_full = compute_segment_mean(fd, start, end);
        double ssgr_full = compute_ssgr(fd, start, end, mean_full);
        *stat = (ssgr_full - min_ssgr) / seg_length; // 标准化
        free_segment_mean(mean_full, fd->n_grid);
    }

    return best_pos;
}

/**
 * 计算segment的SSGR (Segmented Sum of Generalized Residuals)
 * SSGR = Σ_t ∫ |X_t(s) - mean(s)|^2 W(s)ds
 */
double compute_ssgr(FunctionalData *fd, int start, int end, double **segment_mean)
{

    double ssgr = 0.0;
    double ds = (fd->s_max - fd->s_min) / (fd->n_grid - 1.0);

#pragma omp parallel for reduction(+ : ssgr)
    for (int t = start; t <= end; t++)
    {
        for (int i = 0; i < fd->n_grid; i++)
        {
            double residual = fd->data[t][i] - segment_mean[i][0]; // mean is n_grid x 1
            ssgr += residual * residual * ds;                      // 简化权重W(s)=1
        }
    }

    return ssgr;
}

/**
 * 计算segment的平均函数 mean(s) = (1/n) Σ_t X_t(s)
 */
double **compute_segment_mean(FunctionalData *fd, int start, int end)
{

    double **mean = (double **)malloc(fd->n_grid * sizeof(double *));
    for (int i = 0; i < fd->n_grid; i++)
    {
        mean[i] = (double *)calloc(1, sizeof(double));
    }

    int n_obs = end - start + 1;

    for (int t = start; t <= end; t++)
    {
        for (int i = 0; i < fd->n_grid; i++)
        {
            mean[i][0] += fd->data[t][i];
        }
    }

    for (int i = 0; i < fd->n_grid; i++)
    {
        mean[i][0] /= n_obs;
    }

    return mean;
}

void free_segment_mean(double **mean, int n_grid)
{
    for (int i = 0; i < n_grid; i++)
    {
        free(mean[i]);
    }
    free(mean);
}

/**
 * 计算BIC = T * log(SSGR/T) + k * log(T)
 */
double compute_bic(double ssgr, int n_obs, int n_params)
{
    return n_obs * log(ssgr / n_obs) + n_params * log(n_obs);
}

/**
 * 保存断点结果
 */
void save_breakpoints(const char *output_file, BreakpointResult *result)
{
    FILE *fp = fopen(output_file, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "错误: 无法创建输出文件\n");
        return;
    }

    fprintf(fp, "n_breaks=%d\n", result->n_breaks);
    fprintf(fp, "position,test_stat,p_value\n");

    for (int i = 0; i < result->n_breaks; i++)
    {
        fprintf(fp, "%d,%.6f,%.6f\n",
                result->breaks[i].position,
                result->breaks[i].test_stat,
                result->breaks[i].p_value);
    }

    fclose(fp);
}

void free_breakpoint_result(BreakpointResult *result)
{
    free(result->breaks);
    free(result);
}

/* ========== 辅助函数 (与characteristic_function.c相同) ========== */

FunctionalData *load_data(const char *data_file, const char *meta_file)
{
    // [实现与之前相同，这里省略]
    // ...
    return NULL; // placeholder
}

void free_data(FunctionalData *fd)
{
    // [实现与之前相同]
}
