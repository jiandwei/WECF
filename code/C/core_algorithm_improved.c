/*
 * ============================================
 * 改进版断点检测：带BIC和假设检验
 * ============================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#define MAX_LINE_LENGTH 100000
#define PI 3.14159265358979323846
#define MIN_SEGMENT_SIZE 30
#define N_PERMUTATIONS 199 // 置换检验次数

// ============================================
// 数据结构（同前）
// ============================================

typedef struct
{
    double **X;
    int T;
    int grid_size;
    double *t_grid;
} FunctionalData;

typedef struct
{
    int *positions;
    int n_breaks;
    double objective;
    double bic;       // 新增：BIC值
    double *p_values; // 新增：每个断点的p值
} BreakResult;

// ============================================
// 工具函数（省略，与前面相同）
// ============================================

double **read_csv_matrix(const char *filename, int *rows, int *cols);
double *read_csv_vector(const char *filename, int *length);

// ============================================
// BIC计算
// ============================================

/**
 * 计算BIC = T * log(SSGR) + k * log(T) * penalty_factor
 * 其中 k = (M+1) * n_u 是参数个数
 * penalty_factor 控制对复杂度的惩罚程度
 */
double compute_bic(double ssgr, int T, int n_breaks, int n_u,
                   double penalty_factor)
{
    int n_regimes = n_breaks + 1;
    int k = n_regimes * n_u; // 每个regime有n_u个参数（ECF的实部和虚部）

    double log_likelihood = -0.5 * T * log(2 * PI * ssgr) - 0.5 * T;
    double penalty = 0.5 * k * log(T) * penalty_factor;

    double bic = -2 * log_likelihood + 2 * penalty;

    return bic;
}

// ============================================
// 置换检验：评估断点显著性
// ============================================

/**
 * 对给定候选断点进行置换检验
 *
 * @param data 原始数据
 * @param candidate_break 候选断点位置
 * @param u_grid 频率网格
 * @param n_u 频率点数
 * @param weights 权重函数
 * @param n_perm 置换次数
 * @return p值
 */
double permutation_test(FunctionalData *data, int candidate_break,
                        double *u_grid, int n_u, double *weights,
                        int n_perm)
{

    // 计算观测统计量
    int single_break[] = {candidate_break};
    double obs_ssgr = compute_ssgr(data, single_break, 1, u_grid, n_u, weights);
    double null_ssgr = compute_ssgr(data, NULL, 0, u_grid, n_u, weights);
    double obs_statistic = null_ssgr - obs_ssgr;

    // 置换检验
    int n_extreme = 0;

#pragma omp parallel for reduction(+ : n_extreme)
    for (int perm = 0; perm < n_perm; perm++)
    {
        // 创建置换样本（打乱时间顺序但保持曲线完整性）
        int *perm_indices = (int *)malloc(data->T * sizeof(int));
        for (int i = 0; i < data->T; i++)
            perm_indices[i] = i;

        // Fisher-Yates shuffle
        unsigned int seed = time(NULL) ^ (perm << 16);
        for (int i = data->T - 1; i > 0; i--)
        {
            int j = rand_r(&seed) % (i + 1);
            int temp = perm_indices[i];
            perm_indices[i] = perm_indices[j];
            perm_indices[j] = temp;
        }

        // 创建置换数据
        FunctionalData perm_data;
        perm_data.T = data->T;
        perm_data.grid_size = data->grid_size;
        perm_data.t_grid = data->t_grid;
        perm_data.X = (double **)malloc(data->T * sizeof(double *));
        for (int i = 0; i < data->T; i++)
        {
            perm_data.X[i] = data->X[perm_indices[i]];
        }

        // 计算置换统计量
        double perm_ssgr = compute_ssgr(&perm_data, single_break, 1,
                                        u_grid, n_u, weights);
        double perm_statistic = null_ssgr - perm_ssgr;

        if (perm_statistic >= obs_statistic)
        {
            n_extreme++;
        }

        free(perm_data.X);
        free(perm_indices);
    }

    double p_value = (n_extreme + 1.0) / (n_perm + 1.0);
    return p_value;
}

// ============================================
// 改进的二元分割算法
// ============================================

/**
 * 带BIC和假设检验的二元分割
 *
 * @param data 函数型数据
 * @param max_breaks 最大断点数
 * @param significance_level 显著性水平
 * @param bic_penalty BIC惩罚因子
 * @param u_grid 频率网格
 * @param n_u 频率点数
 * @param weights 权重函数
 * @return 检测结果
 */
BreakResult *binary_segmentation_improved(FunctionalData *data,
                                          int max_breaks,
                                          double significance_level,
                                          double bic_penalty,
                                          double *u_grid,
                                          int n_u,
                                          double *weights)
{

    BreakResult *result = (BreakResult *)malloc(sizeof(BreakResult));
    result->positions = (int *)malloc(max_breaks * sizeof(int));
    result->p_values = (double *)malloc(max_breaks * sizeof(double));
    result->n_breaks = 0;

    // 候选segments
    int *candidates = (int *)malloc((max_breaks + 2) * sizeof(int));
    int n_candidates = 2;
    candidates[0] = 0;
    candidates[1] = data->T;

    // 初始状态
    double prev_ssgr = compute_ssgr(data, NULL, 0, u_grid, n_u, weights);
    double prev_bic = compute_bic(prev_ssgr, data->T, 0, n_u, bic_penalty);

    result->objective = prev_ssgr;
    result->bic = prev_bic;

    printf("Initial state:\n");
    printf("  SSGR = %.6f\n", prev_ssgr);
    printf("  BIC = %.6f\n\n", prev_bic);

    // 迭代寻找断点
    for (int iter = 0; iter < max_breaks; iter++)
    {
        printf("Iteration %d:\n", iter + 1);
        printf("  Current # breaks: %d\n", result->n_breaks);

        double best_new_ssgr = INFINITY;
        double best_new_bic = INFINITY;
        int best_break = -1;
        int best_segment = -1;

        // 遍历每个segment
        for (int seg = 0; seg < n_candidates - 1; seg++)
        {
            int seg_start = candidates[seg];
            int seg_end = candidates[seg + 1];
            int seg_length = seg_end - seg_start;

            if (seg_length < 2 * MIN_SEGMENT_SIZE)
                continue;

            int search_start = seg_start + MIN_SEGMENT_SIZE;
            int search_end = seg_end - MIN_SEGMENT_SIZE;

// 在该segment内搜索
#pragma omp parallel for
            for (int k = search_start; k < search_end; k += 3)
            {
                // 构造试探配置
                int *trial_breaks = (int *)malloc((result->n_breaks + 1) * sizeof(int));
                int n_trial = 0;

                for (int i = 0; i < result->n_breaks; i++)
                {
                    trial_breaks[n_trial++] = result->positions[i];
                }
                trial_breaks[n_trial++] = k;

                // 排序
                for (int i = 0; i < n_trial - 1; i++)
                {
                    for (int j = i + 1; j < n_trial; j++)
                    {
                        if (trial_breaks[i] > trial_breaks[j])
                        {
                            int temp = trial_breaks[i];
                            trial_breaks[i] = trial_breaks[j];
                            trial_breaks[j] = temp;
                        }
                    }
                }

                // 计算新SSGR和BIC
                double new_ssgr = compute_ssgr(data, trial_breaks, n_trial,
                                               u_grid, n_u, weights);
                double new_bic = compute_bic(new_ssgr, data->T, n_trial,
                                             n_u, bic_penalty);

#pragma omp critical
                {
                    // 使用BIC作为选择标准
                    if (new_bic < best_new_bic)
                    {
                        best_new_bic = new_bic;
                        best_new_ssgr = new_ssgr;
                        best_break = k;
                        best_segment = seg;
                    }
                }

                free(trial_breaks);
            }
        }

        // 检查1：BIC是否改进
        if (best_break == -1)
        {
            printf("  No valid candidate found.\n");
            break;
        }

        double bic_improvement = prev_bic - best_new_bic;
        printf("  Best candidate: position %d\n", best_break);
        printf("    SSGR: %.6f -> %.6f (Δ = %.6f)\n",
               prev_ssgr, best_new_ssgr, prev_ssgr - best_new_ssgr);
        printf("    BIC: %.6f -> %.6f (Δ = %.6f)\n",
               prev_bic, best_new_bic, bic_improvement);

        if (bic_improvement < 0)
        {
            printf("  BIC did not improve. Stopping.\n");
            break;
        }

        // 检查2：置换检验
        printf("  Running permutation test (%d permutations)...\n", N_PERMUTATIONS);
        double p_value = permutation_test(data, best_break, u_grid, n_u,
                                          weights, N_PERMUTATIONS);
        printf("    p-value = %.4f\n", p_value);

        if (p_value > significance_level)
        {
            printf("  Not significant at α = %.3f. Stopping.\n", significance_level);
            break;
        }

        // 通过所有检验：接受该断点
        printf("  ✓ Break accepted!\n\n");

        result->positions[result->n_breaks] = best_break;
        result->p_values[result->n_breaks] = p_value;
        result->n_breaks++;

        // 排序
        for (int i = 0; i < result->n_breaks - 1; i++)
        {
            for (int j = i + 1; j < result->n_breaks; j++)
            {
                if (result->positions[i] > result->positions[j])
                {
                    int temp_pos = result->positions[i];
                    result->positions[i] = result->positions[j];
                    result->positions[j] = temp_pos;

                    double temp_p = result->p_values[i];
                    result->p_values[i] = result->p_values[j];
                    result->p_values[j] = temp_p;
                }
            }
        }

        // 更新candidates
        for (int i = n_candidates; i > best_segment + 1; i--)
        {
            candidates[i] = candidates[i - 1];
        }
        candidates[best_segment + 1] = best_break;
        n_candidates++;

        // 更新状态
        prev_ssgr = best_new_ssgr;
        prev_bic = best_new_bic;
        result->objective = prev_ssgr;
        result->bic = prev_bic;
    }

    free(candidates);
    return result;
}

// ============================================
// 主函数
// ============================================

int main(int argc, char *argv[])
{
    printf("============================================\n");
    printf("Improved Breakpoint Detection with BIC\n");
    printf("============================================\n\n");

    // 参数解析
    char *input_dir = (argc > 1) ? argv[1] : "data/mean_break";
    char *output_dir = (argc > 2) ? argv[2] : "results";
    double significance_level = (argc > 3) ? atof(argv[3]) : 0.05;
    double bic_penalty = (argc > 4) ? atof(argv[4]) : 1.0;

    printf("Parameters:\n");
    printf("  Input dir: %s\n", input_dir);
    printf("  Output dir: %s\n", output_dir);
    printf("  Significance level: %.3f\n", significance_level);
    printf("  BIC penalty factor: %.2f\n\n", bic_penalty);

    // 读取数据
    char data_path[256], grid_path[256];
    sprintf(data_path, "%s/functional_data.csv", input_dir);
    sprintf(grid_path, "%s/time_grid.csv", input_dir);

    FunctionalData data;
    data.X = read_csv_matrix(data_path, &data.T, &data.grid_size);
    data.t_grid = read_csv_vector(grid_path, &data.grid_size);

    if (!data.X || !data.t_grid)
    {
        fprintf(stderr, "Error reading data!\n");
        return 1;
    }

    printf("Data loaded: T=%d, grid_size=%d\n\n", data.T, data.grid_size);

    // 设置u网格和权重
    int n_u = 50;
    double u_max = 5.0;
    double *u_grid = (double *)malloc(n_u * sizeof(double));
    double *weights = (double *)malloc(n_u * sizeof(double));

    for (int i = 0; i < n_u; i++)
    {
        u_grid[i] = -u_max + i * (2 * u_max) / (n_u - 1);
        weights[i] = exp(-u_grid[i] * u_grid[i] / 2.0) / sqrt(2 * PI);
    }

    // 运行改进算法
    BreakResult *result = binary_segmentation_improved(
        &data, 5, significance_level, bic_penalty, u_grid, n_u, weights);

    // 输出结果
    printf("\n============================================\n");
    printf("Final Results:\n");
    printf("============================================\n");
    printf("Number of breaks: %d\n", result->n_breaks);

    if (result->n_breaks > 0)
    {
        printf("\nBreak Details:\n");
        printf("%-10s %-15s %-10s\n", "Position", "p-value", "Fraction");
        printf("%-10s %-15s %-10s\n", "--------", "-------", "--------");
        for (int i = 0; i < result->n_breaks; i++)
        {
            printf("%-10d %-15.4f %-10.3f\n",
                   result->positions[i],
                   result->p_values[i],
                   (double)result->positions[i] / data.T);
        }
    }

    printf("\nModel Selection:\n");
    printf("  Final SSGR: %.6f\n", result->objective);
    printf("  Final BIC: %.6f\n", result->bic);

    // 保存结果
    char output_path[256];
    sprintf(output_path, "%s/detected_breaks_improved.csv", output_dir);
    FILE *fp = fopen(output_path, "w");
    fprintf(fp, "position,p_value,fraction\n");
    for (int i = 0; i < result->n_breaks; i++)
    {
        fprintf(fp, "%d,%.6f,%.6f\n",
                result->positions[i],
                result->p_values[i],
                (double)result->positions[i] / data.T);
    }
    fclose(fp);

    printf("\nResults saved to %s\n", output_path);

    // 清理
    for (int i = 0; i < data.T; i++)
        free(data.X[i]);
    free(data.X);
    free(data.t_grid);
    free(u_grid);
    free(weights);
    free(result->positions);
    free(result->p_values);
    free(result);

    return 0;
}
