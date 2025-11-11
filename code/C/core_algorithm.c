/*
 * ============================================
 * 函数型数据结构断点检测 - 核心算法
 * ============================================
 * 编译命令（MinGW）:
 * gcc -o breakpoint_detector core_algorithm.c -lm -fopenmp -O3
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define MAX_LINE_LENGTH 100000
#define PI 3.14159265358979323846
#define MIN_SEGMENT_SIZE 30 // 最小segment大小

// ============================================
// 数据结构定义
// ============================================

typedef struct
{
    double **X;     // T × grid_size 数据矩阵
    int T;          // 样本量
    int grid_size;  // 网格点数
    double *t_grid; // 时间网格
} FunctionalData;

typedef struct
{
    int *positions;   // 断点位置
    int n_breaks;     // 断点个数
    double objective; // 目标函数值
} BreakResult;

// ============================================
// 工具函数
// ============================================

// 读取CSV文件到矩阵
double **read_csv_matrix(const char *filename, int *rows, int *cols)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Error: Cannot open %s\n", filename);
        return NULL;
    }

    // 第一遍：计算行列数
    char line[MAX_LINE_LENGTH];
    *rows = 0;
    *cols = 0;

    while (fgets(line, MAX_LINE_LENGTH, fp))
    {
        if (*rows == 0)
        {
            // 计算列数
            char *token = strtok(line, ",");
            while (token)
            {
                (*cols)++;
                token = strtok(NULL, ",");
            }
        }
        (*rows)++;
    }

    // 分配内存
    double **matrix = (double **)malloc(*rows * sizeof(double *));
    for (int i = 0; i < *rows; i++)
    {
        matrix[i] = (double *)malloc(*cols * sizeof(double));
    }

    // 第二遍：读取数据
    rewind(fp);
    int row = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) && row < *rows)
    {
        int col = 0;
        char *token = strtok(line, ",");
        while (token && col < *cols)
        {
            matrix[row][col] = atof(token);
            col++;
            token = strtok(NULL, ",");
        }
        row++;
    }

    fclose(fp);
    return matrix;
}

// 读取向量
double *read_csv_vector(const char *filename, int *length)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Error: Cannot open %s\n", filename);
        return NULL;
    }

    *length = 0;
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, fp))
    {
        (*length)++;
    }

    double *vector = (double *)malloc(*length * sizeof(double));
    rewind(fp);

    int i = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) && i < *length)
    {
        vector[i] = atof(line);
        i++;
    }

    fclose(fp);
    return vector;
}

// ============================================
// ECF计算函数
// ============================================

// 计算某段数据的经验特征函数
void compute_segment_ecf(FunctionalData *data, int start, int end,
                         double *u_grid, int n_u, double *ecf_real, double *ecf_imag)
{
    int n = end - start + 1;

// 对每个u值计算ECF
#pragma omp parallel for
    for (int k = 0; k < n_u; k++)
    {
        double u = u_grid[k];
        double sum_real = 0.0;
        double sum_imag = 0.0;

        for (int i = start; i <= end; i++)
        {
            // 计算积分：\int X_i(t) dt （梯形法则）
            double integral = 0.0;
            for (int j = 0; j < data->grid_size - 1; j++)
            {
                double dt = data->t_grid[j + 1] - data->t_grid[j];
                integral += 0.5 * (data->X[i][j] + data->X[i][j + 1]) * dt;
            }

            // 累加 exp(i * u * integral)
            sum_real += cos(u * integral);
            sum_imag += sin(u * integral);
        }

        ecf_real[k] = sum_real / n;
        ecf_imag[k] = sum_imag / n;
    }
}

// 计算SSGR目标函数
double compute_ssgr(FunctionalData *data, int *break_positions, int n_breaks,
                    double *u_grid, int n_u, double *weights)
{
    double ssgr = 0.0;

    // 定义segments
    int *segment_starts = (int *)malloc((n_breaks + 2) * sizeof(int));
    int *segment_ends = (int *)malloc((n_breaks + 2) * sizeof(int));

    segment_starts[0] = 0;
    for (int i = 0; i < n_breaks; i++)
    {
        segment_ends[i] = break_positions[i] - 1;
        segment_starts[i + 1] = break_positions[i];
    }
    segment_ends[n_breaks] = data->T - 1;

    int n_segments = n_breaks + 1;

    // 对每个segment计算SSGR贡献
    for (int seg = 0; seg < n_segments; seg++)
    {
        int start = segment_starts[seg];
        int end = segment_ends[seg];
        int n = end - start + 1;

        if (n < MIN_SEGMENT_SIZE)
            continue;

        // 计算该segment的ECF
        double *ecf_real = (double *)malloc(n_u * sizeof(double));
        double *ecf_imag = (double *)malloc(n_u * sizeof(double));
        compute_segment_ecf(data, start, end, u_grid, n_u, ecf_real, ecf_imag);

// 计算每个观测到该segment ECF的平方距离
#pragma omp parallel for reduction(+ : ssgr)
        for (int i = start; i <= end; i++)
        {
            for (int k = 0; k < n_u; k++)
            {
                // 计算 \int X_i(t) dt
                double integral = 0.0;
                for (int j = 0; j < data->grid_size - 1; j++)
                {
                    double dt = data->t_grid[j + 1] - data->t_grid[j];
                    integral += 0.5 * (data->X[i][j] + data->X[i][j + 1]) * dt;
                }

                // 计算 exp(i*u*integral) - ecf
                double real_diff = cos(u_grid[k] * integral) - ecf_real[k];
                double imag_diff = sin(u_grid[k] * integral) - ecf_imag[k];

                // 累加加权平方距离
                ssgr += (real_diff * real_diff + imag_diff * imag_diff) * weights[k];
            }
        }

        free(ecf_real);
        free(ecf_imag);
    }

    free(segment_starts);
    free(segment_ends);

    return ssgr / data->T; // 归一化
}

// ============================================
// 二元分割算法
// ============================================

BreakResult *binary_segmentation(FunctionalData *data, int max_breaks,
                                 double epsilon, double *u_grid, int n_u,
                                 double *weights)
{
    BreakResult *result = (BreakResult *)malloc(sizeof(BreakResult));
    result->positions = (int *)malloc(max_breaks * sizeof(int));
    result->n_breaks = 0;

    // 候选断点列表（当前segments的边界）
    int *candidates = (int *)malloc((max_breaks + 2) * sizeof(int));
    int n_candidates = 2;
    candidates[0] = 0;
    candidates[1] = data->T;

    // 初始目标函数值（无断点）
    int no_breaks[] = {};
    double prev_ssgr = compute_ssgr(data, no_breaks, 0, u_grid, n_u, weights);
    result->objective = prev_ssgr;

    printf("Initial SSGR (no breaks): %.6f\n", prev_ssgr);

    // 迭代寻找断点
    for (int iter = 0; iter < max_breaks; iter++)
    {
        double best_improvement = 0.0;
        int best_break = -1;
        int best_segment = -1;

        printf("\nIteration %d: Searching for break point...\n", iter + 1);

        // 遍历每个segment
        for (int seg = 0; seg < n_candidates - 1; seg++)
        {
            int seg_start = candidates[seg];
            int seg_end = candidates[seg + 1];
            int seg_length = seg_end - seg_start;

            if (seg_length < 2 * MIN_SEGMENT_SIZE)
                continue;

            printf("  Segment [%d, %d): length=%d\n", seg_start, seg_end, seg_length);

            // 在该segment内搜索最优分割点
            int search_start = seg_start + MIN_SEGMENT_SIZE;
            int search_end = seg_end - MIN_SEGMENT_SIZE;

#pragma omp parallel for
            for (int k = search_start; k < search_end; k += 5)
            { // 步长5加速
                // 构造试探断点配置
                int *trial_breaks = (int *)malloc((result->n_breaks + 1) * sizeof(int));
                int n_trial = 0;

                for (int i = 0; i < result->n_breaks; i++)
                {
                    trial_breaks[n_trial++] = result->positions[i];
                }
                trial_breaks[n_trial++] = k;

                // 对trial_breaks排序
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

                // 计算新的SSGR
                double new_ssgr = compute_ssgr(data, trial_breaks, n_trial,
                                               u_grid, n_u, weights);
                double improvement = prev_ssgr - new_ssgr;

#pragma omp critical
                {
                    if (improvement > best_improvement)
                    {
                        best_improvement = improvement;
                        best_break = k;
                        best_segment = seg;
                    }
                }

                free(trial_breaks);
            }
        }

        // 检查是否找到有效改进
        if (best_break == -1 || best_improvement < 1e-6)
        {
            printf("No significant improvement found. Stopping.\n");
            break;
        }

        printf("  Best break at position %d, improvement: %.6f\n",
               best_break, best_improvement);

        // 更新断点列表
        result->positions[result->n_breaks++] = best_break;

        // 对断点排序
        for (int i = 0; i < result->n_breaks - 1; i++)
        {
            for (int j = i + 1; j < result->n_breaks; j++)
            {
                if (result->positions[i] > result->positions[j])
                {
                    int temp = result->positions[i];
                    result->positions[i] = result->positions[j];
                    result->positions[j] = temp;
                }
            }
        }

        // 更新candidates列表
        for (int i = n_candidates; i > best_segment + 1; i--)
        {
            candidates[i] = candidates[i - 1];
        }
        candidates[best_segment + 1] = best_break;
        n_candidates++;

        // 更新目标函数值
        prev_ssgr -= best_improvement;
        result->objective = prev_ssgr;

        printf("  Updated SSGR: %.6f\n", prev_ssgr);
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
    printf("Functional Data Structural Break Detection\n");
    printf("============================================\n\n");

    // 设置输入输出路径
    char *input_dir = (argc > 1) ? argv[1] : "data/mean_break";
    char *output_dir = (argc > 2) ? argv[2] : "results";

    char data_path[256], grid_path[256];
    sprintf(data_path, "%s/functional_data.csv", input_dir);
    sprintf(grid_path, "%s/time_grid.csv", input_dir);

    // 读取数据
    printf("Reading data from %s...\n", input_dir);
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
        // 正态权重: W(u) = exp(-u^2/2) / sqrt(2*pi)
        weights[i] = exp(-u_grid[i] * u_grid[i] / 2.0) / sqrt(2 * PI);
    }

    // 运行二元分割
    printf("Running binary segmentation...\n");
    printf("Parameters: max_breaks=5, epsilon=0.1\n\n");

    BreakResult *result = binary_segmentation(&data, 5, 0.1, u_grid, n_u, weights);

    // 输出结果
    printf("\n============================================\n");
    printf("Detection Results:\n");
    printf("============================================\n");
    printf("Number of breaks detected: %d\n", result->n_breaks);
    printf("Break positions: ");
    for (int i = 0; i < result->n_breaks; i++)
    {
        printf("%d ", result->positions[i]);
    }
    printf("\nFinal SSGR: %.6f\n", result->objective);

    // 保存结果
    char output_path[256];
    sprintf(output_path, "%s/detected_breaks.csv", output_dir);
    FILE *fp = fopen(output_path, "w");
    fprintf(fp, "position\n");
    for (int i = 0; i < result->n_breaks; i++)
    {
        fprintf(fp, "%d\n", result->positions[i]);
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
    free(result);

    return 0;
}
