/*
 * ============================================================================
 * 文件名: characteristic_function.c
 * 功能: 计算函数型数据的特征函数及SSGR (Segmented Sum of Generalized Residuals)
 * 编译: gcc -O3 -fopenmp -o char_func characteristic_function.c -lm
 * 使用: ./char_func <data_file> <meta_file> <output_file> <u_min> <u_max> <n_u>
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define MAX_LINE 1024
#define PI 3.14159265358979323846

/* ========== 数据结构 ========== */

typedef struct
{
    double **data;  // T x n_grid 数据矩阵
    int T;          // 观测数量
    int n_grid;     // 函数网格点数
    double s_min;   // 定义域下界
    double s_max;   // 定义域上界
    double *s_grid; // 网格坐标
} FunctionalData;

typedef struct
{
    double *real; // 实部
    double *imag; // 虚部
    int length;
} ComplexArray;

/* ========== 函数声明 ========== */

FunctionalData *load_data(const char *data_file, const char *meta_file);
void free_data(FunctionalData *fd);
ComplexArray *compute_characteristic_function(FunctionalData *fd, double *u_grid, int n_u);
void free_complex_array(ComplexArray *ca);
double compute_ssgr_segment(FunctionalData *fd, int t_start, int t_end,
                            double *u_grid, int n_u);
void save_results(const char *output_file, double *ssgr, int n_breakpoints);

/* ========== 主函数 ========== */

int main(int argc, char *argv[])
{

    if (argc != 7)
    {
        fprintf(stderr, "用法: %s <data_file> <meta_file> <output_file> <u_min> <u_max> <n_u>\n", argv[0]);
        return 1;
    }

    // 解析命令行参数
    const char *data_file = argv[1];
    const char *meta_file = argv[2];
    const char *output_file = argv[3];
    double u_min = atof(argv[4]);
    double u_max = atof(argv[5]);
    int n_u = atoi(argv[6]);

    printf("================================================\n");
    printf("特征函数计算模块\n");
    printf("================================================\n");
    printf("数据文件: %s\n", data_file);
    printf("元数据: %s\n", meta_file);
    printf("u范围: [%.2f, %.2f], 网格数: %d\n", u_min, u_max, n_u);

    // 加载数据
    printf("\n[1/3] 加载数据...\n");
    FunctionalData *fd = load_data(data_file, meta_file);
    if (fd == NULL)
    {
        fprintf(stderr, "错误: 无法加载数据\n");
        return 1;
    }
    printf("  √ 成功加载 %d 个观测, 每个有 %d 个网格点\n", fd->T, fd->n_grid);

    // 构建u网格
    double *u_grid = (double *)malloc(n_u * sizeof(double));
    for (int i = 0; i < n_u; i++)
    {
        u_grid[i] = u_min + (u_max - u_min) * i / (n_u - 1.0);
    }

    // 计算特征函数
    printf("\n[2/3] 计算特征函数 (使用 %d 个OpenMP线程)...\n", omp_get_max_threads());

    double start_time = omp_get_wtime();
    ComplexArray *char_func = compute_characteristic_function(fd, u_grid, n_u);
    double end_time = omp_get_wtime();

    printf("  √ 完成! 用时: %.2f 秒\n", end_time - start_time);

    // 计算SSGR (这里简化为全样本的残差平方和)
    printf("\n[3/3] 计算SSGR...\n");
    double ssgr_total = 0.0;

#pragma omp parallel for reduction(+ : ssgr_total)
    for (int t = 0; t < fd->T; t++)
    {
        for (int i = 0; i < n_u; i++)
        {
            double residual_real = char_func->real[t * n_u + i] - char_func->real[i]; // 简化版
            double residual_imag = char_func->imag[t * n_u + i] - char_func->imag[i];
            ssgr_total += residual_real * residual_real + residual_imag * residual_imag;
        }
    }

    printf("  √ 总SSGR = %.6e\n", ssgr_total);

    // 保存结果
    double results[1] = {ssgr_total};
    save_results(output_file, results, 1);

    // 清理
    free(u_grid);
    free_complex_array(char_func);
    free_data(fd);

    printf("\n================================================\n");
    printf("计算完成! 结果已保存到: %s\n", output_file);
    printf("================================================\n");

    return 0;
}

/* ========== 函数实现 ========== */

/**
 * 加载二进制数据文件
 */
FunctionalData *load_data(const char *data_file, const char *meta_file)
{

    FunctionalData *fd = (FunctionalData *)malloc(sizeof(FunctionalData));

    // 读取元数据
    FILE *meta_fp = fopen(meta_file, "r");
    if (meta_fp == NULL)
    {
        fprintf(stderr, "错误: 无法打开元数据文件 %s\n", meta_file);
        free(fd);
        return NULL;
    }

    char line[MAX_LINE];
    while (fgets(line, MAX_LINE, meta_fp))
    {
        if (sscanf(line, "T=%d", &fd->T) == 1)
            continue;
        if (sscanf(line, "n_grid=%d", &fd->n_grid) == 1)
            continue;
        if (sscanf(line, "s_min=%lf", &fd->s_min) == 1)
            continue;
        if (sscanf(line, "s_max=%lf", &fd->s_max) == 1)
            continue;
    }
    fclose(meta_fp);

    // 分配内存
    fd->data = (double **)malloc(fd->T * sizeof(double *));
    for (int t = 0; t < fd->T; t++)
    {
        fd->data[t] = (double *)malloc(fd->n_grid * sizeof(double));
    }

    // 读取二进制数据 (假设按行主序存储)
    FILE *data_fp = fopen(data_file, "rb");
    if (data_fp == NULL)
    {
        fprintf(stderr, "错误: 无法打开数据文件 %s\n", data_file);
        free_data(fd);
        return NULL;
    }

    for (int t = 0; t < fd->T; t++)
    {
        size_t read_count = fread(fd->data[t], sizeof(double), fd->n_grid, data_fp);
        if (read_count != (size_t)fd->n_grid)
        {
            fprintf(stderr, "错误: 数据读取不完整 (第 %d 行)\n", t);
            fclose(data_fp);
            free_data(fd);
            return NULL;
        }
    }
    fclose(data_fp);

    // 生成s网格
    fd->s_grid = (double *)malloc(fd->n_grid * sizeof(double));
    for (int i = 0; i < fd->n_grid; i++)
    {
        fd->s_grid[i] = fd->s_min + (fd->s_max - fd->s_min) * i / (fd->n_grid - 1.0);
    }

    return fd;
}

/**
 * 释放数据结构
 */
void free_data(FunctionalData *fd)
{
    if (fd == NULL)
        return;

    if (fd->data != NULL)
    {
        for (int t = 0; t < fd->T; t++)
        {
            free(fd->data[t]);
        }
        free(fd->data);
    }

    free(fd->s_grid);
    free(fd);
}

/**
 * 计算特征函数 φ_t(u) = E[exp(iu ∫ X_t(s)ds)]
 * 使用梯形法则进行数值积分
 */
ComplexArray *compute_characteristic_function(FunctionalData *fd, double *u_grid, int n_u)
{

    ComplexArray *ca = (ComplexArray *)malloc(sizeof(ComplexArray));
    ca->length = fd->T * n_u;
    ca->real = (double *)calloc(ca->length, sizeof(double));
    ca->imag = (double *)calloc(ca->length, sizeof(double));

    double ds = (fd->s_max - fd->s_min) / (fd->n_grid - 1.0);

#pragma omp parallel for collapse(2)
    for (int t = 0; t < fd->T; t++)
    {
        for (int j = 0; j < n_u; j++)
        {

            double u = u_grid[j];

            // 计算积分 ∫ X_t(s)ds 使用梯形法则
            double integral = 0.0;
            for (int i = 0; i < fd->n_grid - 1; i++)
            {
                integral += 0.5 * (fd->data[t][i] + fd->data[t][i + 1]) * ds;
            }

            // 计算 exp(iu * integral)
            int idx = t * n_u + j;
            ca->real[idx] = cos(u * integral);
            ca->imag[idx] = sin(u * integral);
        }
    }

    return ca;
}

/**
 * 释放复数数组
 */
void free_complex_array(ComplexArray *ca)
{
    if (ca == NULL)
        return;
    free(ca->real);
    free(ca->imag);
    free(ca);
}

/**
 * 保存结果
 */
void save_results(const char *output_file, double *ssgr, int n_breakpoints)
{
    FILE *fp = fopen(output_file, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "错误: 无法创建输出文件 %s\n", output_file);
        return;
    }

    fprintf(fp, "SSGR\n");
    for (int i = 0; i < n_breakpoints; i++)
    {
        fprintf(fp, "%.10e\n", ssgr[i]);
    }

    fclose(fp);
}
