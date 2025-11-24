/*
 * ============================================================================
 * 文件名: utils.c
 * 功能: 共享函数的实现
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "utils.h"

#define MAX_LINE 1024

/* ========== 函数实现 ========== */

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

    // 读取二进制数据
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

double compute_ssgr_weighted(FunctionalData *fd, int start, int end,
                             double **segment_mean, double *s_weights)
{
    double ssgr = 0.0;
    double ds = (fd->s_max - fd->s_min) / (fd->n_grid - 1.0);

#pragma omp parallel for reduction(+ : ssgr)
    for (int t = start; t <= end; t++)
    {
        for (int i = 0; i < fd->n_grid; i++)
        {
            double residual = fd->data[t][i] - segment_mean[i][0];

            // 使用自定义权重函数
            double weight = (s_weights != NULL) ? s_weights[i] : 1.0;

            ssgr += residual * residual * weight * ds;
        }
    }
    return ssgr;
}

double compute_ssgr(FunctionalData *fd, int start, int end, double **segment_mean)
{

    double ssgr = 0.0;
    double ds = (fd->s_max - fd->s_min) / (fd->n_grid - 1.0);

#pragma omp parallel for reduction(+ : ssgr)
    for (int t = start; t <= end; t++)
    {
        for (int i = 0; i < fd->n_grid; i++)
        {
            double residual = fd->data[t][i] - segment_mean[i][0];
            ssgr += residual * residual * ds;
        }
    }

    return ssgr;
}
