/*
 * ============================================================================
 * 文件名: utils.h
 * 功能: 共享函数的头文件声明
 * ============================================================================
 */

#ifndef UTILS_H
#define UTILS_H

/* ========== 数据结构定义 ========== */

typedef struct
{
    double **data;
    int T;
    int n_grid;
    double s_min;
    double s_max;
    double *s_grid;
} FunctionalData;

/* ========== 函数声明 ========== */

/**
 * 加载数据
 */
FunctionalData *load_data(const char *data_file, const char *meta_file);

/**
 * 释放数据
 */
void free_data(FunctionalData *fd);

/**
 * 计算segment的平均函数
 */
double **compute_segment_mean(FunctionalData *fd, int start, int end);

/**
 * 释放segment均值
 */
void free_segment_mean(double **mean, int n_grid);

/**
 * 计算segment的SSGR
 */
double compute_ssgr(FunctionalData *fd, int start, int end, double **segment_mean);

#endif /* UTILS_H */
