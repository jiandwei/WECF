# 算法目录

本目录用于存放算法相关的文件和代码。

## 用途

您可以在此目录下存放：
- 算法伪代码
- 算法实现代码
- 算法说明文档

## 使用方法

如果您创建了算法相关的文件，可以在 LaTeX 文档中引用：

```latex
% 引入算法伪代码
\input{算法/algorithm1.tex}
```

## 推荐的文件格式

- **TEX**: LaTeX 格式的算法伪代码
- **PY**: Python 实现
- **R**: R 语言实现
- **TXT**: 纯文本说明

## 命名规范

使用描述性的文件名：
- `algorithm1_optimization.tex` - 优化算法
- `algorithm2_estimation.py` - 估计算法的 Python 实现
- `sorting_algorithm.tex` - 排序算法

## 示例

### LaTeX 中使用算法包

```latex
\usepackage{algorithm}
\usepackage{algorithmic}

\begin{algorithm}
    \caption{Your Algorithm Name}
    \label{alg:your_algorithm}
    \begin{algorithmic}[1]
        \STATE Initialize variables
        \FOR{each iteration}
            \STATE Perform computation
        \ENDFOR
        \RETURN result
    \end{algorithmic}
\end{algorithm}
```
