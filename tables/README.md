# Tables Directory

This directory is for storing standalone table files if you prefer to separate them from the main text.

## Usage

You can create separate `.tex` files for complex tables and input them:

```latex
\input{tables/table1_descriptive.tex}
```

## Example Table File

Create a file `table1_descriptive.tex`:

```latex
\begin{table}[H]
    \centering
    \caption{Descriptive Statistics}
    \label{tab:descriptive}
    \begin{tabular}{lcccc}
        \toprule
        Variable & Mean & Std. Dev. & Min & Max \\
        \midrule
        Variable 1 & 0.00 & 0.00 & 0.00 & 0.00 \\
        \bottomrule
    \end{tabular}
\end{table}
```
