# Packages Directory

This directory is for custom LaTeX packages or style files.

## Usage

If you create custom style files, place them here and load them in `main.tex`:

```latex
\usepackage{packages/mystyle}
```

## Example Custom Package

You might create files like:
- `econometrics.sty` - Custom commands for econometrics
- `theorems.sty` - Custom theorem environments
- `figures.sty` - Custom figure settings

## Note

For standard packages, simply use `\usepackage{packagename}` in the preamble. Only place custom or local packages here.
