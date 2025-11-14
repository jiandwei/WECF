# ============================================================================
# Makefile for Functional Breakpoint Detection
# 编译器: MinGW gcc
# OpenMP: 支持
# ============================================================================

# 编译器设置
CC = gcc
CFLAGS = -O3 -Wall -std=c11 -fopenmp
LDFLAGS = -lm -lgomp

# 目录
SRC_DIR = C
BIN_DIR = .

# 目标文件
TARGETS = $(BIN_DIR)/char_func $(BIN_DIR)/break_detect $(BIN_DIR)/bootstrap

# 默认目标
all: $(TARGETS)
	@echo "=========================================="
	@echo "所有程序编译完成!"
	@echo "=========================================="

# 特征函数计算
$(BIN_DIR)/char_func: $(SRC_DIR)/characteristic_function.c
	@echo "编译 characteristic_function.c..."
	$(CC) $(CFLAGS) -o $@.exe $< $(LDFLAGS)
	@echo "  √ 完成: $@.exe"

# 断点检测
$(BIN_DIR)/break_detect: $(SRC_DIR)/breakpoint_detection.c
	@echo "编译 breakpoint_detection.c..."
	$(CC) $(CFLAGS) -o $@.exe $< $(LDFLAGS)
	@echo "  √ 完成: $@.exe"

# Bootstrap
$(BIN_DIR)/bootstrap: $(SRC_DIR)/bootstrap.c
	@echo "编译 bootstrap.c..."
	$(CC) $(CFLAGS) -o $@.exe $< $(LDFLAGS)
	@echo "  √ 完成: $@.exe"

# 清理
clean:
	@echo "清理编译文件..."
	rm -f $(BIN_DIR)/*.exe
	rm -f $(SRC_DIR)/*.o
	rm -f data/temp*
	@echo "  √ 清理完成"

# 创建必要目录
dirs:
	mkdir -p data results

# 运行测试
test: all dirs
	@echo "=========================================="
	@echo "运行测试..."
	@echo "=========================================="
	Rscript -e "source('R/03_main_analysis.R'); complete_analysis(T=300, breakpoints=c(0.5), run_bootstrap=FALSE)"

# 帮助
help:
	@echo "可用目标:"
	@echo "  make all      - 编译所有程序"
	@echo "  make clean    - 清理编译文件"
	@echo "  make dirs     - 创建必要目录"
	@echo "  make test     - 运行测试"
	@echo "  make help     - 显示此帮助"

.PHONY: all clean dirs test help
