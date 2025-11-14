# ============================================================================
# Makefile for Functional Breakpoint Detection
# ============================================================================

CC = gcc
CFLAGS = -O3 -Wall -std=c11 -fopenmp
LDFLAGS = -lm -lgomp

SRC_DIR = C
BIN_DIR = .

# 对象文件
UTILS_OBJ = $(SRC_DIR)/utils.o

# 目标文件
TARGETS = $(BIN_DIR)/char_func $(BIN_DIR)/break_detect $(BIN_DIR)/bootstrap

all: $(TARGETS)
	@echo "=========================================="
	@echo "所有程序编译完成!"
	@echo "=========================================="

# 编译utils对象文件
$(SRC_DIR)/utils.o: $(SRC_DIR)/utils.c $(SRC_DIR)/utils.h
	@echo "编译 utils.c..."
	$(CC) $(CFLAGS) -c -o $@ $<
	@echo "  √ 完成: utils.o"

# 特征函数计算
$(BIN_DIR)/char_func: $(SRC_DIR)/characteristic_function.c
	@echo "编译 characteristic_function.c..."
	$(CC) $(CFLAGS) -o $@.exe $< $(LDFLAGS)
	@echo "  √ 完成: $@.exe"

# 断点检测
$(BIN_DIR)/break_detect: $(SRC_DIR)/breakpoint_detection.c $(UTILS_OBJ)
	@echo "编译 breakpoint_detection.c..."
	$(CC) $(CFLAGS) -o $@.exe $^ $(LDFLAGS)
	@echo "  √ 完成: $@.exe"

# Bootstrap
$(BIN_DIR)/bootstrap: $(SRC_DIR)/bootstrap.c $(UTILS_OBJ)
	@echo "编译 bootstrap.c..."
	$(CC) $(CFLAGS) -o $@.exe $^ $(LDFLAGS)
	@echo "  √ 完成: $@.exe"

clean:
	@echo "清理编译文件..."
	rm -f $(BIN_DIR)/*.exe
	rm -f $(SRC_DIR)/*.o
	rm -f data/temp*
	@echo "  √ 清理完成"

dirs:
	mkdir -p data results

test: all dirs
	@echo "=========================================="
	@echo "运行测试..."
	@echo "=========================================="
	Rscript -e "source('R/03_main_analysis.R'); complete_analysis(T=300, breakpoints=c(0.5), run_bootstrap=FALSE)"

help:
	@echo "可用目标:"
	@echo "  make all      - 编译所有程序"
	@echo "  make clean    - 清理编译文件"
	@echo "  make dirs     - 创建必要目录"
	@echo "  make test     - 运行测试"
	@echo "  make help     - 显示此帮助"

.PHONY: all clean dirs test help
