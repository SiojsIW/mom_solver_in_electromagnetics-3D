# 编译器
FC = gfortran
# 编译选项：检查所有警告、启用浮点陷阱
FFLAGS = -Wall -Wextra -fcheck=all

# 源文件列表（按依赖顺序写更保险，但 gfortran 一般能自动排）
# 添加新文件：在 MAIN.F90 之前插入一行，行尾加反斜杠 \，最后一行不要加\

SRC = EM_TYPES.f90 \
      MESH_GEOMETRY.f90 \
      NUMERICAL_INTEGRATION.f90 \
      GREEN_FUNCTIONS.f90 \
      RWG_BASIS_BUILD.f90 \
      Z_MATRIX.f90 \
      RHS.f90 \
      MAIN.f90

# 目标可执行文件
TARGET = mom_solver

# 默认规则：编译所有
all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(TARGET)

# 编译并运行（一键执行）
run: $(TARGET)
	./$(TARGET)

# 清理编译垃圾
clean:
	rm -f $(TARGET) *.o *.mod

.PHONY: all run clean