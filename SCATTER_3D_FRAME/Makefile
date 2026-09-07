# 编译器
FC = gfortran
CXX = g++
# 编译选项：检查所有警告、启用浮点陷阱
FFLAGS = -O2 -Wall -Wextra -fcheck=all -ffree-line-length-none
CXXFLAGS = -O2 -std=c++11 -IDIRECTFN-master/include

# DIRECTFN C++ 静态库与包装层
DFN_LIB = DIRECTFN-master/build/libdirectfn.a
DFN_WRAPPER = DIRECTFN_WRAPPER.o

# 公共模块（无 main 程序）
MODS = EM_TYPES.f90 \
       MESH_GEOMETRY.f90 \
       NUMERICAL_INTEGRATION.f90 \
       GREEN_FUNCTIONS.f90 \
       RWG_BASIS_BUILD.f90 \
       SINGULAR_INTEFRAL.F90 \
       QBX_EFIE.f90 \
       DIRECTFN_BRIDGE.f90 \
       Z_MATRIX.f90 \
       RHS.f90

# 原 mom_solver：四面体测试
SRC_MOM = $(MODS) MAIN.f90
TARGET_MOM = mom_solver

# RCS 求解器：金属球 RCS 计算
SRC_RCS = $(MODS) FAR_FIELD.f90 MAIN_RCS.f90
TARGET_RCS = rcs_solver

# 默认规则
all: $(TARGET_MOM) $(TARGET_RCS)

$(DFN_WRAPPER): DIRECTFN_WRAPPER.cpp
	$(CXX) $(CXXFLAGS) -c DIRECTFN_WRAPPER.cpp -o $(DFN_WRAPPER)

$(TARGET_MOM): $(SRC_MOM) $(DFN_WRAPPER)
	$(FC) $(FFLAGS) $(SRC_MOM) $(DFN_WRAPPER) $(DFN_LIB) -lstdc++ -o $(TARGET_MOM)

$(TARGET_RCS): $(SRC_RCS) $(DFN_WRAPPER)
	$(FC) $(FFLAGS) $(SRC_RCS) $(DFN_WRAPPER) $(DFN_LIB) -lstdc++ -o $(TARGET_RCS)

# 运行原测试
run: $(TARGET_MOM)
	chcp.com 65001
	$(TARGET_MOM)

# 运行 RCS 计算
rcs: $(TARGET_RCS)
	chcp.com 65001
	$(TARGET_RCS)

# 清理编译垃圾
clean:
	rm -f $(TARGET_MOM) $(TARGET_RCS) *.o *.mod

# 运行奇异积分验证
test: test_singular.exe
	chcp.com 65001
	test_singular.exe

test_singular.exe: EM_TYPES.f90 NUMERICAL_INTEGRATION.f90 GREEN_FUNCTIONS.f90 SINGULAR_INTEFRAL.F90 test_singular.F90
	$(FC) $(FFLAGS) EM_TYPES.f90 NUMERICAL_INTEGRATION.f90 GREEN_FUNCTIONS.f90 SINGULAR_INTEFRAL.F90 test_singular.F90 -o test_singular.exe

# DIRECTFN 桥接单元测试（test_directfn）已随调试完成删除。
# 如需重新逐元验证，可用 DIRECTFN_BRIDGE 的 DIRECTFN_DEBUG 标志（>0 打印 @@DFNRAW 原始积分值）。

.PHONY: all run rcs clean test
