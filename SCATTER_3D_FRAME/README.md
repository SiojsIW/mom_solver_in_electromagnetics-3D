# 3D MoM 电磁求解器

从零开始构建的矩量法（MoM）电磁散射求解器。

## 当前进度

- [x] RWG 基函数构建
- [x] 网格几何处理（法向、边提取）
- [x] 高斯数值积分
- [x] 格林函数计算
- [x] EFIE/MFIE 非对角阻抗矩阵
- [x] CFIE 组合与右端项
- [ ] 矩阵直接求解（LU分解）
- [ ] 矩阵迭代求解（GMRES）
- [ ] 快速算法（MLFMM/ACA）
- [ ] 并行计算（OpenMP/MPI）

## 编译方法

```bash
cd SCATTER_3D_FRAME
mingw32-make