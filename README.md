# 3D MoM Electromagnetic Solver

从零开始构建的矩量法（MoM）电磁散射求解器。

## 项目结构

- `SCATTER_3D_FRAME/` — 核心求解器源码（Fortran）
  - 8 个模块：EM_TYPES, MESH_GEOMETRY, NUMERICAL_INTEGRATION, GREEN_FUNCTIONS, RWG_BASIS_BUILD, Z_MATRIX, RHS, MAIN
  - Makefile 一键编译

## 当前进度

详见 [SCATTER_3D_FRAME/README.md](SCATTER_3D_FRAME/README.md)

## 构建

```bash
cd SCATTER_3D_FRAME
mingw32-make