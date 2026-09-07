// ============================================================================
// DIRECTFN_WRAPPER.cpp
//
// Fortran 可调用的 DIRECTFN 三角形奇异积分接口（extern "C"）。
//
// 供 SCATTER_3D_FRAME（Fortran MoM 求解器）通过 ISO_C_BINDING 调用，
// 作为 EFIE 近场积分的第三种处理方法（DIRECTFN 全数值法）。
//
// 入口函数：dfn_tri_iss
//
//   输入：
//     adj     - 邻接类型：3 = 同一三角形(ST, 3点)
//                          2 = 共边(EA, 4点)
//                          1 = 共顶点(VA, 5点)
//     pts     - contour 顶点坐标，DIRECTFN 约定顺序，每个点 3 个 double：
//               ST: t1 t2 t3
//               EA: 公共边端点 e1, e2, 观察三角形第三顶点, 源三角形第三顶点
//               VA: 公共顶点, 观察三角形其余两顶点, 源三角形其余两顶点
//     k0      - 波数
//     n_gauss - 一维高斯-勒让德阶数（4 维积分各维相同）
//
//   输出（复数均以 (re, im) 连续存放，与 Fortran COMPLEX(C_DOUBLE_COMPLEX) 兼容）：
//     iss_rwg   - 9 个值，线性布局为 (源槽位, 观察槽位) 行优先：
//                 iss_rwg[(a-1) + 3*(b-1)] =
//                 ∬ f_b^obs · f_a^src * G dS_obs dS_src，
//                 其中 f 为标准 RWG 基函数 l/(2A)*(r - v_free)，
//                 a = DIRECTFN 源三角形局部顶点号(1..3)，b = 观察三角形(1..3)
//                 （依据：directfn_kernel_tri.cpp 的 pmf_rwg_val_ 分量顺序
//                   f1g1,f1g2,f1g3,f2g1,...，f 来自 rp=源侧、g 来自 rq=观察侧；
//                   已用独立闭式参考对 VA/EA 全 9 分量验证至 1e-11）
//                 Fortran 侧取 (观察=IO, 源=JO) 分量：ISS_RWG(JO + 3*(IO-1))
//     iss_const - 1 个值，∬ G dS_obs dS_src，
//                 G = exp(-j*k*R) / (4*pi*R)
//
//   返回 0 表示成功，非 0 表示邻接类型非法。
// ============================================================================

#include <cstring>

#include "directfn_algorithm_st.h"
#include "directfn_algorithm_ea.h"
#include "directfn_algorithm_va.h"
#include "directfn_contour.h"
#include "directfn_kernel_tri.h"

using namespace Directfn;

namespace {

// 复用的算法对象（性能：避免每次调用都构造对象并重新生成 Gauss 节点；
// 库的 set(contour) 只转发给 kernel，不清除已设的高斯阶数，可安全重用）。
// 每个对象独立跟踪 (k0, n)：按邻接类型使用不同阶数时互不复位。
struct AlgState { double k0; int n; };

template <typename Alg>
void sync_one(Alg &a, AlgState &s, const double k0, const int n) {
    if (s.k0 == k0 && s.n == n) return;
    s.k0 = k0; s.n = n;
    a.set_wavenumber(k0);
    a.set_Gaussian_orders_4(n, n, n, n);
}

struct DfnCtx {
    Triangular_ST<TriangularKernel_RWG_WS>      st_rwg;
    Triangular_EA<TriangularKernel_RWG_WS>      ea_rwg;
    Triangular_VA<TriangularKernel_RWG_WS>      va_rwg;
    Triangular_ST<TriangularKernel_Constant_ST> st_cst;
    Triangular_EA<TriangularKernel_Constant_EA> ea_cst;
    Triangular_VA<TriangularKernel_Constant_VA> va_cst;
    AlgState s_st_rwg, s_ea_rwg, s_va_rwg, s_st_cst, s_ea_cst, s_va_cst;
    DfnCtx() : s_st_rwg{-1.0, -1}, s_ea_rwg{-1.0, -1}, s_va_rwg{-1.0, -1},
               s_st_cst{-1.0, -1}, s_ea_cst{-1.0, -1}, s_va_cst{-1.0, -1} {}
};

DfnCtx &ctx() {
    static DfnCtx c;
    return c;
}

// 按邻接类型计算 RWG_WS 核的 Iss（9 个分量）
void calc_rwg_ws(const int adj, const SingularContour3xn &cntr,
                 const double k0, const size_t N, dcomplex *out) {
    DfnCtx &c = ctx();
    const int n = static_cast<int>(N);
    switch (adj) {
    case 3:
        sync_one(c.st_rwg, c.s_st_rwg, k0, n);
        c.st_rwg.set(cntr);
        c.st_rwg.calc_Iss();
        std::memcpy(out, c.st_rwg.Iss(), 9 * sizeof(dcomplex));
        break;
    case 2:
        sync_one(c.ea_rwg, c.s_ea_rwg, k0, n);
        c.ea_rwg.set(cntr);
        c.ea_rwg.calc_Iss();
        std::memcpy(out, c.ea_rwg.Iss(), 9 * sizeof(dcomplex));
        break;
    default:
        sync_one(c.va_rwg, c.s_va_rwg, k0, n);
        c.va_rwg.set(cntr);
        c.va_rwg.calc_Iss();
        std::memcpy(out, c.va_rwg.Iss(), 9 * sizeof(dcomplex));
        break;
    }
}

// 按邻接类型计算 Constant 核的 Iss（1 个分量）
void calc_constant(const int adj, const SingularContour3xn &cntr,
                   const double k0, const size_t N, dcomplex *out) {
    DfnCtx &c = ctx();
    const int n = static_cast<int>(N);
    switch (adj) {
    case 3:
        sync_one(c.st_cst, c.s_st_cst, k0, n);
        c.st_cst.set(cntr);
        c.st_cst.calc_Iss();
        out[0] = c.st_cst.Iss()[0];
        break;
    case 2:
        sync_one(c.ea_cst, c.s_ea_cst, k0, n);
        c.ea_cst.set(cntr);
        c.ea_cst.calc_Iss();
        out[0] = c.ea_cst.Iss()[0];
        break;
    default:
        sync_one(c.va_cst, c.s_va_cst, k0, n);
        c.va_cst.set(cntr);
        c.va_cst.calc_Iss();
        out[0] = c.va_cst.Iss()[0];
        break;
    }
}

} // anonymous namespace

extern "C" int dfn_tri_iss(const int adj, const double *pts, const double k0,
                           const int n_gauss, double *iss_rwg, double *iss_const) {
    if (adj < 1 || adj > 3) return 1;

    SingularContour3xn cntr;
    if (adj == 3) {
        cntr.set_points(pts,     pts + 3, pts + 6);
    } else if (adj == 2) {
        cntr.set_points(pts,     pts + 3, pts + 6, pts + 9);
    } else {
        cntr.set_points(pts,     pts + 3, pts + 6, pts + 9, pts + 12);
    }

    const size_t N = static_cast<size_t>(n_gauss);
    calc_rwg_ws(adj, cntr, k0, N, reinterpret_cast<dcomplex *>(iss_rwg));
    calc_constant(adj, cntr, k0, N, reinterpret_cast<dcomplex *>(iss_const));
    return 0;
}
