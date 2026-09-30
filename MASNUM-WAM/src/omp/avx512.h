#ifndef __ZAVX_H__
#define __ZAVX_H__
#include <immintrin.h>

#define DOUBLE double
#define VDOUBLE __m512d 
#define FLOAT float
#define VFLOAT __m512
#define VINT __m512i

#define ADD(A,B) svadd_f64_m(pm,A,B)
#define SUB(A,B) svsub_f64_m(pm,A,B)
#define MLA(A,B,C) svmla_f64_m(pm,A,B,C)
#define MUL(A,B) svmul_f64_m(pm,A,B)
#define VVD(P) *(svfloat64_t*)((DOUBLE*)(P))
#define VSUM(v) svaddv_f64(pm,v)
#define VLD(P) svld1_f64(pm,P)
#define VLFD(P) svcvt_f64_f32_z(pm,svld1_f32(pm,P))
#define VST(P,v) svst1_f64(pm,P,v)
#define VFABS(v) svabs_f64_x(pm,v)
#define VDUP(A) svdup_f64(A)
#define VEXP(v) svml_exp_f64((v))
#define VPOW(v,e) svml_pow_f64((v),(e))

#define LDZA(tile, slice, em, pm) svld1_hor_za64(tile, slice, pm, em)
#define STZA(tile, slice, em, pm) svst1_hor_za64(tile, slice, pm, em)
#define MOPA(tile,pg_m,pg_n,v0,v1) svmopa_za64_m(tile, pg_m, pg_n, v0, v1)
#define MUL_N(v, c) svmul_n_f64_m(pm, v, c)
#define VTRUE svptrue_b64()
#endif

