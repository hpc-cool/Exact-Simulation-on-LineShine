#ifndef __ZARM_H__

#define __ZARM_H__
#ifdef __ARM_ACLE
#include <arm_acle.h>
#endif /* __ARM_ACLE */
#ifdef __ARM_FEATURE_FP16_SCALAR_ARITHMETIC
#include <arm_fp16.h>
#endif /* __ARM_FEATURE_FP16_SCALAR_ARITHMETIC */
#ifdef __ARM_FEATURE_BF16
#include <arm_bf16.h>
#endif /* __ARM_FEATURE_BF16 */

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif /* __ARM_NEON */
#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */
#ifdef __ARM_NEON_SVE_BRIDGE
#include <arm_neon_sve_bridge.h>
#endif /* __ARM_NEON_SVE_BRIDGE */
#if (__ARM_FEATURE_MVE & 3) == 3
#include <arm_mve.h>
/* MVE integer and floating point intrinsics are now available to use. */
#elif __ARM_FEATURE_MVE & 1
#include <arm_mve.h>
/* MVE integer intrinsics are now available to use. */
#endif
#ifdef __ARM_FEATURE_SME
//#include <arm_sme_draft_spec_subject_to_change.h>
#include <arm_sme.h>
#endif
#define DOUBLE float64_t
#define VDOUBLE svfloat64_t
#define FLOAT float32_t
#define VBOOL svbool_t
#define VINT svint32_t

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
