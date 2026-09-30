#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#if defined(__ARM_FEATURE_SME)
	#include <arm_sve.h>
  #if (defined(__clang_major__) && (__clang_major__ >= 17))
    #include<arm_sme.h>
    #include<sys/time.h>
  #else
    #include <arm_sme_draft_spec_subject_to_change.h>
  #endif
#endif

#if defined(SME_VSCALE_CHECK_AND_SET_ENABLE) && (SME_VSCALE_CHECK_AND_SET_ENABLE == 1)
    #include <sys/prctl.h>
#endif
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define ZB0 0
#define ZB1 1
#define ZB2 2
#define ZB3 3

void open_sme()
{
	__asm__ volatile("SMSTART":::"z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10", "z11", "z12", 
    "z13", "z14", "z15", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", 
    "z29", "z30", "z31");
	__asm__ volatile("ISB");
#if defined(SME_VSCALE_CHECK_AND_SET_ENABLE) && (SME_VSCALE_CHECK_AND_SET_ENABLE == 1)
    int vscale;
    __asm__ volatile("CNTD %x0":"=r"(vscale)::); // 这里必须保证在smstart之后，探测SME下的向量长度，防止编译器优化到开启SME前
    if (vscale != 8) { // 检查vscale为8，保证Streaming SVE下的SVE长度512
        __asm__ volatile("SMSTOP" :::"z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10", "z11", "z12", 
        "z13", "z14", "z15", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", 
        "z29", "z30", "z31");
        __asm__ volatile("ISB");
        prctl(PR_SME_SET_VL, 512); // 实测小规模下prctl开销很高，仅长度异常时配置
        __asm__ volatile("SMSTART":::"z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10", "z11", "z12", 
        "z13", "z14", "z15", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", 
        "z29", "z30", "z31");
        __asm__ volatile("ISB");
    }
#endif
}

void close_sme()
{
    __asm__ volatile("SMSTOP":::"z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10", "z11", "z12", 
    "z13", "z14", "z15", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", 
    "z29", "z30", "z31");
    __asm__ volatile("ISB");
}
void stencil_calculation_0_(int *j_s, int *j_e, int *k_s, int *k_e,
    int *data_ld, int *data_hd,
    int *result_ld, int *result_hd,
    double *A, double *B,  double *coef_0, double *coef_1, double *coef_2, int *data_no_l_jk, int *result_no_l_jk)
{
    open_sme();
    int data_ld_size = *data_ld;
    int result_ld_size = *result_ld;
    int data_hd_size = *data_hd;
    int result_hd_size = *result_hd;
    int j_start = *j_s;
    int j_end = *j_e;
    int k_start = *k_s;
    int k_end = *k_e;
    int data_no_loop_jk = *data_no_l_jk;
    int result_no_loop_jk = *result_no_l_jk;
    int t, i, j, k;
    svbool_t pm, pred;
    svfloat64_t g_0,g_1,g_2;
    svfloat64_t s0_0,s0_1,s0_2;
    svfloat64_t s1_0,s1_1,s1_2;
    svfloat64_t s2_0,s2_1,s2_2;
    svfloat64_t s3_0,s3_1,s3_2;
    int h_point = 3;
    int v_point = 3;
    int ld_row = v_point + 7;
    pm = svwhilelt_b64_u64(0, 8);
    long ld_offset[8];
    long st_offset[8];
    for (t = 0; t < 8; ++t) {
        st_offset[t] = t * 1;
    }
    svint64_t sv_st_offset = svld1_s64(pm, &st_offset[0]);
    for (j = j_start; j < j_end - 7; j += 8) {
        for (k = k_start; k < k_end - 31; k += 32) {
            svzero_za();
            g_0 = svld1_f64(pm, &coef_0[0 * 8]);
            g_1 = svld1_f64(pm, &coef_1[0 * 8]);
            g_2 = svld1_f64(pm, &coef_2[0 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[1 * 8]);
            g_1 = svld1_f64(pm, &coef_1[1 * 8]);
            g_2 = svld1_f64(pm, &coef_2[1 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[2 * 8]);
            g_1 = svld1_f64(pm, &coef_1[2 * 8]);
            g_2 = svld1_f64(pm, &coef_2[2 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[3 * 8]);
            g_1 = svld1_f64(pm, &coef_1[3 * 8]);
            g_2 = svld1_f64(pm, &coef_2[3 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[4 * 8]);
            g_1 = svld1_f64(pm, &coef_1[4 * 8]);
            g_2 = svld1_f64(pm, &coef_2[4 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[5 * 8]);
            g_1 = svld1_f64(pm, &coef_1[5 * 8]);
            g_2 = svld1_f64(pm, &coef_2[5 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[6 * 8]);
            g_1 = svld1_f64(pm, &coef_1[6 * 8]);
            g_2 = svld1_f64(pm, &coef_2[6 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[7 * 8]);
            g_1 = svld1_f64(pm, &coef_1[7 * 8]);
            g_2 = svld1_f64(pm, &coef_2[7 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[8 * 8]);
            g_1 = svld1_f64(pm, &coef_1[8 * 8]);
            g_2 = svld1_f64(pm, &coef_2[8 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[9 * 8]);
            g_1 = svld1_f64(pm, &coef_1[9 * 8]);
            g_2 = svld1_f64(pm, &coef_2[9 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            for (int rp = 0; rp < 8; ++rp) {
                svst1_hor_za64(ZB0, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0) * 1]);
                svst1_hor_za64(ZB1, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 8) * 1]);
                svst1_hor_za64(ZB2, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 16) * 1]);
                svst1_hor_za64(ZB3, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 24) * 1]);
            }
        }
        for (; k < k_end - 15; k += 16) {
            svzero_za();
            g_0 = svld1_f64(pm, &coef_0[0 * 8]);
            g_1 = svld1_f64(pm, &coef_1[0 * 8]);
            g_2 = svld1_f64(pm, &coef_2[0 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[1 * 8]);
            g_1 = svld1_f64(pm, &coef_1[1 * 8]);
            g_2 = svld1_f64(pm, &coef_2[1 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[2 * 8]);
            g_1 = svld1_f64(pm, &coef_1[2 * 8]);
            g_2 = svld1_f64(pm, &coef_2[2 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[3 * 8]);
            g_1 = svld1_f64(pm, &coef_1[3 * 8]);
            g_2 = svld1_f64(pm, &coef_2[3 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[4 * 8]);
            g_1 = svld1_f64(pm, &coef_1[4 * 8]);
            g_2 = svld1_f64(pm, &coef_2[4 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[5 * 8]);
            g_1 = svld1_f64(pm, &coef_1[5 * 8]);
            g_2 = svld1_f64(pm, &coef_2[5 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[6 * 8]);
            g_1 = svld1_f64(pm, &coef_1[6 * 8]);
            g_2 = svld1_f64(pm, &coef_2[6 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[7 * 8]);
            g_1 = svld1_f64(pm, &coef_1[7 * 8]);
            g_2 = svld1_f64(pm, &coef_2[7 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[8 * 8]);
            g_1 = svld1_f64(pm, &coef_1[8 * 8]);
            g_2 = svld1_f64(pm, &coef_2[8 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[9 * 8]);
            g_1 = svld1_f64(pm, &coef_1[9 * 8]);
            g_2 = svld1_f64(pm, &coef_2[9 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            for (int rp = 0; rp < 8; ++rp) {
                svst1_hor_za64(ZB0, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0) * 1]);
                svst1_hor_za64(ZB1, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 8) * 1]);
            }
        }
        for (; k < k_end; k += 8) {
            pred = svwhilelt_b64_u64(0, MIN(8, k_end - k));
            svzero_za();
            g_0 = svld1_f64(pm, &coef_0[0 * 8]);
            g_1 = svld1_f64(pm, &coef_1[0 * 8]);
            g_2 = svld1_f64(pm, &coef_2[0 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[1 * 8]);
            g_1 = svld1_f64(pm, &coef_1[1 * 8]);
            g_2 = svld1_f64(pm, &coef_2[1 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[2 * 8]);
            g_1 = svld1_f64(pm, &coef_1[2 * 8]);
            g_2 = svld1_f64(pm, &coef_2[2 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[3 * 8]);
            g_1 = svld1_f64(pm, &coef_1[3 * 8]);
            g_2 = svld1_f64(pm, &coef_2[3 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[4 * 8]);
            g_1 = svld1_f64(pm, &coef_1[4 * 8]);
            g_2 = svld1_f64(pm, &coef_2[4 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[5 * 8]);
            g_1 = svld1_f64(pm, &coef_1[5 * 8]);
            g_2 = svld1_f64(pm, &coef_2[5 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[6 * 8]);
            g_1 = svld1_f64(pm, &coef_1[6 * 8]);
            g_2 = svld1_f64(pm, &coef_2[6 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[7 * 8]);
            g_1 = svld1_f64(pm, &coef_1[7 * 8]);
            g_2 = svld1_f64(pm, &coef_2[7 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[8 * 8]);
            g_1 = svld1_f64(pm, &coef_1[8 * 8]);
            g_2 = svld1_f64(pm, &coef_2[8 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[9 * 8]);
            g_1 = svld1_f64(pm, &coef_1[9 * 8]);
            g_2 = svld1_f64(pm, &coef_2[9 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            for (int rp = 0; rp < 8; ++rp) {
                svst1_hor_za64(ZB0, rp, pred, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0) * 1]);
            }
        }
    }
    for (; j < j_end; ++j) {
        for (k = k_start; k < k_end; k += 8) {
            pred = svwhilelt_b64_u64(0, MIN(8, (k_end - k)));
            s1_0 = svdup_f64(0.);
            g_0 = svdup_f64(coef_0[0 * 8]);
            g_1 = svdup_f64(coef_1[0 * 8]);
            g_2 = svdup_f64(coef_2[0 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            s0_0 = svmul_f64_m(pred, s0_0, g_0);
            s0_1 = svmul_f64_m(pred, s0_1, g_1);
            s0_2 = svmul_f64_m(pred, s0_2, g_2);
            s1_0 = svadd_f64_m(pred, s1_0, s0_0);
            s1_0 = svadd_f64_m(pred, s1_0, s0_1);
            s1_0 = svadd_f64_m(pred, s1_0, s0_2);

            g_0 = svdup_f64(coef_0[1 * 8]);
            g_1 = svdup_f64(coef_1[1 * 8]);
            g_2 = svdup_f64(coef_2[1 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            s0_0 = svmul_f64_m(pred, s0_0, g_0);
            s0_1 = svmul_f64_m(pred, s0_1, g_1);
            s0_2 = svmul_f64_m(pred, s0_2, g_2);
            s1_0 = svadd_f64_m(pred, s1_0, s0_0);
            s1_0 = svadd_f64_m(pred, s1_0, s0_1);
            s1_0 = svadd_f64_m(pred, s1_0, s0_2);

            g_0 = svdup_f64(coef_0[2 * 8]);
            g_1 = svdup_f64(coef_1[2 * 8]);
            g_2 = svdup_f64(coef_2[2 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            s0_0 = svmul_f64_m(pred, s0_0, g_0);
            s0_1 = svmul_f64_m(pred, s0_1, g_1);
            s0_2 = svmul_f64_m(pred, s0_2, g_2);
            s1_0 = svadd_f64_m(pred, s1_0, s0_0);
            s1_0 = svadd_f64_m(pred, s1_0, s0_1);
            s1_0 = svadd_f64_m(pred, s1_0, s0_2);

            svst1_f64(pred, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0) * result_ld_size + (k + 0) * 1], s1_0);
        }
    }
    close_sme();
}

void stencil_calculation_1_(int *j_s, int *j_e, int *k_s, int *k_e,
    int *data_ld, int *data_hd,
    int *result_ld, int *result_hd,
    double *A, double *B,  double *coef_0, double *coef_1, double *coef_2, int *data_no_l_jk, int *result_no_l_jk)
{
    open_sme();
    int data_ld_size = *data_ld;
    int result_ld_size = *result_ld;
    int data_hd_size = *data_hd;
    int result_hd_size = *result_hd;
    int j_start = *j_s;
    int j_end = *j_e;
    int k_start = *k_s;
    int k_end = *k_e;
    int data_no_loop_jk = *data_no_l_jk;
    int result_no_loop_jk = *result_no_l_jk;
    int t, i, j, k;
    svbool_t pm, pred;
    svfloat64_t g_0,g_1,g_2;
    svfloat64_t s0_0,s0_1,s0_2;
    svfloat64_t s1_0,s1_1,s1_2;
    svfloat64_t s2_0,s2_1,s2_2;
    svfloat64_t s3_0,s3_1,s3_2;
    int h_point = 3;
    int v_point = 3;
    int ld_row = v_point + 7;
    pm = svwhilelt_b64_u64(0, 8);
    long ld_offset[8];
    long st_offset[8];
    for (t = 0; t < 8; ++t) {
        st_offset[t] = t * 1;
    }
    svint64_t sv_st_offset = svld1_s64(pm, &st_offset[0]);
    for (j = j_start; j < j_end - 7; j += 8) {
        for (k = k_start; k < k_end - 31; k += 32) {
            svzero_za();
            g_0 = svld1_f64(pm, &coef_0[0 * 8]);
            g_1 = svld1_f64(pm, &coef_1[0 * 8]);
            g_2 = svld1_f64(pm, &coef_2[0 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[1 * 8]);
            g_1 = svld1_f64(pm, &coef_1[1 * 8]);
            g_2 = svld1_f64(pm, &coef_2[1 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[2 * 8]);
            g_1 = svld1_f64(pm, &coef_1[2 * 8]);
            g_2 = svld1_f64(pm, &coef_2[2 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[3 * 8]);
            g_1 = svld1_f64(pm, &coef_1[3 * 8]);
            g_2 = svld1_f64(pm, &coef_2[3 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[4 * 8]);
            g_1 = svld1_f64(pm, &coef_1[4 * 8]);
            g_2 = svld1_f64(pm, &coef_2[4 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[5 * 8]);
            g_1 = svld1_f64(pm, &coef_1[5 * 8]);
            g_2 = svld1_f64(pm, &coef_2[5 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[6 * 8]);
            g_1 = svld1_f64(pm, &coef_1[6 * 8]);
            g_2 = svld1_f64(pm, &coef_2[6 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[7 * 8]);
            g_1 = svld1_f64(pm, &coef_1[7 * 8]);
            g_2 = svld1_f64(pm, &coef_2[7 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[8 * 8]);
            g_1 = svld1_f64(pm, &coef_1[8 * 8]);
            g_2 = svld1_f64(pm, &coef_2[8 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            g_0 = svld1_f64(pm, &coef_0[9 * 8]);
            g_1 = svld1_f64(pm, &coef_1[9 * 8]);
            g_2 = svld1_f64(pm, &coef_2[9 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 0]);
            s2_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 16 + 0]);
            s3_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 24 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 1]);
            s2_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 16 + 1]);
            s3_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 24 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 2]);
            s2_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 16 + 2]);
            s3_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 24 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB2, pm, pm, g_0, s2_0);
            svmopa_za64_f64_m(ZB3, pm, pm, g_0, s3_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB2, pm, pm, g_1, s2_1);
            svmopa_za64_f64_m(ZB3, pm, pm, g_1, s3_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);
            svmopa_za64_f64_m(ZB2, pm, pm, g_2, s2_2);
            svmopa_za64_f64_m(ZB3, pm, pm, g_2, s3_2);

            for (int rp = 0; rp < 8; ++rp) {
                svst1_hor_za64(ZB0, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0) * 1]);
                svst1_hor_za64(ZB1, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 8) * 1]);
                svst1_hor_za64(ZB2, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 16) * 1]);
                svst1_hor_za64(ZB3, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 24) * 1]);
            }
        }
        for (; k < k_end - 15; k += 16) {
            svzero_za();
            g_0 = svld1_f64(pm, &coef_0[0 * 8]);
            g_1 = svld1_f64(pm, &coef_1[0 * 8]);
            g_2 = svld1_f64(pm, &coef_2[0 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[1 * 8]);
            g_1 = svld1_f64(pm, &coef_1[1 * 8]);
            g_2 = svld1_f64(pm, &coef_2[1 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[2 * 8]);
            g_1 = svld1_f64(pm, &coef_1[2 * 8]);
            g_2 = svld1_f64(pm, &coef_2[2 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[3 * 8]);
            g_1 = svld1_f64(pm, &coef_1[3 * 8]);
            g_2 = svld1_f64(pm, &coef_2[3 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[4 * 8]);
            g_1 = svld1_f64(pm, &coef_1[4 * 8]);
            g_2 = svld1_f64(pm, &coef_2[4 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[5 * 8]);
            g_1 = svld1_f64(pm, &coef_1[5 * 8]);
            g_2 = svld1_f64(pm, &coef_2[5 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[6 * 8]);
            g_1 = svld1_f64(pm, &coef_1[6 * 8]);
            g_2 = svld1_f64(pm, &coef_2[6 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[7 * 8]);
            g_1 = svld1_f64(pm, &coef_1[7 * 8]);
            g_2 = svld1_f64(pm, &coef_2[7 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[8 * 8]);
            g_1 = svld1_f64(pm, &coef_1[8 * 8]);
            g_2 = svld1_f64(pm, &coef_2[8 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            g_0 = svld1_f64(pm, &coef_0[9 * 8]);
            g_1 = svld1_f64(pm, &coef_1[9 * 8]);
            g_2 = svld1_f64(pm, &coef_2[9 * 8]);
            s0_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 0]);
            s1_0 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 0]);
            s0_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 1]);
            s1_1 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 1]);
            s0_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 2]);
            s1_2 = svld1_f64(pm, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 8 + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB1, pm, pm, g_0, s1_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB1, pm, pm, g_1, s1_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);
            svmopa_za64_f64_m(ZB1, pm, pm, g_2, s1_2);

            for (int rp = 0; rp < 8; ++rp) {
                svst1_hor_za64(ZB0, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0) * 1]);
                svst1_hor_za64(ZB1, rp, pm, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0 + 8) * 1]);
            }
        }
        for (; k < k_end; k += 8) {
            pred = svwhilelt_b64_u64(0, MIN(8, k_end - k));
            svzero_za();
            g_0 = svld1_f64(pm, &coef_0[0 * 8]);
            g_1 = svld1_f64(pm, &coef_1[0 * 8]);
            g_2 = svld1_f64(pm, &coef_2[0 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[1 * 8]);
            g_1 = svld1_f64(pm, &coef_1[1 * 8]);
            g_2 = svld1_f64(pm, &coef_2[1 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[2 * 8]);
            g_1 = svld1_f64(pm, &coef_1[2 * 8]);
            g_2 = svld1_f64(pm, &coef_2[2 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[3 * 8]);
            g_1 = svld1_f64(pm, &coef_1[3 * 8]);
            g_2 = svld1_f64(pm, &coef_2[3 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 3) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[4 * 8]);
            g_1 = svld1_f64(pm, &coef_1[4 * 8]);
            g_2 = svld1_f64(pm, &coef_2[4 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 4) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[5 * 8]);
            g_1 = svld1_f64(pm, &coef_1[5 * 8]);
            g_2 = svld1_f64(pm, &coef_2[5 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 5) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[6 * 8]);
            g_1 = svld1_f64(pm, &coef_1[6 * 8]);
            g_2 = svld1_f64(pm, &coef_2[6 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 6) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[7 * 8]);
            g_1 = svld1_f64(pm, &coef_1[7 * 8]);
            g_2 = svld1_f64(pm, &coef_2[7 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 7) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[8 * 8]);
            g_1 = svld1_f64(pm, &coef_1[8 * 8]);
            g_2 = svld1_f64(pm, &coef_2[8 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 8) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            g_0 = svld1_f64(pm, &coef_0[9 * 8]);
            g_1 = svld1_f64(pm, &coef_1[9 * 8]);
            g_2 = svld1_f64(pm, &coef_2[9 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 9) * data_ld_size + k + (-1) + 2]);
            svmopa_za64_f64_m(ZB0, pm, pm, g_0, s0_0);
            svmopa_za64_f64_m(ZB0, pm, pm, g_1, s0_1);
            svmopa_za64_f64_m(ZB0, pm, pm, g_2, s0_2);

            for (int rp = 0; rp < 8; ++rp) {
                svst1_hor_za64(ZB0, rp, pred, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0 + rp) * result_ld_size + (k + 0) * 1]);
            }
        }
    }
    for (; j < j_end; ++j) {
        for (k = k_start; k < k_end; k += 8) {
            pred = svwhilelt_b64_u64(0, MIN(8, (k_end - k)));
            s1_0 = svdup_f64(0.);
            g_0 = svdup_f64(coef_0[0 * 8]);
            g_1 = svdup_f64(coef_1[0 * 8]);
            g_2 = svdup_f64(coef_2[0 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 0) * data_ld_size + k + (-1) + 2]);
            s0_0 = svmul_f64_m(pred, s0_0, g_0);
            s0_1 = svmul_f64_m(pred, s0_1, g_1);
            s0_2 = svmul_f64_m(pred, s0_2, g_2);
            s1_0 = svadd_f64_m(pred, s1_0, s0_0);
            s1_0 = svadd_f64_m(pred, s1_0, s0_1);
            s1_0 = svadd_f64_m(pred, s1_0, s0_2);

            g_0 = svdup_f64(coef_0[1 * 8]);
            g_1 = svdup_f64(coef_1[1 * 8]);
            g_2 = svdup_f64(coef_2[1 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 1) * data_ld_size + k + (-1) + 2]);
            s0_0 = svmul_f64_m(pred, s0_0, g_0);
            s0_1 = svmul_f64_m(pred, s0_1, g_1);
            s0_2 = svmul_f64_m(pred, s0_2, g_2);
            s1_0 = svadd_f64_m(pred, s1_0, s0_0);
            s1_0 = svadd_f64_m(pred, s1_0, s0_1);
            s1_0 = svadd_f64_m(pred, s1_0, s0_2);

            g_0 = svdup_f64(coef_0[2 * 8]);
            g_1 = svdup_f64(coef_1[2 * 8]);
            g_2 = svdup_f64(coef_2[2 * 8]);
            s0_0 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 0]);
            s0_1 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 1]);
            s0_2 = svld1_f64(pred, &A[(data_no_loop_jk+0) * data_ld_size * data_hd_size + (j + (-1) + 2) * data_ld_size + k + (-1) + 2]);
            s0_0 = svmul_f64_m(pred, s0_0, g_0);
            s0_1 = svmul_f64_m(pred, s0_1, g_1);
            s0_2 = svmul_f64_m(pred, s0_2, g_2);
            s1_0 = svadd_f64_m(pred, s1_0, s0_0);
            s1_0 = svadd_f64_m(pred, s1_0, s0_1);
            s1_0 = svadd_f64_m(pred, s1_0, s0_2);

            svst1_f64(pred, &B[(result_no_loop_jk+0) * result_ld_size * result_hd_size + (j + 0) * result_ld_size + (k + 0) * 1], s1_0);
        }
    }
    close_sme();
}
