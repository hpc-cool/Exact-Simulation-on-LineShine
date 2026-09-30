#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>

#define NEWN(T,N) (T*)malloc(sizeof(T)*(N))
#define NEWPN(p,N) *(void**)&p=malloc(sizeof(*p)*(N))
#define NEWZN(p,N)   *(void**)&p=malloc(sizeof(*p)*(N));memset(p,0,sizeof(*p)*(N))
#define VFREE(p) if(p)free(p)
/*###################################################################################################
   *            Copyright (C) 2025  Wei Zhao
   *            MODULE NAME : irrp_split
   * --- NOTE for describing of subroutine / function :
   *  B. The describe for the parameters of subroutine / function, started with:
   *   * It means input prameter;
   *   # It means output prameter;
   *   @ It means input and output prameter(it will be changed inside).
   *
   *-------------------------------------------------------------------------------------------------
   * ***                                 INTERFACE DESCRIBE                                       ***
   *
   *  1. sub. irrp_set_pemask : Do irregular partition using MASK and return the results in PEMASK.
   *
   *   irrp_set_pemask(mask_, pemask_, im_, jm_, npe_)
   *
   *   * int  mask_   = Mask for horizonal space/points in which 0 means land/useless points.
   *   # int  pemask_ = Arrary to save partition reaults each point with the value equal to
   *                          the ID of the PE that this point is belonged. This ID of PEs is
   *                          started from 1 which is different with the way for MPI.
   *   * int  im_     = The first dimension size of mask_ and pemask_.
   *   * int  jm_     = The second dimension size of mask_ and pemask_.
   *   * int  npe_    = The total number of PEs used for this partition.
   *    The shape for each partition will be like the following.
   *    !-----------------------------------------------------------------------------------!
   *    !                                                                                   !
   *    !                      Shape of points in each PE.                                  !
   *    !                                                                                   !
   *    !    (i1, j2) ------------------------(i4, j2)                                      !
   *    !        +    -----------------------------------------------(i6, j6)               !
   *    !        +    -----------------------------------------------    +                  !
   *    !        +    -----------------------------------------------    +                  !
   *    !        +    -----------------------------------------------    +                  !
   *    !    (i1, j3) -----------------------------------------------    +                  !
   *    !             +    ------------------------------------------    +                  !
   *    !             +    ------------------------------------------    +                  !
   *    !             +    -------------------------------------------------- (i2, j4)      !
   *    !             +    ---------------------------------------------------    +         !
   *    !             +    ---------------------------------------------------    +         !
   *    !             +    ---------------------------------------------------    +         !
   *    !             +    ---------------------------------------------------    +         !
   *    !          (i5, j5) --------------------------------------------------    +         !
   *    !                               (i3, j1) ----------------------------- (i2, j1)     !
   *    !                                                                                   !
   *    !-----------------------------------------------------------------------------------!*/
/*
   irrp_set_pemask:
     set_pemask_iteration:
       loop:
         set_pemask_one_cycle:
           compute_ratio:
           set_pemask_columns:
 * */
struct pi_pos_type {
  int i ;
  int j ;
} ;
typedef struct irrp_split{
  int pid  ;
  //MPI_Comm *comm;
  int im      ;  // First dimension size of MASK & PEMASK.
  int jm      ;  // Second dimension size of MASK & PEMASK.
  int npe     ;  // Number of PEs used in this partition.
  int sumnp   ;  // Sum of the computer points in MASK.
  int ncols   ;  // Number of columns during partition
  double avenp;  //  Averaged number of computing points for each PE.
  double savenp; //  sqrt(avenp);
  struct pi_pos_type *s2ij;       // Teporary arrary used to sequentializing the 2
                 // dimensional points. This will be help for the
                 // action of partition. The sequentializing is
                 // performed by connect all the column
  int* mask  ;   // The mask for 2 dimensional space in which 0 means
                 // land or the point is not need to be computed.
                 // The shape if this arrary is (im, jm).
  int* pemask;   //It has the same shape with mask and it is used to
                 // record the partition results for each points with
                 // the value equal to the ID of the PE these points
                 // are belonged.
                 // NOTE: This ID is started from 1 but not from zero.
                 //       It is different with the arranging method of
                 //       PIDs in MPI.
  float *iflag, *jflag, *niflag, *njflag;
}Irrp_Split;
/*-----------------------------------------------------------------------------------------------
 * sub. irrp_set_pemask : Do irregular partition using MASK and return the results in PEMASK.
 *    * int  mask_   = Mask for horizonal space/points: 0 means land/useless points.
 *    # int  pemask_ = Arrary to save partition reaults each point with the value equal to
 *                           the ID of the PE that this point is belonged. This ID of PEs is
 *                           started from 1 which is different with the way for MPI.
 *    * int  im_     = The first dimension size of mask_ and pemask_
 *    * int  jm_     = The second dimension size of mask_ and pemask_
 *    * int  npe_    = The total number of PEs used for this partition.
 *-----------------------------------------------------------------------------------------------
*/
void irrp_set_pemask(int *mask_,int*pemask_, int im_, int jm_, int npe_,int pid_);
void irrp_set_pemask_(int *mask_, int*pemask_,int *im_, int *jm_, int *npe_,int *pid_){
  irrp_set_pemask(mask_,pemask_,*im_,*jm_,*npe_,*pid_);
}
static void set_pemask_iteration(Irrp_Split*pis, int flag, double delta_, double *ratio_);
void irrp_set_pemask(int *mask_,int*pemask_, int im_, int jm_, int npe_,int pid_){
  Irrp_Split is={0};
  is.im = im_; is.jm = jm_; is.npe = npe_;// Record the dimension size and number of PEs in the
  is.pid=pid_;
  is.mask = mask_; 
  is.pemask = pemask_;       // Using pointer to refer to the arrary of mask and
  memset(is.pemask,0,is.im*is.jm);
  NEWZN(is.s2ij,is.im*is.jm+1);
  int nn=0;
  struct pi_pos_type *pij=is.s2ij;
  for(int j=0;j<is.jm;j++){
    for(int i=0;i<is.im;i++){
      if(is.mask[i+j*is.im]){
        nn++;pij++;
        pij->i=i;pij->j=j;
      }
    }
  }
  is.sumnp=nn;
  is.avenp = nn / (float)(is.npe); // Averaged computing points for each PE.
  is.savenp=sqrt(is.avenp);
  set_pemask_iteration(&is,1, 0.,NULL); // Do partition by iteration based on the gridient of
                                 // the global maximum aspect ratio of each PE respect
                                 // to parameter of delta.
  VFREE(is.s2ij);
}
/*------------------------------------------------------------------------------------------------
* sub. set_pemask_columns: Set pemask from n1 to n2 as a column of PEs.
*   set_pemask_columns(n1, n2)
*   * int  n1 = The start ID of PEs needed to be set for PEMASK.
*   * int  n2 = The end ID of PEs needed to be set for PEMASK.
*-------------------------------------------------------------------------------------------------
*/
static int set_pemask_columns(Irrp_Split*pis,int n1, int n2){
  int i, nwp, s1, s2, i1, i2, j1, j2, idpe, j, ib, ie, nwpe;
  /* set i1, i2, j1, j2 for points: s1 --> s2, pe: n1 --> n2
   * The start point which should belong to PE with ID of n1.
   */
  s1 = pis->avenp * (n1-1) +0.5+ 1.e-5+1  ;// the end of last column
  /* The end point which should belong to PE with ID of n2.
   * For last PE, s2 should same as sum number of computing points.
   */
  s2 = pis->avenp *  n2 +0.5+ 1.e-5;
  if(n2 == pis->npe)s2 = pis->sumnp;
  /* * Set i1 & j1: the index in matrix for the start point of s1. */
  if(s1 <= 1){
    i1 = 0; j1 = 0;
  }else{
    i1 = pis->s2ij[s1].i; j1 = pis->s2ij[s1].j;
    if(j1 >=pis->jm){
      i1 = i1 + 1; j1 = 0;
    }
  }
  /*-----------------------------------------------------------------------------------*
   *                                                                                   *
   *                      Shape of points in each PE.                                  *
   *                                                                                   *
   *    (i1, j2) ------------------------(i4, j2)                                      *
   *        +    -----------------------------------------------(i6, j6)               *
   *        +    -----------------------------------------------    +                  *
   *        +    -----------------------------------------------    +                  *
   *        +    -----------------------------------------------    +                  *
   *    (i1, j3) -----------------------------------------------    +                  *
   *             +    ------------------------------------------    +                  *
   *             +    ------------------------------------------    +                  *
   *             +    -------------------------------------------------- (i2, j4)      *
   *             +    ---------------------------------------------------    +         *
   *             +    ---------------------------------------------------    +         *
   *             +    ---------------------------------------------------    +         *
   *             +    ---------------------------------------------------    +         *
   *          (i5, j5) --------------------------------------------------    +         *
   *                               (i3, j1) ----------------------------- (i2, j1)     *
   *                                                                                   *
   *-------------------------------------------------------------- --------------------*/

  i2 = pis->s2ij[s2].i; j2 = pis->s2ij[s2].j;     // Set i2 & j2: the index in matrix for
                                                  // the end point of s2.
  idpe = n1 - 1;
  nwp = (pis->avenp * (n1-1) +0.5+ 1.e-5);      // The number of points with ID of 1 to (n1-1).
  nwpe=(pis->avenp * idpe +0.5+ 1.e-5)   ;      // The number of points with ID of 1 to idpe.
  for(j=0;j<pis->jm;j++){
    ib = i1; if(j<j1)ib = i1 + 1;          // Start i for current j.
    ie = i2; if(j>j2)ie = i2 - 1;          // End i for current j.
                                           // They are determined by the PE shape.
    int in=j*pis->im;
    for(i=ib;i<=ie;i++){
      if(pis->mask[i+in]){
        nwp = nwp + 1;                     // New part Begin
        if(nwp > nwpe){
          idpe = idpe + 1;                 // The new ID of PE.
          nwpe = (pis->avenp*idpe +0.5+ 1.e-5); // The number of points with ID of 1 to idpe.
          if(idpe >= n2){                  // For the PE of n2, the nunber of points
                                           // with ID of 1 to idpe should be s2.
            idpe = n2; nwpe = s2;
          }
        }
        pis->pemask[i+in] = idpe;         // Set the point (i, j) belong to PE with id of idpe.
      }
    }
  }
}
/*----------------------------------------------------------------------------------------------
 * sub. set_pemask_one_cycle: do the partition using a given delta.
 *   set_pemask_one_cycle(delta, ratio, dratio, flag)
 *   * double  delta  = The parameter used for computing aspect ratio of each PE.
 *   # double  ratio  = Maximum aspect ratio without thinking of delta.
 *   # double  dratio = Maximum aspect ratio with thinking of delta.
 *   * int  flag   = The flag for update PEMASK. 0: not update PEMASK; else: update PEMASK.
 *---------------------------------------------------------------------------------------------*/
static void compute_ratio(Irrp_Split*pis,int s1,int s2,int n1,int n2,double rt){
  int  s, i, j, sum_niflag, sum_njflag;
  double  sqmp, sum_ij, sum_ji;
  memset(pis->niflag,0,sizeof(int)*pis->im);// To record the number of none-empty-columns.
  memset(pis-> iflag,0,sizeof(int)*pis->im);// To record weights of each point-columns.
  memset(pis->njflag,0,sizeof(int)*pis->jm);// To record the number of none-empty-lines.
  memset(pis-> jflag,0,sizeof(int)*pis->jm);// To record weights of point-lines.
  for(s = s1;s<= s2;s++){
    i = pis->s2ij[s].i; pis->niflag[i] = 1; pis->iflag[i] ++;
    j = pis->s2ij[s].j; pis->njflag[j] = 1; pis->jflag[j] ++;
  }
  sum_niflag=0; sum_njflag=0;    
  sum_ij = 0; sum_ji = 0;
  for(int i=0;i<pis->im;i++){sum_niflag+=pis->niflag[i];} // Total number of none-empty-columns.
  for(int i=0;i<pis->jm;i++){sum_njflag+=pis->njflag[i];} // Total number of none-empty-lines.
  sum_niflag=1./sum_niflag; 
  sum_njflag=1./sum_njflag;    

  // Here, a power of 0.3 is used to emphasizing the effective of land points. 
  // If use 1, it will be no emphasizing.
  for(int i=0;i<pis->im;i++)sum_ij+= pow(pis->iflag[i]*sum_njflag,0.3);  // Weighted number of point-columns.
  for(int i=0;i<pis->jm;i++)sum_ji+= pow(pis->jflag[i]*sum_niflag,0.3);  // Weighted number of point-lines.
  rt = sum_ij *(double)(n2 - n1 + 1)/sum_ji;// aspect ratio of this PE column.
}
static void set_pemask_one_cycle(Irrp_Split*pis,double delta, double *ratio_, double *dratio_, int flag){
  double   ratio, dratio;
  int  n1, n2, s1, s2, m, i, minline, ml, mml, onemore, mn;
  double  r1, r1t, rt, r, rt0, rt1;
  n2 = 0;                                 // initializing n2: end PE of last PE column.
  pis->ncols = 0;                              // initializing number of columns.
  ratio = 0.;                             // initializing global maximum ratio without delta.
  dratio = 0.;                            // initializing global maximum ratio with delta.
  mml = pis->jm * 0.7;                         // estimated minimum PEs for a point-column.
  minline = sqrt(pis->npe * 1.0 * pis->jm / pis->im) / 3;// One third of the expected minimum of number of PEs
                                                   // for one PE-column.
  mn = 0.8 * pis->jm / pis->savenp;            // 0.8 times of the number of PEs for one ponit-column.
  for(;n2 != pis->npe;){                   // Find s2, n2
    n1 = n2 + 1; s1 = pis->avenp * n2 + 1.e-5+0.5 + 1;
    pis->ncols = pis->ncols + 1; r1 = 1000000.;
    m = n1 + mn - 1; if(m >= pis->npe) m = pis->npe - 1;
    /*-----------------------------------------------------------------------------------------
     * Find a best n2 which will lead to a smaller aspect ratio with a given delta.
     * The maximum ratio and dratio will be recorded during this process.
     *-----------------------------------------------------------------------------------------*/
    for(;m < pis->npe;){
      m = m + 1;                          // add one PE for each search continue.
      s2 = pis->avenp * m + 1.e-5+0.5;       // end point of current PE with id of m.
                                           // (pis->s2ij[s2].i - pis->s2ij[s1].i): number of point-columns
                                           // (m - n1 + 1): number of PEs for this PE-column.
      ml = (pis->s2ij[s2].i - pis->s2ij[s1].i) * (m - n1 + 1);
      // ml / mml : gross aspect ratio of PEs from n1 to m.
      if(ml < mml)continue;                  // Speed up by adding more PEs.
                                             // If the rest points of this point-colomn is greater
                                             // than the points of one PE, add more PEs for these
                                             // PE-column and the computing of aspect ratio is not
                                             // needed.
      if(pis->avenp < pis->jm/2 && s2 < pis->sumnp){
        onemore = pis->avenp * (m + 1) + 1.e-5+0,5;
        // End point when add one more PE.
        if(pis->s2ij[onemore].i == pis->s2ij[s2].i)continue;
        // Further speed up by adding more PEs.
        // After one more PE is added, if it still in
        // the same point-column,  do not compute
        // the aspect ratio.
      }
      compute_ratio(pis,s1, s2, n1, m, rt);
      // Compute the aspect ratio for PEs from n1 to m.
      // rt: the pure aspect ratio around to 1.
      rt0 = 1 - rt;                       // The aspect ratio compared with 1.
      rt1 = 1 + delta - rt;               // The aspect ratio compared with (1 + delta).
      if(abs(rt1) < abs(r1)){
        r1 = rt1; r1t = rt0; n2 = m;      // If this ratio is accetable, set n2 eqqual to m,
                                          // and record the aspect ratios.
      }else if(abs(rt1) > 1.){
        break;                             // If this ratio is too big, finish this searching.
      }
    }
    /*------------------------------------------------------------------------------------------
     * Patch the last PE column.
     * n1 > n2 means last n2 is not fund.
     * n2 + minline >= npe means if use n2, the rest PEs after n2 is too less.
     * n2 != npe: the last PE-colum is incorrect.
     *------------------------------------------------------------------------------------------
     */
    if((n1 > n2 || n2 + minline >= pis->npe) && n2 != pis->npe){
      n2 = pis->npe; s2 = pis->sumnp;
      compute_ratio(pis,s1, s2, n1, n2,  rt);
      // Compute the aspect ratio for PEs from n1 to n2 (npe).
      rt0 = 1 - rt; rt1 = 1 + delta - rt;
      if(abs(rt1) < abs(r1)){
        r1 = rt1; r1t = rt0;              // Record the aspect ratio for last PE-column.
      }
    }
    if(abs(ratio)  < abs(r1 ))ratio  = r1;// Record the ratio and dratio for this column of PEs.
    if(abs(dratio) < abs(r1t))dratio = r1t;
    if(flag!=0)set_pemask_columns(pis,n1, n2); // Set pemask for this PE column.
  }
  if(ratio_)*ratio_=ratio; 
  if(dratio_)*dratio_=dratio;
}
/*----------------------------------------------------------------------------------------------
 * sub. set_pemask_iteration : Do partition using iteration based on the gridient of the global
 *                             maximum aspect ratio of each PE respect to parameter of delta.
 *
 *   set_pemask_iteration(flag[, delta_, ratio_])
 *
 *   @ int  flag   = The flag for update PEMASK. 0: not update PEMASK; else: update PEMASK.
 *   * double  delta_ = Parameter used to compute the global maximum aspect ratio of each PE.
 *   * double  ratio_ = Aspect ratio for a given delta.
 *---------------------------------------------------------------------------------------------*/
static void set_pemask_iteration(Irrp_Split*pis,int flag, double delta_, double *ratio_){
  double  delta;
  int  ncols0, i0, i1, n, i, j;
  double  ratio, dratio;
  double  grad, lamb, grad0;
  double  minr1, rm0, mindelta, delta1, dstep, lastr1, r1all0, dt0, dt1;
  double  mmr1, minlamb, small;

  NEWPN(pis->iflag ,pis->im); NEWPN(pis->niflag,pis->im);
  NEWPN(pis->jflag ,pis->jm); NEWPN(pis->njflag,pis->jm);
  delta = 0.;                             // set delta as zero for first guess.
  delta = delta_;      // Use the given delta.

  //DBGB0,"setpemaskI 0",delta,pis->avenp
  set_pemask_one_cycle(pis,delta, &ratio, &dratio, flag);
  //DBGB0,"setpemaskI 1",delta, ratio
                                          // Using the given value of delta to obtain the first
                                          // guess of the maximum aspect ratio.
  if(flag != 0){
    if(ratio_)*ratio_ = ratio;            // Can be used to do parallized partition.
    return;
  }
  rm0 = 1;  lamb = 0.;                    // Set the initial value of minr1, rm0 and lamb.
  minr1 = ratio; mindelta = delta;        // Record the minimum aspect ratio during the iteration.
  ncols0 = pis->ncols;                      // Record the number of columns at the first guess.
  dstep   = 5. / ncols0;                  // Set the step size for the deepest method using the
  small   = dstep * 0.0001;               // number of PE-columns.
  minlamb = dstep * 0.001;                // Limit for minimum lamb.
  mmr1    = 0.1 * ncols0 / pis->im;         // Limit for minimum aspect ratio.
  i0      = 0;                            // To record the iteration cycle.
  grad0   = 0;                            // Set the first gridient to zero for the first
                                          // iteration cycle.
  i1 = 0;                                 // To record the times of calling set_pemask_one_cycle
  dt0 = -0.5; dt1 = 0.5;                  // Bound limit for delta.
  // Deepest Method:
  for(;;){                                // find best delta.
    lastr1 = ratio;                       // Record the last ratio.
    set_pemask_one_cycle(pis,delta+small, &ratio, &dratio, 0);
    //DBGB0,"setpemaskI 2",delta+small, ratio
                                          // Using a small value added to delta to obtain
                                          // the second guess of the maximum aspect ratio.
    i1 = i1 + 1;                          // Record the times calling the set_pemask_one_cycle
                                          // during the iteration.
    grad = (ratio*ratio - lastr1*lastr1) / small;
                                          // Compute the gradient using the first and second
                                          // guess of the maximum aspect ratio.
    if(pis->ncols != ncols0){               // Once the number of PE collums is different with
                                          // former iteration cycle, the small value added
                                          // to delta is set to be a negitive one.
      set_pemask_one_cycle(pis,delta-small, &ratio, &dratio, 0);
      //DBGB0,"setpemaskI 3",delta-small, ratio
      i1 = i1 + 1;
      if(pis->ncols != ncols0)break;            // If the number of PE columns is still different
                                          // with former iteration cycle, break the iteration.
      grad = (ratio*ratio - lastr1*lastr1) / small;
                                          // Compute the gradient using based on the negitive
                                          // small value added to delta.
    }
    if(abs(minr1) > abs(ratio)){
      minr1 = ratio; mindelta = delta;    // Record the optimal delta during the iteration.
    }
    if(abs(grad) < small)break;            // If the gradient is smaller than the added small value,
                                          // finish the iteration.
    if(grad0 * grad > 0){
      break;                               // If the current gradient and the former gradient are
                                          // in the same direction, finish this iteration.
    }else{
      lamb = grad * dstep;                // Set the step size for steepest iteration method.
    }

    grad0 = grad;                         // Record the current gradient as the former one in
                                          // next iteration cycle.
    r1all0 = 100; i1 = 0;
    /*------------------------------------------------------------------------------------------
    * Along the direction of the current gradient,
    * search for a delta which will cause a minimum aspect ratio.
    *-----------------------------------------------------------------------------------------*/
    for(;;){                        // find best lamb
      dratio = 0; r1all0 = ratio;
      if(lamb > 0){
        for(;delta - lamb < dt0;){
           lamb = lamb / 2;
        }
      }else{
        for(;delta - lamb > dt1;){
          lamb = lamb / 2;
        }
      }
      delta = delta - lamb;               // Attempt new delta for minimum aspect ratio.
      set_pemask_one_cycle(pis,delta, &ratio, &dratio, 0);
      //DBGB0,"setpemaskI 4",delta, ratio
      i1 = i1 + 1;
      for(;pis->ncols != ncols0;){          // Keep the number of PE columns same as former one.
        if(abs(lamb) < minlamb)break;
        delta = delta + lamb; lamb = lamb / 5; delta = delta - lamb;
        set_pemask_one_cycle(pis,delta, &ratio, &dratio, 0);
        //DBGB0,"setpemaskI 5",delta, ratio
        i1 = i1 + 1;
      }
      for (;r1all0*ratio<0;){           // If the aspect ratio has a different sign with former,
                                          // search in a different direction.
        if(abs(lamb) < minlamb)break;      // if the step size is smaller enough, finish this search.
        delta = delta + lamb; lamb = lamb / 3; delta = delta - lamb;
        set_pemask_one_cycle(pis,delta, &ratio, &dratio, 0);
        //DBGB0,"setpemaskI 6",delta, ratio
        i1 = i1 + 1;                      // Attempt a different step size.
      }
      if(abs(minr1) > abs(ratio)){
        minr1 = ratio; mindelta = delta;  // record the minimum aspect ratio and delta.
      }
      /*---------------------------------------------------------------------------------------
      * if the step size is smaller enough or the minimum aspect ratio is smaller enough,
      * finish this searching.
      *---------------------------------------------------------------------------------------*/
      if(abs(lamb)<minlamb || abs(minr1)<mmr1)break;

      if(abs(r1all0) < abs(ratio)){    // adjust for next iteration cycle: dt0, dt1 and dstep.
        if(lamb < 0){
          if(dt0 < delta + lamb)dt0 = delta + lamb;
          if(dt1 > delta)dt1 = delta;
        }else{
          if(dt1 > delta + lamb)dt1 = delta + lamb;
          if(dt0 < delta)dt0 = delta;
        }
        dstep = dstep / 2;
        break;
      }
    }                                 // find best lamb
    i0 = i0 + 1;
    /*------------------------------------------------------------------------------------------
    * if the step size is smaller enough or the minimum aspect ratio is smaller enough,
    * break this iteration.
    *------------------------------------------------------------------------------------------*/
    if(abs(lamb) < minlamb || abs(minr1) < mmr1)break;
  }                                   // find best delta.
  delta = mindelta;                       // Using the optimal delta to update the PEMASK.
  //DBGB0,"setpemaskI 7",delta
  set_pemask_one_cycle(pis,delta, &ratio, &dratio, 1);
  //DBGB0,"setpemaskI 8",delta, ratio
  VFREE(pis->iflag ); VFREE(pis->jflag ); VFREE(pis->niflag); VFREE(pis->njflag);
}
