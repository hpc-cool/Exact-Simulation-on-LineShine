
__BEGIN_DECLS
double dclktime();
void  Resettimer(int iid);
void endstarttimer(int iide,int iids);
void starttimer(int iid);
void endtimer(int iid);
double difftimer(int iid);
double gettimer(int iid);
void  resettimer_(int *iid);
void endstarttimer_(int *iide,int *iids);
void starttimer_(int *iid);
void endtimer_(int *iid);
double difftimer_(int *iid);
double gettimer_(int *iid);
void checkenv_();
void msleep(int ms);
void msleep_(int *ms);
void setmask1_(int *p,int *nn);
void setgrid_(int *p,int *nn,float*pvb,float *pv0,float *pv1);
void countrect_(int*p,int *pm,int *pim,int*pjm,int *prect,int*pnc);
int IntVolatile(void*p);
void print_call_stack(int lock) ;
void register_segfault_handler_() ;
int is_ptr_valid(const void *ptr, size_t size) ;
int is_ptr_valid_(const void *ptr, int *size) ;
void zunlock(int *p);
void zlock(int *p);
void zunlock_();
void zlock_();
//#define MID md.si.selfid
int ccheckee(double *ee,int n);
void ccheckee_(double *ee,int *n,int *err);
__END_DECLS
