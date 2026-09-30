#include "svars.h"
#include "ctools.h"
#include <stdlib.h>
#define USE_SMOOTH8
#undef ARM_VEC
#ifdef ARM_VEC
typedef struct vspc_interg{
  float pp0[8],pp1[8],pp2[8],pp3[8];//32*4     128
  short  kpos[4];                    // 2*4=8   128+64+8
  int flag;                                     
}vspc_interg;
typedef struct vgeo_interg{
#ifdef VECN
  float q[8],r[8];                  //32*2     128+64
#else
  float ps0[8],ps1[8],ps2[8],ps3[8];//32*2     256
  int    ipos[3];                    // 4*3=12  128+64+24
#endif
  short  iquad,flag;                 // 2*2=4   128+64+12
  //int    offset[8];                // 8*4=32
}vgeo_interg;
#endif
void vinitPropgats();
#if 0
extern "C" void perr();
__thread struct {
  int ind;
  int iac,iacb,iace;
  void *pb,*pe;
  void *pi,*pc;
}teinf={0};
void perr(){
  char*pb,*pe;
  long size=mkj/8L*(gd->nwpc+1);
  int nl=0;
  if(teinf.ind==1){
    nl=sizeof(vspc_interg );
    pb=(char*)gd->psvis;
  } else {
    nl=sizeof(vgeo_interg);
    pb=(char*)gd->pgvis;
  }
  pe=pb+size*nl;
  printf("EI:%d %d %p %p %p %p:%p %p %d %d %d \n",
         teinf.ind,teinf.iac,
         teinf.pi,(char*)teinf.pi+nl,pb,pe,
         gd->pb,gd->pe,
         teinf.iacb,teinf.iace,gd->nwpc+1);
}
#endif
void InitCPropgats(){
#ifndef MMTHREAD
#ifdef ARM_VEC
  if(gd->psvis)return;
  int iac,iacb=1,iace=gd->nwpc;
  long size=mkj/8L*(gd->nwpc+1);
  long sa=256+size*(sizeof(vspc_interg)+ sizeof(vgeo_interg));
  char *p;
  p=HNEWN(char ,sa);
  gd->psvis=(vspc_interg*)p;p+=(sizeof(vspc_interg)*size+63)&(-64);
  gd->pgvis=(vgeo_interg*)p;
  for(iac=iacb;iac<=iace;iac++){
    size_t it=(iac*mkj);
    int kj,*ii,iquad=-1;
    int *l_pos=pd.ipos12+iac*12;
    spc_interg *psi=pd.psis+it;  // (mkj,0:nwpc)
    geo_interg *pgi=pd.pgis+it;  // (mkj,0:nwpc)
    vspc_interg*psv=gd->psvis+it/8;
    vgeo_interg*pgv=gd->pgvis+it/8;
    if(1){
      for(kj=0;kj<mkj/8;kj++,psi+=8,pgi+=8){
        int i,k;
        int idf=0,igf=0;
        if(0){
          iquad=pgi->iquad;
          int *ii=l_pos+pgi->iquad*3;
          int i,j,it[4],ie[4];
          ie[0]=iac;
          for(i=0;i<4;i++){
            it[i]=psi[0].kpos[i];
            int *ii=l_pos+pgi[0].iquad*3;
            if(i)ie[i]=ii[i-1];
            for(j=1;j<8;j++){
              int *ii=l_pos+pgi[j].iquad*3;
              int ig=0, id=psi[j].kpos[i]-it[i]-j;
              if(i) ig=ii[i-1]-ie[i];
              if(id){
                if(psi[j].pab[i]<1e-5)id=0;
              }
              idf|=id;igf|=ig;
            }
          }
          pgv[kj].flag= (idf||igf);
          psv[kj].flag= (idf||igf);
        }else {
          pgv[kj].flag= 0;
          psv[kj].flag= 0;
        }
        pgv[kj].iquad  =pgi[0].iquad;
        for(i=0;i<8;i++){
          psv[kj].pp0[i]=psi[i].pab [0];
          psv[kj].pp1[i]=psi[i].pab [1];
          psv[kj].pp2[i]=psi[i].pab [2];
          psv[kj].pp3[i]=psi[i].pab [3];
          psv[kj].kpos[0]=psi[i].kpos[0];
          psv[kj].kpos[1]=psi[i].kpos[1];
          psv[kj].kpos[2]=psi[i].kpos[2];
          psv[kj].kpos[3]=psi[i].kpos[3];
          pgv[kj].ps0[i]=pgi[i].pabp[0];
          pgv[kj].ps1[i]=pgi[i].pabp[1];
          pgv[kj].ps2[i]=pgi[i].pabp[2];
          pgv[kj].ps3[i]=pgi[i].pabp[3];
          if(idf||igf){//fixme
             //pgv[kj].offset[]  =0;
          }
        }
      }
    }
  }
#endif
#endif
}
void initPropgats();
void g_initPropgats(){
#ifdef ARM_VEC
  vinitPropgats();//fixme mmthread && no arm_vec,
#else
#ifdef MMTHREAD
  initPropgats();//fixme mmthread && no arm_vec,
#endif
#endif
}
#ifndef ARM_VEC
#ifdef MMTHREAD
void initPropgats(){
  if(gd->psis)return;
  int iacb=1,iace=gd->nwpc;
  long size=mkj*(gd->nwpc+1);
  long sa=64*5+size*(sizeof(spc_interg)+ sizeof(geo_interg)+(12+8)*sizeof(int));
  char *p=HNEWN(char ,sa);
  gd->ipos12=(int*)p;p+=12*sizeof(int);
  gd->ipos8 =(propinf*)p;p+=8*sizeof(int);
  gd->psis=(spc_interg*)p;p+=(sizeof(spc_interg)*size+63)&(-64);
  gd->pgis=(geo_interg*)p;
  for(size_t iac=iacb;iac<=iace;iac++){
    size_t giac=gd->l2gind[iac];
    spc_interg *psi=pd.psis+giac*mkj;  // (mkj,0:nwpc)
    geo_interg *pgi=pd.pgis+giac*mkj;  // (mkj,0:nwpc)
    int *gipos12=pd.ipos12+giac*12;
    int *ipos12=gd->ipos12+iac*12;
    int*gipos8=pd .ipos8[giac].i8;
    int* ipos8=gd->ipos8[ iac].i8;
    for(int n=0;n<8;n++){
      int mt=ipos8[n]=gd->gg2l[gipos8[n]];// assume part is ok
      if(mt<0||mt>gd->nwpa){
        ipos12[n]=0;//fixme
        printf("EEEEEEEEBA %d %d %5d %5d %d %d\n",mpi_id,gi->igrp,
               mt,gd->nwpa,n,gipos12[n]);
      }
    }
    for(int n=0;n<12;n++){
      int mt=ipos12[n]=gd->gg2l[gipos12[n]];// assume part is ok
      if(mt<0||mt>gd->nwpa){
        ipos12[n]=0;//fixme
        printf("EEEEEEEEBA %d %d %5d %5d %d %d\n",mpi_id,gi->igrp,
               mt,gd->nwpa,n,gipos12[n]);
      }
    }
    spc_interg*ps=gd->psis+iac*mkj;
    geo_interg*pg=gd->pgis+iac*mkj;
    for(int im=0;im<12;im++){
      int mt=ipos12[im];
      if(mt<0||mt>gd->nwpa){
        printf("EEEEEEEEAA %d %d %5d %5d\n",mpi_id,gi->igrp,mt,gd->nwpa);
        exit(0);
      }
    }
    for(int kj=0;kj<mkj;kj++,psi++,pgi++){
      int idf=0,igf=0;
      int iquad=pgi->iquad;
      pg[kj].iquad  =iquad;
      ps[kj].kpos[0]=psi->kpos[0];
      ps[kj].kpos[1]=psi->kpos[1];
      ps[kj].kpos[2]=psi->kpos[2];
      ps[kj].kpos[3]=psi->kpos[3];
      ps[kj].pab [0]=psi->pab [0];
      ps[kj].pab [1]=psi->pab [1];
      ps[kj].pab [2]=psi->pab [2];
      ps[kj].pab [3]=psi->pab [3];
      pg[kj].pabp[0]=pgi->pabp[0];
      pg[kj].pabp[1]=pgi->pabp[1];
      pg[kj].pabp[2]=pgi->pabp[2];
      pg[kj].pabp[3]=pgi->pabp[3];
    }
  }
}
#endif
#endif
// Main Core Code
extern "C" void offsetinit();
//static double*ppab,*ppabp,*ppabt=NULL;
#ifdef ARM_VEC
/*
 *svfloat64_t{p0,p1,p2,p3} ...
 * not surport load from float32_t to svfloat64_t ,
 * speed for cache
 * */
void vinitPropgats(){
  if(gd->psvis)return;
  int iacb=1,iace=gd->nwpc;
  long size=mkj/8L*(gd->nwpc+1);
  long sa=256+size*(sizeof(vspc_interg)+ sizeof(vgeo_interg));
  char *p=HNEWN(char ,sa);
  //gd->pb=p; gd->pe=p+sa;
  gd->psvis=(vspc_interg*)p;p+=(sizeof(vspc_interg)*size+63)&(-64);
  gd->pgvis=(vgeo_interg*)p;
  for(int iac=iacb;iac<=iace;iac++){
#ifdef MMTHREAD
    int giac=gd->l2gind[iac];
    spc_interg *psi=pd.psis+giac*mkj;  // (mkj,0:nwpc)
    geo_interg *pgi=pd.pgis+giac*mkj;  // (mkj,0:nwpc)
    int *gl_pos=pd.ipos12+giac*12;
    int l_pos[12]={0};
    for(int n=0;n<12;n++){
      int mt=l_pos[n]=gd->gg2l[gl_pos[n]];// assume part is ok
      if(mt<0||mt>gd->nwpa){
        l_pos[n]=0;//fixme
        printf("EEEEEEEEBA %d %d %5d %5d %d %d\n",mpi_id,gi->igrp,
               mt,gd->nwpa,n,gl_pos[n]);
      }
    }
#else
    spc_interg *psi=pd.psis+iac*mkj;  // (mkj,0:nwpc)
    geo_interg *pgi=pd.pgis+iac*mkj;  // (mkj,0:nwpc)
    int *l_pos=pd.ipos12+iac*12;
#endif
    vspc_interg*psv=gd->psvis+iac*mkj/8;
    vgeo_interg*pgv=gd->pgvis+iac*mkj/8;
    for(int im=0;im<12;im++){
      //int mt=pgv[kj].ipos[im];
      int mt=l_pos[im];
      if(mt<0||mt>gd->nwpa){
#ifdef MMTHREAD
        printf("EEEEEEEEAA %d %d %5d %5d\n",mpi_id,gi->igrp,mt,gd->nwpa);
#else
        printf("EEEEEEEEAB %d %d %5d %5d\n",mpi_id,0,mt,gd->nwpa);
#endif
        exit(0);
      }
    }
    for(int kj=0;kj<mkj/8;kj++,psi+=8,pgi+=8){
      int idf=0,igf=0;
      int iquad=pgi->iquad;
      if(1){
        int *ii=l_pos+pgi->iquad*3;
        int it[4],ie[4];
        ie[0]=iac;
        for(int i=0;i<4;i++){
          it[i]=psi[0].kpos[i];
          if(i)ie[i]=ii[i-1];
          for(int j=1;j<8;j++){
            int *ii=l_pos+pgi[j].iquad*3;
            int ig=0, id=psi[j].kpos[i]-it[i]-j;
            if(i) ig=ii[i-1]-ie[i];
            if(id){
              if(psi[j].pab[i]<1e-5)id=0;
            }
            idf|=id;igf|=ig;
          }
        }
        pgv[kj].flag= (idf||igf);
      }else pgv[kj].flag= 0;
      pgv[kj].iquad  =iquad;
      pgv[kj].ipos[0]=l_pos[iquad*3+0];
      pgv[kj].ipos[1]=l_pos[iquad*3+1];
      pgv[kj].ipos[2]=l_pos[iquad*3+2];
     {
        for(int im=0;im<3;im++){
          //int mt=pgv[kj].ipos[im];
          int mt=l_pos[iquad*3+im];
          if(mt<0||mt>gd->nwpa){
#ifdef MMTHREAD
            printf("EEEEEEEEEA %d %d %5d %5d\n",mpi_id,gi->igrp,mt,gd->nwpa);
#else
            printf("EEEEEEEEEB %d %d %5d %5d\n",mpi_id,0,mt,gd->nwpa);
#endif
            exit(0);
          }
        }
      }
      psv[kj].kpos[0]=psi[0].kpos[0];
      psv[kj].kpos[1]=psi[0].kpos[1];
      psv[kj].kpos[2]=psi[0].kpos[2];
      psv[kj].kpos[3]=psi[0].kpos[3];
      for(int i=0;i<8;i++){
        psv[kj].pp0[i]=psi[i].pab [0];
        psv[kj].pp1[i]=psi[i].pab [1];
        psv[kj].pp2[i]=psi[i].pab [2];
        psv[kj].pp3[i]=psi[i].pab [3];
#ifdef VECN
        pgv[kj].r[i]=pgi[i].pabp[3]+pgi[i].pabp[2];
        pgv[kj].q[i]=pgi[i].pabp[3]+pgi[i].pabp[1];
#else
        pgv[kj].ps0[i]=pgi[i].pabp[0];
        pgv[kj].ps1[i]=pgi[i].pabp[1];
        pgv[kj].ps2[i]=pgi[i].pabp[2];
        pgv[kj].ps3[i]=pgi[i].pabp[3];
#endif
        if(idf||igf){//fixme
          //pgv[kj].offset[]  =0;
        }
      }
    }
  }
}
/*
 *svfloat64_t{p0,p1,p2,p3} ...
 * not surport load from float32_t to svfloat64_t ,
 * speed for cache
 * */
//CALA 36*NK*NJ
void  propagats_spec1(double *ed,double *es,int iac){
  es--; // kpos is index from 1
  int it=(iac*mkj), kj;
  vspc_interg *pi=gd->psvis+it/8;  // (mkj,0:nwpc)
  double *pe0=es+it;
  double *pe=ed+it;
  svbool_t pm=svptrue_b64();
  for(kj=0;kj<mkj;kj+=8,pi++){
    //teinf.pi=pi;
    if(!pi->flag){
      short *kpos=&pi->kpos[0];
#define LDO(l) VLD(pe0+kpos[l])
      VDOUBLE  p0,p1,p2,p3;
      p0=VLFD(pi->pp0 ); 
      p1=VLFD(pi->pp1 ); 
      p2=VLFD(pi->pp2 ); 
      p3=VLFD(pi->pp3 );
      VVD(pe+kj) =MLA(MLA(MLA(MUL(p0,LDO(0)),p1,LDO(1)),p2,LDO(2)),p3,LDO(3));
#undef LDO
    }else{
      // fixme...
    }
  }
}
void  propagats_geo(double *ed,double *es,int iacb,int iace){
  svbool_t pm=svptrue_b64();
  //teinf.ind=2;
  for(int iac=iacb;iac<=iace;iac++){
    if(gd->nsp[iac]!=1)continue;
    //teinf.iac=iac;
    int it=(iac*mkj);
    int iquad=-1;
    vgeo_interg *pi=gd->pgvis+it/8;  // (mkj,0:nwpc)
    double *pe,*pe0,*pe1=NULL,*pe2=NULL,*pe3=NULL;
    pe0=es+it;
    pe=ed+it;
    for(int kj=0;kj<mkj;kj+=8,pi++){
      //teinf.pi=pi;
      if(!pi->flag){
        if(iquad!=pi->iquad){
          iquad=pi->iquad;
          int *ii=pi->ipos;
          pe1=es+ii[0]*mkj;
          pe2=es+ii[1]*mkj;
          pe3=es+ii[2]*mkj;
        }
#ifdef VECN
        VDOUBLE ones=VDUP(1.0);
        VDOUBLE vr,vq,osr,osq,p0,p1,p2,p3;
        vr=VLFD(pi->r);vq=VLFD(pi->q);
        osr=SUB(ones,vr);osq=SUB(ones,vq);
        p0=MUL(osq,osr); p1=MUL(vq,osr);
        p2=MUL(osq, vr); p3=MUL(vq, vr);
        VVD(pe+kj) =MLA(MLA(MLA(MUL(
               p0,VLD(pe0+kj)),p1,VLD(pe1+kj)),
               p2,VLD(pe2+kj)),p3,VLD(pe3+kj));
#else
        VVD(pe+kj) =MLA(MLA(MLA(MUL(
            VLFD(pi->ps0 ),VLD(pe0+kj)), VLFD(pi->ps1 ),VLD(pe1+kj)),
            VLFD(pi->ps2 ),VLD(pe2+kj)), VLFD(pi->ps3 ),VLD(pe3+kj));
#endif
      }
    }
  }
}
#else
//CALA 36*NK*NJ
void vinitPropgats(){}
extern "C" void pppsis_(int *id,spc_interg *pp);
void pppsis_(int *id,spc_interg *pp){
  spc_interg *pi;  // (mkj,0:nwpc)
  char fn[256];                   
#ifdef MTHREAD
  sprintf(fn,"m/pis%2.2d.txt",*id);
#else
  sprintf(fn,"v/pis%2.2d.txt",*id);
#endif
  FILE*fo=fopen(fn,"wt");
  pi=pd.psis+mkj;  // (mkj,0:nwpc)
  for(int kj=0;kj<mkj;kj++,pi++){
    short *kpos=pi->kpos;
    float *pab=pi->pab;
    fprintf(fo,"%d %8.5f %8.5f %8.5f %8.5f :%5d %5d %5d %5d\n",kj
            ,pab[0],pab[1],pab[2],pab[3]
            ,kpos[0],kpos[0],kpos[0],kpos[0]);

  }
  fclose(fo);
}
void  propagats_spec1(double *ed,double *es,int iac){
  es--; // kpos is index from 1
  int it=(iac*mkj), kj;
  spc_interg *pi=gd->psis+it;  // (mkj,0:nwpc)
  double *pe0=es+it;
  double *pe=ed+it;
  for(int kj=0;kj<mkj;kj++,pi++){
    short *kpos=pi->kpos;
    float *pab=pi->pab;
    double vvt;
#define PE(L) pab[L]*pe0[kpos[L]]
    pe[kj]=PE(0)+PE(1)+PE(2)+PE(3);
#undef PE
  }
}
void  propagats_geo(double *ed,double *es,int iacb,int iace){
  for(int iac=iacb;iac<=iace;iac++){
    if(gd->nsp[iac]!=1)continue;
    int it=(iac*mkj);
    int iquad=-1;
    geo_interg *pi=gd->pgis+it;  // (mkj,0:nwpc)
    int *l_pos=gd->ipos12+iac*12;
    double *pe0,*pe1,*pe2,*pe3,*pe;
    pe0=es+it;
    pe=ed+it;
    for(int kj=0;kj<mkj;kj+=8,pi++){
      if(iquad!=pi->iquad){
        int iqt=(iquad=pi->iquad)*3;
        int *ii=l_pos+iqt;
        pe1=es+ii[0]*mkj;
        pe2=es+ii[1]*mkj;
        pe3=es+ii[2]*mkj;
      }
     {
        float *pabp=pi->pabp;
        double et;
#define PE(L) pabp[L]*pe##L[kj]
        pe[kj]=PE(0)+PE(1)+PE(2)+PE(3);
#undef PE
      }
    }
  }
}
#endif
void  propagats_spec(double *ed,double *es,int iacb,int iace){
  for(int iac=iacb;iac<=iace;iac++){
    if(gd->nsp[iac]!=1)continue;
    propagats_spec1(ed,es,iac) ;
  }
}
void propagats_spec_(double *ed,double *es,int *iacb,int *iace){
  propagats_spec(ed ,es,*iacb,*iace) ;
}
void propagats_geo_(double *ed,double *es,int *iacb,int *iace){
  propagats_geo (ed,es ,*iacb,*iace) ;
}
#ifdef USEARM_SME
void init_SSmooth(){
  SSmooth *smt=&td->smt;
  double val= 0.1;
  double *stc=smt->stc,*stv=smt->stcv;
  stc[0]=val;stc[1]=val    ;stc[2]=val;
  stc[3]=val;stc[4]=1-val*8;stc[5]=val;
  stc[6]=val;stc[7]=val    ;stc[8]=val;
  stv[0]=stc[0];stv[1]=stc[1];stv[2]=stc[2];
  stv[3]=stc[3];             ;stv[4]=stc[5];
  stv[5]=stc[6];stv[6]=stc[7];stv[7]=stc[8];
  smt->clcy=0;smt->inited_neb=false;
  int *ind=smt->ind;
  ind[0*3+0]=0;ind[0*3+1]=1;ind[0*3+2]=2;
  ind[1*3+0]=1;ind[1*3+1]=2;ind[1*3+2]=0;
  ind[2*3+0]=2;ind[2*3+1]=0;ind[2*3+2]=1;
  smt->inited_upem= false;
}
void evcpy(double*pd,double*ps){
  const int VL=svcntd();
  VDOUBLE v0,v1,v2,v3;
  //assume gt 4
  v0=VLD(ps+0*VL);v1=VLD(ps+1*VL);
  v2=VLD(ps+2*VL);v3=VLD(ps+3*VL);ps+=4*VL;
  int j=mkj/VL-7;
  for(;j;j-=4,ps+=4*VL,pd+=4*VL){
    VST(pd+0*VL,v0);v0=VLD(ps+0*VL);
    VST(pd+1*VL,v1);v1=VLD(ps+1*VL);
    VST(pd+2*VL,v2);v2=VLD(ps+2*VL);
    VST(pd+3*VL,v3);v3=VLD(ps+3*VL);
  }
  VST(pd+0*VL,v0);VST(pd+1*VL,v1);
  VST(pd+2*VL,v2);VST(pd+3*VL,v3);pd+=4*VL;
  if(j+=3){
    v0=VLD(ps);ps+=VL;
    for(;j;j--){
      VST(pd,v0);pd+=VL;v0=VLD(ps);ps+=VL;
    }
    VST(pd,v0);pd+=VL;
  }
}
void upem(double *em, double *et, int *neb, bool reverse, bool up_neb_all){
  const int VL=svcntd();
  int clcy=td->smt.clcy;
  int *old_neb=td->smt.old_neb;
  int *ind=td->smt.ind;
  svbool_t pm =VTRUE;
  if(reverse){
    for(int n=0;n<9;n++){
      evcpy(et+old_neb[n],em+n*mkj);
    }
  } else{
    if(td->smt.inited_upem){
      if(up_neb_all){
        for(int i=0;i<3;i++){
          evcpy(et+old_neb[ind[clcy*3+i]+0],em+(ind[clcy*3+i]+0)*mkj);
          evcpy(et+old_neb[ind[clcy*3+i]+3],em+(ind[clcy*3+i]+3)*mkj);
          evcpy(et+old_neb[ind[clcy*3+i]+6],em+(ind[clcy*3+i]+6)*mkj);
        }
      }else{
        evcpy(et+old_neb[ind[clcy*3+2]+0],em+(ind[clcy*3+2]+0)*mkj);
        evcpy(et+old_neb[ind[clcy*3+2]+3],em+(ind[clcy*3+2]+3)*mkj);
        evcpy(et+old_neb[ind[clcy*3+2]+6],em+(ind[clcy*3+2]+6)*mkj);
      }
    }
    if(up_neb_all){
      for(int i=0;i<3;i++){
        evcpy(em+(i+0)*mkj,et+old_neb[i+0]);
        evcpy(em+(i+3)*mkj,et+old_neb[i+3]);
        evcpy(em+(i+6)*mkj,et+old_neb[i+6]);
      }
      svbool_t pg=svwhilelt_b32(0,9);
      VINT vac=svld1_s32(pg, neb);
      svst1_s32(pg,old_neb, vac);
    }else{
      evcpy(em+(ind[clcy*3+2]+0)*mkj,et+neb[2]);
      evcpy(em+(ind[clcy*3+2]+3)*mkj,et+neb[5]);
      evcpy(em+(ind[clcy*3+2]+6)*mkj,et+neb[8]);
      old_neb[ind[clcy*3+2]+0]=neb[2];
      old_neb[ind[clcy*3+2]+3]=neb[5];
      old_neb[ind[clcy*3+2]+6]=neb[8];
    }
    td->smt.inited_upem=true;
  }
}
bool getneb(int *neb_pos, int pos){
  int*ind=td->smt.ind;
  int*ipos8=((gd->ipos8)+pos)->i8;
  int clcy=td->smt.clcy;
  int *all=td->smt.all;
  bool up_neb_all;//是否更新全部的邻居节点
  /* 4---3---2
   * |   |   |
   * 5---+---1
   * |   |   |
   * 6---7---8  */
  if(td->smt.inited_neb && all[5]== pos *mkj){
    clcy=(clcy+1)%3;
    all[ ind[clcy*3+2]+ 0]= ipos8[1]*mkj;
    all[ ind[clcy*3+2]+ 3]= ipos8[0]*mkj;
    all[ ind[clcy*3+2]+ 6]= ipos8[7]*mkj;
    up_neb_all=0;
  } else{
    clcy =0;
    all[0]=ipos8[3]*mkj;all[1]=ipos8[2]*mkj; all[2]=ipos8[1]*mkj;
    all[3]=ipos8[4]*mkj;all[4]= pos    *mkj; all[5]=ipos8[0]*mkj;
    all[6]=ipos8[5]*mkj;all[7]=ipos8[6]*mkj; all[8]=ipos8[7]*mkj;
    up_neb_all=1;
  }
  if(!td->smt.inited_neb) td->smt.inited_neb=1;
#define AIND(n,m) all[ind[clcy*3+n]+m]
  neb_pos[0]=AIND(0,0);neb_pos[1]=AIND(1,0);neb_pos[2]=AIND(2,0);
  neb_pos[3]=AIND(0,3);neb_pos[4]=AIND(1,3);neb_pos[5]=AIND(2,3);
  neb_pos[6]=AIND(0,6);neb_pos[7]=AIND(1,6);neb_pos[8]=AIND(2,6);
#undef AIND
  return up_neb_all;
}
#define LDZAIM(it,n,m) LDZA(tile,it,em+(ind[clcy*3+n]+m)*mkj,pm)
#define STZAIM(it,n,m) STZA(tile,it,em+(ind[clcy*3+n]+m)*mkj,pm)
#define LOAD_ZA(tile,em, pm)\
    LDZAIM(0,0,0);LDZAIM(1,1,0);LDZAIM(2,2,0);LDZAIM(3,0,3);\
    LDZAIM(4,2,3);LDZAIM(5,0,6);LDZAIM(6,1,6);LDZAIM(7,2,6)
#define STORE_EM(tile,em, pm) \
    STZAIM(0,0,0);STZAIM(1,1,0);STZAIM(2,2,0);STZAIM(3,0,3);\
    STZAIM(4,2,3);STZAIM(5,0,6);STZAIM(6,1,6);STZAIM(7,2,6);
__arm_new("za")
void diffract_smooth(double *et, double *ee, int iacb, int iace) ARM_STREAMING{
  svbool_t pm=VTRUE;
  int *ind=td->smt.ind;
  int VL=svcntd();
  int clcy=td->smt.clcy;
  VDOUBLE vacA, vacB, vacC, vacD;
  int neb_pos[9];
  DOUBLE c=td->smt.stc[4];
  vacA=VLD(ee);
  int i;
  for(i =iacb*mkj+VL;i<=iace*mkj;i+=VL*2){
    vacB=VLD(ee+i);
    vacC=MUL_N(vacA,c);
    VST(et+i-VL,vacC);
    vacA=VLD(ee+i+VL);
    vacD=MUL_N(vacB, c);
    VST(et+i, vacD);
  }
  vacC=MUL_N(vacA,c);
  VST(et+i-VL,vacC);
  VDOUBLE vac=VLD(td->smt.stcv);
  VDOUBLE vac0,vac1,vac2,vac3;
  svbool_t pg_m=VTRUE, pg_n=VTRUE;
  int vl4 =VL*4,vl3=VL*3,vl2=VL*2;
  bool up_neb_all=false;//是否更新全部的邻居节点
  double *em=td->smt.em;
  for(int iac =iacb;iac<=iace;iac++){
    if(gd->nsp[iac]!=1)continue;
    up_neb_all=getneb(neb_pos, iac);
    int kj,pos=iac*mkj;
    upem(em,et,neb_pos,false, up_neb_all);
    double *pos_em=em,*pos_ee=ee+iac*mkj;
    LOAD_ZA(0,em,pm);
    vac0=VLD(pos_ee);
    for(kj=0;kj< mkj; kj += vl2){
      pos_em+= VL; LOAD_ZA(1,pos_em, pm);
      pos_ee+= VL; vac1=VLD(pos_ee);
      MOPA(0,pg_m,pg_n, vac, vac0);
      STORE_EM(0,pos_em-VL,pm);
      pos_em +=VL;LOAD_ZA(0,pos_em,pm);
      pos_ee +=VL;vac0=VLD(pos_ee);
      MOPA(1,pg_m, pg_n, vac, vac1);
      STORE_EM(1,pos_em-VL, pm);
    }
    MOPA(0,pg_m,pg_n,vac, vac0);
    STORE_EM(0,pos_em,pm);
  }
  upem(em, et, neb_pos, true, up_neb_all);
}
#else
/*
 * et: dst var,size=mkj*nwpa,
 * ee: src var,size=mkj*nwpa,
 * iacb,iace:range of cal points ,0<iacb<=iace<=nwpc<=nwpa
 * ipos:inds of neighbor points,size=8*nwpa,8 is 8 neighbor
 * */
void diffract_smooth(double *ed,double *es, int iacb, int iace){
  double p1=1./16/8;
  for(int iac=iacb;iac<=iace;iac++){
    int it=iac*mkj;
    double sps[8];
    int *ip=(gd->ipos8+iac)->i8;
#ifdef USE_SMOOTH8
    for(int k=0;k<8;k++){
      sps[k]=ip[k]>0?p1:0;
    }
#else
    for(int k=0;k<8;k+=2){
      sps[k]=ip[k]>0?p1:0;
    }
#endif
    for(int kj=0;kj<mkj;kj++){
      double ev=es[kj+it];
      double et=0;
#ifdef USE_SMOOTH8
      for(int k=0;k<8;k++){
        et+=(es[kj+ip[k]*mkj]-ev)*sps[k];
      }
#else
      for(int k=0;k<8;k+=2){
        et+=(es[kj+ip[k]*mkj]-ev)*sps[k];
      }
#endif
      ed[kj+it]=ev+et;
    }
  }
}
#endif
void diffract_smooth_(double *et, double *ee,int *iacb, int *iace){
  diffract_smooth(et,ee, *iacb, *iace);
}
void init_SSmooth(){
}
