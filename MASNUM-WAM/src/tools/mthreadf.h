#ifndef MTHREADF_H_INCLUDED
#define PTI HTHREADINFO ti
#define PTIM HTHREADINFO ti,
#define VTI ti
#define VTIM ti,
#ifdef MTHREAD
#define MLINE ,__LINE__
#ifdef MMTHREAD
!MG
#define MSS(RFB)  call msetgrps  (RFB MLINE)
#define NMWS(RFB,nt)  call mwaitgrps (RFB,nt MLINE)
#define NMWSR(RFB,nt) call mwaitgrpsr(RFB,nt MLINE)
!GM        
#define SSM(RFB)  call gsetmain  (RFB MLINE) !set grp ok     
#define NSWM(RFB,nt)  call gwaitmain (RFB,nt MLINE) !wait mmthread  
#define NSWMR(RFB,nt) call gwaitmainr(RFB,nt MLINE) !wait mmthread  
!SG       
#define SSG(RFB)  call ssetgrp   (RFB MLINE) !set sub  ready 
#define NSWG(RFB,nt)  call swaitgrp  (RFB,nt MLINE) !wait grp ready 
#define NSWGR(RFB,nt) call swaitgrpr (RFB,nt MLINE) !wait grp ready 
!GS      
#define GSS(RFB)  call gsetsubs  (RFB MLINE) !set gmthread ok
#define NGWS(RFB,nt)  call gwaitsubs (RFB,nt MLINE) !wait sub thread
#define NGWSR(RFB,nt) call gwaitsubsr(RFB,nt MLINE) !wait sub thread
#else   
!MS    
#define MSS(RFB)  call msetsubs  (RFB MLINE)
#define NMWS(RFB,nt)  call mwaitsubs (RFB,nt MLINE)
#define NMWSR(RFB,nt) call mwaitsubsr(RFB,nt MLINE)
!SM   
#define SSM(RFB)  call  ssetmain (RFB MLINE) !set sub  ready 
#define NSWM(RFB,nt)  call swaitmain (RFB,nt MLINE) !wait grp ready 
#define NSWMR(RFB,nt) call swaitmainr(RFB,nt MLINE) !wait grp ready 
#endif
#define MWS(RFB) NMWS(RFB,1)  
#define MWSR(RFB) NMWSR(RFB,1) 
#define SWM(RFB) NSWM(RFB,1)  
#define SWMR(RFB) NSWMR(RFB,1) 
#define SWG(RFB) NSWG(RFB,1)  
#define SWGR(RFB) NSWGR(RFB,1) 
#define GWS(RFB) NGWS(RFB,1)  
#define GWSR(RFB) NGWSR(RFB,1) 

#endif
!#define DTIMES
#ifdef DTIMES
#define TBT(i)  call tscb(i)
#define TET(i)  call tsce(i)
#define TEBT(i) call tsceb(i)
#define TPT()   call prtsc()
#else
#define TBT(i)
#define TET(i)
#define TEBT(i)
#define TPT()
#endif
#endif
