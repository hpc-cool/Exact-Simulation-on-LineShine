#define DBGB0(msg)		if(pid==0)then; write(6,'(f8.3,"s ",i7," ",a)') Difftimer(),iwalltime(),msg;call flush(6);endif

!#################################################################################################
!-------------------------------------------------------------------------------------------------
!*************************************************************************************************
!-------------------------------------------------------------------------------------------------
!                                                     Copyright (C) 2013 Xunqiang Yin & Wei Zhao
!                                                     MODULE NAME : irrp_pemask_mod
!                                                     PRESENT VERSION : 2013-04-30
!
! --- USAGE : To convert time among different types: datestr, datevec, datenum.
! --- DEPEND: None
!
! --- NOTE for describing of subroutine / function :
!  A. The parameters bracketed with [], means optional parameter.
!  B. The describe for the parameters of subroutine / function, started with:
!   * It means input prameter;
!   # It means output prameter;
!   @ It means input and output prameter(it will be changed inside).
!
!-------------------------------------------------------------------------------------------------
! ***                                 INTERFACE DESCRIBE                                       ***
!-------------------------------------------------------------------------------------------------
!
!  1. sub. irrp_set_pemask : Do irregular partition using MASK and return the results in PEMASK.
!
!   irrp_set_pemask(mask_, pemask_, im_, jm_, npe_)
!
!   * integer :: mask_   = Mask for horizonal space/points in which 0 means land/useless points.
!   # integer :: pemask_ = Arrary to save partition reaults each point with the value equal to
!                          the ID of the PE that this point is belonged. This ID of PEs is
!                          started from 1 which is different with the way for MPI.
!   * integer :: im_     = The first dimension size of mask_ and pemask_.
!   * integer :: jm_     = The second dimension size of mask_ and pemask_.
!   * integer :: npe_    = The total number of PEs used for this partition.
!
!    The shape for each partition will be like the following.
!
!    !-----------------------------------------------------------------------------------!
!    !                                                                                   !
!    !                      Shape of points in each PE.                                  !
!    !                                                                                   !
!    !    (i1, j2) ------------------------(i4, j2)                                      !
!    !        +    -----------------------------------------------(i6, j6)               !
!    !        +    -----------------------------------------------    +                  !
!    !        +    -----------------------------------------------    +                  !
!    !        +    -----------------------------------------------    +                  !
!    !    (i1, j3) -----------------------------------------------    +                  !
!    !             +    ------------------------------------------    +                  !
!    !             +    ------------------------------------------    +                  !
!    !             +    -------------------------------------------------- (i2, j4)      !
!    !             +    ---------------------------------------------------    +         !
!    !             +    ---------------------------------------------------    +         !
!    !             +    ---------------------------------------------------    +         !
!    !             +    ---------------------------------------------------    +         !
!    !          (i5, j5) --------------------------------------------------    +         !
!    !                               (i3, j1) ----------------------------- (i2, j1)     !
!    !                                                                                   !
!    !-----------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------------------------
!*************************************************************************************************
!-------------------------------------------------------------------------------------------------
  module irrp_split_mod
  implicit none
  public :: irrp_set_pemask,pi_pos_type
  private
!-------------------------------------------------------------------------------------------------
  type pi_pos_type                          ! Self-defined type to record 2 dimensional index.
    !                                       ! This will be fast than 2 dimensional arrary.
    integer :: i                            ! Index of the first dimension
    integer :: j                            ! Index of the second dimension
  end type pi_pos_type
!-------------------------------------------------------------------------------------------------
  integer :: im                             ! First dimension size of MASK & PEMASK.
  integer :: jm                             ! Second dimension size of MASK & PEMASK.
  integer :: npe                            ! Number of PEs used in this partition.
  integer :: pid
  integer :: comm
  integer :: sumnp                          ! Sum of the computer points in MASK.
  integer :: ncols                          ! Number of columns during partition
  real(8) :: avenp                          ! Averaged number of computing points for each PE.

!-------------------------------------------------------------------------------------------------

  integer, pointer :: ncwp(:)               ! number of work points at every col
  integer, pointer :: ncwps(:)               ! number of work points to every cols
  integer, pointer :: mask(:, :)            ! The mask for 2 dimensional space in which 0 means
                                            ! land or the point is not need to be computed.
                                            ! The shape if this arrary is (im, jm).
  integer, pointer :: pemask(:, :)          ! It has the same shape with mask and it is used to
                                            ! record the partition results for each points with
                                            ! the value equal to the ID of the PE these points
                                            ! are belonged.
                                            ! NOTE: This ID is started from 1 but not from zero.
                                            !       It is different with the arranging method of
                                            !       PIDs in MPI.
!-------------------------------------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------------------------------
  ! sub. irrp_set_pemask : Do irregular partition using MASK and return the results in PEMASK.
  !    * integer :: mask_   = Mask for horizonal space/points: 0 means land/useless points.
  !    # integer :: pemask_ = Arrary to save partition reaults each point with the value equal to
  !                           the ID of the PE that this point is belonged. This ID of PEs is
  !                           started from 1 which is different with the way for MPI.
  !    * integer :: im_     = The first dimension size of mask_ and pemask_
  !    * integer :: jm_     = The second dimension size of mask_ and pemask_
  !    * integer :: npe_    = The total number of PEs used for this partition.
  !-----------------------------------------------------------------------------------------------

  subroutine irrp_set_pemask(mask_, im_, jm_, npe_,pid_,comm_)
    integer, intent(in) :: im_, jm_, npe_,comm_
    integer, target, intent(inout) :: mask_(im_, jm_)
    !integer, target, intent(out) :: pemask_(im_, jm_)
    integer  pid_
    integer,allocatable::ncwpt(:)
    integer ::ib,ie,ierr
    integer :: nn                           ! Temp index for sequentializing.
    integer :: i, j                         ! Temp index of the first & second dimension.
    real(8) :: delta                        ! Temp value used to determine the best partition.

    im = im_; jm = jm_; npe = npe_          ! Record the dimension size and number of PEs in the
    !                                       ! global value of this module.
    pid=pid_;comm=comm_
    mask => mask_; !pemask => pemask_        ! Using pointer to refer to the arrary of mask and
                                            ! pemask. This will avoid the copy of arraries which
                                            ! caused the computing speed is increased and the
                                            ! inquirement of memory is reduced.
    allocate(ncwps(0:im))
		DBGB0("setpemask 0")
		call setmask1(mask,im*jm);
		DBGB0("setpemask 0a")
#if 0
    nn=(im+npe)/npe
    if(nn<1)nn=1
    ib=pid*nn+1
    ie=ib+nn-1
    if(ie>im)ie=im
    allocate(ncwpt(0:nn-1))
    allocate(ncwp(nn*npe))
    ncwpt=0;ncwp=0
    do i = ib, ie                            ! Sequentializing computing points.
      do j = 1, jm                          ! This sequentializing is performed by connect
        if(mask(i, j) > 0)then              ! all the point-columns. The direction is from
          ncwpt(i-ib)=ncwpt(i-ib)+mask(i, j)
        endif
      enddo
    enddo
    !print'(100i8)',pid,im,ib,ie,nn,sum(ncwpt),ncwpt;
    call MPI_GATHER(ncwpt,nn,MPI_INTEGER4,ncwp ,nn,MPI_INTEGER4,0,comm,ierr)
    call MPI_BCAST(ncwp,im,MPI_INTEGER4,0,comm,ierr)
    !if(pid==0)print*,pid,nn,sum(ncwp),ncwp;
    deallocate(ncwpt);
#else
    allocate(ncwp(1:im))                    !
    nn = 0;    ncwp=0
    do j = 1, jm                          ! This sequentializing is performed by connect
    	do i = 1, im                            ! Sequentializing computing points
      	                                    ! all the point-columns. The direction is f
      	                                    ! down to up and from left to right.
         ncwp(i)=ncwp(i)+mask(i, j)
     	enddo
   	enddo
#endif
    ncwps(0)=0;nn=0
    do i = 1, im                            !
    	nn=nn+ncwp(i)
    	ncwps(i)=ncwps(i-1)+ncwp(i)
    enddo
    !if(pid==1)print*,pid,nn,sum(ncwp),ncwp;
		DBGB0("setpemask 1")
    sumnp = nn;                             ! Set the sum number of computing points
    avenp = sumnp / float(npe)              ! Averaged computing points for each PE.

    call set_pemask_iteration(1, 0.d0)      ! Do partition by iteration based on the gridient of
                                            ! the global maximum aspect ratio of each PE respect
                                            ! to parameter of delta.
    nullify(mask, pemask)                   ! Finalizing for partition: release temporary arraries.
                                            ! Pointers need to be nullified.
		deallocate(ncwp,ncwps);
		DBGB0("setpemask E")
  end subroutine irrp_set_pemask

  !-----------------------------------------------------------------------------------------------
  ! sub. set_pemask_iteration : Do partition using iteration based on the gridient of the global
  !                             maximum aspect ratio of each PE respect to parameter of delta.
  !
  !   set_pemask_iteration(flag[, delta_, ratio_])
  !
  !   @ integer :: flag   = The flag for update PEMASK. 0: not update PEMASK; else: update PEMASK.
  !   * real(8) :: delta_ = Parameter used to compute the global maximum aspect ratio of each PE.
  !   * real(8) :: ratio_ = Aspect ratio for a given delta.
  !-----------------------------------------------------------------------------------------------

  subroutine set_pemask_iteration(flag, delta_, ratio_)
    integer, intent(in)            :: flag
    real(8), intent(in),  optional :: delta_
    real(8), intent(out), optional :: ratio_

    real(8) :: delta
    integer :: ncols0, i0, i1, n, i, j
    real(8) :: ratio, dratio
    real(8) :: grad, lamb, grad0
    real(8) :: minr1, rm0, mindelta, delta1, dstep, lastr1, r1all0, dt0, dt1
    real(8) :: mmr1, minlamb, small

    delta = 0.                              ! set delta as zero for first guess.
    if(present(delta_))delta = delta_       ! Use the given delta.

 		DBGB0("setpemaskI 0")
    call set_pemask_one_cycle(delta, ratio, dratio, flag)
 		DBGB0("setpemaskI 1")
                                            ! Using the given value of delta to obtain the first
                                            ! guess of the maximum aspect ratio.
    if(flag /= 0)then
      if(present(ratio_))ratio_ = ratio     ! Can be used to do parallized partition.
      return
    endif
    rm0 = 1;  lamb = 0.                     ! Set the initial value of minr1, rm0 and lamb.
    minr1 = ratio; mindelta = delta         ! Record the minimum aspect ratio during the iteration.
    ncols0 = ncols                          ! Record the number of columns at the first guess.
    dstep   = 5. / ncols0                   ! Set the step size for the deepest method using the
    small   = dstep * 0.0001                ! number of PE-columns.
    minlamb = dstep * 0.001;                ! Limit for minimum lamb.
    mmr1    = 0.1 * ncols0 / im             ! Limit for minimum aspect ratio.
    i0      = 0                             ! To record the iteration cycle.
    grad0   = 0                             ! Set the first gridient to zero for the first
                                            ! iteration cycle.
    i1 = 0                                  ! To record the times of calling set_pemask_one_cycle
    dt0 = -0.5; dt1 = 0.5                   ! Bound limit for delta.
    !---------------------------------------------------------------------------------------------
    ! Deepest Method:
    !
    !---------------------------------------------------------------------------------------------
    do while(.true.)                             !{ find best delta.
      lastr1 = ratio                        ! Record the last ratio.
      call set_pemask_one_cycle(delta+small, ratio, dratio, 0)
 		DBGB0("setpemaskI 2")
                                            ! Using a small value added to delta to obtain
                                            ! the second guess of the maximum aspect ratio.
      i1 = i1 + 1                           ! Record the times calling the set_pemask_one_cycle
                                            ! during the iteration.
      grad = (ratio*ratio - lastr1*lastr1) / small
                                            ! Compute the gradient using the first and second
                                            ! guess of the maximum aspect ratio.
      if(ncols /= ncols0)then               ! Once the number of PE collums is different with
                                            ! former iteration cycle, the small value added
      	                                    ! to delta is set to be a negitive one.
        call set_pemask_one_cycle(delta-small, ratio, dratio, 0)
 		DBGB0("setpemaskI 3")
        i1 = i1 + 1
        if(ncols /= ncols0)exit             ! If the number of PE columns is still different
                                            ! with former iteration cycle, exit the iteration.
        grad = (ratio*ratio - lastr1*lastr1) / small
                                            ! Compute the gradient using based on the negitive
                                            ! small value added to delta.
      endif
      if(abs(minr1) > abs(ratio))then
        minr1 = ratio; mindelta = delta     ! Record the optimal delta during the iteration.
      endif
      if(abs(grad) < small)exit             ! If the gradient is smaller than the added small value,
                                            ! finish the iteration.
      if(grad0 * grad > 0)then
      	exit                                ! If the current gradient and the former gradient are
      	                                    ! in the same direction, finish this iteration.
      else
        lamb = grad * dstep                 ! Set the step size for steepest iteration method.
      endif

      grad0 = grad                          ! Record the current gradient as the former one in
                                            ! next iteration cycle.
      r1all0 = 100; i1 = 0
      !-------------------------------------------------------------------------------------------
      ! Along the direction of the current gradient,
      ! search for a delta which will cause a minimum aspect ratio.
      !-------------------------------------------------------------------------------------------
      do while (.true.) ! { find best lamb
        dratio = 0; r1all0 = ratio
        if(lamb > 0)then
          do while(delta - lamb < dt0)
             lamb = lamb / 2
          enddo
        else
          do while(delta - lamb > dt1)
            lamb = lamb / 2
          enddo
        endif
        delta = delta - lamb                ! Attempt new delta for minimum aspect ratio.
        call set_pemask_one_cycle(delta, ratio, dratio, 0)
 		DBGB0("setpemaskI 4")
        i1 = i1 + 1
        do while(ncols /= ncols0)           ! Keep the number of PE columns same as former one.
          if(abs(lamb) < minlamb)exit
          delta = delta + lamb; lamb = lamb / 5; delta = delta - lamb
          call set_pemask_one_cycle(delta, ratio, dratio, 0)
 		DBGB0("setpemaskI 5")
          i1 = i1 + 1
        enddo
        do while(r1all0*ratio<0)            ! If the aspect ratio has a different sign with former,
                                            ! search in a different direction.
          if(abs(lamb) < minlamb)exit       ! if the step size is smaller enough, finish this search.
          delta = delta + lamb; lamb = lamb / 3; delta = delta - lamb
          call set_pemask_one_cycle(delta, ratio, dratio, 0)
 		DBGB0("setpemaskI 6")
          i1 = i1 + 1                       ! Attempt a different step size.
        enddo
        if(abs(minr1) > abs(ratio))then
          minr1 = ratio; mindelta = delta   ! record the minimum aspect ratio and delta.
        endif
        !-----------------------------------------------------------------------------------------
        ! if the step size is smaller enough or the minimum aspect ratio is smaller enough,
        ! finish this searching.
        !-----------------------------------------------------------------------------------------
        if(abs(lamb)<minlamb .or. abs(minr1)<mmr1)exit

        if(abs(r1all0) < abs(ratio))then    ! adjust for next iteration cycle: dt0, dt1 and dstep.
          if(lamb < 0)then
            if(dt0 < delta + lamb)dt0 = delta + lamb
            if(dt1 > delta)dt1 = delta
          else
            if(dt1 > delta + lamb)dt1 = delta + lamb
            if(dt0 < delta)dt0 = delta
          endif
          dstep = dstep / 2
          exit
        endif
      enddo                                 !} find best lamb
      i0 = i0 + 1
      !-------------------------------------------------------------------------------------------
      ! if the step size is smaller enough or the minimum aspect ratio is smaller enough,
      ! exit this iteration.
      !-------------------------------------------------------------------------------------------
      if(abs(lamb) < minlamb .or. abs(minr1) < mmr1)exit
    enddo                                   !} while(.true.) ! find best delta.
    delta = mindelta                        ! Using the optimal delta to update the PEMASK.
 		DBGB0("setpemaskI 7")
    call set_pemask_one_cycle(delta, ratio, dratio, 1)
 		DBGB0("setpemaskI 8")
  end subroutine set_pemask_iteration

  !-----------------------------------------------------------------------------------------------
  ! sub. set_pemask_one_cycle: do the partition using a given delta.
  !
  !   set_pemask_one_cycle(delta, ratio, dratio, flag)
  !
  !   * real(8) :: delta  = The parameter used for computing aspect ratio of each PE.
  !   # real(8) :: ratio  = Maximum aspect ratio without thinking of delta.
  !   # real(8) :: dratio = Maximum aspect ratio with thinking of delta.
  !   * integer :: flag   = The flag for update PEMASK. 0: not update PEMASK; else: update PEMASK.
  !-----------------------------------------------------------------------------------------------

  subroutine set_pemask_one_cycle(delta, ratio, dratio, flag)
    real(8), intent(in)  :: delta
    real(8), intent(out) :: ratio, dratio
    integer, intent(in)  :: flag

    integer :: n1, n2, s1, s2, m, i, minline, ml, mml, onemore, mn
    real(8) :: r1, r1t, rt, r, rt0, rt1
    integer :: i1,i2,mi
    integer, save, allocatable :: iflag(:), jflag(:), niflag(:), njflag(:)
    integer, save :: imt = 0, jmt = 0
    if(imt /= im .or. jmt /= jm)then        ! Initialize the temporary arraries.
      if(imt /= 0)then
        deallocate(iflag, jflag, niflag, njflag)
      endif
      imt = im; jmt = jm
      allocate(iflag(imt), jflag(jmt), niflag(imt), njflag(jmt))
    endif
    !if(flag /= 0)pemask = 0                 ! initializing pemask.
    n2 = 0                                  ! initializing n2: end PE of last PE column.
    ncols = 0                               ! initializing number of columns.
    ratio = 0.                              ! initializing global maximum ratio without delta.
    dratio = 0.                             ! initializing global maximum ratio with delta.
    mml = jm * 0.7                          ! estimated minimum PEs for a point-column.
    minline = sqrt(npe * 1.0 * jm / im)/3    ! One third of the expected minimum of number of PEs
                                            ! for one PE-column.		
		if(minline <1)minline =1
		mi=sqrt(avenp)*0.7		
    s2 = 1;		i2 = 1;
    do while (n2 /= npe)                    !{ Find s2, n2
      s1=s2;i1=i2
      n1 = n2 + 1;
      m=n1
      ncols = ncols + 1; r1 = 1000000.
      if(m >= npe) m = npe - 1
      !-------------------------------------------------------------------------------------------
      ! Find a best n2 which will lead to a smaller aspect ratio with a given delta.
      ! The maximum ratio and dratio will be recorded during this process.
      !-------------------------------------------------------------------------------------------
			i2=i1+mi-1			
      do while(m < npe)
      	i2=i2+1
      	s2=avenp*(m+1);
      	if(ncwps(i2)<s2)cycle
      	m=(ncwps(i2)+0.001)/avenp      	
      	if(i2<im)then
	        ml = (i2-i1) * (m - n1 + 1)
	                                            ! ml / mml : gross aspect ratio of PEs from n1 to m.
	        if(ml < mml)cycle                   ! Speed up by adding more PEs.
	                                            ! If the rest points of this point-colomn is greater
	                                            ! than the points of one PE, add more PEs for these
	                                            ! PE-column and the computing of aspect ratio is not
	                                            ! needed.
        endif
        s2 = nint(avenp * m + 1.e-5)        ! end point of current PE with id of m.

      	if((m + minline >= npe) )then
      		m=npe;i2=im;s2=sumnp
      	endif
        call compute_ratio(s1, s2, n1, m, rt,i1,i2)
                                            ! Compute the aspect ratio for PEs from n1 to m.
                                            ! rt: the pure aspect ratio around to 1.
        rt0 = 1 - rt                        ! The aspect ratio compared with 1.
        rt1 = 1 + delta - rt                ! The aspect ratio compared with (1 + delta).
        if(abs(rt1) < abs(r1))then
          r1 = rt1; r1t = rt0; n2 = m       ! If this ratio is accetable, set n2 eqqual to m,                                             ! and record the aspect ratios.
        endif
	      if(abs(rt1) > abs(r1)*2)then
	      	!if(n2<n1)n2=m
	        exit                              ! If this ratio is too big, finish this searching.
	      endif
      enddo
      if(n1>n2)then
      	n2=npe
      endif
      !-------------------------------------------------------------------------------------------
      ! Patch the last PE column.
      ! n1 > n2 means last n2 is not fund.
      ! n2 + minline >= npe means if use n2, the rest PEs after n2 is too less.
      ! n2 /= npe: the last PE-colum is incorrect.
      !-------------------------------------------------------------------------------------------

      if(abs(ratio)  < abs(r1 ))ratio  = r1 ! Record the ratio and dratio for this column of PEs.
      if(abs(dratio) < abs(r1t))dratio = r1t
      if(flag/=0)call set_pemask_columns(n1, n2)
                                            ! Set pemask for this PE column.
    enddo                                   !} (n2 /= npe) ! find n2

    if(flag /= 0)then
      imt = 0; jmt = 0                      ! Finalizing temporary arraries when the setting
                                            ! of pemask is finished.
      deallocate(iflag, jflag, niflag, njflag)
    endif
    !----------------------------------------------------------------
    contains
    !----------------------------------------------------------------
    ! inner subroutine to compute the aspect ratio of one PE column.
    !----------------------------------------------------------------
    subroutine compute_ratio(s1, s2, n1, n2,  rt,i1,i2)
      integer, intent(in)  :: s1, s2, n1, n2,i1,i2
      real(8), intent(out) :: rt
      integer :: s, i, j,j1,j2,is,ie,mt
      real(8) :: sqmp, sum_niflag, sum_njflag, sum_ij, sum_ji

      niflag = 0                            ! To record the number of none-empty-columns.
      njflag = 0                            ! To record the number of none-empty-lines.
      iflag  = 0                            ! To record weights of each point-columns.
      jflag  = 0                            ! To record weights of point-lines.
  	j1=ssi2j(s1,i1)
  	j2=ssi2j(s2,i2)
    do j = 1, jm
      is = i1; if(j<j1)is = i1 + 1          ! Start i for current j.
      ie = i2; if(j>j2)ie = i2 - 1          ! End i for current j.                                             ! They are determined by the PE shape.
      do i = is, ie
      	mt=mask(i, j)
        !if(mt /= 0)then
  	 	    niflag(i)=ior(niflag(i),mt); iflag(i) = iflag(i) + mt
          njflag(j)=ior(njflag(j),mt); jflag(j) = jflag(j) + mt        	
        !endif
      enddo
    enddo

      sum_niflag = sum(niflag)              ! Total number of none-empty-columns.
      sum_njflag = sum(njflag)              ! Total number of none-empty-lines.
      sum_ij = sum((iflag/sum_njflag)**0.3) ! Weighted number of point-columns.
      sum_ji = sum((jflag/sum_niflag)**0.3) ! Weighted number of point-lines.
                                            ! Here, a power of 0.3 is used to emphasizing
                                            ! the effective of land points. If use 1, it will be
                                            ! no emphasizing.
      rt = sum_ij *dble(n2 - n1 + 1)/sum_ji ! aspect ratio of this PE column.
    end subroutine compute_ratio
  end subroutine set_pemask_one_cycle
	
  integer function ssi2j(s,i)
  	integer,intent(in)::  s,i
  	integer st,j
    st=ncwps(i-1)
    ssi2j=0
    do j = 1, jm                          ! This sequentializing is performed by connect
      if(mask(i, j) > 0)then              ! all the point-columns. The direction is from
        st = st + 1                       ! down to up and from left to right.
        if(st==s)then
        	ssi2j=j
        	exit
        endif
      endif
    enddo
	end function ssi2j
  subroutine ss2ij(s,i,j)
  	integer,intent(in)::s
  	integer,intent(out)::i,j
  	integer st
		integer ,save::iold
    if(s <= 1)then
    	iold=1;
    endif
    do i=iold,im
    	if(s<=ncwps(i))exit
    enddo
    st=ncwps(i-1)
    do j = 1, jm                          ! This sequentializing is performed by connect
      if(mask(i, j) > 0)then              ! all the point-columns. The direction is from
        st = st + 1                       ! down to up and from left to right.
        if(st==s)exit
      endif
    enddo
    iold=i;
	end subroutine ss2ij	
!-------------------------------------------------------------------------------------------------
! sub. set_pemask_columns: Set pemask from n1 to n2 as a column of PEs.
!
!   set_pemask_columns(n1, n2)
!
!   * integer :: n1 = The start ID of PEs needed to be set for PEMASK.
!   * integer :: n2 = The end ID of PEs needed to be set for PEMASK.
!-------------------------------------------------------------------------------------------------
	
  subroutine set_pemask_columns(n1, n2)
    integer, intent(in) :: n1, n2
    integer :: i, nwp, s1, s2, i1, i2, j1, j2, idpe, j, is, ie, nwpe
    !-----------------------------------------------------------------
    ! --- set i1, i2, j1, j2 for points: s1 --> s2, pe: n1 --> n2
    ! The start point which should belong to PE with ID of n1.
    !-----------------------------------------------------------------
    s1 = nint(avenp * (n1-1) + 1.e-5) + 1  ! the end of last column
    !-----------------------------------------------------------------
    ! The end point which should belong to PE with ID of n2.
    ! For last PE, s2 should same as sum number of computing points.
    !-----------------------------------------------------------------
    s2 = nint(avenp *  n2    + 1.e-5); if(n2 == npe)s2 = sumnp
    !-----------------------------------------------------------------
    ! Set (i1 & j1),(i2,j2): the index in matrix for the start point of s1&s2.
    !-----------------------------------------------------------------
    call ss2ij(s1,i1,j1)
    call ss2ij(s2,i2,j2)
		!if(pid==0)print'("pemask", 8i8)',n1,n2,s1,s2,i1,i2,j1,j2,sumnp
    !-----------------------------------------------------------------------------------!
    !                                                                                   !
    !                      Shape of points in each PE.                                  !
    !                                                                                   !
    !    (i1, j2) ------------------------(i4, j2)                                      !
    !        +    -----------------------------------------------(i6, j6)               !
    !        +    -----------------------------------------------    +                  !
    !        +    -----------------------------------------------    +                  !
    !        +    -----------------------------------------------    +                  !
    !    (i1, j3) -----------------------------------------------    +                  !
    !             +    ------------------------------------------    +                  !
    !             +    ------------------------------------------    +                  !
    !             +    -------------------------------------------------- (i2, j4)      !
    !             +    ---------------------------------------------------    +         !
    !             +    ---------------------------------------------------    +         !
    !             +    ---------------------------------------------------    +         !
    !             +    ---------------------------------------------------    +         !
    !          (i5, j5) --------------------------------------------------    +         !
    !                               (i3, j1) ----------------------------- (i2, j1)     !
    !                                                                                   !
    !-----------------------------------------------------------------------------------!

    idpe = n1 - 1
    nwp = nint(avenp * (n1-1) + 1.e-5)      ! The number of points with ID of 1 to (n1-1).
    nwpe=nint(avenp * idpe + 1.e-5)         ! The number of points with ID of 1 to idpe.
    do j = 1, jm
      is = i1; if(j<j1)is = i1 + 1          ! Start i for current j.
      ie = i2; if(j>j2)ie = i2 - 1          ! End i for current j.
                                            ! They are determined by the PE shape.

      do i = is, ie
        if(mask(i, j) /= 0)then
          nwp = nwp + 1                     ! New part Begin
          if(nwp > nwpe)then
            idpe = idpe + 1                 ! The new ID of PE.
            nwpe = nint(avenp*idpe + 1.e-5) ! The number of points with ID of 1 to idpe.
            if(idpe >= n2)then              ! For the PE of n2, the nunber of points                                             ! with ID of 1 to idpe should be s2.
              idpe = n2; nwpe = s2
            endif
          endif
          mask(i, j) = idpe                 ! Set the point (i, j) belong to PE with id of idpe.
        endif
      enddo
    enddo
  end subroutine set_pemask_columns

!-------------------------------------------------------------------------------------------------

  end module irrp_split_mod

!-------------------------------------------------------------------------------------------------
!#################################################################################################
