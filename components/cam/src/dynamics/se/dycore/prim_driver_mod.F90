#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"file: ",__FILE__," line: ",__LINE__," ithr: ",hybrid%ithr
#define _DBG_
module prim_driver_mod
  use kinds, only : real_kind, longdouble_kind
  use cam_logfile, only : iulog
  use dimensions_mod, only : np, nlev, nelem, nelemd, nelemdmax, GlobalUniqueCols, qsize, nc,nhc
  use hybrid_mod, only : hybrid_t, config_thread_region, PrintHybrid
  use filter_mod, only : filter_t
  use derivative_mod, only : derivative_t
  use fvm_mod, only : fvm_init1,fvm_init2, fvm_init3, physgrid_init2
  use fvm_control_volume_mod, only : fvm_struct

  use element_mod, only : element_t, timelevels,  allocate_element_desc
  use thread_mod, only :  omp_get_num_threads
  use thread_mod , only : max_num_threads, horz_num_threads, vert_num_threads, tracer_num_threads
  use perf_mod, only : t_startf, t_stopf
  use prim_init, only : gp, fvm_corners, fvm_points

  implicit none
  private
  public :: prim_init2, prim_run_subcycle, prim_finalize
  public :: smooth_topo_datasets, prim_set_mass

  type (filter_t)       :: flt             ! Filter struct for v and p grid
  type (filter_t)       :: flt_advection   ! Filter struct for v grid for advection only
  real*8  :: tot_iter

contains

!=============================================================================!

  subroutine prim_init2(elem, fvm, hybrid, nets, nete, tl, hvcoord)
    use dimensions_mod, only: irecons_tracer, qsize_condensate_loading,ldry_mass_vertical_coordinates
    use dimensions_mod, only: fv_nphys, ntrac
    use parallel_mod, only : parallel_t, haltmp, syncmp, abortmp, iam 
    use time_mod, only : timelevel_t, tstep, phys_tscale, timelevel_init, nendstep, smooth, nsplit, TimeLevel_Qdp
    use prim_state_mod, only : prim_printstate
    use filter_mod, only : filter_t, fm_filter_create, taylor_filter_create, &
         fm_transfer, bv_transfer
    use control_mod, only : runtype, integration, filter_mu, filter_mu_advection, test_case, &
         debug_level, vfile_int, filter_freq, filter_freq_advection, &
         transfer_type, vform, vfile_mid, filter_type, kcut_fm, wght_fm, p_bv, &
         s_bv, topology,columnpackage, rsplit, qsplit, rk_stage_user,&
         sub_case, &
         limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
         hypervis_subcycle_q
    use fvm_control_volume_mod, only : fvm_supercycling
    use physical_constants, only : p0
    use derivative_mod, only : subcell_integration

    use thread_mod, only : max_num_threads, omp_get_thread_num
    use derivative_mod, only : derivinit, interpolate_gll2fvm_points, v2pinit
    use global_norms_mod, only : test_global_integral, print_cfl
    use hybvcoord_mod, only : hvcoord_t
    use prim_advection_mod, only: prim_advec_init2, deriv
    use prim_advance_mod, only: prim_advance_init

    use fvm_control_volume_mod, only : n0_fvm, np1_fvm
    use control_mod, only : tracer_transport_type


    type (element_t), intent(inout) :: elem(:)
    type (fvm_struct), intent(inout)    :: fvm(:)
    type (hybrid_t), intent(in) :: hybrid

    type (TimeLevel_t), intent(inout)    :: tl              ! time level struct
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)


    ! ==================================
    ! Local variables
    ! ==================================

    real (kind=real_kind) :: dt              ! "timestep dependent" timestep
!   variables used to calculate CFL
    real (kind=real_kind) :: dtnu            ! timestep*viscosity parameter
    real (kind=real_kind) :: dt_dyn_vis      ! viscosity timestep used in dynamics
    real (kind=real_kind) :: dt_tracer_vis      ! viscosity timestep used in tracers

    real (kind=real_kind) :: dp

    real (kind=real_kind) :: ps(np,np)

    character(len=80)     :: fname
    character(len=8)      :: njusn
    character(len=4)      :: charnum

    real (kind=real_kind) :: Tp(np)     ! transfer function

    logical :: ltrac(qsize)

    integer :: simday
    integer :: i,j,k,ie,iptr,t,q
    integer :: ierr
    integer :: nfrc
    integer :: n0_qdp


    ! ==========================
    ! begin executable code
    ! ==========================
    call prim_advance_init(hybrid%par,elem,integration)

    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)
    end if


    ! compute most restrictive dt*nu for use by variable res viscosity:
    if (tstep_type == 0) then
       ! LF case: no tracers, timestep seen by viscosity is 2*tstep
       dt_tracer_vis = 0
       dt_dyn_vis = 2*tstep
       dtnu = 2.0d0*tstep*max(nu,nu_div)
    else
       ! compute timestep seen by viscosity operator:
       dt_dyn_vis = tstep
       if (qsplit>1 .and. tstep_type == 1) then
          ! tstep_type==1: RK2 followed by LF.  internal LF stages apply viscosity at 2*dt
          dt_dyn_vis = 2*tstep
       endif
       dt_tracer_vis=tstep*qsplit

       ! compute most restrictive condition:
       ! note: dtnu ignores subcycling
       dtnu=max(dt_dyn_vis*max(nu,nu_div), dt_tracer_vis*nu_q)
       ! compute actual viscosity timesteps with subcycling
       dt_tracer_vis = dt_tracer_vis/hypervis_subcycle_q
       dt_dyn_vis = dt_dyn_vis/hypervis_subcycle
    endif


    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call Prim_Advec_Init2(hybrid, fvm_corners, fvm_points)

    ! ================================================
    ! fvm initialization
    ! ================================================
    if (fv_nphys>0) then
      call fvm_init2(elem,fvm,hybrid,nets,nete)
      call physgrid_init2(elem,fvm,hybrid,nets,nete,irecons_tracer)
    endif
    ! ====================================
    ! In the semi-implicit case:
    ! initialize vertical structure and
    ! related matrices..
    ! ====================================
    ! ==========================================
    ! Initialize pressure and velocity grid
    ! filter matrix...
    ! ==========================================
    if (transfer_type == "bv") then
       Tp    = bv_transfer(p_bv,s_bv,np)
    else if (transfer_type == "fm") then
       Tp    = fm_transfer(kcut_fm,wght_fm,np)
    end if
    if (filter_type == "taylor") then
       flt           = taylor_filter_create(Tp, filter_mu,gp)
       flt_advection = taylor_filter_create(Tp, filter_mu_advection,gp)
    else if (filter_type == "fischer") then
       flt           = fm_filter_create(Tp, filter_mu, gp)
       flt_advection = fm_filter_create(Tp, filter_mu_advection, gp)
    end if



    if (hybrid%masterthread) then
       if (filter_freq>0 .or. filter_freq_advection>0) then
          write(iulog,*) "transfer function type in preq=",transfer_type
          write(iulog,*) "filter type            in preq=",filter_type
          write(*,'(a,99f10.6)') "dynamics: I-mu + mu*Tp(:) = ",&
               (1-filter_mu)+filter_mu*Tp(:)
          write(*,'(a,99f10.6)') "advection: I-mu + mu*Tp(:) = ",&
               (1-filter_mu_advection)+filter_mu_advection*Tp(:)
       endif
    endif

    !$OMP BARRIER
    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if
    !$OMP BARRIER

    if (topology /= "cube") then
       call abortmp('Error: only cube topology supported for primaitve equations')
    endif

     !
     ! set diagnostic pswet
     !     
     do ie=1,nelemd
       elem(ie)%state%pswet(:,:) = elem(ie)%state%ps(:,:,1)
       if (ldry_mass_vertical_coordinates) then
         !
         ! %ps is dry - add water variable(s)
         !
         do k=1,nlev
           do j = 1,np
             do i = 1,np   
               elem(ie)%state%pswet(i,j) = elem(ie)%state%pswet(i,j)+sum(elem(ie)%state%Qdp(i,j,k,1:qsize_condensate_loading,1))
             end do
           end do
         end do
       else
         do j=1,np
           do i=1,np
             elem(ie)%state%pswet(i,j) = elem(ie)%state%ps(i,j,1)
           end do
         end do
       end if
     end do


     if (fv_nphys>0) then
       call fvm_init3(elem,fvm,hybrid,nets,nete,irecons_tracer)
     endif

! for restart runs, we read in Qdp for exact restart, and rederive Q

    !
    ! phl - still necessary
    !
    if (runtype==1) then
       call TimeLevel_Qdp( tl, qsplit, n0_qdp)
       do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, t, q, i, j, dp)
#endif
          do k=1,nlev    !  Loop inversion (AAM)
             do t=tl%n0,tl%n0
                do q=1,qsize
                   do i=1,np
                      do j=1,np
                         dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(i,j,t)
                         elem(ie)%state%Q(i,j,k,q)=elem(ie)%state%Qdp(i,j,k,q, n0_qdp)/dp
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


    ! timesteps to use for advective stability:  tstep*qsplit and tstep
    call print_cfl(elem,hybrid,nets,nete,dtnu)

    if (hybrid%masterthread) then
       ! CAM has set tstep based on dtime before calling prim_init2(),
       ! so only now does HOMME learn the timstep.  print them out:
       write(iulog,'(a,2f9.2)') "dt_remap: (0=disabled)   ",tstep*qsplit*rsplit

       if (ntrac>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (fvm)          ",tstep*qsplit*fvm_supercycling
       end if
       if (qsize>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (SE), per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       end if
       write(iulog,'(a,2f9.2)')    "dt_dyn:                  ",tstep
       write(iulog,'(a,2f9.2)')    "dt_dyn (viscosity):      ",dt_dyn_vis
       write(iulog,'(a,2f9.2)')    "dt_tracer (viscosity):   ",dt_tracer_vis


       if (phys_tscale/=0) then
          write(iulog,'(a,2f9.2)') "CAM physics timescale:       ",phys_tscale
       endif
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*nsplit*qsplit*max(rsplit,1)

       if (ldry_mass_vertical_coordinates) then
          write(iulog,*) "Using dry-mass vertical coordinates"
       else
          write(iulog,*) "Using wet-mass vertical coordinates"
       end if
    end if


    if (hybrid%masterthread) write(iulog,*) "initial state:"
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)

  end subroutine prim_init2

!=======================================================================================================!


  subroutine prim_run_subcycle(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)
!
!   advance all variables (u,v,T,ps,Q,C) from time t to t + dt_q
!
! for the RK schemes: 
!   input:
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
! for the implicit schemes: 
!
!   input:
!       tl%nm1   variables at t-1 level are stored fro BDF2 scheme 
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!       generally dt_q = t for BDF2, so its t+1
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, timelevel_qdp, nsplit
    use control_mod, only: statefreq,&
          qsplit, rsplit, test_cfldep, disable_diagnostics
    use prim_advance_mod, only : applycamforcing
    use prim_advance_mod, only : calc_tot_energy_dynamics
    use prim_state_mod, only : prim_printstate!phlxxx , prim_diag_scalars, prim_energy_halftimes
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use prim_advection_mod, only : vertical_remap
    use fvm_control_volume_mod, only : n0_fvm
    use thread_mod, only : omp_get_thread_num
    use dimensions_mod        , only : qsize_condensate_loading, ldry_mass_vertical_coordinates
    use perf_mod   , only : t_startf, t_stopf
    use fvm_mod    , only : fill_halo_fvm
    use dimensions_mod, only: ntrac,fv_nphys

    type (element_t) , intent(inout)        :: elem(:)

    type(fvm_struct), intent(inout) :: fvm(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                     :: nsubstep  ! nsubstep = 1 .. nsplit
    real(kind=real_kind) :: st, st1, dp, dt_q, dt_remap
    integer :: ie, t, q,k,i,j,n
    integer :: n0_qdp,np1_qdp,r, nstep_end
    integer :: ithr,region_num_threads

    real (kind=real_kind)                          :: maxcflx, maxcfly
    real (kind=real_kind) :: dp_np1(np,np),dp_np1_inv(np,np)
    logical :: compute_diagnostics
    type (hybrid_t) :: vybrid

    ! ===================================
    ! Main timestepping loop
    ! ===================================
    dt_q = dt*qsplit
    dt_remap = dt_q
    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       dt_remap=dt_q*rsplit   ! rsplit=0 means use eulerian code, not vert. lagrange
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif



    ! compute diagnostics and energy for STDOUT
    ! compute energy if we are using an energy fixer
    compute_diagnostics=.false.

    if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
       compute_diagnostics=.true.
    endif

    if(disable_diagnostics) compute_diagnostics=.false.


    call TimeLevel_Qdp( tl, qsplit, n0_qdp)

    call calc_tot_energy_dynamics(elem,nets,nete,tl%n0,n0_qdp,'dAF')
    call ApplyCAMForcing(elem, fvm, hvcoord,tl%n0,n0_qdp, dt_remap,nets,nete,nsubstep,hybrid)



    ! loop over rsplit vertically lagrangian timesteps
    call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,1)
    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,r)
    enddo
    ! defer final timelevel update until after remap and diagnostics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  apply vertical remap
    !  always for tracers
    !  if rsplit>0:  also remap dynamics and compute reference level ps
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    ! note: time level update for fvm tracers takes place in fvm_mod

    call calc_tot_energy_dynamics(elem,nets,nete,tl%np1,np1_qdp,'dAD')

    call t_startf('vertical_remap')
    call vertical_remap(hybrid,elem,fvm,hvcoord,dt_remap,tl%np1,np1_qdp,n0_fvm,nets,nete)
    call t_stopf('vertical_remap')



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! lnps (we should get rid of this)
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps(:,:,tl%np1))!phl - necessary?
#if (defined COLUMN_OPENMP)
       !$omp parallel do default(shared), private(k,q,j,i,dp_np1,dp_np1_inv)
#endif
       elem(ie)%state%pswet(:,:)=hvcoord%hyai(1)*hvcoord%ps0 !ptop
       do k=1,nlev    !  Loop inversion (AAM)
           !if (k == 1) then
           !  write(*,*) "In prim run there are ", omp_get_num_threads(), " in the current team in parallel region"
           !endif
          dp_np1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(:,:,tl%np1)
          elem(ie)%state%dp3d(:,:,k,tl%np1) = dp_np1(:,:)
!         dp_np1(:,:) = elem(ie)%state%dp3d(:,:,k,tl%np1) !bug in previous version - dp3d not updated for dp_coupling
!          dp_np1_inv(:,:) = 1.0_real_kind/dp_np1(:,:)
          do q=1,qsize
            !             elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*dp_np1_inv(:,:)
            elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp_np1(:,:)
          enddo         
          !
          ! update full surface pressure (diagnostic)
          !          
          if (ldry_mass_vertical_coordinates) then
             do j=1,np
                do i=1,np
                   elem(ie)%state%pswet(i,j)=elem(ie)%state%pswet(i,j)+&
                        (1.0D0+SUM(elem(ie)%state%q(i,j,k,1:qsize_condensate_loading)))*dp_np1(i,j)
                end do
             end do
          end if
        enddo
        if (.not.ldry_mass_vertical_coordinates) then
          elem(ie)%state%pswet(:,:)=elem(ie)%state%ps(:,:,tl%np1)
        end if
    enddo

    call calc_tot_energy_dynamics(elem,nets,nete,tl%np1,np1_qdp,'dAR')


    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap



    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")
    ! note: time level update for fvm tracers takes place in fvm_mod

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined


    ! ============================================================
    ! Print some diagnostic information
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)
    end if

    if (ntrac>0.and.nsubstep==nsplit.and.nc.ne.fv_nphys) then
      !
      ! fill the fvm halo for mapping in d_p_coupling if
      ! physics grid resolution is different than fvm resolution
      !
      call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm,nhc,1,nlev)
    end if

  end subroutine prim_run_subcycle






  subroutine prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord, rstep)
!
!   Take qsplit dynamics steps and one tracer step
!   for vertically lagrangian option, this subroutine does only the horizontal step
!
!   input:
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, nsplit
    use control_mod, only: statefreq, integration, qsplit, nu_p, test_cfldep, rsplit
    use control_mod, only : tracer_grid_type, TRACER_GRIDTYPE_GLL, tracer_transport_type
    use thread_mod, only : omp_get_thread_num
    use prim_advance_mod, only : prim_advance_exp
    use prim_advection_mod, only : prim_advec_tracers_remap, prim_advec_tracers_fvm, deriv
    use derivative_mod, only : subcell_integration
    use parallel_mod, only : abortmp, iam 
    use reduction_mod, only : parallelmax
    use time_mod,    only : time_at
    use fvm_control_volume_mod, only : fvm_supercycling
    use hybrid_mod, only : set_region_num_threads, config_thread_region
    use dimensions_mod, only: ntrac
  
    type (element_t) , intent(inout)        :: elem(:)

     type(fvm_struct), intent(inout) :: fvm(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                     :: rstep ! vertical remap subcycling step

    type (hybrid_t) :: hybridnew
    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie,t,q,k,i,j,n, n_Q
    integer :: ithr
    integer :: region_num_threads

    real (kind=real_kind)                          :: maxcflx, maxcfly

    real (kind=real_kind) ::  tempdp3d(np,np), x
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)

    real (kind=real_kind) :: dp_np1(np,np)

    dt_q = dt*qsplit
    if (ntrac>0.and.rstep==1) then
      do ie=nets,nete
        elem(ie)%sub_elem_mass_flux=0
      end do
    end if
 
    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0     ! mean vertical mass flux
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      elem(ie)%derived%omega=0
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif

      if (rsplit==0) then
        ! save dp at time t for use in tracers
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k)
#endif
         do k=1,nlev
            elem(ie)%derived%dp(:,:,k)=&
                 ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(:,:,tl%n0)
         enddo
      else
         ! dp at time t:  use floating lagrangian levels:
         elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
      endif
    enddo

    ! ===============
    ! Dynamical Step
    ! ===============
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t.  need to save this
                 ! FV tracers still carry 3 timelevels
                 ! SE tracers only carry 2 timelevels

    call t_startf('prim_advance_exp')
!    ithr   = 0 ! omp_get_thread_num()
!    vybrid = hybrid_create(hybrid%par,ithr)

    call prim_advance_exp(elem, deriv, hvcoord,   &
         hybrid, dt, tl, nets, nete)

    call t_stopf('prim_advance_exp')

    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")

       call t_startf('prim_advance_exp')

       call prim_advance_exp(elem, deriv, hvcoord,   &
            hybrid, dt, tl, nets, nete)

    call t_stopf('prim_advance_exp')

       ! defer final timelevel update until after Q update.
    enddo
#ifdef HOMME_TEST_SUB_ELEMENT_MASS_FLUX
    if (ntrac>0.and.rstep==1) then
      do ie=nets,nete
      do k=1,nlev
        tempdp3d = elem(ie)%state%dp3d(:,:,k,tl%np1) - &
                   elem(ie)%derived%dp(:,:,k) 
        tempmass = subcell_integration(tempdp3d, np, nc, elem(ie)%metdet)
        tempflux = dt_q*elem(ie)%sub_elem_mass_flux(:,:,:,k)
        do i=1,nc
        do j=1,nc
          x = SUM(tempflux(i,j,:))
          if (ABS(tempmass(i,j)).lt.1e-11_real_kind .and. 1e-11_real_kind.lt.ABS(x)) then
            print *,__FILE__,__LINE__,"**********",ie,k,i,j,tempmass(i,j),x
          elseif (1e-5_real_kind.lt.ABS((tempmass(i,j)-x)/tempmass(i,j))) then
            print *,__FILE__,__LINE__,"**********",ie,k,i,j,tempmass(i,j),x,&
                   ABS((tempmass(i,j)-x)/tempmass(i,j))
          endif
        end do
        end do
      end do
      end do
    end if
#endif

    ! current dynamics state variables:
    !    derived%dp              =  dp at start of timestep
    !    derived%vstar           =  velocity at start of tracer timestep
    !    derived%vn0             =  mean horiz. flux:   U*dp
    ! rsplit=0
    !        state%v(:,:,:,np1)      = velocity on reference levels
    !        state%ps(:,:,:,np1)   = ps
    ! rsplit>0
    !        state%v(:,:,:,np1)      = velocity on lagrangian levels 
    !        state%dp3d(:,:,:,np1)   = dp3d
    !


    ! ===============
    ! Tracer Advection.  
    ! in addition, this routine will apply the DSS to:
    !        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
    !        derived%omega           =
    ! Tracers are always vertically lagrangian.  
    ! For rsplit=0: 
    !   if tracer scheme needs v on lagrangian levels it has to vertically interpolate
    !   if tracer scheme needs dp3d, it needs to derive it from ps
    ! ===============
    ! Advect tracers if their count is > 0.  
    ! special case in CAM: if CSLAM tracers are turned on , qsize=1 but this tracer should 
    ! not be advected.  This will be cleaned up when the physgrid is merged into CAM trunk
    ! Currently advecting all species
    if (qsize > 0) then

      call t_startf('prim_advec_tracers_remap')
      region_num_threads = tracer_num_threads
#ifdef _OPENMP
      call omp_set_nested(.true.)
#endif
!JMD      !$OMP PARALLEL NUM_THREADS(region_num_threads), DEFAULT(SHARED), PRIVATE(hybridnew)
!JMD     hybridnew = config_thread_region(hybrid,'tracer')
      call Prim_Advec_Tracers_remap(elem, deriv,hvcoord,flt_advection,hybrid,&
           dt_q,tl,nets,nete)
!JMD      !$OMP END PARALLEL
#ifdef _OPENMP
      call omp_set_nested(.false.)
#endif
      call t_stopf('prim_advec_tracers_remap')


    end if
    !
    ! only run fvm transport every fvm_supercycling rstep
    !
    if (ntrac>0 .and. (mod(rstep,fvm_supercycling) == 0)) then
       !
       ! FVM transport
       !
      call Prim_Advec_Tracers_fvm(elem, fvm, deriv,hvcoord,hybrid,&
           dt_q,tl,nets,nete)

       if(test_cfldep) then
         maxcflx=0.0D0
         maxcfly=0.0D0
         do k=1, nlev

!            maxcflx = parallelmax(fvm(:)%maxcfl(1,k),hybrid)
!            maxcfly = parallelmax(fvm(:)%maxcfl(2,k),hybrid)
           maxcflx = max(maxcflx,parallelmax(fvm(:)%maxcfl(1,k),hybrid))
           maxcfly = max(maxcfly,parallelmax(fvm(:)%maxcfl(2,k),hybrid))
          end do

           if(hybrid%masterthread) then
             write(*,*) "nstep",tl%nstep,"dt_fvm=", dt_q*fvm_supercycling, "maximum over all Level"
             write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly
             write(*,*) "WARNING: maxcfl currently not set if running consistent_se_cslam!"
             print *
           endif
       endif

       !overwrite SE density by fvm(ie)%psc
!        call overwrite_SEdensity(elem,fvm,dt_q,hybrid,nets,nete,tl%np1)
    endif

  end subroutine prim_step


!=======================================================================================================!


  subroutine prim_finalize(hybrid)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize





    subroutine smooth_topo_datasets(phis,sghdyn,sgh30dyn,elem,hybrid,nets,nete)
    use control_mod, only : smooth_phis_numcycle,smooth_sgh_numcycle
    use hybrid_mod, only : hybrid_t
    use bndry_mod, only : bndry_exchangev
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use viscosity_mod, only : biharmonic_wk
    use prim_advance_mod, only : smooth_phis
    use prim_advection_mod, only: deriv
    implicit none

    integer , intent(in) :: nets,nete
    real (kind=real_kind), intent(inout)   :: phis(np,np,nets:nete)
    real (kind=real_kind), intent(inout)   :: sghdyn(np,np,nets:nete)
    real (kind=real_kind), intent(inout)   :: sgh30dyn(np,np,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    ! local
    integer :: ie
    real (kind=real_kind) :: minf

    minf=-9.e9_real_kind
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to PHIS"
    call smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,smooth_phis_numcycle)

    minf=0
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH"
    call smooth_phis(sghdyn,elem,hybrid,deriv,nets,nete,minf,smooth_sgh_numcycle)
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH30"
    call smooth_phis(sgh30dyn,elem,hybrid,deriv,nets,nete,minf,smooth_sgh_numcycle)

    end subroutine smooth_topo_datasets


  subroutine prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)
  use kinds, only : real_kind
  use control_mod, only : initial_total_mass
  use physical_constants, only : g
  use element_mod, only : element_t
  use time_mod, only : timelevel_t 
  use hybvcoord_mod, only : hvcoord_t 
  use hybrid_mod, only : hybrid_t
  use dimensions_mod, only : np
  use global_norms_mod, only : global_integral 

  type (element_t), intent(inout) :: elem(:)
  type (TimeLevel_t), target, intent(in) :: tl
  type (hybrid_t),intent(in)     :: hybrid
  type (hvcoord_t), intent(in)   :: hvcoord
  integer,intent(in)             :: nets,nete
  
  ! local 
  real (kind=real_kind)  :: tmp(np,np,nets:nete)
  real (kind=real_kind)  :: scale,mass0
  integer :: n0,nm1,np1,ie

  if (initial_total_mass == 0) return;
  
  n0=tl%n0
  nm1=tl%nm1
  np1=tl%np1
  
  scale=1/g                                  ! assume code is using Pa
  if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
  ! after scaling, Energy is in J/m**2,  Mass kg/m**2
  
  do ie=nets,nete
     tmp(:,:,ie)=elem(ie)%state%ps(:,:,n0)
  enddo
  mass0 = global_integral(elem, tmp(:,:,nets:nete),hybrid,np,nets,nete)
  mass0 = mass0*scale;  
  
  do ie=nets,nete
     elem(ie)%state%ps(:,:,n0)=elem(ie)%state%ps(:,:,n0)*(initial_total_mass/mass0)
     elem(ie)%state%ps(:,:,np1)=elem(ie)%state%ps(:,:,n0)
     elem(ie)%state%ps(:,:,nm1)=elem(ie)%state%ps(:,:,n0)
  enddo
  if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
     write (*,'(a,e24.15)') "Initializing Total Mass (kg/m^2) = ",initial_total_mass
  endif
  end subroutine prim_set_mass


end module prim_driver_mod
