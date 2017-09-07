#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
  !MODULE FVM_RECONSTRUCTION_MOD--------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
  ! This module contains everything  to do (ONLY) a CUBIC (3rd order) reconstruction  !
  !                                                                                   !
  ! IMPORTANT: the implementation is done for a ncfl > 1, which is not working        !
  !            but it works for ncfl=1                                                !
  !
  ! This module has been recoded for multi-tracer efficiency (May, 2014)
  !
  !-----------------------------------------------------------------------------------!
  module fvm_reconstruction_mod

    use kinds, only                  : int_kind, real_kind
    use coordinate_systems_mod, only : cartesian2D_t,cartesian3D_t
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use parallel_mod, only: abortmp
    implicit none
    private
!    integer, parameter, private:: nh = nhr+(nhe-1) ! = 2 (nhr=2; nhe=1)
    ! = 3 (nhr=2; nhe=2)
  public :: reconstruction, recons_val_cart, extend_panel_interpolate
!reconstruction_gradient,
contains
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE RECONSTRUCTION------------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
  ! DESCRIPTION: controls the cubic (3rd order) reconstructions:                      !
  !                                                                                   !
  ! CALLS: fillhalo_cubic, reconstruction_cubic                                       !
  ! INPUT: fcube    ...  tracer values incl. the halo zone                            !
  !        fvm    ...  structure incl. tracer values aso                            !                                   !
  ! OUTPUT:recons   ...  has the reconstruction coefficients (5) for the 3rd order    !
  !                      reconstruction: dx, dy, dx^2, dy^2, dxdy                     !
  !-----------------------------------------------------------------------------------!
  subroutine reconstruction(fcube,recons,irecons,llimiter,ntrac_in,&
       nc,nhe,nhr,nhc,nht,ns,nh,&
       jx_min,jx_max,jy_min,jy_max,&
       cubeboundary,halo_interp_weight,ibase,&
       spherecentroid,&
       recons_metrics,recons_metrics_integral,&
       rot_matrix,centroid_stretch,&
       vertex_recons_weights,vtx_cart&
       )
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    implicit none
    !
    ! dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc)
    !
    integer, intent(in) :: irecons
    integer, intent(in) :: ntrac_in,nc,nhe,nhr,nhc,nht,ns,nh,cubeboundary
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,ntrac_in), intent(inout) :: fcube
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons,ntrac_in), intent(out)  :: recons
    integer, intent(in) :: jx_min(3), jx_max(3), jy_min(3), jy_max(3)
    integer              , intent(in):: ibase(1-nh:nc+nh,1:nhr,2)
    real (kind=real_kind), intent(in):: halo_interp_weight(1:ns,1-nh:nc+nh,1:nhr,2)
    real (kind=real_kind), intent(in):: spherecentroid(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons-1)
    real (kind=real_kind), intent(in):: recons_metrics(1-nhe:nc+nhe,1-nhe:nc+nhe,3)
    real (kind=real_kind), intent(in):: recons_metrics_integral(1-nhe:nc+nhe,1-nhe:nc+nhe,3)
    integer              , intent(in):: rot_matrix(2,2,1-nhc:nc+nhc,1-nhc:nc+nhc)
    real (kind=real_kind), intent(in):: centroid_stretch(1-nhe:nc+nhe,1-nhe:nc+nhe,1:7)
    real (kind=real_kind), intent(in):: vertex_recons_weights(4,1:irecons-1,1-nhe:nc+nhe,1-nhe:nc+nhe)
    real (kind=real_kind), intent(in):: vtx_cart(1-nhc:nc+nhc,1-nhc:nc+nhc,2,4)

    logical, intent(in) :: llimiter(ntrac_in)

    real (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,3) :: f

    integer :: i,j,ir,in,h,itr
    integer,               dimension(2,3)                              :: jx,jy

    jx(1,1)=jx_min(1); jx(2,1)=jx_max(1)-1
    jx(1,2)=jx_min(2); jx(2,2)=jx_max(2)-1
    jx(1,3)=jx_min(3); jx(2,3)=jx_max(3)-1

    jy(1,1)=jy_min(1); jy(2,1)=jy_max(1)-1
    jy(1,2)=jy_min(2); jy(2,2)=jy_max(2)-1
    jy(1,3)=jy_min(3); jy(2,3)=jy_max(3)-1

    recons=0.0D0
    if (nhe>0) then
      do itr=1,ntrac_in
        call extend_panel_interpolate(nc,nhc,nhr,nht,ns,nh,&
             fcube(:,:,itr),cubeboundary,halo_interp_weight,ibase,f(:,:,1),f(:,:,2:3))
        if (irecons>1) call get_gradients(f(:,:,:),jx,jy,irecons,recons(:,:,:,itr),&
             rot_matrix,centroid_stretch,nc,nht,nhe,nhc)
      end do
    else
      do itr=1,ntrac_in
        call extend_panel_interpolate(nc,nhc,nhr,nht,ns,nh,&
             fcube(:,:,itr),cubeboundary,halo_interp_weight,ibase,f(:,:,1))
        if (irecons>1) call get_gradients(f(:,:,:),jx,jy,irecons,recons(:,:,:,itr),&
             rot_matrix,centroid_stretch,nc,nht,nhe,nhc)
      end do
    end if

    !
    ! fill in non-existent (in physical space) corner values to simplify
    ! logic in limiter code (min/max operation)
    !
    do itr=1,ntrac_in
      if (llimiter(itr)) then
        if (cubeboundary>4) then
          select case(cubeboundary)
          case (nwest)
            do h=1,nhe+1
              fcube(0,nc+h  ,itr) = fcube(1-h,nc  ,itr)
              fcube(1-h,nc+1,itr) = fcube(1  ,nc+h,itr)
            end do
          case (swest)
            do h=1,nhe+1
              fcube(1-h,0,itr) = fcube(1,1-h,itr)
              fcube(0,1-h,itr) = fcube(1-h,1,itr)
            end do
          case (seast)
            do h=1,nhe+1
              fcube(nc+h,0  ,itr) = fcube(nc,1-h,itr)
              fcube(nc+1,1-h,itr) = fcube(nc+h,1,itr)
            end do
          case (neast)
            do h=1,nhe+1
              fcube(nc+h,nc+1,itr) = fcube(nc,nc+h,itr)
              fcube(nc+1,nc+h,itr) = fcube(nc+h,nc,itr)
            end do
          end select
        end if
        call slope_limiter(nhe,nc,nhc,fcube(:,:,itr),jx,jy,irecons,recons(:,:,:,itr),&
             spherecentroid(1-nhe:nc+nhe,1-nhe:nc+nhe,:),&
             recons_metrics,vertex_recons_weights,vtx_cart )
!        call slope_limiter(fvm,fcube(:,:,itr),jx,jy,irecons,recons(:,:,:,itr))
      end if
    end do

    select case (irecons)
    case(1)
      recons(:,:,1,1:ntrac_in) = fcube(1-nhe:nc+nhe,1-nhe:nc+nhe,1:ntrac_in)
    case(3)
      do j=1-nhe,nc+nhe
        do i=1-nhe,nc+nhe
          recons(i,j,1,1:ntrac_in)  = fcube(i,j,1:ntrac_in) &
               - recons(i,j,2,1:ntrac_in)*spherecentroid(i,j,1) &
               - recons(i,j,3,1:ntrac_in)*spherecentroid(i,j,2)
          recons(i,j,2,1:ntrac_in) = recons(i,j,2,1:ntrac_in)
          recons(i,j,3,1:ntrac_in) = recons(i,j,3,1:ntrac_in)
        end do
      end do
    case(6)
      do itr=1,ntrac_in
        do j=1-nhe,nc+nhe
          do i=1-nhe,nc+nhe
!             recons(i,j,1,itr)  = fcube(i,j,itr)!hack first-order
!             recons(i,j,2:6,itr)  = 0.0D0       !hack first-order

!            write(*,*) "yyy ",i,j,itr,"fcube ",fcube(i,j,itr),&
!                 "recons ",recons(i,j,:,itr)," centroid ",spherecentroid(i,j,:), &
!                 "rmetri ",recons_metrics_integral(i,j,:)

            recons(i,j,1,itr)  = fcube(i,j,itr) &
                 - recons(i,j,2,itr)*spherecentroid(i,j,1) &
                 - recons(i,j,3,itr)*spherecentroid(i,j,2) &
                 + recons(i,j,4,itr)*recons_metrics_integral(i,j,1) &
                 + recons(i,j,5,itr)*recons_metrics_integral(i,j,2) &
                 + recons(i,j,6,itr)*recons_metrics_integral(i,j,3)
            recons(i,j,2,itr) = recons(i,j,2,itr)                 &
                 - recons(i,j,4,itr)*2.0D0*spherecentroid(i,j,1) &
                 - recons(i,j,6,itr)      *spherecentroid(i,j,2)
            recons(i,j,3,itr) = recons(i,j,3,itr)                &
                 - recons(i,j,5,itr)*2.0D0*spherecentroid(i,j,2) &
                 - recons(i,j,6,itr)*spherecentroid(i,j,1)
          !
          ! recons(i,j,4:6) already set in get_gradients
          !
          end do
      end do
    end do
    case default
      write(*,*) "irecons out of range in get_ceof", irecons
    end select

    !          recons(a,b,3) * (centroid(a,b,1)**2 - centroid(a,b,3)) + &
    !          recons(a,b,4) * (centroid(a,b,2)**2 - centroid(a,b,4)) + &
    !          recons(a,b,5) * (centroid(a,b,1) * centroid(a,b,2) - centroid(a,b,5)) + &


    !  call debug_halo(fvm,fcubenew,fpanel)
    !  call debug_halo_recons(fvm,recons,recons_trunk)
    !  call print_which_case(fvm)
    !
    !  call debug_halo_neighbor       (fvm,fotherface,fotherpanel)
    !  call debug_halo_neighbor_recons(fvm,recons,recons_trunk)
  end subroutine reconstruction
  !END SUBROUTINE RECONSTRUCTION--------------------------------------------CE-for FVM!

  subroutine get_gradients(f,jx,jy,irecons,gradient,rot_matrix,centroid_stretch,nc,nht,nhe,nhc)
    implicit none
    integer,                                                               intent(in)  :: irecons,nc,nht,nhe,nhc
    real (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,3)        , intent(in)  :: f
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons  ), intent(out) :: gradient
    integer,               dimension(2,3)                                , intent(in)  :: jx,jy
    integer              , dimension(2,2,1-nhc:nc+nhc,1-nhc:nc+nhc)      , intent(in)  :: rot_matrix
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,1:7)      , intent(in)  :: centroid_stretch

    integer :: i,j,in
    real (kind=real_kind), dimension(2) :: g
    real (kind=real_kind) :: tmp,sign

    select case (irecons)
    case(3)
      in=1
      do j=jy(1,in),jy(2,in)
        do i=jx(1,in),jx(2,in)
          !
          ! df/dx: 4-th-order finite difference: (-f(i+2)+8f(i+1)-8f(i-1)+f(i-2))/12dx
          !
          gradient(i,j,2) = -f(i+2,j  ,in)+8.0D0*f(i+1,j  ,in)-8.0D0*f(i-1,j  ,in)+f(i-2,j  ,in)
          gradient(i,j,3) = -f(i  ,j+2,in)+8.0D0*f(i  ,j+1,in)-8.0D0*f(i  ,j-1,in)+f(i  ,j-2,in)
        end do
      end do
      do in=2,3
        do j=jy(1,in),jy(2,in)
          do i=jx(1,in),jx(2,in)
            g(1) = -f(i+2,j  ,in)+8.0D0*f(i+1,j  ,in)-8.0D0*f(i-1,j  ,in)+f(i-2,j  ,in)
            g(2) = -f(i  ,j+2,in)+8.0D0*f(i  ,j+1,in)-8.0D0*f(i  ,j-1,in)+f(i  ,j-2,in)
            gradient(i,j,2:3) = MATMUL(rot_matrix(:,:,i,j),g(:))
          end do
        end do
      end do
      gradient(:,:,2) = centroid_stretch(:,:,1)*gradient(:,:,2)
      gradient(:,:,3) = centroid_stretch(:,:,2)*gradient(:,:,3)
    case (6)
      in=1
      do j=jy(1,in),jy(2,in)
        do i=jx(1,in),jx(2,in)
          !
          ! df/dx: 4-th-order finite difference: (-f(i+2)+8f(i+1)-8f(i-1)+f(i-2))/12dx
          !
          gradient(i,j,2) = -f(i+2,j  ,in)+ 8.0D0*f(i+1,j  ,in)                 - 8.0D0*f(i-1,j  ,in)+f(i-2,j  ,in)
          gradient(i,j,3) = -f(i  ,j+2,in)+ 8.0D0*f(i  ,j+1,in)                 - 8.0D0*f(i  ,j-1,in)+f(i  ,j-2,in)
          !
          ! d2f/dx2:
          !
          gradient(i,j,4) = -f(i+2,j  ,in)+16.0D0*f(i+1,j  ,in)-30.0D0*f(i,j,in)+16.0D0*f(i-1,j  ,in)-f(i-2,j  ,in)
          gradient(i,j,5) = -f(i  ,j+2,in)+16.0D0*f(i  ,j+1,in)-30.0D0*f(i,j,in)+16.0D0*f(i  ,j-1,in)-f(i  ,j-2,in)

          gradient(i,j,6) =  f(i+1,j+1,in)-       f(i+1,j-1,in)                 -       f(i-1,j+1,in)+f(i-1,j-1,in)
        end do
      end do
      do in=2,3
        if (SUM(rot_matrix(:,:,jx(1,in),jy(1,in)))==0) then
          sign=-1
        else
          sign=1
        end if
        do j=jy(1,in),jy(2,in)
          do i=jx(1,in),jx(2,in)
            g(1) = -f(i+2,j  ,in)+8.0D0*f(i+1,j  ,in)-8.0D0*f(i-1,j  ,in)+f(i-2,j  ,in)
            g(2) = -f(i  ,j+2,in)+8.0D0*f(i  ,j+1,in)-8.0D0*f(i  ,j-1,in)+f(i  ,j-2,in)
            gradient(i,j,2:3) = MATMUL(rot_matrix(:,:,i,j),g(:))
            g(1) = -f(i+2,j  ,in)+16.0D0*f(i+1,j  ,in)-30.0D0*f(i,j,in)+16.0D0*f(i-1,j  ,in)-f(i-2,j  ,in)
            g(2) = -f(i  ,j+2,in)+16.0D0*f(i  ,j+1,in)-30.0D0*f(i,j,in)+16.0D0*f(i  ,j-1,in)-f(i  ,j-2,in)
            gradient(i,j,4:5) = MATMUL(ABS(rot_matrix(:,:,i,j)),g(:))

            gradient(i,j,6) =  sign*(f(i+1,j+1,in)-       f(i+1,j-1,in)                 -       f(i-1,j+1,in)+f(i-1,j-1,in))
          end do
        end do
      end do
      gradient(:,:,2) = centroid_stretch(:,:,1)*gradient(:,:,2)
      gradient(:,:,3) = centroid_stretch(:,:,2)*gradient(:,:,3)

      gradient(:,:,4) = centroid_stretch(:,:,3)*gradient(:,:,4)+centroid_stretch(:,:,6)*gradient(:,:,2)
      gradient(:,:,5) = centroid_stretch(:,:,4)*gradient(:,:,5)+centroid_stretch(:,:,7)*gradient(:,:,3)

      gradient(:,:,6) = centroid_stretch(:,:,5)*gradient(:,:,6)
    case default
      call abortmp('ERROR: irecons out of range in fvm_reconstruction_mod')
    end select
  end subroutine get_gradients


  subroutine slope_limiter(nhe,nc,nhc,fcube,jx,jy,irecons,recons,spherecentroid,recons_metrics,&
       vertex_recons_weights,vtx_cart)
    implicit none
    integer                                                                  , intent(in) :: irecons,nhe,nc,nhc
    real (kind=real_kind), dimension(1-nhc:, 1-nhc:),                      intent(inout)  :: fcube
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons)    , intent(inout):: recons
    integer,               dimension(2,3)                                    , intent(in) :: jx,jy
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons-1)    , intent(in) :: spherecentroid
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,3)            , intent(in) :: recons_metrics
    real (kind=real_kind), dimension(4,1:irecons-1,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in) :: vertex_recons_weights
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,2,4)          , intent(in) :: vtx_cart

    real (kind=real_kind)     :: minval_patch,maxval_patch
    real (kind=real_kind)  :: min_phi, minx, maxx, miny,maxy, phi, min_val, max_val

    real (kind=real_kind) :: disc
    real (kind=real_kind) :: min_x,max_x,min_y,max_y
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons-1) :: recons_weights

    real (kind=real_kind)            :: value,max_diff,min_diff,phi_max,phi_min
    real (kind=real_kind)            :: extrema(2), xminmax(2),xmax,yminmax(2),extrema_value(13)

    real(kind=real_kind) :: invtmp  ! temporary to pre-compute inverses
    integer :: itmp1, itmp2         ! temporary index

    integer :: i,j,in,vertex,ir,n

!    real (kind=real_kind), dimension(-1:5) :: diff_value
    real (kind=real_kind), parameter :: threshold = 0.0D0
    select case (irecons)
      !
      ! PLM limiter
      !

    case(3)
      do in=1,3
        do j=jy(1,in),jy(2,in)
          do i=jx(1,in),jx(2,in)
            !     do j=1-nhe,nc+nhe
            !        do i=1-nhe,nc+nhe
            !           if (mask(i,j)) then

             !rck combined min/max and unrolled inner loop
             !minval_patch = MINVAL(fcube(i-1:i+1,j-1:j+1))
             !maxval_patch = MAXVAL(fcube(i-1:i+1,j-1:j+1))
             minval_patch = fcube(i-1,j-1)
             maxval_patch = fcube(i-1,j-1)
             do itmp2=j-1,j+1
                minval_patch = min(minval_patch,fcube(i-1,itmp2),fcube(i,itmp2),fcube(i+1,itmp2))
                maxval_patch = max(maxval_patch,fcube(i-1,itmp2),fcube(i,itmp2),fcube(i+1,itmp2))
             enddo
             min_phi=1.0_real_kind
             !rck restructured loop
            do vertex=1,4
              extrema_value(vertex) = &
                   SUM(recons(i,j,2:irecons)*vertex_recons_weights(vertex,1:irecons-1,i,j))+fcube(i,j)
              call slopelimiter_val(extrema_value(vertex), fcube(i,j),minval_patch, maxval_patch, min_phi)
            end do
           max_val = MAXVAL(extrema_value(1:4))
           min_val = MINVAL(extrema_value(1:4))
!            if (ABS(min_val-fcube(i,j))<1.0D-16.or.ABS(max_val-fcube(i,j))<1.0D-16) then
!               min_phi=0.0D0
!            else
               if (max_val>maxval_patch) then
                  phi = (maxval_patch-fcube(i,j))/(max_val-fcube(i,j))
                  if (phi<min_phi) min_phi=phi
               end if
               if (min_val<minval_patch) then
                  phi = (minval_patch-fcube(i,j))/(min_val-fcube(i,j))
                  if (phi<min_phi) min_phi=phi
               end if
            ! Apply monotone limiter to all reconstruction coefficients
            recons(i,j,2:irecons)=min_phi*recons(i,j,2:irecons)
          end do
        end do
      end do
      !
      ! PPM limiter
      !
    case(6)
       !
       ! default branch
       !
       do in=1,3
          do j=jy(1,in),jy(2,in)
             do i=jx(1,in),jx(2,in)
                !rck combined min/max and unrolled inner loop
                !minval_patch = MINVAL(fcube(i-1:i+1,j-1:j+1))
                !maxval_patch = MAXVAL(fcube(i-1:i+1,j-1:j+1))
                minval_patch = fcube(i-1,j-1)
                maxval_patch = fcube(i-1,j-1)
                do itmp2=j-1,j+1
                   minval_patch = min(minval_patch,fcube(i-1,itmp2),fcube(i,itmp2),fcube(i+1,itmp2))
                   maxval_patch = max(maxval_patch,fcube(i-1,itmp2),fcube(i,itmp2),fcube(i+1,itmp2))
                enddo
                min_phi=1.0_real_kind
                !rck restructured loop

                extrema_value(1:4) = fcube(i,j)
                do itmp1=1,irecons-1
                   do vertex=1,4
                      extrema_value(vertex) = extrema_value(vertex) + &
                           recons(i,j,itmp1+1)*vertex_recons_weights(vertex,itmp1,i,j)
                   enddo
                enddo
                extrema_value(5:13) = extrema_value(1)
                !
                ! coordinate bounds (could be pre-computed!)
                !
                xminmax(1) = MINVAL(vtx_cart(i,j,1,:)); xminmax(2) = MAXVAL(vtx_cart(i,j,1,:));
                yminmax(1) = MINVAL(vtx_cart(i,j,2,:)); yminmax(2) = MAXVAL(vtx_cart(i,j,2,:));

                ! Check if the quadratic is minimized within the element
                ! Extrema in the interior of the element (there might be just one candidate)
                ! DO NOT NEED ABS here, if disc<0 we have a saddle point (no maximum or minimum)
                disc =  4.0D0 * recons(i,j,4) * recons(i,j,5) - recons(i,j,6)**2
                if (abs(disc) > threshold) then
                   extrema(1) = recons(i,j,6) * recons(i,j,3) - 2.0D0 * recons(i,j,5) * recons(i,j,2)
                   extrema(2) = recons(i,j,6) * recons(i,j,2) - 2.0D0 * recons(i,j,4) * recons(i,j,3)

                   disc=1.0D0/disc
                   extrema(1) = extrema(1) * disc + spherecentroid(i,j,1)
                   extrema(2) = extrema(2) * disc + spherecentroid(i,j,2)
                   if ( (extrema(1) - xminmax(1) > -threshold) .and. &    !xmin
                        (extrema(1) - xminmax(2) <  threshold) .and. &    !xmax
                        (extrema(2) - yminmax(1) > -threshold) .and. &    !ymin
                        (extrema(2) - yminmax(2) <  threshold)) then      !ymax
                      call recons_val_cart(fcube(i,j), extrema(1), extrema(2), spherecentroid(i,j,:), &
                           recons_metrics(i,j,:), recons(i,j,:), extrema_value(5))
                   endif
                endif
                !
                ! Check all potential minimizer points along element boundaries
                !
                if (abs(recons(i,j,6)) > threshold) then
                   invtmp = 1.0d0 / (recons(i,j,6) + spherecentroid(i,j,2))
                   do n=1,2
                      ! Left edge, intercept with du/dx = 0
                      extrema(2) = invtmp * (-recons(i,j,2) - 2.0D0 * recons(i,j,4) * (xminmax(n) - spherecentroid(i,j,1)))
                      if ((extrema(2) > yminmax(1)-threshold) .and. (extrema(2) < yminmax(2)+threshold)) then
                         call recons_val_cart(fcube(i,j), xminmax(n), extrema(2), spherecentroid(i,j,:), &
                              recons_metrics(i,j,:), recons(i,j,:), extrema_value(5+n))
                      endif
                   enddo
                   ! Top/bottom edge, intercept with du/dy = 0
                   invtmp = 1.0d0 / recons(i,j,6) + spherecentroid(i,j,1)
                   do n = 1,2
                      extrema(1) = invtmp * (-recons(i,j,3) - 2.0D0 * recons(i,j,5) * (yminmax(n) - spherecentroid(i,j,2)))
                      if ((extrema(1) > xminmax(1)-threshold) .and. (extrema(1) < xminmax(2)+threshold)) then
                         call recons_val_cart(fcube(i,j), extrema(1), yminmax(n),spherecentroid(i,j,:), &
                              recons_metrics(i,j,:), recons(i,j,:), extrema_value(7+n))
                      endif
                   enddo
                endif

                ! Top/bottom edge, y=const., du/dx=0
                if (abs(recons(i,j,4)) > threshold) then
                   invtmp = 1.0d0 / (2.0D0 * recons(i,j,4))! + spherecentroid(i,j,1)
                   do n = 1,2
                      extrema(1) = spherecentroid(i,j,1)+&
                           invtmp * (-recons(i,j,2) - recons(i,j,6) * (yminmax(n) - spherecentroid(i,j,2)))

                      if ((extrema(1) > xminmax(1)-threshold) .and. (extrema(1) < xminmax(2)+threshold)) then
                         call recons_val_cart(fcube(i,j), extrema(1), yminmax(n), spherecentroid(i,j,:),&
                              recons_metrics(i,j,:),recons(i,j,:), extrema_value(9+n))
                      endif
                   enddo
                endif
                ! Left/right edge, x=const., du/dy=0
                if (abs(recons(i,j,5)) > threshold) then
                   invtmp = 1.0d0 / (2.0D0 * recons(i,j,5))
                   do n = 1,2
                      extrema(2) = spherecentroid(i,j,2)+&
                           invtmp * (-recons(i,j,3) - recons(i,j,6) * (xminmax(n) - spherecentroid(i,j,1)))

                      if ((extrema(2)>yminmax(1)-threshold) .and. (extrema(2) < yminmax(2)+threshold)) then
                         call recons_val_cart(fcube(i,j), xminmax(n), extrema(2), spherecentroid(i,j,:), &
                              recons_metrics(i,j,:), recons(i,j,:), extrema_value(11+n))
                      endif
                   enddo
                endif
                !rck - combined min/max calculation and unrolled
                !           max_val = MAXVAL(extrema_value)
                !           min_val = MINVAL(extrema_value)
                max_val = extrema_value(13)
                min_val = extrema_value(13)
                do itmp1 = 1,12,4
                   max_val = max(max_val, extrema_value(itmp1),extrema_value(itmp1+1),extrema_value(itmp1+2),extrema_value(itmp1+3))
                   min_val = min(min_val, extrema_value(itmp1),extrema_value(itmp1+1),extrema_value(itmp1+2),extrema_value(itmp1+3))
                enddo
                !rck

                if (max_val>maxval_patch) then
                   phi = (maxval_patch-fcube(i,j))/(max_val-fcube(i,j))
                   if (phi<min_phi) min_phi=phi
                end if
                if (min_val<minval_patch) then
                   phi = (minval_patch-fcube(i,j))/(min_val-fcube(i,j))
                   if (phi<min_phi) min_phi=phi
                end if
                recons(i,j,2:6)=min_phi*recons(i,j,2:6)
             end do
          end do
       end do
    case default
      call abortmp('ERROR: irecons out of range in fvm_reconstruction_mod')
    end select

  end subroutine slope_limiter



  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE RECONS_VAL_CART-----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
  ! DESCRIPTION: returns the value from the reconstruction (3rd order Taylor polynom) !
  !              at the point (cartx,carty) -> in cube CARTESIAN coordinates          !
  !                                                                                   !
  ! INPUT: fcube  ...  tracer values incl. the halo zone                              !
  !        cartx ...  x cartesian coordinate of the evaluation point                  !
  !        carty ...  y cartesian coordinate of the evaluation point                  !
  !        centroid..  x,y,x^2,y^2,xy                                                 !
  !        recons ...  array of reconstructed coefficients                            !
  ! OUTPUT: value ... evaluation at a given point                                     !
  !-----------------------------------------------------------------------------------!
  SUBROUTINE recons_val_cart(fcube, cartx, carty, centroid, pre_computed_metrics, recons, value)
    IMPLICIT NONE
    REAL(KIND=real_kind), intent(in) :: fcube
    REAL(KIND=real_kind), intent(in) :: cartx, carty
    REAL(KIND=real_kind), dimension(:), intent(in) :: centroid
    REAL(KIND=real_kind), dimension(3), intent(in) :: pre_computed_metrics
    REAL(KIND=real_kind), dimension(:), intent(in) :: recons
    REAL(KIND=real_kind), intent(out) :: value
    real(kind=real_kind) :: dx, dy
    dx = cartx - centroid(1)
    dy = carty - centroid(2)
    ! Evaluate constant order terms
    value = fcube + &
         ! Evaluate linear order terms
         recons(2) * dx + &
         recons(3) * dy + &
         ! Evaluate second order terms
         recons(4) * (pre_computed_metrics(1) + dx*dx) + &
         recons(5) * (pre_computed_metrics(2) + dy*dy) + &
         recons(6) * (pre_computed_metrics(3) + dx*dy)
END SUBROUTINE recons_val_cart


  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE SLOPELIMITER_VAL----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
  ! DESCRIPTION: returns the value from the reconstruction (3rd order Taylor polynom) !
  !              at the point (cartx,carty) -> in cube CARTESIAN coordinates          !
  !                                                                                   !
  ! INPUT: value  ...  point value (calculated here by recons_val_cart)               !
  !        cell_value ...  tracer value (in the cell center) of the cell              !
  !        local_min ...  minmal value in the patch                                   !
  !        local_max ...  maximal value in the patch                                  !
  ! INPUT/OUTPUT: min_phi ... slope limiter, inout because we go through any possible !
  !                           extrema on the cell                                     !
  !-----------------------------------------------------------------------------------!
  subroutine slopelimiter_val(value, cell_value, local_min, local_max, min_phi)
    implicit none
    real (kind=real_kind), intent(in)    :: value, cell_value
    real (kind=real_kind), intent(in)    :: local_min, local_max
    real (kind=real_kind), intent(inout) :: min_phi
    real (kind=real_kind) :: phi

    phi= 0.0D0
    ! Check against the minimum bound on the reconstruction
    if (value - cell_value > 1.0e-12_real_kind * value) then
      phi = (local_max - cell_value) / (value - cell_value)
      if (phi < min_phi) then
        min_phi = phi
      endif
      ! Check against the maximum bound on the reconstruction
    elseif (value - cell_value < -1.0e-12_real_kind * value) then
      phi = (local_min - cell_value) / (value - cell_value)
      if(phi < min_phi) then
        min_phi = phi
      endif
    endif
  end subroutine slopelimiter_val
  !END SUBROUTINE SLOPELIMITER_VAL------------------------------------------CE-for FVM!

  function matmul_w(w,f,ns)
    implicit none
    real (kind=real_kind)                          :: matmul_w
    real (kind=real_kind),dimension(:), intent(in) :: w,f      !dimension(ns)
    integer                           , intent(in) :: ns
    real (kind=real_kind)                          :: tmp
    integer                                        :: k
    matmul_w = 0.0D0
    do k=1,ns
      matmul_w = matmul_w+w(k)*f(k)
    end do
  end function matmul_w

  ! special hard-coded version of the function where ns=3
  ! for performance optimization
!  function matmul_w(w, f)
!    IMPLICIT NONE
!    REAL(KIND=real_kind), dimension(3), intent(in) :: w
!    REAL(KIND=real_kind), dimension(3), intent(in) :: f
!    REAL(KIND=real_kind) :: matmul_w
!    matmul_w = w(1)*f(1) + w(2)*f(2) + w(3)*f(3)
!  end function matmul_w

  subroutine extend_panel_interpolate(nc,nhc,nhr,nht,ns,nh,fcube,cubeboundary,halo_interp_weight,ibase,&
       fpanel,fotherpanel)
    implicit none
    integer, intent(in) :: cubeboundary,nc,nhr,nht,nh,nhc,ns
    real (kind=real_kind),   &
         dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in)          :: fcube

    real (kind=real_kind), intent(in) :: halo_interp_weight(1:ns,1-nh:nc+nh,1:nhr,2)
    integer              , intent(in) :: ibase(1-nh:nc+nh,1:nhr,2)

    real (kind=real_kind)  , dimension(1-nht:nc+nht, 1-nht:nc+nht ), intent(out)           :: fpanel
    real   (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,2), intent(out), optional :: fotherpanel

    integer (kind=int_kind)                                         :: i, j, halo,ibaseref
    real (kind=real_kind), dimension(1:ns,1-nh:nc+nh,1:nhr) :: w
    !
    !  fpanel = 1.0E19 !dbg
    !
    !
    ! Stencil for reconstruction is:
    !
    !     ---------------------
    !     |   |   | i |   |   |
    !     ---------------------
    !     |   | i | i | i |   |
    !     ---------------------
    !     | i | i | R | i | i |
    !     ---------------------
    !     |   | i | i | i |   |
    !     ---------------------
    !     |   |   | i |   |   |
    !     ---------------------
    !
    ! where
    !
    !   "R" is cell for which we whish to do the reconstruction
    !   "i" is the stencil
    !
    !
    ! If one or more point in the stencil is on another panel(s) then we need to interpolate
    ! to a stencil that is an extension of the panel on which R is located
    ! (this is done using one dimensional cubic Lagrange interpolation along panel side)
    !
    ! Example: say that southern most "s" on Figure above is on another panels projection then the stencil becomes
    !
    !
    !   ---------------------------------
    !   |   |   |   |   |   | i |   |   |
    !   ----------------|----------------
    !   |   |   |   |   | i | i | i |   |
    !   ----------------|----------------
    !   |   |   |   | i | i | R | i | i |
    !   ----------------|----------------
    !   |   |   |   |   | i | i | i |   |
    !   ---------------------------------
    !   /   /   /   /   / S /S&i/ S / S /
    !  /---/---/---/---/---/---/---/---/
    ! /   /   /   /   /   /   /   /   /
    !/---/---/---/---/---/---/---/---/
    !
    !
    ! where "S" are the cell average values used for the cubic interpolation (located on the South panel)
    !
    !
    if (cubeboundary==0) then
      fpanel(1-nht:nc+nht,1-nht:nc+nht)=fcube(1-nht:nc+nht,1-nht:nc+nht)
    else if (cubeboundary==west) then
      !                                                       !
      !                                                       ! Case shown below: nhr=2, nhe=1, nht=nhr+nhe
      !                                                       ! (nhr = reconstruction width along x and y)
      !                                                       ! (nhe = max. Courant number)
      !                                                       !
      !                                                       !
      !    Figure below shows the element in question         ! In terms of data structure:
      !    (center element) and the surrounding elements      !
      !    on the element in question's projection            !     * "H" is on same panel average value
      !                                                       !     * "w" is west panel average values that need
      !    Notation: "0" marks the element boundaries         !       to be interpolated to main element
      !                                                       !       projection
      !    Elements to the west are on a different projection !     * "i" is extra halo required by the cubic
      !                                                       !       interpolation
      !     0                                                 !
      !     |0000                                             !
      !     |   |00000                                        !
      !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---x---x---x---x---x---x---x---x---x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   |   | i |   |   |   |   |   |   |   |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | i | H | H | H | H | H |   |   |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | w | H | H | H | H | H | H |   |   |
      !     |0000   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     |   |0000   0   |   |   |   0   |   |   |   0     !     |   | w | w | r | r | r | r | r | H | H |   |
      !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---00000000000000000---x---x---x---x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |\--|   \---0---------------0---------------0     !    -------------0---------------0---------------x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |\--|   \---0---------------0---------------0     !    -------------0---------------0---------------x
      !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |0000   |\--0---------------0---------------0     !    -------------0---------------0---------------x
      !     |   |0000   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
      !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---00000000000000000---x---x---x---x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w | r | r | r | r | r | H | H |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | w | H | H | H | H | H | H |   |   |
      !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | i | H | H | H | H | H |   |   |   |
      !      0000   |\--0---------------0---------------0     !    -------------x---------------x---------------x
      !          0000   0   |   |   |   0   |   |   |   0     !     |   |   | i |   |   |   |   |   |   |   |   |
      !              000000000000000000000000000000000000     !    -x---x---x---x---x---x---x---x---x---x---x---x
      !
      !
      !      -2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8
      !
      !
      !
      ! fill in values (incl. halo) that are on the "main" panels projection
      !
      fpanel(1:nc+nht,1-nht:nc+nht)=fcube(1:nc+nht,1-nht:nc+nht)
      !
      ! fill in values that are on the west panels projection
      !
      w = halo_interp_weight(:,:,:,1)
      do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref=ibase(i,halo,1)
          !           ibaseref = ibase(i,halo,1)
          fpanel(1-halo ,i) = matmul_w(w(:,i,halo),fcube(1-halo ,ibaseref:ibaseref+ns-1),ns)
        end do
      end do

      if (present(fotherpanel)) then
        !
        ! fill in values that are on the west panels projection
        !
        fotherpanel (1-nht:0,1-nht:nc+nht,1)=fcube(1-nht:0,1-nht:nc+nht)
        !
        do halo=1,nhr
          do i=halo-nh,nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)
            !
            ! Exploit symmetry in interpolation weights
            !
            fotherpanel(halo,i,1)     = matmul_w(w(:,i,halo),fcube(halo   ,ibaseref:ibaseref+ns-1),ns)
          end do
        end do
      end if
    else if (cubeboundary==east) then
      !
      ! north part is on different panel
      !
      ! stencil
      !
      ! CN<1 case                                             !
      !                                                       !
      !
      !
      !                                                 0     !
      !                                             0000|     !
      !                                         0000|   |     !
      !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---x---x---x---x---x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   |   |   |   |   |   |   | i |   |   |
      !     0---------------0---------------0--/    |--/|     !     x---------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   |   | H | H | H | H | H | i | i |   |
      !     0---------------0---------------0--/|   |--/|     !     x---------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   |   | H | H | H | H | H | H | e | i |   |
      !     0---------------0---------------0--/|   0000|     !     x---------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   0000|   |     !     |   | H | H | r | r | r | r | r | e | e |   |
      !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---00000000000000000---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     0---------------0---------------0--/|   |--/|     !     x---------------0---------------0---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     0---------------0---------------0--/|   |--/|     !     x---------------0---------------0---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     0---------------0---------------0--/|   0000|     !     x---------------0---------------0---x---x---x-
      !     0   |   |   |   0   |   |   |   0   0000|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
      !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---00000000000000000---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r | r | r | r | r | e | e |   |
      !     0---------------0---------------0--/|   |--/|     !     ----------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   | H | H | H | H | H | H | e | i |   |
      !     0---------------0---------------0--/|   |--/|     !     ----------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   |   |   | H | H | H | H | H | i | i |   |
      !     0---------------0---------------0--/|   0000      !     ----------------x---------------x---x---x---x-
      !     0   |   |   |   0   |   |   |   0   0000          !     |   |   |   |   |   |   |   |   | i |   |   |
      !     000000000000000000000000000000000000              !     x---x---x---x---x---x---x---x---x---x---x---x-
      !
      !
      !      -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7
      !
      fpanel      (1-nht:nc     ,1-nht:nc+nht  )=fcube(1-nht:nc     ,1-nht:nc+nht)
      w = halo_interp_weight(:,:,:,1)
      do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref = ibase(i,halo,1)
          fpanel      (nc+halo   ,i  ) = matmul_w(w(:,i,halo),fcube(nc  +halo,ibaseref:ibaseref+ns-1),ns)
        end do
      end do

      if (present(fotherpanel)) then
        fotherpanel (nc+1 :nc+nht ,1-nht:nc+nht,1)=fcube(nc+1 :nc+nht ,1-nht:nc+nht) !
        do halo=1,nhr
          do i=halo-nh,nc+nh-(halo-1)
            !           ibaseref=ibase(i,halo,1 )
            ibaseref = ibase(i,halo,1)
            fotherpanel (nc+1-halo ,i,1) = matmul_w(w(:,i,halo),fcube(nc+1-halo,ibaseref:ibaseref+ns-1),ns)
          end do
        end do
      end if

    else if (cubeboundary==north) then
      !
      ! north part is on different panel
      !
      ! stencil
      !
      ! CN<1 case
      !                                                      !   x---------------x---------------x---------------x
      !                                                      !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !0---\---\---\---0---\---\---\---0---\---\---\---0     !   x---------------x---------------x---------------x
      ! 0   \   \   \   0   \   \   \   0   \   \   \   0    !   |   | i | i | n | n | n | n | n | n | i | i |   |
      !  0---\---\---\---0---\---\---\---0---\---\---\---0   !   x---------------x---------------x---------------x
      !   0   \   \   \   0   \   \   \   0   \   \   \   0  !   | i | i | n | n | n | n | n | n | n | n | i | i |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r | r | r | r | r | r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   | H | H | H | H | H | H | H | H |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   | H | H | H | H | H | H |   |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---x---x---x---x---x---x---x---x---x
      !
      !    -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8    !    -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8
      !
      ! fill in values that are on the same projection as "main" element
      fpanel      (1-nht:nc+nht ,1-nht:nc)=fcube(1-nht:nc+nht ,1-nht:nc)
      w = halo_interp_weight(:,:,:,1)
      do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref = ibase(i,halo,1)
          fpanel      (i,nc+halo    ) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+halo  ),ns) !north
        end do
      end do
      if (present(fotherpanel)) then
        ! fill in halo for north element
        fotherpanel (1-nht:nc+nht ,nc+1:nc+nht,1)=fcube(1-nht:nc+nht ,nc+1:nc+nht)
        !
        do halo=1,nhr
          do i=halo-nh,nc+nh-(halo-1)
            ibaseref = ibase(i,halo,1)
            fotherpanel (i,nc+1-halo,1) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
          end do
        end do
      end if


    else if (cubeboundary==south) then
      !
      ! south part is on different panel
      !
      ! stencil
      !
      !                                                      !
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---x---x---x---x---x---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   | H | H | H | H | H | H |   |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   | H | H | H | H | H | H | H | H |   |   |
      !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r | r | r | r | r | r | H | H |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
      !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
      !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
      !   0   /   /   /   0   /   /   /   0   /   /   /   0  !   | i | i | s | s | s | s | s | s | s | s | i | i |
      !  0---/---/---/---0---/---/---/---0---/---/---/---0   !   x---------------x---------------x---------------x
      ! 0   /   /   /   0   /   /   /   0   /   /   /   0    !   |   | i | i | s | s | s | s | s | s | i | i |   |
      !0---/---/---/---0---/---/---/---0---/---/---/---0     !   x---------------x---------------x---------------x
      !                                                      !   |   |   |   |   |   |   |   |   |   |   |   |   |
      !
      !     0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9           !     0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
      !
      ! fill in values that are on the same projection as "main" element (marked with "i" in Figure above)
      !
      fpanel      (1-nht:nc+nht,1:nc+nht  )=fcube(1-nht:nc+nht,1:nc+nht)
      w = halo_interp_weight(:,:,:,1)
      do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
          ibaseref=ibase(i,halo,1)!ibase(i,halo,2)
          fpanel      (i,1-halo ) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,1-halo),ns)  !south
        end do
      end do
      if (present(fotherpanel)) then
        fotherpanel (1-nht:nc+nht,1-nht:0 ,1)=fcube(1-nht:nc+nht,1-nht:0 )
        do halo=1,nhr
          do i=halo-nh,nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)!ibase(i,halo,2)
            fotherpanel (i,  halo,1) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,  halo),ns)
          end do
        end do
      end if
    else if (cubeboundary==swest) then
      !
      ! south and west neighboring cells are on different panel
      !
      ! stencil
      !
      !
      ! CN<1 case
      !
      !
      !
      !     |000000000000000000000000000000000000   !   x---x---x---x---x---x---x---x---x---x---x---x---x
      !  0000   0   |   |   |   0   |   |   |   0   !   |   |   |   |   |   |   |   |   |   |   |   |   |
      ! 0   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   |   | w | H | H | H | H | H |   |   |   |
      ! |   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w | H | H | H | H | H | H |   |   |
      ! |   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w | r | r | r | r | r | H | H |   |
      ! |   |000000000000000000000000000000000000   !   x---x---x---x---00000000000000000---x---x---x---x
      ! |0000   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! 0   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
      ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
      ! |   |   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! | -/|   000000000000000000000000000000000   !   x---x---x---x---00000000000000000---x---x---x---x
      ! |/  | 0    /   /   /   0   /   /   /   0    !   |   |   | w |   | s | s | s | s | s | s |   |   |
      ! |   0-----/---/---/---0---/---/---/---0     !   x---------------x---------------x---------------x
      ! | 0      /   /   /   0   /   /   /   0      !   |   |   |   | s | s | s | s | s | s |   |   |   |
      ! 0-------/---/---/---0---/---/---/---0       !   x---------------x---------------x---------------x
      !
      !
      !  -1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |       !   |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
      !
      ! fill in values that are on the same projection as "main" element (marked with "i" in Figure above)
      !
      fpanel(1:nc+nht,1:nc+nht)=fcube(1:nc+nht,1:nc+nht)
      !
      ! fill in west part (marked with "w" on Figure above) and south part (marked with "s")
      !
      w = halo_interp_weight(:,:,:,1)
      do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
          ibaseref=ibase(i,halo,1)!ibase(i,halo,1)
          fpanel(1-halo ,i) = matmul_w(w(:,i,halo),fcube(1-halo ,ibaseref:ibaseref+ns-1),ns) !west
          fpanel(i,1-halo ) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,1-halo) ,ns)  !south
        end do
      end do
      !
      ! corner value
      !
      fpanel(0,0)  =0.25D0*(fpanel(0,1)+fpanel(1,0)+fpanel(-1,0)+fpanel(0,-1))
      !
      ! ****************************************************************
      !
      ! fill halo for reconstruction on south neighbor panel projection
      !
      ! ****************************************************************
      !
      ! On the south panel projection the neighbors are arragened as follow (nwest case):
      !
      !
      ! \
      !  \    p
      !   \
      !    \-----
      !    |
      !  w |  s
      !    |
      !
      !
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   | p 0 p | p | p | p 0 p |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   | w | wp0 p | p | p | p 0 p | p |   |   |
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   |   | w | w | r | r | r | r | r | i | i |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   | w | i | i | i | i | i | i |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   | i | i | i | i | i |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   |   |   |   |   |   |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---
      !
      !
      ! fill values on same panel projection ("r" and "i" on Figure above)
      !
      if (present(fotherpanel)) then
        fotherpanel(1:nc+nht,1-nht:0,1)  = fcube(1:nc+nht,1-nht:0)
        !
        ! compute interpolated cell average values in "p" cells on Figure on above
        !
        w = halo_interp_weight(:,:,:,1)
        do halo=1,nhr
          do i=max(halo-nh,0),nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)
            !
            ! use same weights as interpolation south from main panel (symmetric)
            !
            fotherpanel(i,halo,1)  = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,halo),ns)
          end do
        end do
        !
        ! compute interpolated cell average values in "w" cells on Figure on above
        !
        w = halo_interp_weight(:,:,:,2)
        do halo=1,nhr
          do i=nc+halo-nhr,nc+1
            ibaseref=ibase(i,halo,2)-nc
            !
            ! fotherpanel indexing follows main panel indexing
            ! fcube indexing most be "rotated":
            !
            ! ===============================
            ! |              |              |
            ! |  W      ^    |   S          |
            ! |         |    |              |
            ! |       x |    |              |
            ! |         |    |              |
            ! !              |              |
            ! !   <-----     |              |
            ! !      y       |              |
            ! !              |              |
            ! ===============================
            !
            fotherpanel(1-halo,i-nc,1)  = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,halo),ns)
          end do
        end do
        fotherpanel(0,1,1) = 0.25D0*(fotherpanel(-1,1,1)+fotherpanel(1,1,1)+fotherpanel(0,2,1)+fotherpanel(0,0,1))
        !
        ! ****************************************************************
        !
        ! fill halo for reconstruction on west neighbor panel projection
        !
        ! ****************************************************************
        !
        ! On the west panel projection the neighbors are arragened as follow (seast case):
        !
        !   --------
        !   |      |
        !   |  w   |    p
        !   |      |
        !   -------\
        !           \
        !       s    \
        !
        !
        !
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   | i |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   | i | i | e |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   | i | i | r | e | e |   |   |   |   |   |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   |   | s | s | se| e |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   | s | s |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---
        !
        !
        ! fill values on same panel projection ("r" and "i" on Figure above)
        !
        fotherpanel(1-nht:nc,1:nc+nht,2)  = fcube(1-nht:nc,1:nc+nht)
        !
        ! compute interpolated cell average values in "p" cells on Figure on above
        !
        w = halo_interp_weight(:,:,:,1) ! symmetry
        do halo=1,nhr
          do i=max(halo-nh,0),nc+nh-(halo-1)
            ibaseref=ibase(i,halo,1)
            !
            ! use same weights as interpolation south from main panel (symmetric)
            !
            fotherpanel(halo,i,2)  = matmul_w(w(:,i,halo),fcube(halo,ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        !
        ! compute interpolated cell average values in "s" cells on Figure on above
        !
        w = halo_interp_weight(:,:,:,2)
        do halo=1,nhr
          do i=nc+halo-nhr,nc+1
            ibaseref=ibase(i,halo,2)-nc
            !
            ! fotherpanel indexing follows main panel indexing
            ! fcube indexing most be "rotated":
            !
            ! ===============================
            ! |              |              |
            ! |  W      ^    |   S          |
            ! |         |    |              |
            ! |       x |    |              |
            ! |         |    |              |
            ! !              |              |
            ! !   <-----     |              |
            ! !      y       |              |
            ! !              |              |
            ! ===============================
            !
            fotherpanel(i-nc,1-halo,2)  = matmul_w(w(:,i,halo),fcube(halo,ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        fotherpanel(1,0,2) = 0.25D0*(fotherpanel(0,0,2)+fotherpanel(2,0,2)+fotherpanel(1,-1,2)+fotherpanel(1,1,2))
      end if
    else if (cubeboundary==seast) then
      !
      ! south and east neighboring cells are on different panel
      !
      !
      !
      ! 000000000000000000000000000000000000|
      ! 0   |   |   |   0   |   |   |   0   0000    !   |   |   |   |   |   |   |   |   |   |   |   |   |
      ! 0---------------0---------------0--\|   0   !   x---------------x---------------x---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   |   |   |   | H | H | H | H |   |   |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------x---------------x---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   |   |   | H | H | H | H | H | e |   |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------x---------------x---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r | r | r | r | r | e | e |   |   |
      ! 000000000000000000000000000000000000|   |   !   x---x---x---x---00000000000000000---x---x---x---x
      ! 0   |   |   |   0   |   |   |   0   0000|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 0---------------0---------------0--\|   0   !   x---------------0---------------0---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------0---------------0---------------x
      ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 0---------------0---------------0--\|   |   !   x---------------0---------------0---------------x
      ! 0   |   |   |   0   |   |   |   0   |   |   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
      ! 000000000000000000000000000000000   |\- |   !   x---x---x---x---00000000000000000---x---x---x---x
      !  0   \   \   \   0   \   \   \    0 |  \|   !   |   |   | s | s | s | s | s | s |s/e| e |   |   |
      !   0---\---\---\---0---\---\---\-----0   |   !   x---------------x---------------x---------------x
      !    0   \   \   \   0   \   \   \      0 |   !   |   |   |   | s | s | s | s | s | s |   |   |   |
      !     0---\---\---\---0---\---\---\-------0   !   x---------------x---------------x---------------x
      !
      !
      fpanel       (1-nht:nc,1:nc+nht)=fcube(1-nht:nc,1:nc+nht)
      !
      ! east
      !
      w = halo_interp_weight(:,:,:,1)
      do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
          ibaseref = ibase(i,halo,1)
          fpanel(nc+halo,i) = matmul_w(w(:,i,halo),fcube(nc  +halo,ibaseref:ibaseref+ns-1),ns)
        end do
      end do
      !
      ! south
      !
      w = halo_interp_weight(:,:,:,2)
      do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
          ibaseref = ibase(i,halo,2)
          fpanel(i,1-halo ) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,1-halo),ns)  !south
        end do
      end do
      fpanel(nc+1,0   )=0.25D0*(&
           fpanel(nc+1,1)+fpanel(nc,0)+fpanel(nc+2,0)+fpanel(nc+1,-1))
      !
      ! ****************************************************************
      !
      ! fill halo for reconstruction on south neighbor panel projection
      !
      ! ****************************************************************
      !
      ! On the south panel projection the neighbors are arragened as follow (neast case):
      !
      !
      !             /
      !       P    /
      !           /
      !    ------/
      !    |     |  E
      !    |  S  |
      !    |     |
      !
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   |   |   |   |   |   |   |   |   |
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   |   0   |   |   |   0   |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   |   | n 0 n | n | n | n 0 n |   |   |   |
      ! x---x---x---x---0---x---x---x---0---x---x---x---x
      ! |   |   | n | n 0 n | n | n | n 0 ne| e |   |   |
      ! x---x---x---x---00000000000000000---x---x---x---x
      ! |   | i | i | r | r | r | r | r | e | e |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   | i | i | i | i | i | i | e |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   | i | i | i | i | i |   |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      ! |   |   |   |   |   |   |   |   |   |   |   |   |
      ! x---x---x---x---x---x---x---x---x---x---x---x---x
      !
      !
      !
      ! fill values on same panel projection ("r" and "i" on Figure above)
      !
      if (present(fotherpanel)) then
        fotherpanel(1-nht:nc,1-nht:0,1)  = fcube(1-nht:nc,1-nht:0)
        !
        w = halo_interp_weight(:,:,:,2)
        !
        ! fill in "n" on Figure above
        !
        do halo=1,nhr
          do i=halo-nh,min(nc+nh-(halo-1),nc+1)
            ibaseref = ibase(i,halo,2)
            fotherpanel (i,halo,1) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,  halo),ns)
          end do
        end do
        !
        ! fill in "e" on Figure above
        !
        w = halo_interp_weight(:,:,:,1)
        do halo=1,nhr
          do i=0,nht-halo!nc+nh-(halo-1)
            ibaseref = ibase(i,halo,1)
            !
            ! fother panel follows indexing on main panel
            !
            ! use symmetry for weights (same weights as East from main panel but for south panel
            ! projection the indecies are rotated)
            !
            fotherpanel (nc+halo ,1-i,1) = matmul_w(w(:,i,halo),fcube(nc+ibaseref:nc+ibaseref+ns-1,halo),ns)
          end do
        end do
        fotherpanel(nc+1,1,1) = 0.25D0*(fotherpanel(nc+2,1,1)+fotherpanel(nc,1,1)&
             +fotherpanel(nc+1,2,1)+fotherpanel(nc+1,0,1))

        !
        ! ****************************************************************
        !
        ! fill halo for reconstruction on east neighbor panel projection
        !
        ! ****************************************************************
        !
        ! On the south panel projection the neighbors are arragened as follow (neast case):
        !
        !
        !             |     |
        !         P   |  E  |
        !             |-----|
        !            /
        !           /    S
        !          /
        !
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   | i |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   | w | i | i |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   | w | w | r | i | i |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---0---x---x---x---0---x---x---x---x
        ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
        ! x---x---x---x---00000000000000000---x---x---x---x
        ! |   |   |   |   |   |   | w | ws| s | s |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   | s | s |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        ! |   |   |   |   |   |   |   |   |   |   |   |   |
        ! x---x---x---x---x---x---x---x---x---x---x---x---x
        !
        !
        !
        ! fill values on same panel projection ("r" and "i" on Figure above)
        !
        fotherpanel(nc+1:nc+nht,1:nc+nht,2)  = fcube(nc+1:nc+nht,1:nc+nht)
        !
        !
        ! fill in "w" on Figure above
        !
        w = halo_interp_weight(:,:,:,1)
        do halo=1,nhr
          do i=0,nc+nh-(halo-1)
            ibaseref = ibase(i,halo,1)
            fotherpanel(nc+1-halo,i,2) = matmul_w(w(:,i,halo),fcube(nc+1-halo,ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        !
        ! fill in "s" on Figure above
        !
        w = halo_interp_weight(:,:,:,2)
        do halo=1,nhr
          do i=nc+1-nht+halo,nc+1
            !
            !
            ! !  P           |  E
            ! !              |
            ! !              |
            ! ================
            ! |              |
            ! |  S      |    |  <----- y
            ! |         |    |           ^
            ! |       x |    |           |
            ! |         v    |           |
            ! !              |           |
            ! !    ----->    |           x
            ! !      y       |
            ! !              |
            ! ================
            !
            !
            ! shift (since we are using south weights from main panel interpolation
            !
            ibaseref = ibase(i,halo,2)-nc
            !
            ! fotherpanel index: reverse
            !
            ! fcube index: due to rotation (see Figure above)
            !
            fotherpanel(nc+(nc+1-i),1-halo,2) = matmul_w(w(:,i,halo),fcube(nc+1-halo,ibaseref:ibaseref+ns-1),ns)
          end do
        end do
        fotherpanel(nc,0,2) = 0.25D0*(fotherpanel(nc+1,0,2)+fotherpanel(nc-1,0,2)&
             +fotherpanel(nc,1,2)+fotherpanel(nc,-1,2))
      end if
    else if (cubeboundary==nwest) then
      !
      !
      ! 0-------\---\---\---0---\---\---\---0       !   --------x---------------x---------------x
      ! | 0      \   \   \   0   \   \   \   0      !   |   | n | n | n | n | n | n |   |   |   |
      ! |   0-----\---\---\---0---\---\---\---0     !   --------x---------------x---------------x
      ! |   | 0    \   \   \   0   \   \   \   0    !   | w | a | n | n | n | n | n | n |   |   |
      ! |\  |   000000000000000000000000000000000   !   --------00000000000000000---------------x
      ! | -\|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |\--0---------------0---------------0   !   --------0---------------0---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |\--0---------------0---------------0   !   --------0---------------0---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! 0   |\--0---------------0---------------0   !   --------0---------------0---------------x
      ! |0000   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
      ! |   |000000000000000000000000000000000000   !   --------00000000000000000---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w | r | r | r | r | r | H | H |   |
      ! |   |\--0---------------0---------------0   !   --------x---------------x---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   |   | w | H | H | H | H | H | H |   |   |
      ! |   |\--0---------------0---------------0   !   --------x---------------x---------------x
      ! |\--|   0   |   |   |   0   |   |   |   0   !   |   |   | H | H | H | H | H |   |   |   |
      ! 0   |\--0---------------0---------------0   !   --------x---------------x---------------x
      !  0000   0   |   |   |   0   |   |   |   0   !   |   |   |   |   |   |   |   |   |   |   |
      !      000000000000000000000000000000000000   !   --------x---------------x---------------x
      !
      !
      !
      fpanel(1:nc+nht,1-nht:nc)=fcube(1:nc+nht,1-nht:nc)
      !
      ! west
      !
      w = halo_interp_weight(:,:,:,1)
      do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
          ibaseref=ibase(i,halo,1)
          fpanel(1-halo ,i) = matmul_w(w(:,i,halo),fcube(1-halo ,ibaseref:ibaseref+ns-1),ns)
        end do
      end do
      !
      ! north
      !
      w = halo_interp_weight(:,:,:,2)
      do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref = ibase(i,halo,2)
           fpanel(i,nc+halo) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+halo  ),ns) !north
         end do
       end do
       fpanel(0   ,nc+1)=0.25D0*(&
            fpanel(0,nc)+fpanel(1,nc+1)+fpanel(-1,nc+1)+fpanel(0,nc+2))
       !
       ! ****************************************************************
       !
       ! fill halo for reconstruction on north neighbor panel projection
       !
       ! ****************************************************************
       !
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   |   |   |   |   |   |   |   |   |   |   |
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   |   |   | i | i | i | i | i |   |   |   |
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   |   | w | i | i | i | i | i | i |   |   |
       !x---x---x---x---x---x---x---x---x---x---x---x---x
       !|   |   | w | w | r | r | r | r | r | i | i |   |
       !x---x---x---x---00000000000000000---x---x---x---x
       !|   |   | w | ws0 s | s | s | s 0 s | s |   |   |
       !x---x---x---x---0---x---x---x---0---x---x---x---x
       !|   |   |   | s 0 s | s | s | s 0 s |   |   |   |
       !x---x---x---x---0---x---x---x---0---x---x---x---x
       !|   |   |   |   0   |   |   |   0   |   |   |   |
       !x---x---x---x---0---x---x---x---0---x---x---x---x
       !|   |   |   |   0   |   |   |   0   |   |   |   |
       !x---x---x---x---00000000000000000---x---x---x---x
       !
       !
       ! fill values on same panel projection ("r" and "i" on Figure above)
       !
       if (present(fotherpanel)) then
         fotherpanel(1:nc+nht,nc+1:nc+nht,1)  = fcube(1:nc+nht,nc+1:nc+nht)
         !
         !
         ! fill in "s" on Figure above
         !
         ! (use code from north above)
         !
         w = halo_interp_weight(:,:,:,2)
         do halo=1,nhr
           do i=max(halo-nh,0),nc+nh-(halo-1)
             ibaseref = ibase(i,halo,2)
             fotherpanel(i,nc+1-halo,1) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+1-halo  ),ns)
           end do
         end do
         !
         ! fill in "w" on Figure above
         !
         ! (use code from west above)
         !
         w = halo_interp_weight(:,:,:,1)
         do halo=1,nhr
           do i=nc+1-nht+halo,nc+1
             ibaseref=ibase(i,halo,1)-nc
             fotherpanel(1-halo,nc-(i-(nc+1)),1) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
           end do
         end do
         fotherpanel(0,nc,1)=0.25D0*(&
              fotherpanel(1,nc,1)+fotherpanel(-1,nc,1)+fotherpanel(0,nc+1,1)+fotherpanel(0,nc-1,1))

         !
         ! ****************************************************************
         !
         ! fill halo for reconstruction on west neighbor panel projection
         !
         ! ****************************************************************
         !
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   | n | n |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   | n | n | ne| e |   |   |   |   |   |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   | i | i | r 0 e | e |   |   0   |   |   |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   | i | i | r | e | e |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   | i | i | e |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   | i |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---
         !
         !
         ! fill values on same panel projection ("r" and "i" on Figure above)
         !
         fotherpanel(1-nht:nc,1-nht:nc,2)  = fcube(1-nht:nc,1-nht:nc)
         !
         !
         ! fill in "e" on Figure above
         !
         ! (use code from west above)
         !
         w = halo_interp_weight(:,:,:,1)
         do halo=1,nhr
           do i=halo-nh,min(nc+nh-(halo-1),nc+1)
             ibaseref=ibase(i,halo,1)
             fotherpanel(halo ,i,2) = matmul_w(w(:,i,halo),fcube(halo ,ibaseref:ibaseref+ns-1),ns)
           end do
         end do
         !
         !
         ! fill in "n" on Figure above
         !
         ! (use code from north above)
         !
         w = halo_interp_weight(:,:,:,2)
         do halo=1,nhr
           do i=0,nht-halo
             ibaseref = ibase(i,halo,2)+nc
             fotherpanel(1-i,nc+halo,2) = matmul_w(w(:,i,halo),fcube(halo,ibaseref:ibaseref+ns-1),ns) !north
           end do
         end do
         fotherpanel(1,nc+1,2)=0.25D0*(&
              fotherpanel(2,nc+1,2)+fotherpanel(0,nc+1,2)+fotherpanel(1,nc+2,2)+fotherpanel(1,nc,2))
       end if

     else if (cubeboundary==neast) then
       !
       !
       !     0---/---/---/---0---/---/---/-------0     !   x---------------x---------------x--------
       !    0   /   /   /   0   /   /   /      0 |     !   |   |   |   |   | n | n | n | n | n |   |
       !   0---/---/---/---0---/---/---/-----0   |     !   x---------------x---------------x--------
       !  0   /   /   /   0   /   /   /    0 |   |     !   |   |   |   | n | n | n | n | n | a | e |
       ! 000000000000000000000000000000000   |   |     !   x---------------00000000000000000--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 0---------------0---------------0--/|   |     !   x---------------0---------------0--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 0---------------0---------------0--/|   |     !   x---------------0---------------0--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 0---------------0---------------0--/|   0     !   x---------------0---------------0--------
       ! 0   |   |   |   0   |   |   |   0   0000|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
       ! 000000000000000000000000000000000000|   |     !   x---------------00000000000000000--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H | r | r | r | r | e | e |
       ! 0---------------0---------------0--/|   |     !   x---------------x---------------x--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   |   | H | H | H | H | H | e |   |
       ! 0---------------0---------------0--/|   |     !   x---------------x---------------x--------
       ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   |   |   | H | H | H | H |   |   |
       ! 0---------------0---------------0--/|   0     !   x---------------x---------------x--------
       ! 0   |   |   |   0   |   |   |   0   0000      !   |   |   |   |   |   |   |   |   |   |   |
       ! 000000000000000000000000000000000000          !   x---------------x---------------x--------
       !
       !
       !
       fpanel(1-nht:nc,1-nht:nc)=fcube(1-nht:nc,1-nht:nc)
       !     fotherpanel (nc+1 :nc+nht ,1-nht:nc+nht)=fcube(nc+1 :nc+nht ,1-nht:nc+nht)
       !
       ! east
       !
       w = halo_interp_weight(:,:,:,1)
       do halo=1,nhr
         do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=ibase(i,halo,1 )
           fpanel(nc+halo,i) = matmul_w(w(:,i,halo),fcube(nc  +halo,ibaseref:ibaseref+ns-1),ns)
         end do
       end do
       !
       ! north
       !
       !     w = halo_interp_weight(:,:,:,1)
       do halo=1,nhr
         do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=ibase(i,halo,1)
           fpanel(i,nc+halo) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+halo  ),ns) !north
         end do
       end do
       fpanel(nc+1,nc+1)=0.25D0*(&
            fpanel(nc,nc+1)+fpanel(nc+1,nc)+fpanel(nc+1,nc+2)+fpanel(nc+2,nc+1))
       !
       ! ****************************************************************
       !
       ! fill halo for reconstruction on north neighbor panel projection
       !
       ! ****************************************************************
       !
       ! On the north panel projection the neighbors are arragened as follow (seast case):
       !
       !
       !             |     |
       !             |  N  | E
       !             |-----|
       !                   \
       !                S   \
       !                     \
       !
       !
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   |   |   |   |   |   |   |   |   |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   |   |   | i | i | i | i | i |   |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   |   | i | i | i | i | i | i | e |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       ! |   | i | i | r | r | r | r | r | e | e |   |   |
       ! x---x---x---x---00000000000000000---x---x---x---x
       ! |   |   | s | s 0 s | s | s | s 0 se| e |   |   |
       ! x---x---x---x---0---x---x---x---0---x---x---x---x
       ! |   |   |   | s 0 s | s | s | s 0 s |   |   |   |
       ! x---x---x---x---0---x---x---x---0---x---x---x---x
       ! |   |   |   |   0   |   |   |   0   |   |   |   |
       ! x---x---x---x---0---x---x---x---0---x---x---x---x
       ! |   |   |   |   0   |   |   |   0   |   |   |   |
       ! x---x---x---x---00000000000000000---x---x---x---x
       ! |   |   |   |   |   |   |   |   |   |   |   |   |
       ! x---x---x---x---x---x---x---x---x---x---x---x---x
       !
       !
       ! fill values on same panel projection ("r" and "i" on Figure above)
       !
       if (present(fotherpanel)) then
         fotherpanel(1-nht:nc,nc+1:nc+nht,1)  = fcube(1-nht:nc,nc+1:nc+nht)
         !
         ! fill in "s" on Figure above
         !
         ! (use north case from above and shift/reverse j-index
         !
         w = halo_interp_weight(:,:,:,1)
         do halo=1,nhr
           do i=halo-nh,min(nc+nh-(halo-1),nc+1)
             ibaseref=ibase(i,halo,1)
             fotherpanel (i,nc+1-halo,1) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
           end do
         end do
         !
         ! fill in "e" on Figure above
         !
         w = halo_interp_weight(:,:,:,2)
         do halo=1,nhr
           do i=max(halo-nh,0),nht-halo
             ibaseref=ibase(i,halo,2) +nc
             !
             ! fotherpanel uses indexing of main panel's projection
             ! fcube: rotated indexing
             !
             fotherpanel (nc+halo,nc+i,1) = matmul_w(w(:,i,halo),fcube(ibaseref:ibaseref+ns-1,nc+1-halo),ns)
           end do
         end do
         fotherpanel(nc+1,nc,1)=0.25D0*(&
              fotherpanel(nc+2,nc,1)+fotherpanel(nc,nc,1)+fotherpanel(nc+1,nc+1,1)+fotherpanel(nc+1,nc-1,1))
         !
         ! ****************************************************************
         !
         ! fill halo for reconstruction on east neighbor panel projection
         !
         ! ****************************************************************
         !
         ! On the north panel projection the neighbors are arragened as follow (seast case):
         !
         !
         !           \    N
         !            \
         !             \------
         !             |     |
         !         P   |  E  |
         !             |     |
         !             -------
         !
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   | n | n |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   | w | wn| n | n |   |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---0---x---x---x---0---x---x---x---x
         !|   |   |   |   0   |   | w | w 0 r | i | i |   |
         !x---x---x---x---00000000000000000---x---x---x---x
         !|   |   |   |   |   |   | w | w | r | i | i |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   | w | i | i |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   | i |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---x
         !|   |   |   |   |   |   |   |   |   |   |   |   |
         !x---x---x---x---x---x---x---x---x---x---x---x---
         !
         !
         !
         ! fill values on same panel projection ("r" and "i" on Figure above)
         !
         fotherpanel(nc+1:nc+nht,1-nht:nc,2)  = fcube(nc+1:nc+nht,1-nht:nc)
         !
         ! fill in "w" on Figure above
         !
         ! (use east case from above and shift/reverse j-index
         !
         w = halo_interp_weight(:,:,:,1)
         do halo=1,nhr
           do i=halo-nh,min(nc+nh-(halo-1),nc+1)
             ibaseref=ibase(i,halo,1 )
             fotherpanel(nc+1-halo,i,2) = matmul_w(w(:,i,halo),fcube(nc+1-halo,ibaseref:ibaseref+ns-1),ns)
           end do
         end do
         !
         ! fill in "n" on Figure above
         !
         w = halo_interp_weight(:,:,:,2)
         do halo=1,nhr
           do i=max(halo-nh,0),nht-halo
             ibaseref=ibase(i,halo,2) +nc
             !
             ! fotherpanel uses indexing of main panel's projection
             ! fcube: rotated indexing
             !
             fotherpanel (nc+i,nc+halo,2) = matmul_w(w(:,i,halo),fcube(nc+1-halo,ibaseref:ibaseref+ns-1),ns)
           end do
         end do
         fotherpanel(nc,nc+1,2)=0.25D0*(&
              fotherpanel(nc+1,nc+1,2)+fotherpanel(nc-1,nc+1,2)+fotherpanel(nc,nc+2,2)+fotherpanel(nc,nc,2))
       end if
     end if
   end subroutine extend_panel_interpolate
end module fvm_reconstruction_mod
