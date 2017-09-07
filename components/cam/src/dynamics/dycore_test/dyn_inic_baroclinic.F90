module dyn_inic_baroclinic
  !-----------------------------------------------------------------------
  !
  ! Purpose: Set idealized initial conditions for the Ullrich, Melvin, 
  !          Jablonowski and Staniforth (QJRMS, 2014) baroclinic 
  !          instability test.
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  use physconst, only : rair, cpair, gravit, rearth, pi, omega
  use hycoef,     only: hyai, hybi, hyam, hybm, ps0

  implicit none
  private
  real(r8), parameter :: deg2rad = pi/180.D0
    
  !=======================================================================
  !    Baroclinic wave test case parameters
  !=======================================================================
  real(r8), parameter, private :: Mvap = 0.608D0    ! Ratio of molar mass dry air/water vapor
  real(r8), parameter, private :: psurf_moist = 100000.0D0 !moist surface pressure

  real(r8), parameter, private ::     &
       T0E        = 310.d0     ,      & ! temperature at equatorial surface (K)
       T0P        = 240.d0     ,      & ! temperature at polar surface (K)
       B          = 2.d0       ,      & ! jet half-width parameter
       KK         = 3.d0       ,      & ! jet width parameter
       lapse      = 0.005d0             ! lapse rate parameter
  
  real(r8), parameter, private ::     &
       pertu0     = 0.5d0      ,      & ! SF Perturbation wind velocity (m/s)
       pertr      = 1.d0/6.d0  ,      & ! SF Perturbation radius (Earth radii)
       pertup     = 1.0d0      ,      & ! Exp. perturbation wind velocity (m/s)
       pertexpr   = 0.1d0      ,      & ! Exp. perturbation radius (Earth radii)
       pertlon    = pi/9.d0    ,      & ! Perturbation longitude
       pertlat    = 2.d0*pi/9.d0,     & ! Perturbation latitude
       pertz      = 15000.d0   ,      & ! Perturbation height cap
       dxepsilon  = 1.d-5               ! Small value for numerical derivatives
  
  real(r8), parameter, private ::     &
       moistqlat  = 2.d0*pi/9.d0,     & ! Humidity latitudinal width
       moistqp    = 34000.d0,         & ! Humidity vertical pressure width
       moistq0    = 0.018d0             ! Maximum specific humidity


  integer, parameter  :: deep  = 0! Deep (1) or Shallow (0) test case
  integer, parameter  :: pertt = 0!! 0: exponential, 1: streamfunction
  real(r8), parameter :: bigx  = 1.0  ! factor for a reduced size earth
  integer, parameter  :: moist = 1 ! moist (1) or dry (0) baroclinic wave


  ! Public interface
  public :: bc_wav_dyn_set_init,  terminator_tendency, test_func

contains

  subroutine bc_wav_dyn_set_init(vcoord,latvals, lonvals, U, V, T, PS, PHIS, &
       Q, m_cnst, mask, verbose)
    use constituents,  only: cnst_name

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set baroclinic wave initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer           , intent(in)    :: vcoord     !vcoord=0 moist pressure-based vertical coordinate
                                                    !vcoord=1 dry pressure-based vertical coordinate
                                                    !vcoord=2 height-based vertical coordinate
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
                                                    ! z_k for vccord 1) 
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(inout) :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical           , intent(in)    :: verbose    ! For internal use
    ! Local variables
    logical, allocatable              :: mask_use(:)
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'BC_WAV_DYN_SET_INIT'
    real(r8)                          :: ztop,ptop
    real(r8)                          :: uk,vk,Tk,qk,rhok,zk,pk !mid-level state
    real(r8)                          :: thetav,surface_geo,psurface,eta
    real(r8)                          :: wvp,zdummy,qdry
    logical                           :: lU, lV, lT, lPS, lPHIS, lQ, l3d_vars
    real(r8), allocatable             :: latrad(:), lonrad(:) !lat-lon in radians

    real(r8), allocatable             :: pdry_half(:), pwet_half(:)

    if (vcoord.eq.0.or.vcoord.eq.1) then
      !
      ! pressure-based vertical coordinate
      !
      ptop = hyai(1)*ps0
      if (ptop>1D5) then
        call endrun('For iterate_z_given_pressure to work ptop must be less than 100hPa')        
      end if
      ztop = iterate_z_given_pressure(ptop,.false.,ptop,0.0D0,0.0D0,-1000D0) !Find height of top pressure surface      
    else if (vcoord.eq.2) then
      !
      ! height-based vertical coordinate
      !
!      ztop=
      call endrun('z-based vertical coordinate not coded yet')
    else
      call endrun('vcoord value out of range')
    end if

    if(masterproc .and. verbose) then
      write(iulog,*) 'Model top (in km) is at z= ',ztop/1000.0D0
    end if

   
    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun('cnst_init_default: input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    ncol = size(latvals, 1)
    nlev = -1
    allocate(latrad(ncol))
    allocate(lonrad(ncol))
    latrad = latvals!*deg2rad !convert to radians
    lonrad = lonvals!*deg2rad !convert to radians
    !
    !*******************************
    !
    ! initialize surface pressure
    !
    !*******************************
    !    
    if (present(PS)) then
      if (vcoord==0) then
        where(mask_use)
          PS = psurf_moist
        end where
      else if(vcoord==1) then
        !
        ! compute dry surface pressure (subtract water vapor in coloumn)
        !
        do i=1,ncol
          if (mask_use(i)) then 
            wvp = weight_of_water_vapor_given_z(0.0D0,latrad(i),lonrad(i),ztop)
            ps(i) = psurf_moist-wvp
          end if
        end do
      endif

      if(masterproc .and. verbose) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
    end if
    !
    !*******************************
    !
    ! Initialize PHIS
    !
    !*******************************
    !
    if (present(PHIS)) then
      where(mask_use)
        PHIS = 0.0_r8
      end where
      if(masterproc .and. verbose) then
        write(iulog,*) '          PHIS initialized by "',subname,'"'
      end if
    end if
    !
    !*******************************
    !
    ! Initialize 3D vars
    !
    !
    !*******************************
    !
    lu = present(U); lv = present(V); lT = present(T); lq = present(Q);    
    l3d_vars = lu.or.lv.or.lt.or.lq
    nlev = -1
    if (l3d_vars) then
      if (lu) nlev = size(U, 2)
      if (lv) nlev = size(V, 2)
      if (lt) nlev = size(T, 2)
      if (lq) nlev = size(Q, 2)
      if (lq.and.vcoord==1) then
        allocate(pdry_half(nlev+1))
        allocate(pwet_half(nlev+1))
      end if
      do i=1,ncol
        if (mask_use(i)) then 
          if (vcoord==0) then
            psurface = psurf_moist
            wvp = -99
          else if (vcoord==1) then
            !
            ! convert surface pressure to dry
            !
            wvp = weight_of_water_vapor_given_z(0.0D0,latrad(i),lonrad(i),ztop)
            psurface = psurf_moist-wvp            
          end if
          
          do k=1,nlev
            pk =  hyam(k)*ps0 + hybm(k)*psurface
            call baroclinic_wave_test(moist,pk,ptop,zk,uk,vk,tk,thetav,&
                 surface_geo,rhok,qk,&
                 vcoord==1,latrad(i),lonrad(i),ztop)
            if (lt) T(i,k)   = tk
            if (lu) U(i,k)   = uk
            if (lv) V(i,k)   = vk
            if (lq) Q(i,k,1) = qk 
          end do
          if (lq.and.vcoord==1.and.moist.ne.0) then
            !
            ! for dry pressure vertical coordinate
            !
            do k=1,nlev+1
              pdry_half(k) =  hyai(k)*ps0 + hybi(k)*psurf_moist!psurface
              !Find height of pressure surface
              zdummy = iterate_z_given_pressure(pdry_half(k),.true.,ptop,latrad(i),lonrad(i),ztop) 
              pwet_half(k) = moist_pressure_given_z(zdummy,latrad(i),lonrad(i))
            end do
            do k=1,nlev
              qdry =((pwet_half(k+1)-pwet_half(k))/(pdry_half(k+1)-pdry_half(k)))-1.0D0
              !
              ! CAM expects water vapor mixing ratio to be wet - convert from dry to wet:
              !
              Q(i,k,1) = qdry/(1D0+qdry)
              Q(i,k,1) = max(1.0D-12,Q(i,k,1))
            end do
          end if
        end if
      end do
      if(lu.and.masterproc.and. verbose)  write(iulog,*) '          U initialized by "',subname,'"'
      if(lv.and.masterproc.and. verbose)  write(iulog,*) '          V initialized by "',subname,'"'
      if(lt.and.masterproc.and. verbose)  write(iulog,*) '          T initialized by "',subname,'"'
      if(lq.and.masterproc.and. verbose)  write(iulog,*) &
           '          ', trim(cnst_name(m_cnst(1))), ' initialized by "',subname,'"'
    end if
    
    if (lq) then
      ncnst = size(m_cnst, 1)
      if (vcoord==0.or.vcoord==1) then
        do m = 2, ncnst
          do k = 1, nlev
            do i=1,ncol
              if (mask_use(i)) then              
                Q(i,k,m_cnst(m)) = test_func(latrad(i),lonrad(i), k, m)
              end if
            end do
          end do
          if(masterproc .and. verbose) then
            write(iulog,*) '          ', trim(cnst_name(m_cnst(m))), ' initialized by "',subname,'"'
          end if
        end do
      end if
    end if

    deallocate(mask_use)

  end subroutine bc_wav_dyn_set_init

  !-----------------------------------------------------------------------
  !  SUBROUTINE baroclinic_wave_sample(
  !    deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)
  !
  !  Options:
  !     deep    deep atmosphere (1 = yes or 0 = no)
  !    moist    include moisture (1 = yes or 0 = no)
  !    pertt    type of perturbation (0 = exponential, 1 = stream function)
  !        X    Earth scaling factor1
  !
  !  Given a point specified by: 
  !      lon    longitude (radians) 
  !      lat    latitude (radians) 
  !      p/z    pressure (Pa) / height (m)
  !  zcoords    1 if z is specified, 0 if p is specified
  !
  !  the functions will return:
  !        p    pressure if z is specified and zcoords = 1 (Pa)
  !        u    zonal wind (m s^-1)
  !        v    meridional wind (m s^-1)
  !        t    temperature (K)
  !   thetav    virtual potential temperature (K)
  !     phis    surface geopotential (m^2 s^-2)
  !       ps    surface pressure (Pa)
  !      rho    density (kj m^-3)
  !        q    water vapor mixing ratio (kg/kg)
  !
  !
  !  Author: Paul Ullrich
  !          University of California, Davis
  !          Email: paullrich@ucdavis.edu
  !
  !-----------------------------------------------------------------------


  SUBROUTINE baroclinic_wave_test(moist,p,ptop,z,u,v,temp,thetav,phis,rho,q,&
       ldry_mass_vertical_coordinates,lat,lon,ztop)
    IMPLICIT NONE
    
    !-----------------------------------------------------------------------
    !     input/output params parameters at given location
    !-----------------------------------------------------------------------
    integer, INTENT(IN)  :: &
         moist        ! Moist (1) or Dry (0) test case
    
    
    real(r8), INTENT(IN) :: &
         p            ,&! Pressure at the full model level (Pa)
         ptop         ,&!
         lat          ,&! latitude
         lon          ,&! longitude
         ztop           ! model top height

    logical, intent(in) :: ldry_mass_vertical_coordinates
    
    real(r8), INTENT(OUT) :: &
         u,          & ! Zonal wind (m s^-1)
         v,          & ! Meridional wind (m s^-1)
         temp,       & ! Temperature (K)
         thetav,     & ! Virtual potential temperature (K)
         phis,       & ! Surface Geopotential (m^2 s^-2)
!         ps,         & ! Surface Pressure (Pa)
         rho,        & ! density (kg m^-3)
         q,          & ! water vapor mixing ratio (kg/kg)
         z             ! Altitude (m)


    z = iterate_z_given_pressure(p,ldry_mass_vertical_coordinates,ptop,lat,lon,ztop) !Find height of pressure surface
    call uv_given_z(z,u,v,lat,lon)
    temp = Tv_given_z(z,lat,lon)
    phis = 0.d0
    if (moist .eq. 1) then
       q = qv_given_moist_pressure(moist_pressure_given_z(z,lat,lon),lat,lon)       
    else
       q = 0.d0                  ! dry
    end if
    !
    ! Convert virtual temperature to temperature
    !
    temp = temp / (1.d0 + Mvap * q)
    rho = p / (Rair * temp * (1.d0 + 0.61d0 * q))
    thetav = temp * (1.d0 + 0.61d0 * q) * (psurf_moist / p)**(Rair / cpair)
!    if (ldry_mass_vertical_coordinates) then
!       q=q/(1-q)! CAM expects water vapor to be 'wet' mixing ratio so do not convert to dry
!    end if
  END SUBROUTINE baroclinic_wave_test


  real(r8) FUNCTION iterate_z_given_pressure(p,ldry_mass_vertical_coordinates,ptop,lat,lon,ztop)
    implicit none
    real(r8), INTENT(IN)  :: &
         p,              &! Pressure (Pa)
         ptop           ,&! Pressure (Pa)
         lat            ,&! latitude
         lon            ,&! longitude
         ztop

    logical, INTENT(IN)  :: ldry_mass_vertical_coordinates
    
    integer :: ix
    
    real(r8) :: z0, z1, z2
    real(r8) :: p0, p1, p2
    z0 = 0.d0
    z1 = 10000.d0

    if (ldry_mass_vertical_coordinates) then
       p0 = weight_of_dry_air_given_z(z0,ptop,lat,lon,ztop)
       p1 = weight_of_dry_air_given_z(z1,ptop,lat,lon,ztop)
    else
       p0 =  moist_pressure_given_z(z0,lat,lon)
       p1 =  moist_pressure_given_z(z1,lat,lon)
    endif

    DO ix = 1, 100
       z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
       if (ldry_mass_vertical_coordinates) then
          p2 = weight_of_dry_air_given_z(z2,ptop,lat,lon,ztop)
       else
          p2 = moist_pressure_given_z(z2,lat,lon)   
       end if
       
       IF (ABS((p2 - p)/p) .lt. 1.0d-13) THEN
          EXIT
       END IF
       
       z0 = z1
       p0 = p1
       
       z1 = z2
       p1 = p2
    END DO
    if (ix==101) then
      call endrun('iteration did not converge in iterate_z_given_pressure')
    end if    
    iterate_z_given_pressure = z2
  END FUNCTION iterate_z_given_pressure

  real(r8) FUNCTION moist_pressure_given_z(z,lat,lon)
    IMPLICIT NONE
    real(r8), INTENT(IN) :: z,lat,lon
    real(r8) :: aref, omegaref
    real(r8) :: T0, constA, constB, constC, constH, scaledZ
    real(r8) :: tau1, tau2, inttau1, inttau2
    real(r8) :: rratio, inttermT,pwet,wvp
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = rearth / bigX
    omegaref = omega * bigX
    
    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (KK + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rair * T0 / gravit
    
    scaledZ = z / (B * constH)
    
    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    
    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
         + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)
    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if
    
    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**KK &
         - KK / (KK + 2.d0) * (rratio * cos(lat))**(KK + 2.d0)

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    moist_pressure_given_z = psurf_moist * exp(- gravit / Rair * (inttau1 - inttau2 * inttermT))
  END FUNCTION moist_pressure_given_z

  real(r8) FUNCTION Tv_given_z(z,lat,lon)
    IMPLICIT NONE
    real(r8), INTENT(IN) :: z, lat, lon
    real(r8) :: aref, omegaref
    real(r8) :: T0, constA, constB, constC, constH, scaledZ
    real(r8) :: tau1, tau2, inttau1, inttau2
    real(r8) :: rratio, inttermT
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = rearth / bigX
    omegaref = omega * bigX
    
    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (KK + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rair * T0 / gravit
    
    scaledZ = z / (B * constH)
    
    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    
    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
         + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)
    
    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if
    
    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**KK &
         - KK / (KK + 2.d0) * (rratio * cos(lat))**(KK + 2.d0)
    
    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    Tv_given_z = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))
  END FUNCTION Tv_given_z

  SUBROUTINE uv_given_z(z,u,v,lat,lon)
    IMPLICIT NONE
    real(r8), INTENT(IN)  :: z, lat, lon
    real(r8), INTENT(OUT) :: u,v
    real(r8) :: aref, omegaref
    real(r8) :: T0, constH, constC, scaledZ, inttau2, rratio
    real(r8) :: inttermU, bigU, rcoslat, omegarcoslat
    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = rearth / bigx
    omegaref = omega * bigx
    
    T0 = 0.5d0 * (T0E + T0P)
    
    constH = Rair * T0 / gravit
    
    constC = 0.5d0 * (KK + 2.d0) * (T0E - T0P) / (T0E * T0P)
    
    scaledZ = z / (B * constH)
    
    inttau2 = constC * z * exp(- scaledZ**2)
    
    ! radius ratio
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if
    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    inttermU = (rratio * cos(lat))**(KK - 1.d0) - (rratio * cos(lat))**(KK + 1.d0)
    bigU = gravit / aref * KK * inttau2 * inttermU * Tv_given_z(z,lat,lon)
    if (deep .eq. 0) then
       rcoslat = aref * cos(lat)
    else
       rcoslat = (z + aref) * cos(lat)
    end if
    
    omegarcoslat = omegaref * rcoslat
    
    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.d0
    
    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------
!    if (.false.) then !xxxx
    ! Exponential type
    if (pertt .eq. 0) then
       u = u + evaluate_exponential(z,lat,lon)
       
       ! Stream function type
    elseif (pertt .eq. 1) then
       u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
            ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
            - evaluate_streamfunction(lon, lat - dxepsilon, z))
       
       v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
            ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
            - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if
!    endif!xxx
  END SUBROUTINE uv_given_z

  !-----------------------------------------------------------------------
  !    Exponential perturbation function
  !-----------------------------------------------------------------------
  real(r8) FUNCTION evaluate_exponential(z,lat,lon)
    real(r8), INTENT(IN)  :: &
         z             ,&! Altitude (meters)
         lat,lon
    
    real(r8) :: greatcircler, perttaper
    
    ! Great circle distance
    greatcircler = 1.d0 / pertexpr &
         * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))
    
    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
       perttaper = 0.d0
    end if
    
    ! Zonal velocity perturbation
    if (greatcircler < 1.d0) then
       evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
       evaluate_exponential = 0.d0
    end if
    
  END FUNCTION evaluate_exponential
  
  !-----------------------------------------------------------------------
  !    Stream function perturbation function
  !-----------------------------------------------------------------------
  real(r8) FUNCTION evaluate_streamfunction(z,lon_local,lat_local)
    
    real(r8), INTENT(IN)  :: &
         lon_local, lat_local,&
         z             ! Altitude (meters)
    
    real(r8) :: greatcircler, perttaper, cospert
    
    ! Great circle distance
    greatcircler = 1.d0 / pertr &
         * acos(sin(pertlat) * sin(lat_local) + cos(pertlat) * cos(lat_local) * cos(lon_local - pertlon))
    
    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
       perttaper = 0.d0
    end if
    
    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1.d0) then
       cospert = cos(0.5d0 * pi * greatcircler)
    else
       cospert = 0.d0
    end if
    
    evaluate_streamfunction = &
         (- pertu0 * pertr * perttaper * cospert**4)
    
  END FUNCTION evaluate_streamfunction

  real(r8) FUNCTION qv_given_moist_pressure(pwet,lat,lon)
    implicit none
    real(r8), INTENT(IN)  :: pwet, lat, lon

    real(r8)  :: eta
    if (moist==0) then
      qv_given_moist_pressure = 0.0D0
    else
      eta = pwet/psurf_moist 
      if (eta.gt.0.1D0) then  ! intialize q if p > 100 hPa
        qv_given_moist_pressure = moistq0 * exp(- (lat/moistqlat)**4)          & 
             * exp(- ((eta-1.d0)*psurf_moist/moistqp)**2)
      else
        qv_given_moist_pressure = 1.d-12              ! above 100 hPa set q to 1e-12 to avoid supersaturation
      endif
    end if
  END FUNCTION qv_given_moist_pressure

  real(r8) FUNCTION weight_of_water_vapor_given_z(z,lat, lon,ztop)
    implicit none
    real(r8), INTENT(IN)  :: z,lat, lon, ztop
    real (r8)  :: dx,xm,xr,gaussw(10),gaussx(10),integral, tmp1, tmp2
    real(r8)   :: temp, rho, qv, pressure, z1, z2, Tv,pwet, ztmp
    integer   :: jgw
    SAVE gaussw,gaussx    
    DATA gaussw/0.1527533871307258D0,0.1491729864726037D0,0.1420961093183820D0,0.1316886384491766D0,0.1181945319615184D0,&
         0.1019301198172404D0,0.0832767415767048D0,0.0626720483341091D0,0.0406014298003869D0,0.0176140071391521D0/
    DATA gaussx/0.0765265211334973D0,0.2277858511416451D0,0.3737060887154195D0,0.5108670019508271D0,0.6360536807265150D0,&
         0.7463319064601508D0,0.8391169718222188D0,0.9122344282513259D0,0.9639719272779138D0,0.9931285991850949D0/
    
    if (moist==0) then
      !
      ! dry case
      !
      weight_of_water_vapor_given_z = 0.0D0
    else
      z1=z
      z2=ztop
      xm=0.5D0*(z1+z2)
      xr=0.5D0*(z2-z1)
      integral=0 
      do jgw=1,10 
        dx=xr*gaussx(jgw)
        ztmp=xm+dx
        pwet = moist_pressure_given_z(ztmp,lat,lon); qv= qv_given_moist_pressure(pwet,lat,lon);Tv= Tv_given_z(ztmp,lat,lon)
        tmp1=gravit*pwet*qv/(Rair*Tv)
        
        ztmp=xm-dx
        pwet = moist_pressure_given_z(ztmp,lat,lon); qv= qv_given_moist_pressure(pwet,lat,lon);Tv= Tv_given_z(ztmp,lat,lon)
        tmp2=gravit*pwet*qv/(Rair*Tv)
        integral=integral+gaussw(jgw)*(tmp1+tmp2)
      enddo
      integral=xr*integral    ! Scale the answer to the range of integration.    
      
      weight_of_water_vapor_given_z = integral
    end if
  end FUNCTION weight_of_water_vapor_given_z


  real(r8) FUNCTION weight_of_dry_air_given_z(z,ptop,lat,lon,ztop)
    implicit none
    real (r8), INTENT(IN)  :: z,ptop, lat, lon, ztop
    real (r8)  :: dx,xm,xr,gaussw(10),gaussx(10),integral, tmp1, tmp2
    real(r8)   :: temp, rho, qv, pressure, z1, z2, Tv,pwet, ztmp
    integer    :: jgw
    SAVE gaussw,gaussx    
    DATA gaussw/0.1527533871307258D0,0.1491729864726037D0,0.1420961093183820D0,0.1316886384491766D0,0.1181945319615184D0,&
         0.1019301198172404D0,0.0832767415767048D0,0.0626720483341091D0,0.0406014298003869D0,0.0176140071391521D0/
    DATA gaussx/0.0765265211334973D0,0.2277858511416451D0,0.3737060887154195D0,0.5108670019508271D0,0.6360536807265150D0,&
         0.7463319064601508D0,0.8391169718222188D0,0.9122344282513259D0,0.9639719272779138D0,0.9931285991850949D0/
    
    z1=z
    z2=ztop
    xm=0.5*(z1+z2)
    xr=0.5*(z2-z1)
    integral=0 
    do jgw=1,10 
       dx=xr*gaussx(jgw)
       ztmp=xm+dx
       pwet = moist_pressure_given_z(ztmp,lat,lon); qv= qv_given_moist_pressure(pwet,lat,lon);Tv= Tv_given_z(ztmp,lat,lon)
       tmp1=gravit*pwet*(1-qv)/(Rair*Tv)
       
       ztmp=xm-dx
       pwet = moist_pressure_given_z(ztmp,lat,lon); qv= qv_given_moist_pressure(pwet,lat,lon);Tv= Tv_given_z(ztmp,lat,lon)
       tmp2=gravit*pwet*(1-qv)/(Rair*Tv)
       integral=integral+gaussw(jgw)*(tmp1+tmp2)
    enddo
    integral=xr*integral    ! Scale the answer to the range of integration.    
    
    weight_of_dry_air_given_z = integral+ptop
  end FUNCTION weight_of_dry_air_given_z

  ! A simple analytic functions
  ! Test functions
  ! (0) fout = 1D-6
  ! (1) fout = 2 + theta
  ! (2) fout = 2 + cos(theta)
  ! (3) fout = 2 + cos(lambda)
  ! (4) fout = 2 + cos^2(theta)*cos(8*lambda)
  ! (5) fout = 2 + cos^2(theta)*cos(16*lambda)
  ! (6) slotted-cylinder
  ! (7) Cl for terminator chemistry test (Lauritzen et al, 2015)
  ! (8) Cl2 for terminator chemistry test (Lauritzen et al, 2015)
  ! All functions multiplied by factor (default 1.0)
  function test_func(lat, lon, k, funcnum) result(fout)
    use shr_sys_mod,     only: shr_sys_flush
    use physical_constants, only : DD_PI

    real(r8),           intent(in)  :: lon
    real(r8),           intent(in)  :: lat
    integer ,           intent(in)  :: k
    integer,            intent(in)  :: funcnum
    real(r8)                        :: fout
    real(r8)                        :: lon1,lat1,R0,Rg1,Rg2,lon2,lat2,cl,cl2
    real(r8)                        :: eta_c

    real(r8)  :: radius                 = 10.d0 ! radius of the perturbation                                                              
    real(r8)  :: perturb_lon = 20.d0      ! longitudinal position, 20E                                                              
    real(r8)  :: perturb_lat  = 40.d0     ! latitudinal position, 40N
    real(r8)  :: cos_tmp, sin_tmp, eta
    !
    ! For Rossby-Haurwitz velocity field
    !
    real(r8) :: m_rossby_haurwitz,umax

    select case(funcnum)
      !
      ! terminator chemistry initial conditions (Lauritzen et al., 2015, GMD)
      !      
    case(2)
      call initial_value_terminator_chemistry( lat, lon, cl, cl2)
      fout = cl
    case(3)
      call initial_value_terminator_chemistry( lat, lon, cl, cl2)
      fout = cl2
    case(4)     
      !
      !   Non-smooth scalar field (slotted cylinder)
      !
      R0=0.5_r8
      lon1=4.0_r8*DD_PI/5.0_r8
      lat1=0.0_r8
      Rg1 = acos(sin(lat1)*sin(lat)+cos(lat1)*cos(lat)*cos(lon-lon1))
      lon2=6.0_r8*DD_PI/5.0_r8
      lat2=0.0_r8
      Rg2 = acos(sin(lat2)*sin(lat)+cos(lat2)*cos(lat)*cos(lon-lon2))
      
      if ((Rg1 .le. R0) .AND. (abs(lon-lon1).ge. R0/6)) then
        fout = 2.0_r8
      elseif ((Rg2 .le. R0) .AND. (abs(lon-lon2).ge. R0/6)) then
        fout = 2.0_r8
      elseif ((Rg1 .le. R0) .AND. (abs(lon-lon1) < R0/6) &
           .AND. (lat-lat1 < -5.0_r8*R0/12.0_r8)) then
        fout = 2.0_r8
      elseif ((Rg2 .le. R0) .AND. (abs(lon-lon2) < R0/6) &
           .AND. (lat-lat2 > 5.0_r8*R0/12.0_r8)) then
        fout = 2.0_r8
      else
        fout = 1.0_r8
      endif
    case(5)
      !
      ! Smooth Gaussian "ball"
      !
      R0    = 10.d0           ! radius of the perturbation                                                              
      lon1  = 20.d0*deg2rad   ! longitudinal position, 20E                                                              
      lat1  = 40.d0 *deg2rad  ! latitudinal position, 40N
      eta_c = 0.6D0
      sin_tmp = SIN(lat1)*SIN(lat)
      cos_tmp = COS(lat1)*COS(lat)
      Rg1 = ACOS( sin_tmp + cos_tmp*COS(lon-lon1) )    ! great circle distance                                                      
      eta =  (hyam(k)*ps0 + hybm(k)*psurf_moist)/psurf_moist
      fout = EXP(- ((Rg1*R0)**2 + ((eta-eta_c)/0.1D0)**2))
      IF (ABS(fout)<1.0E-8) fout = 0.0D0
    case(6)
      !
      !
      !
      fout = 0.5D0 * ( tanh( 3.D0*abs(lat)-pi ) + 1.D0)      
    case(7)
      fout = 1.0D-8      
    case(8)
      !
      ! approximately Y^2_2 spherical harmonic
      !
      fout = 0.5_r8 + 0.5D0*(cos(lat)*cos(lat)*cos(2.0D0*lon))
    case(9)
      !
      ! approximately Y32_16 spherical harmonic
      !
      fout = 0.5_r8 + 0.5D0*(cos(16*lon)*(sin(2D0*lat)**16))
    case(10)
      fout = 2.0_r8 + lat
    case(11)
      fout = 2.0_r8 + cos(lon)
    case(12)
      !
      ! Idealized static Rossby-Haurwitz vector field: zonal component
      !
      ! D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, P. N. Swarztrauber, A standard test set for
      ! numerical approximations to the shallow water equations in spherical geometry, J. Comput. Phys.
      ! 102 (1) (1992) 211–224.
      !
      !      
      umax=50.0D0
      m_rossby_haurwitz = 8D0
      fout = umax*cos(lat)+umax*cos(lat)**(m_rossby_haurwitz-1)*&
           (m_rossby_haurwitz*sin(lat)**2-cos(lat)**2)*cos(m_rossby_haurwitz*lon)
    case(13)
      !
      ! Idealized static Rossby-Haurwitz vector field: meridional component
      !     
      umax=50.0D0
      m_rossby_haurwitz = 8D0
      fout = -umax*m_rossby_haurwitz*cos(lat)**(m_rossby_haurwitz-1)*sin(lat)*sin(m_rossby_haurwitz*lon)
    case default
      call endrun("Illegal funcnum_arg in test_func")
    end select
  end function test_func
  
  subroutine initial_value_terminator_chemistry( lat, lon, x, x2 )
    real(r8), intent(in)  :: lat, lon  ! latitude and longitude, radians                                                           
    real(r8), intent(out) :: x, x2   ! molar mixing ratio of x and x2    
    real(r8) :: r, det  ! useful algebraic forms                                                                                   
    real(r8) :: k1, k2  ! reaction rates                                                                                           

    real(r8), parameter :: xt_constant = 4.D-6
    real(r8), parameter :: k1_lat_center =   20.d0*pi/180.0D0
    real(r8), parameter :: k1_lon_center =  300.d0*pi/180.0D0
    !
    !  Solar photolysis rate and recombination rate (for terminator chemistry)
    !
    k1 = 1.0D0*max(0.d0,sin(lat)*sin(k1_lat_center) + cos(lat)*cos(k1_lat_center)*cos(lon-k1_lon_center))
    k2 = 1.D0

    r = k1 / (4.D0*k2)
    det = sqrt(r*r + 2.D0*xt_constant*r)

    x  = (det-r)
    x2 = xt_constant/2.D0 - (det-r)/2.D0
  end subroutine initial_value_terminator_chemistry

  subroutine terminator_tendency( lat, lon, x, x2, dt, x_f, x2_f)

    !-----------------------------------------------------------------------
    ! Arguments:
    !-----------------------------------------------------------------------

    real(r8), intent(in)    :: lat, lon  ! latitude and longitude, radians
    real(r8), intent(in)    :: x, x2     ! molar mixing ratio of x and x2
    real(r8), intent(in)    :: dt        ! size of physics time step

    real(r8), intent(out)   :: x_f, x2_f  ! time rate of change of x and x2

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: xt_constant = 4.D-6
    real(r8), parameter :: k1_lat_center =   20.d0*pi/180.0D0
    real(r8), parameter :: k1_lon_center =  300.d0*pi/180.0D0


    real(r8) :: r, det, expdt, el ! useful algebraic quantities used in the computation
    real(r8) :: k1, k2            ! reaction rates
    real(r8) :: xt                ! quantity that should be conserved
    !
    !  Solar photolysis rate and recombination rate (for terminator chemistry)
    !
    k1 = 1.0D0*max(0.d0,sin(lat)*sin(k1_lat_center) + cos(lat)*cos(k1_lat_center)*cos(lon-k1_lon_center))
    k2 = 1.D0


    r = k1 / (4.D0*k2)
    xt = x + 2.D0* x2

    det = sqrt( r*r + 2.D0*r*xt )
    expdt = exp( -4.D0*k2*det*dt )

    if ( abs(det * k2 * dt) .gt. 1D-16 ) then
       el = (1.D0 - expdt) /det /dt
    else
       el = 4.D0*k2
    endif

    x_f  = -el * (x - det + r)*(x + det + r) / (1.D0 + expdt + dt*el*(x + r))
    x2_f = -x_f / 2.D0
  end subroutine terminator_tendency


end module dyn_inic_baroclinic
