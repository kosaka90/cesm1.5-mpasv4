#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
! This module should 'use' only module 'kinds'
!
module physical_constants
  ! ------------------------------
  use kinds, only : real_kind, longdouble_kind
  use physconst, only : pi,                    & ! _EXTERNAL
                        g => gravit,           &
                        rearth,                &
                        omega,                 &
                        Rgas => rair,          &
                        rh2o,                  &
                        Cp => cpair,           &
                        p0 => pstd,            &
                        MWDAIR => mwdry,       &
                        Rwater_vapor => rh2o,  &
                        Cpwater_vapor => cpwv, &
                        Cp_liq => cpliq,       &
                        Cp_ice => cpice,       &
                        kappa => cappa,        &
                        Rd_on_Rv => epsilo,    &
                        Cpd_on_Cpv,            &
                        rrearth => ra
  ! -----------------------------
  implicit none

  private

  real (kind=real_kind), public, parameter :: DD_PI = pi
  real (kind=longdouble_kind), public, parameter :: QQ_PI = 3.141592653589793238462643383279_longdouble_kind
  public                                   :: g              ! m s^-2
  public                                   :: rearth         ! m
  public                                   :: omega          ! s^-1
  public                                   :: Rgas
  public                                   :: Rh2o
  public                                   :: Cp
  public                                   :: p0             ! Pa
  public                                   :: cp_ice
  public                                   :: cp_liq
  public                                   :: MWDAIR
  public                                   :: Rwater_vapor
  public                                   :: Cpwater_vapor
  public                                   :: kappa
  public                                   :: Rd_on_Rv
  public                                   :: Cpd_on_Cpv
  public                                   :: rrearth         ! m

end module physical_constants
