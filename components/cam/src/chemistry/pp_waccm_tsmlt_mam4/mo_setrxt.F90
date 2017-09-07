
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,153) = 0.000258_r8
      rate(:,:,154) = 0.085_r8
      rate(:,:,155) = 1.31e-10_r8
      rate(:,:,156) = 3.5e-11_r8
      rate(:,:,157) = 9e-12_r8
      rate(:,:,158) = 1.2e-10_r8
      rate(:,:,163) = 1.2e-10_r8
      rate(:,:,164) = 1e-20_r8
      rate(:,:,165) = 1.3e-16_r8
      rate(:,:,167) = 4.2e-13_r8
      rate(:,:,169) = 8e-14_r8
      rate(:,:,170) = 3.9e-17_r8
      rate(:,:,177) = 6.9e-12_r8
      rate(:,:,178) = 7.2e-11_r8
      rate(:,:,179) = 1.6e-12_r8
      rate(:,:,185) = 1.8e-12_r8
      rate(:,:,189) = 1.8e-12_r8
      rate(:,:,193) = 7e-13_r8
      rate(:,:,194) = 5e-12_r8
      rate(:,:,203) = 3.5e-12_r8
      rate(:,:,205) = 1e-11_r8
      rate(:,:,206) = 2.2e-11_r8
      rate(:,:,207) = 5e-11_r8
      rate(:,:,242) = 1.7e-13_r8
      rate(:,:,244) = 2.607e-10_r8
      rate(:,:,245) = 9.75e-11_r8
      rate(:,:,246) = 2.07e-10_r8
      rate(:,:,247) = 2.088e-10_r8
      rate(:,:,248) = 1.17e-10_r8
      rate(:,:,249) = 4.644e-11_r8
      rate(:,:,250) = 1.204e-10_r8
      rate(:,:,251) = 9.9e-11_r8
      rate(:,:,252) = 3.3e-12_r8
      rate(:,:,271) = 4.5e-11_r8
      rate(:,:,272) = 4.62e-10_r8
      rate(:,:,273) = 1.2e-10_r8
      rate(:,:,274) = 9e-11_r8
      rate(:,:,275) = 3e-11_r8
      rate(:,:,280) = 2.14e-11_r8
      rate(:,:,281) = 1.9e-10_r8
      rate(:,:,294) = 2.57e-10_r8
      rate(:,:,295) = 1.8e-10_r8
      rate(:,:,296) = 1.794e-10_r8
      rate(:,:,297) = 1.3e-10_r8
      rate(:,:,298) = 7.65e-11_r8
      rate(:,:,312) = 4e-13_r8
      rate(:,:,322) = 6.8e-14_r8
      rate(:,:,323) = 2e-13_r8
      rate(:,:,337) = 7e-13_r8
      rate(:,:,338) = 1e-12_r8
      rate(:,:,342) = 1e-14_r8
      rate(:,:,343) = 1e-11_r8
      rate(:,:,344) = 1.15e-11_r8
      rate(:,:,345) = 4e-14_r8
      rate(:,:,358) = 3e-12_r8
      rate(:,:,359) = 6.7e-13_r8
      rate(:,:,369) = 3.5e-13_r8
      rate(:,:,370) = 5.4e-11_r8
      rate(:,:,373) = 2e-12_r8
      rate(:,:,374) = 1.4e-11_r8
      rate(:,:,377) = 2.4e-12_r8
      rate(:,:,388) = 5e-12_r8
      rate(:,:,398) = 1.6e-12_r8
      rate(:,:,400) = 6.7e-12_r8
      rate(:,:,403) = 3.5e-12_r8
      rate(:,:,406) = 1.3e-11_r8
      rate(:,:,407) = 1.4e-11_r8
      rate(:,:,411) = 2.4e-12_r8
      rate(:,:,412) = 1.4e-11_r8
      rate(:,:,417) = 2.4e-12_r8
      rate(:,:,418) = 4e-11_r8
      rate(:,:,419) = 4e-11_r8
      rate(:,:,421) = 1.4e-11_r8
      rate(:,:,425) = 2.4e-12_r8
      rate(:,:,426) = 4e-11_r8
      rate(:,:,430) = 7e-11_r8
      rate(:,:,431) = 1e-10_r8
      rate(:,:,436) = 2.4e-12_r8
      rate(:,:,451) = 4.7e-11_r8
      rate(:,:,464) = 2.1e-12_r8
      rate(:,:,465) = 2.8e-13_r8
      rate(:,:,473) = 1.7e-11_r8
      rate(:,:,479) = 8.4e-11_r8
      rate(:,:,481) = 1.9e-11_r8
      rate(:,:,482) = 1.2e-14_r8
      rate(:,:,483) = 2e-10_r8
      rate(:,:,490) = 2.4e-12_r8
      rate(:,:,491) = 2e-11_r8
      rate(:,:,495) = 2.3e-11_r8
      rate(:,:,496) = 2e-11_r8
      rate(:,:,500) = 3.3e-11_r8
      rate(:,:,501) = 1e-12_r8
      rate(:,:,502) = 5.7e-11_r8
      rate(:,:,503) = 3.4e-11_r8
      rate(:,:,506) = 2.3e-12_r8
      rate(:,:,507) = 1.2e-11_r8
      rate(:,:,508) = 5.7e-11_r8
      rate(:,:,509) = 2.8e-11_r8
      rate(:,:,510) = 6.6e-11_r8
      rate(:,:,511) = 1.4e-11_r8
      rate(:,:,514) = 1.9e-12_r8
      rate(:,:,525) = 6.34e-08_r8
      rate(:,:,526) = 6.34e-08_r8
      rate(:,:,529) = 1.9e-11_r8
      rate(:,:,530) = 1.2e-14_r8
      rate(:,:,531) = 2e-10_r8
      rate(:,:,536) = 1.34e-11_r8
      rate(:,:,540) = 1.34e-11_r8
      rate(:,:,542) = 1.7e-11_r8
      rate(:,:,563) = 6e-11_r8
      rate(:,:,566) = 1e-12_r8
      rate(:,:,567) = 4e-10_r8
      rate(:,:,568) = 2e-10_r8
      rate(:,:,569) = 1e-10_r8
      rate(:,:,570) = 5e-16_r8
      rate(:,:,571) = 4.4e-10_r8
      rate(:,:,572) = 9e-10_r8
      rate(:,:,575) = 1.2e-14_r8
      rate(:,:,586) = 1.2e-10_r8
      rate(:,:,589) = 2.8e-13_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,159) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      rate(:,:,160) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,161) = 2.64e-11_r8 * exp_fac(:,:)
      rate(:,:,162) = 6.6e-12_r8 * exp_fac(:,:)
      rate(:,:,166) = 3.6e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,168) = 1.8e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,171) = 3.5e-11_r8 * exp( -135._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2060._r8 * itemp(:,:) )
      rate(:,:,172) = 8e-12_r8 * exp_fac(:,:)
      rate(:,:,588) = 8e-12_r8 * exp_fac(:,:)
      rate(:,:,175) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,176) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,427) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,534) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,580) = 1.05e-14_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,181) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,269) = 5.5e-12_r8 * exp_fac(:,:)
      rate(:,:,308) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,327) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,354) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,362) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,366) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,382) = 2.3e-11_r8 * exp_fac(:,:)
      rate(:,:,392) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,402) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,429) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,437) = 1.52e-12_r8 * exp_fac(:,:)
      rate(:,:,443) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,446) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,450) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,466) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,470) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,476) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,480) = 3.8e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -490._r8 * itemp(:,:) )
      rate(:,:,182) = 1e-14_r8 * exp_fac(:,:)
      rate(:,:,578) = 1e-14_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -470._r8 * itemp(:,:) )
      rate(:,:,183) = 1.4e-10_r8 * exp_fac(:,:)
      rate(:,:,579) = 1.4e-10_r8 * exp_fac(:,:)
      rate(:,:,184) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,186) = 4.8e-11_r8 * exp_fac(:,:)
      rate(:,:,267) = 1.7e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,187) = 1.8e-11_r8 * exp_fac(:,:)
      rate(:,:,340) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,353) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,361) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,390) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,410) = 4.4e-12_r8 * exp_fac(:,:)
      rate(:,:,416) = 4.4e-12_r8 * exp_fac(:,:)
      rate(:,:,489) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,494) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,499) = 4.2e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -940._r8 * itemp(:,:) )
      rate(:,:,188) = 1.7e-12_r8 * exp_fac(:,:)
      rate(:,:,587) = 1.7e-12_r8 * exp_fac(:,:)
      rate(:,:,192) = 1.3e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,195) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,196) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,197) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,198) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,199) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,200) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,201) = 1.2e-13_r8 * exp_fac(:,:)
      rate(:,:,227) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,584) = 1.2e-13_r8 * exp_fac(:,:)
      rate(:,:,204) = 1.5e-11_r8 * exp( 170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,208) = 3.3e-12_r8 * exp_fac(:,:)
      rate(:,:,223) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,237) = 7.4e-12_r8 * exp_fac(:,:)
      rate(:,:,336) = 8.1e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,209) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,268) = 5.8e-12_r8 * exp_fac(:,:)
      rate(:,:,585) = 3e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,211) = 7.26e-11_r8 * exp_fac(:,:)
      rate(:,:,212) = 4.64e-11_r8 * exp_fac(:,:)
      rate(:,:,219) = 8.1e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,220) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:,:) )
      rate(:,:,221) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,222) = 1.1e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,224) = 3.6e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,225) = 2.3e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,226) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,228) = 1e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,229) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,230) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,231) = 6.4e-12_r8 * exp_fac(:,:)
      rate(:,:,261) = 4.1e-13_r8 * exp_fac(:,:)
      rate(:,:,439) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,453) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,456) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,459) = 7.5e-12_r8 * exp_fac(:,:)
      rate(:,:,232) = 6.5e-12_r8 * exp( 135._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,234) = 3.6e-12_r8 * exp_fac(:,:)
      rate(:,:,283) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,235) = 1.2e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,236) = 2.8e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,238) = 6e-13_r8 * exp_fac(:,:)
      rate(:,:,258) = 1.5e-12_r8 * exp_fac(:,:)
      rate(:,:,266) = 1.9e-11_r8 * exp_fac(:,:)
      rate(:,:,239) = 1e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,240) = 1.8e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,241) = 3.4e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,243) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,277) = 1.4e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,255) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,282) = 6.3e-12_r8 * exp_fac(:,:)
      rate(:,:,256) = 4.8e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,257) = 1.6e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,259) = 9.5e-13_r8 * exp( 550._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,260) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,263) = 8.8e-12_r8 * exp_fac(:,:)
      rate(:,:,262) = 4.5e-12_r8 * exp( 460._r8 * itemp(:,:) )
      rate(:,:,265) = 1.9e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,270) = 1.2e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,276) = 1.6e-10_r8 * exp( -260._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,278) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,280) = 2.14e-11_r8 * exp_fac(:,:)
      rate(:,:,281) = 1.9e-10_r8 * exp_fac(:,:)
      rate(:,:,294) = 2.57e-10_r8 * exp_fac(:,:)
      rate(:,:,295) = 1.8e-10_r8 * exp_fac(:,:)
      rate(:,:,296) = 1.794e-10_r8 * exp_fac(:,:)
      rate(:,:,297) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,298) = 7.65e-11_r8 * exp_fac(:,:)
      rate(:,:,312) = 4e-13_r8 * exp_fac(:,:)
      rate(:,:,322) = 6.8e-14_r8 * exp_fac(:,:)
      rate(:,:,323) = 2e-13_r8 * exp_fac(:,:)
      rate(:,:,337) = 7e-13_r8 * exp_fac(:,:)
      rate(:,:,338) = 1e-12_r8 * exp_fac(:,:)
      rate(:,:,342) = 1e-14_r8 * exp_fac(:,:)
      rate(:,:,343) = 1e-11_r8 * exp_fac(:,:)
      rate(:,:,344) = 1.15e-11_r8 * exp_fac(:,:)
      rate(:,:,345) = 4e-14_r8 * exp_fac(:,:)
      rate(:,:,358) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,359) = 6.7e-13_r8 * exp_fac(:,:)
      rate(:,:,369) = 3.5e-13_r8 * exp_fac(:,:)
      rate(:,:,370) = 5.4e-11_r8 * exp_fac(:,:)
      rate(:,:,373) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,374) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,377) = 2.4e-12_r8 * exp_fac(:,:)
      rate(:,:,388) = 5e-12_r8 * exp_fac(:,:)
      rate(:,:,398) = 1.6e-12_r8 * exp_fac(:,:)
      rate(:,:,400) = 6.7e-12_r8 * exp_fac(:,:)
      rate(:,:,403) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,406) = 1.3e-11_r8 * exp_fac(:,:)
      rate(:,:,407) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,411) = 2.4e-12_r8 * exp_fac(:,:)
      rate(:,:,412) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,417) = 2.4e-12_r8 * exp_fac(:,:)
      rate(:,:,418) = 4e-11_r8 * exp_fac(:,:)
      rate(:,:,419) = 4e-11_r8 * exp_fac(:,:)
      rate(:,:,421) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,425) = 2.4e-12_r8 * exp_fac(:,:)
      rate(:,:,426) = 4e-11_r8 * exp_fac(:,:)
      rate(:,:,430) = 7e-11_r8 * exp_fac(:,:)
      rate(:,:,431) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,436) = 2.4e-12_r8 * exp_fac(:,:)
      rate(:,:,451) = 4.7e-11_r8 * exp_fac(:,:)
      rate(:,:,464) = 2.1e-12_r8 * exp_fac(:,:)
      rate(:,:,465) = 2.8e-13_r8 * exp_fac(:,:)
      rate(:,:,473) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,479) = 8.4e-11_r8 * exp_fac(:,:)
      rate(:,:,481) = 1.9e-11_r8 * exp_fac(:,:)
      rate(:,:,482) = 1.2e-14_r8 * exp_fac(:,:)
      rate(:,:,483) = 2e-10_r8 * exp_fac(:,:)
      rate(:,:,490) = 2.4e-12_r8 * exp_fac(:,:)
      rate(:,:,491) = 2e-11_r8 * exp_fac(:,:)
      rate(:,:,495) = 2.3e-11_r8 * exp_fac(:,:)
      rate(:,:,496) = 2e-11_r8 * exp_fac(:,:)
      rate(:,:,500) = 3.3e-11_r8 * exp_fac(:,:)
      rate(:,:,501) = 1e-12_r8 * exp_fac(:,:)
      rate(:,:,502) = 5.7e-11_r8 * exp_fac(:,:)
      rate(:,:,503) = 3.4e-11_r8 * exp_fac(:,:)
      rate(:,:,506) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,507) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,508) = 5.7e-11_r8 * exp_fac(:,:)
      rate(:,:,509) = 2.8e-11_r8 * exp_fac(:,:)
      rate(:,:,510) = 6.6e-11_r8 * exp_fac(:,:)
      rate(:,:,511) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,514) = 1.9e-12_r8 * exp_fac(:,:)
      rate(:,:,525) = 6.34e-08_r8 * exp_fac(:,:)
      rate(:,:,526) = 6.34e-08_r8 * exp_fac(:,:)
      rate(:,:,529) = 1.9e-11_r8 * exp_fac(:,:)
      rate(:,:,530) = 1.2e-14_r8 * exp_fac(:,:)
      rate(:,:,531) = 2e-10_r8 * exp_fac(:,:)
      rate(:,:,536) = 1.34e-11_r8 * exp_fac(:,:)
      rate(:,:,540) = 1.34e-11_r8 * exp_fac(:,:)
      rate(:,:,542) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,563) = 6e-11_r8 * exp_fac(:,:)
      rate(:,:,566) = 1e-12_r8 * exp_fac(:,:)
      rate(:,:,567) = 4e-10_r8 * exp_fac(:,:)
      rate(:,:,568) = 2e-10_r8 * exp_fac(:,:)
      rate(:,:,569) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,570) = 5e-16_r8 * exp_fac(:,:)
      rate(:,:,571) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,572) = 9e-10_r8 * exp_fac(:,:)
      rate(:,:,575) = 1.2e-14_r8 * exp_fac(:,:)
      rate(:,:,586) = 1.2e-10_r8 * exp_fac(:,:)
      rate(:,:,589) = 2.8e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,279) = 6e-12_r8 * exp_fac(:,:)
      rate(:,:,375) = 5e-13_r8 * exp_fac(:,:)
      rate(:,:,408) = 5e-13_r8 * exp_fac(:,:)
      rate(:,:,413) = 5e-13_r8 * exp_fac(:,:)
      rate(:,:,422) = 5e-13_r8 * exp_fac(:,:)
      rate(:,:,433) = 5e-13_r8 * exp_fac(:,:)
      rate(:,:,284) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:,:) )
      rate(:,:,285) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1520._r8 * itemp(:,:) )
      rate(:,:,286) = 1.64e-12_r8 * exp_fac(:,:)
      rate(:,:,394) = 8.5e-16_r8 * exp_fac(:,:)
      rate(:,:,583) = 8.5e-16_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1100._r8 * itemp(:,:) )
      rate(:,:,287) = 2.03e-11_r8 * exp_fac(:,:)
      rate(:,:,513) = 3.4e-12_r8 * exp_fac(:,:)
      rate(:,:,288) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:,:) )
      rate(:,:,289) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,290) = 9e-13_r8 * exp( -360._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,291) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,301) = 3.4e-11_r8 * exp_fac(:,:)
      rate(:,:,292) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,293) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:,:) )
      rate(:,:,299) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      rate(:,:,300) = 6e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,302) = 5.5e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,303) = 5e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,304) = 1.9e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,305) = 4.1e-13_r8 * exp( 750._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,306) = 2.8e-12_r8 * exp_fac(:,:)
      rate(:,:,365) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,307) = 2.9e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,309) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,313) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,324) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,339) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,352) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,360) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,364) = 8.6e-13_r8 * exp_fac(:,:)
      rate(:,:,376) = 8e-13_r8 * exp_fac(:,:)
      rate(:,:,389) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,399) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,409) = 8e-13_r8 * exp_fac(:,:)
      rate(:,:,414) = 8e-13_r8 * exp_fac(:,:)
      rate(:,:,423) = 8e-13_r8 * exp_fac(:,:)
      rate(:,:,434) = 8e-13_r8 * exp_fac(:,:)
      rate(:,:,441) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,445) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,448) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,461) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,468) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,474) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,477) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,488) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,493) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,498) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,314) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,315) = 2.6e-12_r8 * exp( 265._r8 * itemp(:,:) )
      rate(:,:,316) = 1.08e-10_r8 * exp( 105._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2630._r8 * itemp(:,:) )
      rate(:,:,321) = 1.2e-14_r8 * exp_fac(:,:)
      rate(:,:,576) = 1.2e-14_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 365._r8 * itemp(:,:) )
      rate(:,:,325) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,442) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,447) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,449) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,462) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,469) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,475) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,478) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,326) = 6.9e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,328) = 7.2e-11_r8 * exp( -70._r8 * itemp(:,:) )
      rate(:,:,329) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,330) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,350) = 6.5e-15_r8 * exp_fac(:,:)
      rate(:,:,577) = 6.5e-15_r8 * exp_fac(:,:)
      rate(:,:,331) = 4.63e-12_r8 * exp( 350._r8 * itemp(:,:) )
      rate(:,:,332) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,333) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,334) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,363) = 7.1e-13_r8 * exp_fac(:,:)
      rate(:,:,384) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,487) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,492) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,497) = 2e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,335) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,385) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,438) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,452) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,455) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,458) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,341) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,349) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,351) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,355) = 8.7e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,356) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:,:) )
      rate(:,:,357) = 8.4e-13_r8 * exp( 830._r8 * itemp(:,:) )
      rate(:,:,371) = 4.8e-12_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:,372) = 5.1e-14_r8 * exp( 693._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,378) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,379) = 1.3e-13_r8 * exp_fac(:,:)
      rate(:,:,381) = 9.6e-12_r8 * exp_fac(:,:)
      rate(:,:,387) = 5.3e-12_r8 * exp_fac(:,:)
      rate(:,:,424) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,435) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,380) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,383) = 4.6e-12_r8 * exp_fac(:,:)
      rate(:,:,386) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,391) = 2.3e-12_r8 * exp( -170._r8 * itemp(:,:) )
      rate(:,:,395) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,401) = 5.4e-14_r8 * exp( 870._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,404) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,405) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,415) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -446._r8 * itemp(:,:) )
      rate(:,:,420) = 3.03e-12_r8 * exp_fac(:,:)
      rate(:,:,533) = 3.03e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 410._r8 * itemp(:,:) )
      rate(:,:,428) = 2.54e-11_r8 * exp_fac(:,:)
      rate(:,:,535) = 2.54e-11_r8 * exp_fac(:,:)
      rate(:,:,432) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -193._r8 * itemp(:,:) )
      rate(:,:,440) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,532) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,444) = 5.9e-12_r8 * exp( 225._r8 * itemp(:,:) )
      rate(:,:,463) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 352._r8 * itemp(:,:) )
      rate(:,:,471) = 1.7e-12_r8 * exp_fac(:,:)
      rate(:,:,541) = 1.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 490._r8 * itemp(:,:) )
      rate(:,:,484) = 1.2e-12_r8 * exp_fac(:,:)
      rate(:,:,537) = 1.2e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -580._r8 * itemp(:,:) )
      rate(:,:,485) = 6.3e-16_r8 * exp_fac(:,:)
      rate(:,:,538) = 6.3e-16_r8 * exp_fac(:,:)
      rate(:,:,582) = 6.3e-16_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 440._r8 * itemp(:,:) )
      rate(:,:,486) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,539) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,504) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:,:) )
      rate(:,:,505) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:,:) )
      rate(:,:,512) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:,:) )
      rate(:,:,515) = 2.7e-11_r8 * exp( 335._r8 * itemp(:,:) )
      rate(:,:,518) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,519) = 9.6e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,520) = 1.7e-12_r8 * exp( -710._r8 * itemp(:,:) )
      rate(:,:,581) = 4.4e-15_r8 * exp( -2500._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.4e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,180), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.9e-31_r8 * itemp(:,:)**1._r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,190), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.5e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,202), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3e-11_r8
      call jpl( rate(1,1,210), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.9e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 4e-12_r8 * itemp(:,:)**0.3_r8
      call jpl( rate(1,1,213), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.4e-30_r8 * itemp(:,:)**3._r8
      kinf(:,:) = 1.6e-12_r8 * itemp(:,:)**(-0.1_r8)
      call jpl( rate(1,1,214), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.8e-30_r8 * itemp(:,:)**3._r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,215), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.8e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,233), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.9e-32_r8 * itemp(:,:)**3.6_r8
      kinf(:,:) = 3.7e-12_r8 * itemp(:,:)**1.6_r8
      call jpl( rate(1,1,253), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.2e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,264), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.9e-33_r8 * itemp(:,:)**1._r8
      kinf(:,:) = 1.1e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,310), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.3e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,311), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 5.2e-30_r8 * itemp(:,:)**2.4_r8
      kinf(:,:) = 2.2e-10_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,318), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.5e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2._r8)
      call jpl( rate(1,1,319), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.6e-29_r8 * itemp(:,:)**3.3_r8
      kinf(:,:) = 3.1e-10_r8 * itemp(:,:)
      call jpl( rate(1,1,320), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.6e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,346), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.7e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.3e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,347), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3e-11_r8
      call jpl( rate(1,1,367), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3e-11_r8
      call jpl( rate(1,1,393), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 9.7e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.3e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,454), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.7e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.3e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,457), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.7e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.3e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,460), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.7e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.3e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,467), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)

      rate(:,:kbot,164) = 1e-20_r8
      rate(:,:kbot,165) = 1.3e-16_r8
      rate(:,:kbot,169) = 8e-14_r8
      rate(:,:kbot,170) = 3.9e-17_r8
      rate(:,:kbot,177) = 6.9e-12_r8
      rate(:,:kbot,193) = 7e-13_r8
      rate(:,:kbot,194) = 5e-12_r8
      rate(:,:kbot,563) = 6e-11_r8
      rate(:,:kbot,566) = 1e-12_r8
      rate(:,:kbot,567) = 4e-10_r8
      rate(:,:kbot,568) = 2e-10_r8
      rate(:,:kbot,569) = 1e-10_r8
      rate(:,:kbot,571) = 4.4e-10_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,160) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,161) = 2.64e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,162) = 6.6e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,166) = 3.6e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,168) = 1.8e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,171) = 3.5e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,172) = 8e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,181) = 3e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,182) = 1e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,183) = 1.4e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,186) = 4.8e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,187) = 1.8e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,188) = 1.7e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,195) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,199) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,200) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:kbot,208) = 3.3e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,209) = 3e-12_r8 * exp( -1500._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.4e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,180) = wrk(:,:)























      end subroutine setrxt_hrates

      end module mo_setrxt
