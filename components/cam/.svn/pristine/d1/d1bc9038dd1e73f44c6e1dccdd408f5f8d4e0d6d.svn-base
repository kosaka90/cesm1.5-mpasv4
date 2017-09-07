#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module perfmodel_mod


    use kinds , only : real_kind, log_kind, int_kind

implicit none
private

    type,  public :: perf_t
       real(kind=real_kind)   :: onnode
       real(kind=real_kind)   :: offnode
        end type perf_t

    type, public :: commtime_t
       real(kind=real_kind)  :: latency
       real(kind=real_kind)  :: bandwidth
    end type commtime_t

    public :: SetPtoPNetworkParams
        public :: SetReduceNetworkParams
    public :: SetSerialParamsExplicit
        public :: SetSerialParamsPrimExplicit
    public :: SetSerialParamsImplicit

contains

    subroutine SetPtoPNetworkParams(offnode,onnode,network,found)

    type (commtime_t),intent(inout)     :: offnode
    type (commtime_t),intent(inout)     :: onnode
        character(len=80)                   :: network
        logical(kind=log_kind),intent(inout)  :: found
        real (kind=real_kind) :: contention

    ! Macine specific characteristics
    ! bluesky Power4 with dual rail 8-way LPAR
        ! bandwidth: microseconds per byte
        ! latency:   microseconds

    found=.FALSE.
        if(network(1:7) == "bluesky") then
           offnode%latency   = 17.9_real_kind
           onnode%latency    = .9_real_kind
           offnode%bandwidth = 0.002777_real_kind
           !offnode%bandwidth = 3._real_kind*0.002777_real_kind
           onnode%bandwidth  = 0.001203_real_kind
           found=.TRUE.
        elseif (network(1:11) == "blackforest") then
           offnode%latency   = 30.0_real_kind
           onnode%latency    = 2._real_kind
           offnode%bandwidth = 0.0125_real_kind
           onnode%bandwidth  = 0.00416_real_kind
           found=.TRUE.
        elseif (network(1:8) == "protoBGL") then
       ! This is the from the description of the prototype BGL from Oct 2003 meeting
           offnode%latency   = 6._real_kind
           onnode%latency    = .1_real_kind
           !contention = 3.0_real_kind
           !contention = 1.0_real_kind
           contention =  9.0_real_kind
           offnode%bandwidth = contention/420._real_kind! Note the value of 420 is 3 times what I get on point-to-point
                          ! messaging... I use the three because they only get 3 receive
                          ! streams going at once.  The other three is due to the fact
                          ! that the corner processors only have 3 network connections
                          ! to simulate contention 8 neighbors/6 links * 2.25 = 3
           onnode%bandwidth  =1.0_real_kind/(3._real_kind*420._real_kind)
           found=.TRUE.
        elseif (network(1:8) == "BGL-fast") then
       ! BGL with network as described in the Paper "???"
           offnode%latency   = 1.5_real_kind
           onnode%latency    = .1_real_kind
           !contention = 2.25_real_kind
           contention = 4.5_real_kind
           offnode%bandwidth = contention/1050._real_kind !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           onnode%bandwidth  =1.0_real_kind/(3._real_kind*1050._real_kind)
           found=.TRUE.
        elseif (network(1:8) == "BGL-slow") then
       ! BGL with network running at 1/2 of described in the Paper "???"
           offnode%latency   = 4.5_real_kind
           onnode%latency    = .1_real_kind
           !contention = 2.25_real_kind
           contention = 4.5_real_kind
           offnode%bandwidth = contention/525._real_kind    !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           onnode%bandwidth  =1.0_real_kind/(3._real_kind*525._real_kind)
           found=.TRUE.
        else
       print *,"SetPtoPNetworkParams:  Undefined network name"
       print *,"SetPtoPNetworkParams:  Please select from the following or create your own:"
           print *,"SetPtoPNetworkParams:  bluesky     == Power4 with dual rail colony and 8-way LPARs"
           print *,"SetPtoPNetworkParams:  blackforest == Power3 with TBMX network"
           print *,"SetPtoPNetworkParams:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
           print *,"SetPtoPNetworkParams:  BGL-fast    == BGL network as described in Paper ??? "
           print *,"SetPtoPNetworkParams:  BGL-slow    == 1/2 perf. of BGL network as described in Paper ???? "
        endif

    end subroutine SetPtoPNetworkParams

    subroutine SetReduceNetworkParams(latency,bandwidth,network,found)

    type (perf_t),intent(inout)         :: latency
    type (perf_t),intent(inout)         :: bandwidth
        character(len=80)                   :: network
        logical(kind=log_kind),intent(inout)  :: found

    ! Macine specific characteristics
    ! bluesky Power4 with dual rail 8-way LPAR
        ! bandwidth: microseconds per byte
        ! latency:   microseconds

    found=.FALSE.
        if(network(1:7) == "bluesky") then
           latency%offnode   = 17.9_real_kind
           latency%onnode    = .9_real_kind
           bandwidth%offnode = 0.002777_real_kind
           bandwidth%onnode  = 0.001203_real_kind
           found=.TRUE.
        elseif (network(1:8) == "protoBGL") then
       ! This is the from the description of the prototype BGL from Oct 2003 meeting
           latency%offnode   = 3._real_kind
           latency%onnode    = .1_real_kind
           bandwidth%offnode = 1._real_kind/350._real_kind    ! Add 3 here because the corner processors only
                      ! have three network connections (just a 3-D mesh) at the moment
           bandwidth%onnode  =1.0_real_kind/(3._real_kind*420._real_kind)
           found=.TRUE.
        elseif (network(1:8) == "BGL-fast") then
       ! BGL with network as described in the Paper "???"
           latency%offnode   = 1.5_real_kind
           latency%onnode    = .1_real_kind
           bandwidth%offnode = 1.0_real_kind/350._real_kind !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           bandwidth%onnode  =1.0_real_kind/(3._real_kind*350._real_kind)
           found=.TRUE.
        elseif (network(1:8) == "BGL-slow") then
       ! BGL with network running at 1/2 of described in the Paper "???"
           latency%offnode   = 1.5_real_kind
           latency%onnode    = .1_real_kind
           bandwidth%offnode = 1.0_real_kind/350._real_kind    !JMD I added 2.25 because I
                                            !JMD wanted simulate contention
                                            !JMD 8 neighbors/6 links * 2.25 = 3
           bandwidth%onnode  =1.0_real_kind/(3._real_kind*350._real_kind)
           found=.TRUE.
        else
       print *,"SetReduceNetworkParams:  Undefined network name"
       print *,"SetReduceNetworkParams:  Please select from the following or create your own:"
           print *,"SetReduceNetworkParams:  bluesky    == Power4 with dual rail colony and 8-way LPARs"
           print *,"SetReduceNetworkParams:  protoBGL   == Prototype BGL machine (November 2003 configuration)"
           print *,"SetReduceNetworkParams:  BGL-fast   == BGL network as described in Paper ??? "
           print *,"SetReduceNetworkParams:  BGL-slow   == 1/2 perf. of BGL network as described in Paper ???? "
        endif

    end subroutine SetReduceNetworkParams

    subroutine SetSerialParamsExplicit(microsec,np,nlev,machine,found)

        real (kind=real_kind),intent(out)   :: microsec
        integer(kind=int_kind),intent(in)   :: np
        integer(kind=int_kind),intent(in)   :: nlev
        character(len=80)                   :: machine
        logical(kind=log_kind),intent(inout)  :: found
        logical(kind=log_kind)                :: m_found,np_found,nlev_found

    m_found    = .FALSE.
    np_found   = .FALSE.
    nlev_found = .FALSE.
        if(machine(1:7) == "bluesky") then
       m_found = .TRUE.
           if(np .eq. 6 ) then
         np_found=.TRUE.
             select case (nlev)
                case(96)
!                   microsec         = 4243.D0 ! SMP numbers (nproc=8)
                    microsec         = 3299.D0 ! single threaded
            nlev_found=.TRUE.
                case(60)
!                 microsec         = 2545.D0 ! SMP number (nproc=8)
!                 microsec         = 2194.D0 ! SMP number (nproc=2)
                  microsec         = 2013.D0 ! single threaded
          nlev_found=.TRUE.
                case(48)
!                 microsec         = 1914.D0 ! SMP number (nproc=8)
                  microsec         = 1702.D0 ! single threaded
          nlev_found=.TRUE.
                case(30)
!                 microsec         = 1142.D0 ! SMP number (nproc=8)
                  microsec         = 1023.D0 ! single threaded
          nlev_found=.TRUE.
                case(16)
                  microsec         =  540.D0
          nlev_found=.TRUE.
                case(2)
                  microsec         =   76.D0
          nlev_found=.TRUE.
        case default
              print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
         end select
          elseif( np .eq. 8) then
         np_found=.TRUE.
             select case (nlev)
                case(16)
                  microsec         =  556.D0
          nlev_found=.TRUE.
                case(2)
                  microsec         =   76.D0
          nlev_found=.FALSE.
        case default
              print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
         end select
      elseif (np .eq. 10) then
         np_found = .TRUE.
             select case (nlev)
                case(96)
          microsec        = 9570.D0 ! single threaded
          nlev_found = .TRUE.
        case default
              print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
         end select
          else
         print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
         print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np
          endif
        elseif(machine(1:11) == "blackforest") then
       m_found = .TRUE.
           if(np .eq. 8 ) then
         np_found=.TRUE.
             select case (nlev)
                case(16)
                  microsec         =  1652.D0  ! single threaded
          nlev_found=.TRUE.
        case default
              print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsExplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
         end select
      elseif (np .eq. 10) then
         np_found = .TRUE.
             select case (nlev)
                case(96)
          !microsec        = 9570.D0 ! Please Fix Me
          nlev_found = .FALSE.
        case default
              print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsExplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
         end select
          else
         print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
         print *,"SetSerialParamsExplicit:  machine: ",machine(1:11)," np := ",np
          endif
        elseif(machine(1:8) == "protoBGL") then
       m_found = .TRUE.
          if (np .eq. 8) then
         np_found = .TRUE.
             select case (nlev)
                case(16)
          microsec        = 2087.D0 ! single threaded
          nlev_found = .TRUE.
        case default
              print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
         end select
      elseif( np .eq. 10 ) then
         np_found = .TRUE.
             select case (nlev)
                case(96)
          microsec        = (900._real_kind/210._real_kind)*9570._real_kind ! This is an estimate
          nlev_found = .TRUE.
        case default
              print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
         end select
          else
         print *,"SetSerialParamsExplicit:  Could not find performance estimate for: "
         print *,"SetSerialParamsExplicit:  machine: ",machine(1:8)," np := ",np
          endif
        elseif(machine(1:8) == "BGL-fast") then
       m_found = .TRUE.
       if (np .eq. 10) then
         np_found = .TRUE.
             select case (nlev)
                case(96)
          ! Note: Bluesky runs the code at 1140 Mflops in 5761.7
              !       Estimate BGL-fast runs at about 750 Mflops
          !       The
          microsec   = (1140.0_real_kind/750.0_real_kind)* 5761.7_real_kind
          nlev_found = .TRUE.
        case default
              print *,"SetSerialParams:  Could not find performance estimate for: "
              print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
         end select
           else
         print *,"SetSerialParams:  Could not find performance estimate for: "
         print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np
           endif
        elseif(machine(1:8) == "BGL-slow") then
       m_found = .TRUE.
       if (np .eq. 10) then
         np_found = .TRUE.
             select case (nlev)
                case(96)
          ! Note: Bluesky runs the code at 1140 Mflops in 5761.7
              !       Estimate BGL-slow runs at about 450 Mflops
          !       The
          microsec   = (1140.0_real_kind/450.0_real_kind)* 5761.7_real_kind
          nlev_found = .TRUE.
        case default
              print *,"SetSerialParams:  Could not find performance estimate for: "
              print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
          end select
            else
          print *,"SetSerialParams:  Could not find performance estimate for: "
          print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np
        endif
         else
       print *,"SetSerialParamsExplicit:   Undefined machine name"
       print *,"SetSerialParamsExplicit:  Please select from the following or create your own:"
           print *,"SetSerialParamsExplicit:  bluesky     == Power4 with dual rail colony and 8-way LPARs"
           print *,"SetSerialParamsExplicit:  blackforest == Power3 with TBMX network"
           print *,"SetSerialParamsExplicit:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
           print *,"SetSerialParamsExplicit:  BGL-fast    == BGL processor clock scaled based on Power3 "
           print *,"SetSerialParamsExplicit:  BGL-slow    == BGL processor clock scaled based on Power4 "
       return
         endif


    if(m_found) then
      if(np_found) then
        if(nlev_found) found=.TRUE.
      endif
    endif


    end subroutine SetSerialParamsExplicit

    subroutine SetSerialParamsPrimExplicit(microsec,np,nlev,machine,found)

        real (kind=real_kind),intent(out)   :: microsec
        integer(kind=int_kind),intent(in)   :: np
        integer(kind=int_kind),intent(in)   :: nlev
        character(len=80)                   :: machine
        logical(kind=log_kind),intent(inout)  :: found
        logical(kind=log_kind)                :: m_found,np_found,nlev_found

    m_found    = .FALSE.
    np_found   = .FALSE.
    nlev_found = .FALSE.
        if(machine(1:8) == "protoBGL") then
       m_found = .TRUE.
          if (np .eq. 8) then
         np_found = .TRUE.
             select case (nlev)
                case(26)
          microsec        = 10706.8_real_kind ! single threaded
          nlev_found = .TRUE.
        case default
              print *,"SetSerialParamsPrimExplicit:  Could not find performance estimate for: "
              print *,"SetSerialParamsPrimExplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
         end select
          else
         print *,"SetSerialParamsPrimExplicit:  Could not find performance estimate for: "
         print *,"SetSerialParamsPrimExplicit:  machine: ",machine(1:8)," np := ",np
          endif
        elseif(machine(1:8) == "BGL-slow") then
       m_found = .TRUE.
       if (np .eq. 8) then
         np_found = .TRUE.
             select case (nlev)
                case(26)
          ! Assume it scales like clock speed from the prototype BGL
          microsec   = (500._real_kind/750._real_kind)* 10706.8_real_kind
          nlev_found = .TRUE.
        case default
              print *,"SetSerialPrimParams:  Could not find performance estimate for: "
              print *,"SetSerialPrimParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
          end select
            else
          print *,"SetSerialPrimParams:  Could not find performance estimate for: "
          print *,"SetSerialPrimParams:  machine: ",machine(1:8)," np := ",np
        endif
         else
       print *,"SetSerialParamsPrimExplicit:   Undefined machine name"
       print *,"SetSerialParamsPrimExplicit:  Please select from the following or create your own:"
           print *,"SetSerialParamsPrimExplicit:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
           print *,"SetSerialParamsPrimExplicit:  BGL-fast    == BGL processor clock scaled based on Power3 "
           print *,"SetSerialParamsPrimExplicit:  BGL-slow    == BGL processor clock scaled based on Power4 "
       return
         endif


    if(m_found) then
      if(np_found) then
        if(nlev_found) found=.TRUE.
      endif
    endif




    end subroutine SetSerialParamsPrimExplicit

    subroutine SetSerialParamsImplicit(microsec,microsec_per_iter,cg_iters,np,nelem,nlev,machine,found)

       real (kind=real_kind),intent(out)   :: microsec
       real (kind=real_kind),intent(out)   :: microsec_per_iter
       real (kind=real_kind),intent(out)   :: cg_iters
       integer(kind=int_kind),intent(in)   :: np
       integer(kind=int_kind),intent(in)   :: nelem
       integer(kind=int_kind),intent(in)   :: nlev
       character(len=80)                   :: machine
       logical(kind=log_kind),intent(inout)  :: found
       logical(kind=log_kind)                :: m_found,np_found,nlev_found,nelem_found

       print *,'SetSerialParamsImplicit:  Note all values of broken... please fix'
       m_found    = .FALSE.
       np_found   = .FALSE.
       nelem_found   = .FALSE.
       nlev_found = .FALSE.
       if(machine(1:7) == "bluesky") then
          m_found = .TRUE.
          if(np .eq. 6 ) then
             np_found=.TRUE.
             select case (nlev)
             case(96)
!                   microsec         = 4243.D0 ! SMP numbers (nproc=8)
                microsec         = 3299.D0 ! single threaded
                nlev_found=.TRUE.
             case(60)
!                 microsec         = 2545.D0 ! SMP number (nproc=8)
!                 microsec         = 2194.D0 ! SMP number (nproc=2)
                microsec         = 2013.D0 ! single threaded
                nlev_found=.TRUE.
             case(48)
!                 microsec         = 1914.D0 ! SMP number (nproc=8)
                microsec         = 1702.D0 ! single threaded
                nlev_found=.TRUE.
             case(30)
!                 microsec         = 1142.D0 ! SMP number (nproc=8)
                microsec         = 1023.D0 ! single threaded
                nlev_found=.TRUE.
             case(16)
                microsec         =  540.D0
                nlev_found=.TRUE.
             case(2)
                microsec         =   76.D0
                nlev_found=.TRUE.
             case default
                print *,"SetSerialParamsImplicit:  Could not find performance estimate for: "
                print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
             end select
          elseif (np .eq. 8) then
             np_found = .TRUE.
             select case (nlev)
             case(16)
                microsec          = 800._real_kind
                microsec_per_iter = 450._real_kind
                nlev_found = .TRUE.
             case default
                print *,"SetSerialParamsImplicit: Could not find performance estimate for: "
                print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
             end select
          elseif (np .eq. 10) then
             np_found = .TRUE.
             select case (nlev)
             case(16)
                microsec          = 1303._real_kind
                microsec_per_iter =  793._real_kind
                nlev_found = .TRUE.
             case(96)
                microsec        = 9570._real_kind ! single threaded
                nlev_found = .FALSE.
             case default
                print *,"SetSerialParamsImplicit:  Could not find performance estimate for: "
                print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
             end select
          else
             print *,"SetSerialParamsImplicit:  Could not find performance estimate for: "
             print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np
          endif
       elseif(machine(1:11) == "blackforest") then
          m_found = .TRUE.
          if (np .eq. 8) then
             np_found = .TRUE.
             select case (nlev)
             case(16)
                microsec          = 2077._real_kind
                microsec_per_iter = 1252._real_kind
                nlev_found = .TRUE.
             case default
                print *,"SetSerialParamsImplicit: Could not find performance estimate for: "
                print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
             end select
          elseif( np .eq. 10 ) then
             np_found = .TRUE.
             select case (nlev)
             case(16)
                microsec          = 3200._real_kind
                microsec_per_iter = 2320._real_kind
                nlev_found = .TRUE.
             case(96)
                nlev_found = .FALSE.
             case default
                print *,"SetSerialParamsImplicit: Could not find performance estimate for: "
                print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np,"nlev: ",nlev
             end select
          else
             print *,"SetSerialParamsImplicit: Could not find performance estimate for: "
             print *,"SetSerialParamsImplicit:  machine: ",machine(1:11)," np := ",np
          endif
       elseif(machine(1:8) == "protoBGL") then
          m_found = .TRUE.
          if (np .eq. 8) then
             np_found = .TRUE.
             select case (nlev)
             case(16)
                microsec        = 9187.4_real_kind ! single threaded
                cg_iters = 3.5528_real_kind
                nlev_found = .TRUE.
             case default
                print *,"SetSerialParamsImplicit: Could not find performance estimate for: "
                print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
             end select
          elseif( np .eq. 10 ) then
             np_found = .TRUE.
             select case (nlev)
             case(96)
                microsec        = (900._real_kind/210._real_kind)*9570._real_kind ! This is an estimate
                nlev_found = .TRUE.
             case default
                print *,"SetSerialParamsImplicit: Could not find performance estimate for: "
                print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
             end select
          else
             print *,"SetSerialParamsImplicit: Could not find performance estimate for: "
             print *,"SetSerialParamsImplicit:  machine: ",machine(1:8)," np := ",np
          endif
       elseif(machine(1:8) == "BGL-fast") then
          m_found = .TRUE.
          if (np .eq. 10) then
             np_found = .TRUE.
             select case (nlev)
             case(96)
                ! This is an estimate
                microsec   = (375._real_kind/700._real_kind) * 24218._real_kind ! single threaded
                nlev_found = .TRUE.
             case default
                print *,"SetSerialParams: Could not find performance estimate for: "
                print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
             end select
          else
             print *,"SetSerialParams: Could not find performance estimate for: "
             print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np
          endif
       elseif(machine(1:8) == "BGL-slow") then
          m_found = .TRUE.
          if (np .eq. 10) then
             np_found = .TRUE.
             select case (nlev)
             case(96)
                ! This is an estimate
                microsec        = (1300._real_kind/700._real_kind) * 9570._real_kind
                nlev_found = .TRUE.
             case default
                print *,"SetSerialParams: Could not find performance estimate for: "
                print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np,"nlev: ",nlev
             end select
          else
             print *,"SetSerialParams: Could not find performance estimate for: "
             print *,"SetSerialParams:  machine: ",machine(1:8)," np := ",np
          endif
       else
          print *,"SetSerialParams:  Undefined machine name"
          print *,"SetSerialParams:  Please select from the following or create your own:"
          print *,"SetSerialParams:  bluesky     == Power4 with dual rail colony and 8-way LPARs"
          print *,"SetSerialParams:  blackforest == Power3 with TBMX network"
          print *,"SetSerialParams:  protoBGL    == Prototype BGL machine (November 2003 configuration)"
          print *,"SetSerialParams:  BGL-fast    == BGL processor clock scaled based on Power3 "
          print *,"SetSerialParams:  BGL-slow    == BGL processor clock scaled based on Power4 "
          return
       endif

       ! =====================================
       !   Set the expected iteration count
       ! =====================================
       select case(np)
       case(8)
          select case (nelem)
          case(24)
             cg_iters = 2.1042_real_kind
             nelem_found = .TRUE.
          case(54)
             cg_iters = 2.9722_real_kind
             nelem_found = .TRUE.
          case(96)
             cg_iters = 3.2778_real_kind
             nelem_found = .TRUE.
          case(384)
             cg_iters = 3.5556_real_kind
             nelem_found = .TRUE.
          case(864)
             cg_iters = 3.9708_real_kind
             nelem_found = .TRUE.
          case default
             print *,"SetSerialParamsImplicit: Could not iteration count estimate for: "
             print *,'SetSerialParamsImplicit: np := ',np, 'nelem := ',nelem
             nelem_found = .FALSE.
          end select
       case(10)
          select case (nelem)
          case(24)
             cg_iters = 2.9583_real_kind
             nelem_found = .TRUE.
          case(54)
             cg_iters = 3.4427_real_kind
             nelem_found = .TRUE.
          case(96)
             cg_iters = 3.6958_real_kind
             nelem_found = .TRUE.
          case(384)
             cg_iters = 3.7979_real_kind
             nelem_found = .TRUE.
          case(1536)
             cg_iters = 3.9708_real_kind
             nelem_found = .TRUE.
          case(55296)
             cg_iters = 5.8511_real_kind
             nelem_found = .TRUE.
          case default
             print *,"SetSerialParamsImplicit: Could not iteration count estimate for: "
             print *,'SetSerialParamsImplicit: np := ',np, 'nelem := ',nelem
             nelem_found = .FALSE.
          end select
       case default
          print *,"SetSerialParamsImplicit: Could not iteration count estimate for: "
          print *,'SetSerialParamsImplicit: np := ',np, 'nelem := ',nelem
       end select


       if(m_found) then
          if(np_found) then
             if(nelem_found) then
                if(nlev_found) found=.TRUE.
             endif
          endif
       endif


    end subroutine SetSerialParamsImplicit

end module perfmodel_mod
