#include "fabm_driver.h"

module ersem_benthic_fauna

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

  type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_fauna
  type (type_state_variable_id) :: id_O2o
  type (type_dependency_id), allocatable,dimension(:) :: id_foodc,id_foodn,id_foodp,id_foods
  type (type_dependency_id) :: id_ETW    
     real(rk) :: qnYIcX,qpYIcX
     real(rk) :: q10YX
     real(rk) :: hO2YX,rlO2YX
     real(rk) :: xclYX,xcsYX,xchYX
     real(rk) :: suYX,luYX,huYX
     real(rk),allocatable :: pueYX(:),pufood(:)
     real(rk) :: pudilX
  contains
     procedure :: initialize
     procedure :: do_bottom
  end type
contains

  subroutine initialize(self,configunit)
    class (type_ersem_benthic_fauna),intent(inout),target :: self
    integer,                                 intent(in)           ::configunit
    integer           :: ifood
    character(len=16) :: index
    logical           :: foodispom
    real(rk)          :: pueYX, pueQX
    ! Register parameters
      call self%get_parameter(self%qnYIcX,  'qnYc',  'mmol N/mg C','Maximum nitrogen to carbon ratio')
      call self%get_parameter(self%qpYIcX,  'qpYc',  'mmol P/mg C','Maximum phosphorus to carbon ratio')
      call self%get_parameter(self%q10YX,   'q10Y',  '-',          'Regulating temperature factor Q10')
      call self%get_parameter(self%hO2YX,   'hO2Y',  'mmol/m^3',   'Michaelis-Menten constant for oxygen limitation')
      call self%get_parameter(self%rlO2YX,  'rlO2Y', 'mmol/m^3',   'Minimum treshold of oxygen required')
      call self%get_parameter(self%xclYX,   'xclY',  'mg C/m^2',   'Lower treshold for crowding effect on food uptake')
      call self%get_parameter(self%xcsYX,   'xcsY',  'mg C/m^2',   'Michaelis-Menten constant for the impact of crowding')
      call self%get_parameter(self%xchYX,   'xchY',  'mg C/m^2',   'Concentration determining asymptotic treshold of shading limitation (-> xchYXi/(1+xchYXi) for Yc-> inf)')
      call self%get_parameter(self%suYX,    'suY',   '1/d',        'Specific maximal uptake at reference temperature')
      call self%get_parameter(self%luYX,    'luY',   'mg C/m^2',   'Michaelis-Menten constant for food species uptake')
      call self%get_parameter(self%huYX,    'huY',   'mg C/m^2',   'Michaelis-Menten constant for gross carbon uptake')
      call self%get_parameter(pueYX,        'pueY',  '-',          'Excreted fraction of food uptake')
      call self%get_parameter(pueQX,        'pueQ',  '-',          'Excreted fraction of POM uptake')
      call self%get_parameter(pudilX,       'pudil', '-',          'Relative nutrient content in the fecal pellets')
    ! Determine number of food sources
      call self%get_parameter(self%nfood,'nfood','','number of food sources',default=0)   

    allocate(self%id_food(self%nfood))
    do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%register_dependency(self%id_foodc(ifood), 'food'//trim(index)//'c','mmol C m-2','Food '//trim(index)//' C') 
         call self%register_dependency(self%id_foodn(ifood), 'food'//trim(index)//'n','mmol C m-2','Food '//trim(index)//' N')
         call self%register_dependency(self%id_foodp(ifood), 'food'//trim(index)//'p','mmol C m-2','Food '//trim(index)//' P') 
    end do
    allocate(self%pueYX(self%nfood))  
    do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(foodispom,'food'//trim(index)//'ispom','','food type '//trim(index)//' is POM',default=.false.)
         if (foodispom) then
            self%pueY(ifood) = pueQX
         else
            self%pueY(ifood) = pueYX
         end if
    end do
    allocate(self%pufood(self%nfood))
    do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(self%pufood(ifood),'pufood'//trim(index),'-','Food preference of benthic faunal group for food source'//trim(index))
    end do

    ! Environmental dependencies
      call self%register_dependency(self%id_ETW,standard_variables%temperature)

    ! Dependencies on state variables of external modules.
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O/m^3','oxygen')
  end subroutine

  subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

     class (type_ersem_benthic_fauna),intent(in) :: self
     _DECLARE_ARGUMENTS_DO_BOTTOM_
     
     real(rk) :: Yc,Yn,Yp,O2o
     real(rk) :: eT,eO,eC,ETW,Y,x
     real(rk) :: rate
     real(rk) :: AQ6c,AQ6n,AQ6p,AQ6c2
     real(rk),dimension(self%nfood) :: foodcP
     real(rk),dimension(self%nfood) :: feed, sflux
     real(rk) :: foodsum,mm
     real(rk),dimension(self%nfood) :: grossfluxc,grossfluxn,grossfluxp
     real(rk),dimension(self%nfood) :: netfluxc,netfluxn,netfluxp
  
     _HORIZONTAL_LOOP_BEGIN_

     _GET_HORIZONTAL_(self%id_c,Yc)

     _GET_(self%id_O2o,O2o)
     _GET_HORIZONTAL_(self%id_AQ6c,AQ6c)
     _GET_HORIZONTAL_(self%id_AQ6n,AQ6n)
     _GET_HORIZONTAL_(self%id_AQ6p,AQ6p)
     _GET_HORIZONTAL_(self%id_AQ6c2,AQ6c2)

     Yn = Yc*self%qnYIcX
     Yp = Yc*self%qpYIcX
     
     ! Calculate temperature limitation factor
     eT = q10YX**((ETW-10._rk)/10._rk) - q10YX**((ETW-32._rk)/3._rk)

     ! Calculate oxygen limitation factor
     eO = (O2o-self%rlO2YX)**3/((O2o-self%rlO2YX)**3+(self%hO2YX-self%rlO2YX)**3)
     
     ! Calculate overcrowding limitation factor
     Y = Yc - self%xclYX
     if ( Y .gt. 0._rk ) then
       x = Y * Y/(Y+self%xcsYX)
       eC = 1._rk - x/(x+self%xchYX)
     else
       eC = 1._rk
     end if

    ! Calculate uptake rate................................................
    rate = self%suYX * Yc * eT * eO * eC
 
    !!!Need to get organic matter from certain horizons

    !Get food concentrations !!!!!Benthic and PELAGIC!
    do ifood=1,self%nfood
       _GET_(self%id_foodc(ifood),foodcP(ifood))
       _GET_(self%id_foodn(ifood),foodnP(ifood))
       _GET_(self%id_foodp(ifood),foodpP(ifood))
    end do

    ! Food Partition
    feed = self%pufood * (self%pufood * foodcP/(self%pufood * foodcP + self%luY))
    foodsum = sum(feed * foodcP)
    mm = foodsum + self%huYX

    if (foodsum .gt. 0._rk) then
     sflux = rate * feed / mm
      else
     sflux = 0._rk
    end if

   ! Food uptake fluxes

   grossfluxc = sflux * foodcP
   grossfluxn = sflux * foodnP
   grossfluxp = sflux * foodpP
   netfluxc = grossfluxc * (1._rk -             self%pueY(ifood))
   netfluxn = grossfluxn * (1._rk - self%pudilX*self%pueY(ifood))
   netfluxp = grossfluxp * (1._rk - self%pudilX*self%pueY(ifood))
     
  
     _HORIZONTAL_LOOP_END_

  end subroutine do_bottom

end module

