#include "fabm_driver.h"

module ersem_benthic_fauna

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base

   implicit none

   private

  type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_fauna
  type (type_state_variable_id)   :: id_O2o
  type (type_bottom_state_variable_id) :: id_Q6c,id_Q6n,id_Q6p,id_Q6s,id_Q6c2
  type (type_bottom_state_variable_id) :: id_G3c,id_G2o,id_K4n,id_K1p,id_K4n2,id_K1p2
  type (type_horizontal_dependency_id), allocatable,dimension(:) :: id_foodc,id_foodn,id_foodp,id_foods
  type (type_dependency_id), allocatable,dimension(:) :: id_foodpelc,id_foodpeln,id_foodpelp,id_foodpels
  type (type_dependency_id) :: id_ETW
  type (type_model_id),allocatable,dimension(:) :: id_food
     integer  :: nfood    
     real(rk) :: qnYIcX,qpYIcX
     real(rk) :: q10YX
     real(rk) :: hO2YX,rlO2YX
     real(rk) :: xclYX,xcsYX,xchYX
     real(rk) :: suYX,luYX,huYX
     real(rk),allocatable :: pueYX(:),pufood(:),pu_anX(:)
     logical,allocatable ::foodispel(:)
     real(rk) :: pudilX
     real(rk) :: sdYX,sdmO2YX,sdcYX,xdcYX
     real(rk) :: srYX,purYX
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
    logical           :: foodispom,food_anaerobic,foodispel
    real(rk)          :: pueYX,pueQX
    self%dt = 86400._rk
    ! Register parameters
      call self%get_parameter(self%qnYIcX, 'qnYc',  'mmol N/mg C','Maximum nitrogen to carbon ratio')
      call self%get_parameter(self%qpYIcX, 'qpYc',  'mmol P/mg C','Maximum phosphorus to carbon ratio')
      call self%get_parameter(self%q10YX,  'q10Y',  '-',          'Regulating temperature factor Q10')
      call self%get_parameter(self%hO2YX,  'hO2Y',  'mmol/m^3',   'Michaelis-Menten constant for oxygen limitation')
      call self%get_parameter(self%rlO2YX, 'rlO2Y', 'mmol/m^3',   'Minimum treshold of oxygen required')
      call self%get_parameter(self%xclYX,  'xclY',  'mg C/m^2',   'Lower treshold for crowding effect on food uptake')
      call self%get_parameter(self%xcsYX,  'xcsY',  'mg C/m^2',   'Michaelis-Menten constant for the impact of crowding')
      call self%get_parameter(self%xchYX,  'xchY',  'mg C/m^2',   'Concentration determining asymptotic treshold of shading limitation (-> xchYXi/(1+xchYXi) for Yc-> inf)')
      call self%get_parameter(self%suYX,   'suY',   '1/d',        'Specific maximal uptake at reference temperature')
      call self%get_parameter(self%luYX,   'luY',   'mg C/m^2',   'Michaelis-Menten constant for food species uptake')
      call self%get_parameter(self%huYX,   'huY',   'mg C/m^2',   'Michaelis-Menten constant for gross carbon uptake')
      call self%get_parameter(pueYX,       'pueY',  '-',          'Excreted fraction of food uptake')
      call self%get_parameter(pueQX,       'pueQ',  '-',          'Excreted fraction of POM uptake')
      call self%get_parameter(self%pudilX, 'pudil', '-',          'Relative nutrient content in the fecal pellets')
      call self%get_parameter(self%sdYX,   'sdY',   '1/d',        'Specific mortality at reference temperature')
      call self%get_parameter(self%sdmO2YX,'sdmO2Y','1/d',        'Specific maximal additional mortality of deposit feeders due to oxygen stress')
      call self%get_parameter(self%sdcYX,  'sdcY',  '1/d',        'Specific maximal additional mortality of deposit feeders induced by cold temperature')
      call self%get_parameter(self%xdcYX,  'xdcY',  '1/degree C', 'e-folding temperature factor of cold mortality response')
    ! Determine number of food sources
      call self%get_parameter(self%nfood,  'nfood', '',           'Number of food sources',default=0)   
      call self%get_parameter(self%srYX,   'srY',   '1/d',        'Specific rest respiration at reference temperature')
      call self%get_parameter(self%purYX,  'purY',  '-',          'Respired fraction of uptake')
   
      call self%add_constituent('c',3000._rk,self%qnYIcX,self%qpYIcX)
  
    allocate(self%id_food(self%nfood))
    allocate(self%foodispel(self%nfood))
    allocate(self%id_foodpelc(self%nfood))
    allocate(self%id_foodpeln(self%nfood))
    allocate(self%id_foodpelp(self%nfood))
    allocate(self%id_foodpels(self%nfood))
    allocate(self%id_foodc(self%nfood))
    allocate(self%id_foodn(self%nfood))
    allocate(self%id_foodp(self%nfood))
    allocate(self%id_foods(self%nfood))

   ! Allocate components of food sources
    do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(self%foodispel(ifood),'food'//trim(index)//'ispel','','food type '//trim(index)//'is pelagic',default=.false.)
         call self%register_model_dependency(self%id_food(ifood),'food'//trim(index))
         if (foodispel) then
         call self%register_dependency(self%id_foodpelc(ifood), 'food'//trim(index)//'c','mmol C m-2','Food '//trim(index)//' C') 
         call self%register_dependency(self%id_foodpeln(ifood), 'food'//trim(index)//'n','mmol C m-2','Food '//trim(index)//' N')
         call self%register_dependency(self%id_foodpelp(ifood), 'food'//trim(index)//'p','mmol C m-2','Food '//trim(index)//' P')
         call self%register_dependency(self%id_foodpels(ifood), 'food'//trim(index)//'s','mmol C m-2','Food '//trim(index)//' Si')
         call self%request_coupling_to_model(self%id_foodpelc(ifood),self%id_food(ifood),standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_foodpeln(ifood),self%id_food(ifood),standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_foodpelp(ifood),self%id_food(ifood),standard_variables%total_phosphorus)
         call self%request_coupling_to_model(self%id_foodpels(ifood),self%id_food(ifood),standard_variables%total_silicate)
           else
         call self%register_dependency(self%id_foodc(ifood),'food'//trim(index)//'c','mmol C m-2','Food '//trim(index)//' C') 
         call self%register_dependency(self%id_foodn(ifood),'food'//trim(index)//'n','mmol C m-2','Food '//trim(index)//' N')
         call self%register_dependency(self%id_foodp(ifood),'food'//trim(index)//'p','mmol C m-2','Food '//trim(index)//' P')
         call self%register_dependency(self%id_foods(ifood),'food'//trim(index)//'s','mmol C m-2','Food '//trim(index)//' Si')
         call self%request_coupling_to_model(self%id_foodc(ifood),self%id_food(ifood),standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_foodn(ifood),self%id_food(ifood),standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_foodp(ifood),self%id_food(ifood),standard_variables%total_phosphorus)
         call self%request_coupling_to_model(self%id_foods(ifood),self%id_food(ifood),standard_variables%total_silicate)
         end if 
   end do

    ! Allocate excretion rates
    allocate(self%pueYX(self%nfood))
    allocate(self%pu_anX(self%nfood))  
    do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(foodispom,'food'//trim(index)//'ispom','','food type '//trim(index)//' is POM',default=.false.)
         if (foodispom) then
            self%pueYX(ifood) = pueQX
         else
            self%pueYX(ifood) = pueYX
         end if
         call self%get_parameter(food_anaerobic,'food'//trim(index)//'anaerobic','','food type '//trim(index)//' is from anaerobic layer',default=.false.)
        if (food_anaerobic) then
           self%pu_anX(ifood) = 1._rk
        else
           self%pu_anX(ifood) = 0._rk
        end if
    end do

    ! Allocate food source preferences
    allocate(self%pufood(self%nfood))
    do ifood=1,self%nfood
         write (index,'(i0)') ifood
         call self%get_parameter(self%pufood(ifood),'pufood'//trim(index),'-','Food preference of benthic faunal group for food source'//trim(index))
    end do

    ! Environmental dependencies
      call self%register_dependency(self%id_ETW,standard_variables%temperature)

    ! Dependencies on state variables of external modules
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O/m^3','oxygen')
      call self%register_state_dependency(self%id_Q6c,'Q6c','mg C/m^2', 'particulate organic carbon')
      call self%register_state_dependency(self%id_Q6n,'Q6n','mmol N/m^2','particulate organic nitrogen')
      call self%register_state_dependency(self%id_Q6p,'Q6p','mmol P/m^2','particulate organic phosphorus')
      call self%register_state_dependency(self%id_Q6s,'Q6s','mmol Si/m^2','particulate organic silicate')
      call self%register_state_dependency(self%id_G3c,'G3c','mmol C m-2','Carbon Dioxide')
      call self%register_state_dependency(self%id_G2o,'G2o','mmol O m-2','Oxygen')
      call self%register_state_dependency(self%id_K1p2,'K1p2','mmol P/m^2','benthic phosphate in 2nd layer')
      call self%register_state_dependency(self%id_K4n2,'K4n2','mmol N/m^2','benthic ammonium in 2nd layer')
      call self%register_state_dependency(self%id_K1p,'K1p','mmol P/m^2','benthic phosphate in 1st layer')
      call self%register_state_dependency(self%id_K4n,'K4n','mmol N/m^2','benthic ammonium in 1st layer')
 
  end subroutine

  subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

     class (type_ersem_benthic_fauna),intent(in) :: self
     _DECLARE_ARGUMENTS_DO_BOTTOM_
     
     real(rk) :: Yc,Yn,Yp,O2o,foodP
     real(rk) :: eT,eO,eC,ETW,Y,x
     real(rk) :: rate
     real(rk) :: SQ6c
     real(rk),dimension(self%nfood) :: foodcP,foodnP,foodpP,foodsP
     real(rk),dimension(self%nfood) :: feed, sflux
     real(rk) :: foodsum,mm
     real(rk),dimension(self%nfood) :: grossfluxc,grossfluxn,grossfluxp,grossfluxs
     real(rk),dimension(self%nfood) :: netfluxc,netfluxn,netfluxp
     real(rk) :: mort_act,mort_ox,mort_cold,mortflux
     integer  :: ifood,istate
     real(rk) :: fBTYc,nfBTYc,fYG3c,p_an,excn,excp
     
     _HORIZONTAL_LOOP_BEGIN_

     _GET_HORIZONTAL_(self%id_c,Yc)
     _GET_(self%id_O2o,O2o)


     Yn = Yc*self%qnYIcX
     Yp = Yc*self%qpYIcX
     
     ! Calculate temperature limitation factor
     eT = self%q10YX**((ETW-10._rk)/10._rk) - self%q10YX**((ETW-32._rk)/3._rk)

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
       if (self%foodispel(ifood)) then
       _GET_(self%id_foodpelc(ifood),foodcP(ifood))
       _GET_(self%id_foodpeln(ifood),foodnP(ifood))
       _GET_(self%id_foodpelp(ifood),foodpP(ifood))
       _GET_(self%id_foodpels(ifood),foodsP(ifood))
       else
       _GET_HORIZONTAL_(self%id_foodc(ifood),foodcP(ifood))
       _GET_HORIZONTAL_(self%id_foodn(ifood),foodnP(ifood))
       _GET_HORIZONTAL_(self%id_foodp(ifood),foodpP(ifood))
       _GET_HORIZONTAL_(self%id_foods(ifood),foodsP(ifood))
       end if
    end do

    ! Food Partition
    feed = self%pufood * (self%pufood * foodcP/(self%pufood * foodcP + self%luYX))
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
   grossfluxs = sflux * foodsP
   netfluxc = grossfluxc * (1._rk -             self%pueYX(ifood))
   netfluxn = grossfluxn * (1._rk - self%pudilX*self%pueYX(ifood))
   netfluxp = grossfluxp * (1._rk - self%pudilX*self%pueYX(ifood))

   do ifood=1,self%nfood
      if (self%foodispel(ifood)) then
       do istate=1,size(self%id_food(ifood)%state)
         _GET_(self%id_food(ifood)%state(istate),foodP)
         _SET_ODE_(self%id_food(ifood)%state(istate),-sflux(ifood)*foodP)   
       end do
      else
       do istate=1,size(self%id_food(ifood)%bottom_state)
         _GET_HORIZONTAL_(self%id_food(ifood)%bottom_state(istate),foodP)
         _SET_BOTTOM_ODE_(self%id_food(ifood)%bottom_state(istate),-sflux(ifood)*foodP)
       end do
      end if 
   end do
   
   fBTYc = sum(grossfluxc)
   nfBTYc= sum(netfluxc)
   
   _SET_BOTTOM_ODE_(self%id_c,nfBTYc)
   _SET_BOTTOM_ODE_(self%id_Q6c,fBTYc - nfBTYc)
   _SET_BOTTOM_ODE_(self%id_Q6n,sum(grossfluxn) - sum(netfluxn))
   _SET_BOTTOM_ODE_(self%id_Q6p,sum(grossfluxp) - sum(netfluxp))
   _SET_BOTTOM_ODE_(self%id_Q6s,sum(grossfluxs))

   ! Respiration fluxes = activity respiration + basal respiration

   fYG3c = self%srYX * Yc * eT + self%purYX * nfBTYc
   
   _SET_BOTTOM_ODE_(self%id_c,-fYG3c)
   _SET_BOTTOM_ODE_(self%id_G3c, fYG3c/CMass)
   _SET_BOTTOM_ODE_(self%id_G2o,-fYG3c/CMass)

   ! Mortality fluxes     
   mort_act  = self%sdYX * eT
   mort_ox   = self%sdmO2YX * (1._rk - eO) * eT
   mort_cold = self%sdcYX * exp(ETW * self%xdcYX)

   mortflux = mort_act + mort_ox + mort_cold
   
   _SET_BOTTOM_ODE_(self%id_c, -mortflux * Yc)
   _SET_BOTTOM_ODE_(self%id_Q6c,mortflux * Yc)
   _SET_BOTTOM_ODE_(self%id_Q6n,mortflux * Yn)
   _SET_BOTTOM_ODE_(self%id_Q6p,mortflux * Yp)

   ! Adjust fixed nutrients
   call Adjust_fixed_nutrients(Yc,Yn,Yp,self%qnYIcX,self%qpYIcX,excn,excp,SQ6c)

    p_an = (sum(self%pu_anX*grossfluxc))/max(fBTYc,1.e-8_rk)

   _SET_BOTTOM_ODE_(self%id_K4n,(1._rk-p_an) * excn)
   _SET_BOTTOM_ODE_(self%id_K4n2,p_an * excn)
   _SET_BOTTOM_ODE_(self%id_K1p,(1._rk-p_an) * excp)
   _SET_BOTTOM_ODE_(self%id_K1p2,p_an * excp)
   _SET_BOTTOM_ODE_(self%id_Q6c,SQ6c)
     _HORIZONTAL_LOOP_END_

  end subroutine do_bottom

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adjust_fixed_nutrients \label{sec:AdjustFixedNutrients}
!
! !DESCRIPTION:
!  TODO - description
!
!  This routine determines the amount of excess c,n or p in the
!  source terms and excretes the appropriate amount(s) to Q6 and
!  nutrients, so that the fixed nutrient ratio is re-established
!
!  IN:  Source terms of fixed quota bio-state.(cnp)......SXx
!       Fixed quota values (np)..........................qx
!
!  INT: Excess nutrient in bio-state (np)................ExcessX
!
!  OUT: Source terms for Predator (SX), Q6 and nutrients
!\\
!\\
! !INTERFACE:
      subroutine adjust_fixed_nutrients ( SXc, SXn, SXp,qn, qp, SKn, SKp, SQc )
!
! !INPUT/OUTPUT PARAMETERS:
       real(rk), intent(inout)  :: SXc, SXn, SXp, SKn, SKp, SQc
       real(rk), intent(in)     :: qn, qp
!
! !LOCAL PARAMETERS:
       real(rk) :: ExcessN, ExcessP, ExcessC
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
       ExcessC = max(max(SXc - SXp/qp,SXc - SXn/qn),0._rk)
!
       if (ExcessC>0.0_rk) then
         SXc = SXc - ExcessC
         SQc = SQc + ExcessC
       end if

       ExcessN = max(SXn - SXc*qn,0._rk)
       ExcessP = max(SXp - SXc*qp,0._rk)

       if (ExcessN>0.0_rk) then
         SXn = SXn - ExcessN
         SKn = SKn + ExcessN
       end if

       if (ExcessP>0.0_rk) then
         SXp = SXp - ExcessP
         SKp = SKp + ExcessP
       end if

       end subroutine adjust_fixed_nutrients
!
!EOC
!-----------------------------------------------------------------------

end module
