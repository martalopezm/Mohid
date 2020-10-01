!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Carbonate System
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Feb 2020
! REVISION      : Marta López Mozos - v1.0
! DESCRIPTION   : Zero-dimensional model for marine carbonate system                  
!                 
!------------------------------------------------------------------------------
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License 
!  version 2, as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


    
!******************************************************************************************************
!                               .DAT FILE EXAMPLE
!******************************************************************************************************
!
!    
!    
!    
!    
!    
!   
!
!    
!
!    
!    
!******************************************************************************************************
!                               EXTENDED DESCRIPTION
!******************************************************************************************************
!         
!      The carbonate system module calculates the alkalinity and the dissolved organic carbon 
!      to compute, along with other provided properties, the pH and the rest of the carbonate 
!      system parameters.
!    
!      Alkalinity and DIC are calculated from other MOHID modules, while the rest of the carbonate 
!      system parameters are calculated by mocsy/CO2sys (insert license), which is coupled to MOHID 
!      inside this module.     
!     
!      Two options to compute alkalinity are allowed.
!      
!      
!      The module follows the others MOHID's modules statements structure: 
!      
!       -  Privaciy of subroutines
!       -  Declaration: parameters, variables, types and global module variables
!       -  Constructor: first to be called. Prepares all to be ready for the calculation;
!                       data files are read and matrixes are created and initialized
!       -  Selector   : allows other modules to get the matrixes values from  this Module
!       -  Modifier   : calculations and calculations subroutines; called each time the time changes
!                       and is resposible for updating the matrixes for the new time step.     
!       -  Destructor : where all the matrixes are released and connections between modules destroyed
!       -  Management    
!      
!      
!        
!*************************************************************************************************
!                                KEYWORDS AND UNITS
!*************************************************************************************************      
!      
!     
!          ALKsal            = Alkalinity from salinity              ( umol / kgSW )
!          ALKbio            = Biological alkalinity                 ( umol / kgSW )
!          DIC               = Dissolved Inorganic Carbon            ( umol  /kgSW ) | (mg C m-3)
!
!*************************************************************************************************
!*************************************************************************************************   
    
      
 Module ModuleCarbonateSystem
    
      use ModuleGlobalData
      use ModuleLUD
      use ModuleEnterData
      !use ModuleFunctions, only: OxygenSaturation, PhytoLightLimitationFactor
      use ModuleTime
     
      implicit none              

      private      
      
      
!******************************************************************************
!                           SUBROUTINES PRIVACY
!******************************************************************************
      
     !Constructor      
      public  :: ConstructCarbonateSystem
      private ::      AllocateInstance
      private ::      ReadData
      private ::          ConstructGlobalVariables
  !    private ::          ! ConstructProducers, Species ? 
      
      
    
     !Selector          
      !public  :: GetDTCarbonateSystem
      !public  :: GetCarbonateSystemPropertyList
      !public  :: GetCarbonateSystemSize
      !public  :: GetCarbonateSystemPropIndex
      !public  :: UnGetCarbonateSystem
    
    !Modifier         
      !public  :: ModifyCarbonateSystem
      !private :: ComputeAlkalinity_salt
      !private :: ComputeAlkalinity_bio
      !private :: ComputeDIC
      !private :: mocsy
      !private ::   mocsy_1
    
    !Destructor        
      !public  :: KillCarbonateSystem   
      !private :: DeAllocateInstance
    
    !Management     
      private :: Ready
      private :: LocateObjCarbonateSystem

 
    
!*******************************************************************************************
!                              DECLARATION STATEMENTS
!*******************************************************************************************    
    
        
 
! ------ Types, variables 

    private :: T_AuxiliarParameters
    type       T_AuxiliarParameters       
        real                             :: C_AtomicMass          = 12                   !mgC/molC
        real                             :: H_AtomicMass          = 1                    !mgH/molH
        real                             :: O_AtomicMass          = 16                   !mgO/molO
        real                             :: N_AtomicMass          = 14                   !mgN/molN
        real                             :: P_AtomicMass          = 31                   !mgP/molP
    end type T_AuxiliarParameters



    private :: T_BioChemParam
    type       T_BioChemParam   
        real                             :: Redfield_NC            = null_real          !Redfield N:C ratio
        real                             :: Redfield_PC            = null_real          !Redfield P:C ratio
        real                             :: Redfield_SiC           = null_real          !Standard Si:C ratio
        real                             :: BioSi_Diss             = null_real          !Biogenic silica dissolution rate
        real                             :: O2C_Conversion         = null_real          !Oxygen to carbon conversion factor
        real                             :: CO2C_Conversion        = null_real          !Carbon Dioxide to carbon conversion factor
    end type T_BioChemParam



    private :: T_ComputeOptions
    type       T_ComputeOptions
        logical                          :: NoBiologicalAlkalinity    = .true.          !Compute alkalinity by salinity
        logical                          :: BiologicalAlkalinity      = .false.         !Compute biological alkalinity
    end type T_ComputeOptions





   private ::  T_ExternalVar
    type       T_ExternalVar
    end type T_ExternalVar





    private :: T_ID
    type       T_ID         
        integer                          :: ID, IDNumber             = null_int
        character(len=StringLength)      :: Name                     = null_str
        character(len=StringLength)      :: Description              = null_str
    end type T_ID



   



    private :: T_PropIndex
      type     T_PropIndex
          
      end type T_PropIndex       
      


      private :: T_CarbonateSystem
      type       T_CarbonateSystem

        integer                                         :: InstanceID
        integer                                         :: ObjEnterData         = 0
        integer, dimension(:), pointer                  :: PropertyList         => null()
        real                                            :: DT
        real                                            :: DT_day
        type(T_AuxiliarParameters)                      :: AuxParam
        type(T_BioChemParam      )                      :: BioChemPar
        type(T_ComputeOptions    )                      :: ComputeOptions
        type(T_ExternalVar       )                      :: ExternalVar
        type(T_Size1D            )                      :: Array
        type(T_Size1D            )                      :: Prop
        type(T_PropIndex         )                      :: PropIndex
        type(T_CarbonateSystem   ), pointer             :: Next

      end type T_CarbonateSystem

  !real   :: ALKsal                = null_real      ! Alkalinity concentration, computed throught a salinity parametrization    ( umol  /kgSW ) 
  !real   :: ALKbio                = null_real      ! Alkalinity concentration, derived from alkalinity variations due to...    ( umol  /kgSW ) 
                                                                                                 ! ...biological processes
  !real   :: DIC                   = null_real      ! Dissolved inorganic carbon concentration                                  ( umol  /kgSW ) | (mg C m-3)                                         
  !   = FillValueReal              !FillValueReal valor d einicializacion de variables, comun en MOHID para ello
                   
 ! ------ Global Module Variables ----------------------   
      
       type (T_CarbonateSystem), pointer                        :: FirstObjCarbonateSystem  => null()
       type (T_CarbonateSystem), pointer                        :: Me                       => null()
      
    contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine ConstructCarbonateSystem(ObjCarbonateSystemID, FileName, STAT)
      
       !Arguments
        integer                                         :: ObjCarbonateSystemID      
        character(len=*)                                :: FileName
        integer, optional, intent(OUT)                  :: STAT

       !External
        integer                                         :: ready_, STAT_CALL

        !Local
        integer                                         :: STAT_

       !-----------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mCarbonateSystem_)) then
          nullify (FirstObjCarbonateSystem)
          call RegisterModule (mCarbonateSystem_)         
        endif
        
        !Call Ready, which is placed on the management block 
        call Ready(ObjCarbonateSystemID, ready_)  

        !If ready_, which is the return of the sbrtine Ready, is equal to OFF_ERR_, call the following 
        ! subroutines, placed a little lower
        
cd0 : if (ready_ .EQ. OFF_ERR_) then

          call AllocateInstance
          call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructCarbonateSystem - ModuleCarbonateSystem - ERROR #1'
          call ReadData
          call PropertyIndexNumber
          call ConstrucPropertyList
          call KillEnterData(Me%ObjEnterData,STAT = STAT_CALL)   
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructCarbonateSystem - ModuleCarbonateSystem - ERROR #2'
        
          !Returns ID 

          ObjCarbonateSystemID     = Me%InstanceID

          STAT_ = SUCCESS_

      else 

          stop 'ModuleCarbonateSystem - ConstructCarbonateSystem - ERROR #1'

      end if cd0

    
      ! If the optional argument STAT is present in this subroutine (so, the function present returns 
      ! true) assign to STAT_ (previously assigned as SUCCES_ as STAT, to return the value

        if (present(STAT)) STAT = STAT_    


  end subroutine ConstructCarbonateSystem

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine AllocateInstance

  !Local-----------------------------------------------------------------
        type (T_CarbonateSystem), pointer                         :: NewObjCarbonateSystem
        type (T_CarbonateSystem), pointer                         :: PreviousObjCarbonateSystem


        !Allocates new instance
        allocate (NewObjCarbonateSystem)
        nullify  (NewObjCarbonateSystem%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjCarbonateSystem)) then
            FirstObjCarbonateSystem         => NewObjCarbonateSystem
            Me                              => NewObjCarbonateSystem
        else
            PreviousObjCarbonateSystem      => FirstObjCarbonateSystem
            Me                              => FirstObjCarbonateSystem%Next
            do while (associated(Me))
                PreviousObjCarbonateSystem  => Me
                Me                          => Me%Next
            enddo
            Me                              => NewObjCarbonateSystem
            PreviousObjCarbonateSystem%Next => NewObjCarbonateSystem
        endif

        Me%InstanceID = RegisterNewInstance (mCarbonateSystem_)

  end subroutine AllocateInstance

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine ReadData
  
  !Local-----------------------------------------
   
   call ConstructGlobalVariables
   call ConstructModelOptions
   call ConstructCO2SYS_mocsy

  end subroutine ReadData


  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine ConstructGlobalVariables


  
  end subroutine ConstructGlobalVariables 

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
















!*******************************************************************************************
!                            MODIFIER :D :D 
!*******************************************************************************************    
      ! subroutine ComputeAlkalinity_salt()
    
      ! Alksal = 2305. + 53.97 * (SSS - 35.) + 2.74 *((SSS - 35.) ** 2.) - 1.16 *(SST - 20.) &
      !          & - 0.040*((SST - 20.)** 2.)

         
         
         

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

    subroutine Ready (ObjCarbonateSystemID, ready_) 

        !Arguments
        integer                                     :: ObjCarbonateSystemID
        integer                                     :: ready_

      !----------------------------------------------------------------------
       nullify (Me)

cd1:    if (ObjCarbonateSystemID > 0) then
            call LocateObjCarbonateSystem (ObjCarbonateSystemID)
            ready_ = VerifyReadLock (mCarbonateSystem_, Me%InstanceID)
        else
            ready_ = OFF_ERR_

        end if cd1
      !----------------------------------------------------------------------

    end subroutine Ready         
         
 !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::        
        
 !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
   subroutine LocateObjCarbonateSystem(ObjCarbonateSystemID)

        !Arguments-------------------------------------------------------------
        integer                                    :: ObjCarbonateSystemID

        !Local-----------------------------------------------------------------

        Me => FirstObjCarbonateSystem
        do while (associated (Me))
            if (Me%InstanceID == ObjCarbonateSystemID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                          &
            stop 'ModuleCarbonateSystem - LocateObjCarbonateSystem - ERROR #1'

    end subroutine LocateObjCarbonateSystem

 !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 





!----------------------------------------------------------------------------------------------------------

end module ModuleCarbonateSystem

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------




         
!*******************************************************************************************
!                            FUNCTIONS AND SUBROUTINES
!*******************************************************************************************          
         
        ! function LOQUESEA (arg1, arg2)

               ! real                :: LOQUESEA
    !Arguments
       ! real                :: arg1
        !real                :: arg2                     
       
    !Calculations
    
       ! LOQUESEA = arg1 + arg2
 
        ! end function LOQUESEA
         
         
         !subroutine
         
        ! end subroutine 
         
