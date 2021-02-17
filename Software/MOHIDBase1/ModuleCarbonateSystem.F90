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
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!    
!******************************************************************************************************
!                             MODULE EXTENDED DESCRIPTION
!******************************************************************************************************
!         
!      The carbonate system module calculates the alkalinity and the dissolved organic carbon 
!      to compute, along with other provided properties, the pH and the rest of the carbonate 
!      system parameters.
!    
!      Alkalinity and DIC variables are calculated by other MOHID modules outputs, while the rest of the carbonate 
!      system parameters are calculated through mocsy package (Orr et al., 2018)(https://github.com/jamesorr/mocsy), 
!      which is coupled to MOHID inside this module. 
!     
!      Two options to compute alkalinity are allowed (marta:explicar!!).
!      
!      
!      The module follows others MOHID's modules statements structure: 
!      
!       -  Privacity of subroutines (list of all subroutines contained in the module)
!       -  Declaration: parameters, types and global module variables
!       -  Constructor: first part of the module to be called. Prepares all to be ready for the calculation;
!                       data files are read and arrays are created, allocated and initialized
!       -  Selector   : allows other modules to get parameters and variables values from  this module
!       -  Modifier   : perform calculations; called each time step  and is resposible for updating the
!                       arrays for the new time step.     
!       -  Destructor : where all the matrixes are deallocated and connections between modules destroyed
!       -  Management    
!      
!      
!        Following the MOHID structure, this biogeochemical module is dimensionless, this is: 
!        When the property value arrives to the biogeochem. module, it has already experience physical phenomena  
!        (as adv-diffusion processes) and only remains the sink and sources component of change, that will be 
!        computed in this biogchm. module. 
!        The interface module (which links between the waterproperties module and the biogeochemical modules)
!        converts the 3D, 2D or 1D properties arrays (from water property module) into 1D arrays for biogeochemical 
!        modules. After that, changes in the property value, due to biogeochemical processes (sink-sources)         
!        are computed in the biogeochem. module. The resulting 1D array is send back to interface, wherein 
!        is converted to 3D, 2D or 1D properties arrays and sent to waterproperties for the next time step.   
!  
!******************************************************************************************************
!      IMPLEMENTATION:       CARBONATESYSTEM.DAT FILE EXAMPLE EXTENDED
!******************************************************************************************************
!
!  !---keyword----------------- units ----------------default value--------------description
!
!   ALK_parametrized      Compute option type         .true.  [1]       !Compute alkalinity by salinity parametrization 
!   ALK_biology           Compute option type         .false. [0]       !Compute changes in alkalinity due biological activity
!   DIC_calc              Compute option type         .false. [0]       !Compute DIC taking into account CaCO3 dissol/precip
!   DIC_no_calc           Compute option type         .true.  [0]       !Compute DIC withou taking into account CaCO3 dissol/preci
!   DT                        seconds                     3600.         !Time step between two module calls
!   PELAGIC_MODEL         Compute option type              -            !Pelagic model used to do calculations.Options: WaterQuality
!   REDFIELD_NC             mol N (mol C -1)             16./122.       !Nitrogen to Carbon Redfield ratio            | or LifeModel
! 
!   
!  Example| Example of CarbonateSystem_1.dat for the user:
!  Example|  
!  Example|      DT: 3600. !Not an option. pre-defined in the WaterQuality.dat 
!  Example|      REDFIELD_NC: 16./122.
!  Example|
!  Example| !Alkalinity compute options. Note that you can't use both methods at the same time
!  Example|      ALK_parametrized          : 0
!  Example|      ALK_biology               : 1
!  Example|    
!  Example| !DIC compute options. Note that you can't use both methods at the same time
!  Example|      DIC_calc                  : 0
!  Example|      DIC_no_calc               : 1   
!  Example|    
!  Example| !Pelagic model compute option. Note that you can only use one. 
!  Example|      PELAGIC_MODEL             : WaterQuality (or LifeModel)
    

!*************************************************************************************************
!     IMPLEMENTATION:      WATERPROPERTIES.DAT FILE WITH CARBONATE SYSTEM VARIABLES
!*************************************************************************************************      
!    On the WaterProperties.dat file (implementation file) you should add properties blocks with part 
!    of the following properties. Remember, for alkalinity only one property placed below can be added
!    (as for DIC). Remember setting by one the keyword CARBSYST  : 1  (look at the example below).  
!    
!     
!     'parametrized alkalinity CS'             |ALK_cs_p  = Alkalinity  ( umol / kgSW )
!     'biological alkalinity CS'               |ALK_cs_b  = Alkalinity  ( umol / kgSW )
!     'dissolved inorganic carbon CS calc'     |DIC_cs_c  = Dissolved Inorganic Carbon concentration, calcification 
!     'dissolved inorganic carbon CS nocalc'   |DIC_cs_nc = Dissolved Inorganic Carbon concentration, no calcification 
!                                                                ( umolC  / kgSW ) | (mg C m-3)
!    
!    
! Example| Example of WaterProperties.dat with a carbonate system variable: !@marta mejorar este ejemplo
!    
! Example| <beginproperty>
! Example| NAME                         : parametrized alkalinity CS!Property name 
! Example| UNITS                        : umol/kgSW                 !Property units
! Example| DESCRIPTION                  : alkalinity                !Small description of the property
! Example| INITIALIZATION_METHOD        : CONSTANT                  !Initialization of concentration values. See module FillMatrix
! Example| PARTICULATE                  : 0                         !Property physical state: 0 - Dissolved ; 1 - Particulate
! Example| OLD                          : 0                         !Initialization from previous run (overrides FillMatrix)
! Example| IS_COEF                      : .001                      !Conversion factor to I.S. units (1e-3 = mg/l)
! Example| DEFAULTVALUE                 : 2350                      !Value assumed by default
! Example| DEFAULTBOUNDARY              : 2350                      !Value assumed in open boundaries by default
! Example| BOUNDARY_INITIALIZATION      : INTERIOR                  !Type of boundary initialization: INTERIOR or EXTERIOR
! Example| ADVECTION_DIFFUSION          : 1                         !Compute advection-diffusion
! Example| 	BOUNDARY_CONDITION          : 2                         !See modules WaterProperties and Hydrodynamics
! Example| 	ADV_METHOD_H                : 1                         !See modules WaterProperties and Hydrodynamics
! Example| 	ADV_METHOD_V                : 1                         !See modules WaterProperties and Hydrodynamics
! Example|     ADVECTION_V_IMP_EXP      : 0                         !See modules WaterProperties and Hydrodynamics
! Example| 	ADVECTION_H_IMP_EXP         : 1                         !See modules WaterProperties and Hydrodynamics  	
! Example| SUBMODEL                     : 0                         !Property is influenced by a father model 
! Example| 	SUBMODEL_INI                : 0                         !Property is initialized as being part of a sub model                         
! Example| PARTITION                    : 0						    !Name of the property (oposite phase) to compute partition
! Example| WATER_QUALITY                : 0                         !Compute water quality processes (OFF = 0)
! Example| CARBSYST                     : 1                         !Compute CARBonate SYSTem processes (ON = 1)
! Example| SURFACE_FLUXES               : 0                         !Compute fluxes from water-air interface (always 0 for alkal.)
! Example| BOTTOM_FLUXES                : 0                         !Compute fluxes from sediment-water interface (0 for alkalinity)
! Example| DISCHARGES                   : 0							!Compute discharges (WWTP, river, etc)
! Example| 	DISCHARGES_TRACKING         : 0							!Monitor discharges with outputing a time serie 
! Example| VERTICAL_MOVEMENT            : 0							!Compute vertical movement due to settling velocity
! Example| DATA_ASSIMILATION            : 1							!Add nudging term
! Example| TIME_SERIE                   : 0                         !Ouputs results in time series
! Example| BOX_TIME_SERIE               : 0                         !Ouputs results in box time series
! Example| OUTPUT_HDF                   : 1                         !Ouputs results in HDF5 format
! Example| OUTPUT_SURFACE_HDF           : 0                         !Ouputs results in HDF5 format
! Example| STATISTICS                   : 0                         !Perform statistics for property concentration
! Example| 	STATISTICS_FILE             : C:\...\STATISTICS_alkalinity.dat
! Example| <endproperty>  
    
!*************************************************************************************************
!*************************************************************************************************   
          
 Module ModuleCarbonateSystem
    
      use ModuleGlobalData
      use ModuleEnterData      
      use ModuleTime


      implicit none            
      private      
            
!******************************************************************************
!                            SUBROUTINES 
!******************************************************************************
      
     !Constructor      
      public  :: ConstructCarbonateSystem              !Model preparation-construction
      private ::      AllocateInstance                 !Instance allocation
      private ::      ReadData                         !Calling of reading subroutines
      private ::          ConstructGlobalVariables     !Parameters, keywords in the .dat file
      private ::          ConstructModelOptions        !Calculation option, keywords in the .dat file
      private ::      PropertyIndexNumber
      private ::      ConstructPropertyList  
      private ::      CS_configuration
  
      
     !Selector          
      public  :: GetDTCarbonateSystem 
      public  :: GetCarbonateSystemPropertyList
      public  :: GetCarbonateSystemSize
      public  :: GetCarbonateSystemPropIndex
      public  :: UnGetCarbonateSystem
    
    !Modifier                                          !Model calculations
      public  :: ModifyCarbonateSystem                 !Associates external variables, call subr to do calculations 
      private ::  Biogeochemical_ratios_and_parameters !If alkbio, stores some pelagicmodule param as glovbal CS var          
      private ::  CS_computations                      !Call, in a concrete order, subroutines starting as Compute 
      private ::   ComputeAlkalinity_param             !Calculates alkalinity through salt and temp, depending on grid area
      private ::   ComputeAlkalinity_bio               !Calculates alkalinity taking into account changes due to biological activity 
      private ::   ComputeAlkalinity_bio_calc          !" " but taking into account calcium carbonate processes too      
      private ::   ComputeDIC_calc                     !Calculates dissolved ingoranic carbon taking into account caco3 processes
      private ::   ComputeDIC_no_calc                  !Calculates dissolved ingoranic carbon withou taking into account caco3 process
      private ::     Compute_biogeoch_rates_nitr_denit  
      private ::     Compute_biogeoch_rates_resp
      private ::     Compute_biogeoch_rates_nitrogen_uptake
      private ::     Compute_biogeoch_rates_caco3
      !private ::  Computemocsy
      !private ::    Computemocsy_1_
    
    !Destructor        
      public  :: KillCarbonateSystem   
      private :: DeAllocateInstance
 
    !Management     
      private :: Ready
      private :: LocateObjCarbonateSystem

 
    
!*******************************************************************************************
!                              DECLARATION STATEMENTS
!*******************************************************************************************     
! ------ Types, variables, parameters 

    private :: T_AuxiliarParameters
    type       T_AuxiliarParameters       
        real                            :: C_AtomicMass          = 12                   !mgC/molC
        !real                            :: H_AtomicMass          = 1                    !mgH/molH
        !real                            :: O_AtomicMass          = 16                   !mgO/molO
        real                            :: N_AtomicMass          = 14                   !mgN/molN
        !real                            :: P_AtomicMass          = 31                   !mgP/molP
    end type T_AuxiliarParameters


    private :: T_BioChemParam
    type       T_BioChemParam   
        real                            :: Redfield_NC            = null_real           !Redfield N:C ratio
    !    real                            :: Redfield_PC            = null_real          !Redfield P:C ratio
    !    real                            :: Redfield_SiC           = null_real          !Standard Si:C ratio
    !    real                            :: BioSi_Diss             = null_real          !Biogenic silica dissolution rate
    !    real                            :: O2C_Conversion         = null_real          !Oxygen to carbon conversion factor
    !    real                            :: CO2C_Conversion        = null_real          !Carbon Dioxide to carbon conversion factor
    end type T_BioChemParam


    private :: T_ComputeOptions
    type       T_ComputeOptions
        logical                         :: NoBiologicalAlkalinity    = .true.          !Compute alkalinity by salinity
        logical                         :: BiologicalAlkalinity      = .false.         !Compute biological alkalinity 
        logical                         :: DIC_calcification         = .false.         !Compute DIC with calcif/dissolution
        logical                         :: DIC_no_calcification      = .true.          !Compute DIC without calcif/dissolu
        character(len=StringLength)     :: PelagicModel              = null_str        !Pelagic model to use outputs
        logical                         :: Diatoms                   = .false.         !Diatoms activated in WQ or Life
        logical                         :: Ciliates                  = .false.         !Ciliates activated in WQ or Life
    end type T_ComputeOptions 
    
    
    private :: T_PropIndex                         
    type     T_PropIndex ! CarbonateSystem properties-module index, NOT global data index!
        integer                         :: ALK_cs_p         = null_int   ! Alkalinity computed by algorithm 
        integer                         :: ALK_cs_b         = null_int   ! Alkalinity computed by biological changes 
        integer                         :: DIC_cs_c         = null_int   ! Dissolved inorganic carbon concentration (calcifi) 
        integer                         :: DIC_cs_nc        = null_int   ! Dissolved inorganic carbon concentration (no calc) 
        integer                         :: AM               = null_int   ! External ammonia concentration
        integer                         :: NI               = null_int   ! External nitrite concentration
        integer                         :: NA               = null_int   ! External nitrate concentration
        integer                         :: Oxygen           = null_int   ! External oxygen concentration
        integer                         :: Phyto            = null_int   ! External phytoplank biomass 
        integer                         :: Zoo              = null_int   ! External zooplank biomass
        integer                         :: Diatoms          = null_int   ! External diatoms  biomass
        integer                         :: Ciliates         = null_int   ! External ciliates biomass        
    end type T_PropIndex       
      
    
    private :: T_ID
    type       T_ID         
        integer                         :: ID, IDNumber             = null_int
        character(len=StringLength)     :: Name                     = null_str
        character(len=StringLength)     :: Description              = null_str
    end type T_ID  


   private ::  T_ExternalVar                !Output variables from other modules needed as input information for this one            
    type       T_ExternalVar               
        real, pointer, dimension(:  )    :: Salinity           !Salinity        1D array. Origin: WaterPropertiesModule
        real, pointer, dimension(:  )    :: Temperature        !Temperature     1D array. Origin: WaterPropertiesModule
        real, pointer, dimension(:  )    :: Thickness          !Cell thickness  1D array. Origin: GeometryModule,original name DWZ
        real, pointer, dimension(:,:)    :: Mass               !Property mass   2D array. Origin: WaterPropertiesModule
        real, pointer, dimension(:  )    :: Latitude           !Latitude        1D array. Origin: HorizontalGridModule
        real, pointer, dimension(:  )    :: Longitude          !Longitude       1D array. Origin: HorizontalGridModule 
        real, pointer, dimension(:  )    :: Ratios             !Biogeoch ratios 1D array. Origin: WaterQualityModule or Life 
        real, pointer, dimension(:  )    :: GrossGrowthRate    !Biogeoch rate   1D array. Origin: WaterQualityModule or Life 
        real, pointer, dimension(:  )    :: GrossGrowthRateDiat!Biogeoch rate   1D array. Origin: WaterQualityModule or Life 
        real, pointer, dimension(:  )    :: VolumenZ           !Cell volume (L) 1D array. Origin: GeometryModule
        real, pointer, dimension(:  )    :: Parameters         !Biogeoch param  1D array. Origin: WaterQualityModule or Life
        real, pointer, dimension(:  )    :: Nitrification1     !Nitrif1         1D array. Origin: WaterQualityModule 
        real, pointer, dimension(:  )    :: Nitrification2     !Nitrif2         1D array. Origin: WaterQualityModule
        real, pointer, dimension(:  )    :: Denitrification    !Denit           1D array. Origin: WaterQualityModule
        integer, pointer, dimension(:  ) :: OpenPoints         !Grid points with water where perform calculations. Origin: MapModule              
    end type T_ExternalVar               !Latitude and latitude are geographic coordenates in decimal format
                                         !Latitude, Longitude, Thickness, OpenPoints and Volumen are gotten by WaterPropertiesModule

    private ::  T_ExternalRatio                          
    type       T_ExternalRatio              
        real            :: NC_phyto             = null_real    
        real            :: PC_phyto             = null_real
        real            :: NC_zoo               = null_real   
        real            :: PC_zoo               = null_real
        real            :: NC_diatm             = null_real   
        real            :: PC_diatm             = null_real
        real            :: NC_cilia             = null_real   
        real            :: PC_cilia             = null_real
        integer         :: Diatm_are_calculated = null_int 
        integer         :: Cilia_are_calculated = null_int 
    end type T_ExternalRatio  
    
    private ::  T_ExternalParameter                         
    type       T_ExternalParameter 
        real            :: pelagic_module_dt       = null_real          
        real            :: pelagic_module_dt_day   = null_real
        real            :: MinOxygen               = null_real    
        real            :: TNitrification          = null_real
        real            :: NitrificationSatConst   = null_real   
        real            :: KNitrificationRateK1    = null_real
        real            :: KNitrificationRateK2    = null_real 
        real            :: KDenitrificationRate    = null_real
        real            :: TDenitrification        = null_real
        real            :: DenitrificationSatConst = null_real
        real            :: NSatConst               = null_real
        real            :: DiaNSatConst            = null_real
        real            :: PhytoEndogRepConst      = null_real
        real            :: PhotorespFactor         = null_real
        real            :: ZK1                     = null_real
        real            :: ZK2                     = null_real
        real            :: ZK3                     = null_real
        real            :: ZK4                     = null_real      
        real            :: Toptzoomin              = null_real
        real            :: Toptzoomax              = null_real
        real            :: Tzoomin                 = null_real
        real            :: Tzoomax                 = null_real
        real            :: ZooReferenceRespRate    = null_real
    end type T_ExternalParameter  


     private :: T_CarbonateSystem
     type       T_CarbonateSystem
        integer                                         :: InstanceID
        integer                                         :: ObjEnterData                  = 0
        integer, dimension(:), pointer                  :: PropertyList                  => null()
        integer                                         :: CSConfiguration               = 0
        integer                                         :: PelagicModuleConfiguration = 0
        !real                                           :: DT
        !real                                           :: DT_day
        type(T_AuxiliarParameters)                      :: AuxParam
        type(T_BioChemParam      )                      :: BioChemParam
        type(T_ComputeOptions    )                      :: ComputeOptions
        type(T_ExternalVar       )                      :: ExternalVar
        type(T_ExternalRatio     )                      :: ExternalRatio
        type(T_ExternalParameter )                      :: ExternalParam
        type(T_Size1D            )                      :: Array
        type(T_Size1D            )                      :: Prop
        type(T_PropIndex         )                      :: PropIndex
        type(T_CarbonateSystem   ), pointer             :: Next      

      end type T_CarbonateSystem
                                                      
 ! ------ Global Module Variables ----------------------   
      
       type (T_CarbonateSystem), pointer        :: FirstObjCarbonateSystem  => null()
       type (T_CarbonateSystem), pointer        :: Me                       => null()
       

        
    contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
    

!>@author: Marta López, Maretec
!>@Brief: Register this module, allocates instance, calls the subroutines 
!>which read the information to construct the module (from the .datfile), 
!>gives an index to properties that will be calculated in the module (maybe 
!>activated by the user or necessary for the module; the index is different 
!>from the one in globaldata) and construct the PropertyList with them.     
!>@param[in]  FileName
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine ConstructCarbonateSystem(ObjCarbonateSystemID, &                                      
                                      FileName,             &
                                      Alk_bio,              &
                                      STAT)
    
     !Arguments---------------------------------------------------------------
       integer                        :: ObjCarbonateSystemID   
       character(len=*) , intent(IN)  :: FileName !path to carbonatesystem.dat file
       logical, optional, intent(OUT) :: Alk_bio
       integer, optional, intent(OUT) :: STAT
     !External----------------------------------------------------------------
       integer                                 :: ready_, STAT_CALL
     !Local------------------------------------------------------------------
       integer                                 :: STAT_
     !------------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mCarbonateSystem_)) then
          nullify (FirstObjCarbonateSystem)
          call RegisterModule (mCarbonateSystem_)         
        endif        
        !Call Ready, which is placed on the management block 
        call Ready(ObjCarbonateSystemID, ready_)
        !If ready_, which is the return of the sbrtine Ready, is equal to OFF_ERR_, call the following 
        !subroutines, placed below        
cd0 : if (ready_ .EQ. OFF_ERR_) then

          call AllocateInstance
          call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructCarbonateSystem - ModuleCarbonateSystem - ERROR #1'
          call ReadData
          call PropertyIndexNumber
          call ConstructPropertyList
          call CS_configuration
             Alk_bio = .false.
             if(Me%ComputeOptions%BiologicalAlkalinity)then
             Alk_bio = .true.
             end if
          call KillEnterData(Me%ObjEnterData,STAT = STAT_CALL)   
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructCarbonateSystem - ModuleCarbonateSystem - ERROR #2'
        
        !Returns ID 
          ObjCarbonateSystemID     = Me%InstanceID
          STAT_ = SUCCESS_
      else 
          stop 'ModuleCarbonateSystem - ConstructCarbonateSystem - ERROR #1'
      end if cd0
    
      ! If the optional argument STAT is present in this subroutine (so, the function present returns 
      ! true) assign to STAT_ (previously assigned as SUCCES_ as STAT to 
        if (present(STAT)) STAT = STAT_   
        
  !---------------------------------------------------------------------------               
  end subroutine ConstructCarbonateSystem
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  
  
  !>@author: Marta López, Maretec
  !>@Brief: Allocates instance 
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   subroutine AllocateInstance
   
  !Local---------------------------------------------------------------------
        type (T_CarbonateSystem), pointer      :: NewObjCarbonateSystem
        type (T_CarbonateSystem), pointer      :: PreviousObjCarbonateSystem
  !---------------------------------------------------------------------------
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
  !---------------------------------------------------------------------------
   end subroutine AllocateInstance
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
   
  
 !>@author: Marta López, Maretec
 !>@Brief: Calls the subroutines that read and set the info to construct the 
 !>module (from the carbonatesystem.dat file reads keywords and values; if not
 !present/necessary, sets the default value).   
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
   subroutine ReadData  
  !---------------------------------------------------------------------------
    call ConstructModelOptions
    call ConstructGlobalVariables    
    !call ConstructCO2SYS_mocsy
  !---------------------------------------------------------------------------
   end subroutine ReadData
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
   

  !>@author: Marta López, Maretec
  !>@Brief: Sets module calculation options and their keywords. Through CS.dat 
  !>file this configuration can be changed by the user.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
   subroutine ConstructModelOptions
  
  !Arguments-------------------------------------------------------------------        
     !External---------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        integer                                         :: iflag1, iflag2
        integer                                         :: iflag3, iflag4
      !Local-----------------------------------------------------------------
        integer                                         :: FromFile 
      !----------------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)
        
        !ALKALINITY
        !Calculation option: compute alkalinity by salinity parametrization        
        call GetData(Me%ComputeOptions%NoBiologicalAlkalinity,        &
                     Me%ObjEnterData, iflag1,                         &
                     SearchType   = FromFile,                         &
                     keyword      = 'ALK_parametrized ',              &
                     Default      = .false.,                          &
                     ClientModule = 'ModuleCarbonateSystem',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #1'        
        if(Me%ComputeOptions%NoBiologicalAlkalinity) then
            write(*,*)'WARNING: Alkalinity parametrization through salinity should be used only for'
            write(*,*)'surface water cases of study. Be careful'
        endif
        
        !Calculation option: compute changes on alkalinity by biological interactions
        call GetData(Me%ComputeOptions%BiologicalAlkalinity,           &
                     Me%ObjEnterData, iflag2,                          &
                     SearchType   = FromFile,                          &
                     keyword      = 'ALK_biology ',                    &
                     Default      = .false.,                           &
                     ClientModule = 'ModuleCarbonateSystem',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #2'        
  
        !Verification if alkalinity calculation option is correctly activated  
        if((iflag1==0) .and. (iflag2==0)) then
            write(*,*)'Please define an alkalinity calculation option. Choose one:ALK_parametrized or'
            write(*,*)'ALK_biology in the carbonatesystem.dat file '
            stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #V1a'
        end if         
        
        if(Me%ComputeOptions%NoBiologicalAlkalinity .and. Me%ComputeOptions%BiologicalAlkalinity) then
          write(*,*)' You cannot activate two alkalinity calculation options.Choose one:ALK_parametrized or'
          write(*,*)' ALK_biology. Please, check it in the carbonatesystem.dat file ' 
            stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #V1b' 
        endif    
        
        !DIC
        !Calculation option: compute DIC taking into account calcification/dissolut processes
        call GetData(Me%ComputeOptions%DIC_calcification,             &
                     Me%ObjEnterData, iflag3,                          &
                     SearchType   = FromFile,                         &
                     keyword      = 'DIC_calc ',                      &
                     Default      = .false.,                          &
                     ClientModule = 'ModuleCarbonateSystem',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #3'
        
        !Calculation option: compute DIC without taking into account calcification/dissolut proc.
        call GetData(Me%ComputeOptions%DIC_no_calcification,          &
                     Me%ObjEnterData, iflag4,                          &
                     SearchType   = FromFile,                         &
                     keyword      = 'DIC_no_calc ',                   &
                     Default      = .false.,                           &
                     ClientModule = 'ModuleCarbonateSystem',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #4'
        
        !Verification if DIC calculation option is correctly activated  
        if((iflag3==0) .and. (iflag4==0)) then
            write(*,*)'Please define a DIC calculation option. Choose one: DIC_calc or'
            write(*,*)'DIC_no_calc in the carbonatesystem.dat file '
            stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #V2a'
        end if  
        
        if(Me%ComputeOptions%NoBiologicalAlkalinity .and. Me%ComputeOptions%BiologicalAlkalinity) then
          write(*,*)' You cannot activate two DIC calculation options. Choose one: DIC_calc or'
          write(*,*)' DIC_no_calc. Please, check it in the carbonatesystem.dat file ' 
            stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #V2b'             
        endif
        
        !PELAGIC MODEL
        !Pelagic model that is going to be used
        call GetData(Me%ComputeOptions%PelagicModel,                   &
                    Me%ObjEnterData,  iflag,                            &
                    SearchType   = FromFile,                           &
                    keyword      = 'PELAGIC_MODEL',                    &
                    ClientModule = 'ModuleCarbonateSystem',            &
                    STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #5'

        !Verification if pelagic model option is correctly activated
        if(iflag==0)then
            write(*,*)'Please define the pelagic model to couple with ModuleCarbonateSystem'
            stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #6'
        end if 
        !PROVISIONAL: module life
        if(Me%ComputeOptions%PelagicModel .eq. 'LifeModel') then
            write(*,*)'The Life model option to couple with ModuleCarbonateSystem is still under construction:'
            write(*,*)'Please, use WaterQualityModel instead. Thank you'
            stop 'PelagicModuleLIFE no possible- ConstructModelOptions - ModuleCarbonateSystem - ERROR #6b'
        end if      
                       
       if((Me%ComputeOptions%PelagicModel .ne. WaterQualityModel .and. Me%ComputeOptions%PelagicModel .ne. 'LifeModel'))then
            write(*,*)'Pelagic model to couple with ModuleCarbonateSystem must be one of the following:'
            write(*,*)trim(WaterQualityModel)
            write(*,*)trim(LifeModel)
            stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #7'
       end if     
       
        call GetData(Me%ComputeOptions%Diatoms,                       &
                     Me%ObjEnterData, iflag1,                         &
                     SearchType   = FromFile,                         &
                     keyword      = 'DIATOMS ',                       &
                     Default      = .false.,                          &
                     ClientModule = 'ModuleCarbonateSystem',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #8'   
        
         call GetData(Me%ComputeOptions%Ciliates,                     &
                     Me%ObjEnterData, iflag1,                         &
                     SearchType   = FromFile,                         &
                     keyword      = 'CILIATES',                       &
                     Default      = .false.,                          &
                     ClientModule = 'ModuleCarbonateSystem',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleCarbonateSystem - ERROR #9'   
       
       !...............
        write(*,*)" ------------ CARBONATE SYSTEM OPTIONS  ------------    "
        write(*,*)"                                        "
        write(*,*)" Alkalinity   - parametrized --------  :", Me%ComputeOptions%NoBiologicalAlkalinity
        write(*,*)" Alkalinity   - biological   --------  :", Me%ComputeOptions%BiologicalAlkalinity
        write(*,*)" DisInorgCarb -  no CaCO3    --------  :", Me%ComputeOptions%DIC_no_calcification
        write(*,*)" DisInorgCarb - yes CaCO3    --------  :", Me%ComputeOptions%DIC_calcification
        write(*,*)" Diatoms in pelagic module   --------  :", Me%ComputeOptions%Diatoms
        write(*,*)" Ciliates in pelagic module  --------  :", Me%ComputeOptions%Ciliates
        write(*,*)" Pelagic module              --------  :", Me%ComputeOptions%PelagicModel

  !----------------------------------------------------------------------    
   end subroutine ConstructModelOptions 
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
   
  
   !>@author: Marta López, Maretec
  !>@Brief: Sets model parameters and their keywords. Through CS.dat 
  !>file this configuration can be changed by the user. It also gets the time 
  !>step of the pelagicModel activated as well as the phyto and diatoms C:N:P ratio
  !>used in those modules. 
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
 
   subroutine ConstructGlobalVariables

     !Arguments----------------------------------------------------------------        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL        
        !Local-----------------------------------------------------------------
        integer                                         :: FromFile 
        !Begin-----------------------------------------------------------------
        
        call GetExtractType    (FromFile = FromFile)     
        
        !Nitrogen to Carbon Redfield ratio ( mol N / mol C ) !marta escoges este valor siguiendo a pisces
        !call GetData(Me%BioChemParam%Redfield_NC,                     &
        !             Me%ObjEnterData, iflag,                          &
         !            SearchType   = FromFile,                         &
          !           keyword      = 'REDFIELD_NC',                    &
           !          Default      = 16./122.,                         &
            !         ClientModule = 'ModuleCarbonateSystem',          &
             !        STAT         = STAT_CALL)
        !if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleCarbonateSystem - ERROR #1'
        
        !Below it is written the code for reading DT from the carbsyst.dat file. Because this module needs
        !outputs from other biogeochemical modules, carbonatesystem dt is going to be equal to the model
        !coupled. This is done in WaterPropertiesModule. However, it would be possible to change in the future
        !(which I do not recommend). Uncomment here, in global properties declaration(dt, dt_day),in subroutine
        !CoupleCarbonateSystem placed in WaterProperties and subroutine GetDTCarbonateSystem placed in this module. 
        !DTSecond, time step, in seconds, between 2 module calls (seconds) from WaterQuality or Life Module
        !call GetData(Me%DT,                                           &
                     !Me%ObjEnterData, iflag,                          &
                     !SearchType   = FromFile,                         &
                     !keyword      = 'DT',                             &
                     !Default      = 3600.,                            &
                     !ClientModule = 'ModuleCarbonateSystem',          &
                     !STAT         = STAT_CALL)
        !if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleCarbonateSystem - ERROR #1'
        !For compatibility with the rest of the program,  DT [sec] converted to DT [day]
        !Me%DT_day      = Me%DT / (3600.*24.)
        
  !---------------------------------------------------------------------------  
   end subroutine ConstructGlobalVariables
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
   
  
  !> @author: Marta López, Maretec
  !> @Brief: Gives an index, different form the one given in globaldata!, to the 
  !> properties that are going to be calculated in the module. The minimun number of 
  !> properties for this module is two, and the maximun is !@martaaaaa. The index number
  !> will depend on the order in which the subroutine is written: alkalinity will always
  !> have an propindex of 1, because just one type of it can be calculated and it is the
  !> first property written below.
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
  
   subroutine PropertyIndexNumber  
  !----------------------------------------------------------------------------
   
        Me%Prop%ILB = 1
        Me%Prop%IUB = 0
 
        !Alkalinity
        if(Me%ComputeOptions%NoBiologicalAlkalinity) then
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%ALK_cs_p            = Me%Prop%IUB    
        endif
        
        if(Me%ComputeOptions%BiologicalAlkalinity)  then
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%ALK_cs_b            = Me%Prop%IUB            
        endif  
        
        !DIC 
        if(Me%ComputeOptions%DIC_no_calcification)    then
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%DIC_cs_nc           = Me%Prop%IUB
        endif     
        
        if(Me%ComputeOptions%DIC_calcification)       then
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%DIC_cs_c            = Me%Prop%IUB
        endif  
        
                      
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%AM                  = Me%Prop%IUB
            
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%NI                  = Me%Prop%IUB
            
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%NA                  = Me%Prop%IUB
            
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%Oxygen              = Me%Prop%IUB   
            
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%Phyto               = Me%Prop%IUB
            
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%Zoo                 = Me%Prop%IUB
            
         if(Me%ComputeOptions%Diatoms)       then   
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%Diatoms             = Me%Prop%IUB
         endif
         
         if(Me%ComputeOptions%Ciliates)       then
            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%Ciliates            = Me%Prop%IUB
         endif
  !---------------------------------------------------------------------- 
   end subroutine PropertyIndexNumber  
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
  
   
   
  !> @author: Marta López, Maretec
  !> @Brief: Allocates PropertyList(1D) with max dimension until the number of 
  !> properties that are going to be calculated in this module and stores the
  !> global data property index (not the index given in the module!)
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
   subroutine ConstructPropertyList

   !Arguments-------------------------------------------------------------     
        
        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))        
                
        if(Me%ComputeOptions%NoBiologicalAlkalinity) then
        Me%PropertyList(Me%PropIndex%ALK_cs_p )   = ALK_cs_p_
        endif      
        
        if(Me%ComputeOptions%BiologicalAlkalinity)   then
        Me%PropertyList(Me%PropIndex%ALK_cs_b )   = ALK_cs_b_ 
        endif
        
        if(Me%ComputeOptions%DIC_no_calcification)   then
        Me%PropertyList(Me%PropIndex%DIC_cs_nc)   = DIC_cs_nc_
        endif
        
        if(Me%ComputeOptions%DIC_calcification)      then
        Me%PropertyList(Me%PropIndex%DIC_cs_c )   = DIC_cs_c_
        endif
        
        Me%PropertyList(Me%PropIndex%AM       )   = Ammonia_
        Me%PropertyList(Me%PropIndex%NI       )   = Nitrite_
        Me%PropertyList(Me%PropIndex%NA       )   = Nitrate_
        Me%PropertyList(Me%PropIndex%Oxygen   )   = Oxygen_
        Me%PropertyList(Me%PropIndex%Phyto    )   = Phytoplankton_
        Me%PropertyList(Me%PropIndex%Zoo      )   = Zooplankton_ 
        
        if(Me%ComputeOptions%Diatoms          ) then
        Me%PropertyList(Me%PropIndex%Diatoms  )   = Diatoms_
        endif
        if(Me%ComputeOptions%Ciliates         ) then
        Me%PropertyList(Me%PropIndex%Ciliates  )   = Ciliate_
        endif
        
        
  !---------------------------------------------------------------------- 
   end subroutine ConstructPropertyList  
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
    
   
  !> @author: Marta López, Maretec
  !> @Brief: Depending on the CS.dat configuration, stores an integer number in
  !> Me%CSConfiguration to let modifier know the subroutines to call 
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
   
   subroutine CS_configuration
   
   !--------------------------------------------------------------------------- 
   
  i1:  if(Me%ComputeOptions%NoBiologicalAlkalinity) then
            
            if (Me%ComputeOptions%DIC_no_calcification) then
                
                Me%CSConfiguration = 1
            else
                Me%CSConfiguration = 2
                
            endif    
            
        else i1
      
            if (Me%ComputeOptions%DIC_no_calcification) then
                
                Me%CSConfiguration = 3
            else
                Me%CSConfiguration = 4
                
            endif    
        
        endif i1 
  !---------------------------------------------------------------------------- 
   end subroutine  CS_configuration
  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
   
   
   

   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
   !> @author: Marta López, Maretec
   !> @Brief: Sends the CarbonateSysteMmodule time step and daily time step 
   !> to the module that need it (InterfaceModule). This subroutine is build 
   !> because of the MOHID structure, but actually the value is irrelevant 
   !> because once DT is acquired by WProperties, there is equalized to WQDT
   !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
  
   subroutine GetDTCarbonateSystem(CarbonateSystem_ID, DTDay, DT, STAT)
   
        !Arguments-------------------------------------------------------------
        integer                             :: CarbonateSystem_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DT
        integer, optional, intent(OUT)      :: STAT
        !External--------------------------------------------------------------
        integer                             :: ready_    
        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(CarbonateSystem_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            if (present(DTDay   )) DTDay    = 0.02  !DT_day
            if (present(DT)) DT             = 1800. !DT
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_
    !----------------------------------------------------------------------
   end subroutine GetDTCarbonateSystem
   !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   
   !> @author: Marta López, Maretec
   !> @Brief: Sends the number of properties that are going to be calculated in
   !> this module, during an specific simulation, to the module that need it
   !> (which is the interface module). 
   !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
   subroutine GetCarbonateSystemSize(CarbonateSystem_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: CarbonateSystem_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT
        !External--------------------------------------------------------------
        integer                             :: ready_             
        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(CarbonateSystem_ID, ready_)            
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

   end subroutine GetCarbonateSystemSize
   !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
   
   !> @author: Marta López, Maretec
   !> @Brief: Sends the global data index of properties that are going to be calculated in
   !> this module, during an specific simulation, to the module that need it
   !> (interface module). There, interface will compare if the properties activated in the 
   !> waterpropertiesmodule (for being calculated in CSmodule with CARBSYST keyword on) are 
   !> the same as in the carbonatesystempropertylist.
   !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine GetCarbonateSystemPropertyList(CarbonateSystem_ID, PropertyList, STAT)

        !Arguments-------------------------------------------------------------
        integer                                            :: CarbonateSystem_ID
        integer, dimension(:), pointer, intent(OUT)        :: PropertyList
        integer, optional, intent(OUT)                     :: STAT

        !External--------------------------------------------------------------
        integer                                            :: ready_            
        !Local-----------------------------------------------------------------
        integer                                            :: STAT_
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(CarbonateSystem_ID, ready_)            
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mCarbonateSystem_, Me%InstanceID)
            PropertyList =>  Me%PropertyList
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_
    !----------------------------------------------------------------------
    end subroutine GetCarbonateSystemPropertyList
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    !> @author: Marta López, Maretec
    !> @Brief: Nullify the carbonatesystemproplist in the interface module, after it 
    !> compares if the properties activated in the waterpropertiesmodule 
    !> (for being calculated in CSmodule with CARBSYST keyword on) are 
    !> the same as in the carbonatesystempropertylist.  
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    subroutine UnGetCarbonateSystem(ObjCarbonateSystemID, Array, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: ObjCarbonateSystemID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(ObjCarbonateSystemID, ready_)
        if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)
            call Read_Unlock(mCarbonateSystem_, Me%InstanceID, "UnGetCarbonateSystem3D_I")
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_        
    !----------------------------------------------------------------------
    end subroutine UnGetCarbonateSystem
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    !> @author: Marta López, Maretec
    !> @Brief: Compares the global property ID of all properties activated in 
    !> the water properties module with the global property ID of the PropertyList,
    !> this is, the properties that are going to be calculated in this module. 
    !> This subroutine is called by module interface,first to fill variable mass of 
    !> the properties this module will calculate, and 2nd to check consistence between
    !> module interface and carbonatesystem.
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    
    subroutine GetCarbonateSystemPropIndex (CarbonateSystem_ID, PropIDNumber, PropertyIndex, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: CarbonateSystem_ID
        integer,           intent(IN )      :: PropIDNumber
        integer,           intent(OUT)      :: PropertyIndex
        integer, optional, intent(OUT)      :: STAT
        !External--------------------------------------------------------------
        integer                             :: ready_            
        !Local-----------------------------------------------------------------
        integer                             :: STAT_, CurrentIndex        
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(CarbonateSystem_ID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            do CurrentIndex = Me%Prop%ILB, Me%Prop%IUB
                if (PropIDNumber == Me%PropertyList(CurrentIndex))then
                    PropertyIndex = CurrentIndex
                    STAT_ = SUCCESS_
                    exit
                else
                    STAT_ = NOT_FOUND_ERR_
                end if
            end do
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_
    !--------------------------------------------------------------------------
    end subroutine GetCarbonateSystemPropIndex
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      
     
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !>@author Marta López Maretec
    !>@Brief:
    !> Associates external variables needed in this module, already calculated in other 
    !> modules. Checks the grid points wherein perform calculations and calls
    !> the subroutine which contains the ones that do computations. This subroutine is 
    !> called each time step.
    !>@param[in] ObjCarbonateSystemID,Salinity,Temperature;thickness,Mass,OpenPoints,
    !>@ArraySize,Latitude,Longitude   
    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
    subroutine ModifyCarbonateSystem(ObjCarbonateSystemID,                       &
                                     Salinity,                                   & 
                                     Temperature,                                &
                                     Thickness,                                  &
                                     Mass,                                       &
                                     OpenPoints,                                 &
                                     ArraySize,                                  &
                                     Latitude,                                   &
                                     Longitude,                                  &
                                     Ratios,                                     &
                                     GrossGrowthRate,                            &
                                     GrossGrowthRateDiat,                        &
                                     VolumenZ,                                   &                            
                                     STAT)
    !Arguments-------------------------------------------------------------------------    
      integer                                                   :: ObjCarbonateSystemID
      real,               pointer, intent(IN), dimension(:  )   :: Salinity
      real,               pointer, intent(IN), dimension(:  )   :: Temperature
      real,               pointer, intent(IN), dimension(:  )   :: Thickness
      real,               pointer, intent(IN), dimension(:,:)   :: Mass
      integer, optional,  pointer, intent(IN), dimension(:  )   :: OpenPoints
      real,               pointer, intent(IN), dimension(:  )   :: Latitude
      real,               pointer, intent(IN), dimension(:  )   :: Longitude
      real,               pointer, intent(IN), dimension(:  )   :: Ratios
      real,               pointer, intent(IN), dimension(:  )   :: GrossGrowthRate
      real,   optional,   pointer, intent(IN), dimension(:  )   :: GrossGrowthRateDiat
      real,               pointer, intent(IN), dimension(:  )   :: VolumenZ
      type(T_Size1D)             , intent(IN)                   :: ArraySize
      integer, optional,           intent(OUT)                  :: STAT           
    !Local------------------------------------------------------------------------------      
      integer                                     :: STAT_, ready_
      logical                                     :: CalcPoint    
      integer                                     :: index        !Element of the 1Darray
    !----------------------------------------------------------------------
        STAT_ = UNKNOWN_
        
        call Ready(ObjCarbonateSystemID, ready_)
        if (ready_ .EQ. IDLE_ERR_) then

         !Associates external variables (salinity) to variables for this module (Me%ExternalVar%Salinity)
            Me%ExternalVar%Salinity                   => Salinity
            if (.NOT. associated(Me%ExternalVar%Salinity))         &
                stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR01' 
            Me%ExternalVar%Temperature                => Temperature
            if (.NOT. associated(Me%ExternalVar%Temperature))        &
                stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR02'       
            Me%ExternalVar%Thickness                  => Thickness
            if (.NOT. associated(Me%ExternalVar%Thickness))          &
                stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR03'  
            Me%ExternalVar%Mass                       => Mass
            if (.NOT. associated(Me%ExternalVar%Mass))               &
                stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR04' 
            Me%ExternalVar%OpenPoints                 => OpenPoints
            if (.NOT. associated(Me%ExternalVar%OpenPoints))         &
               stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR05' 
            Me%ExternalVar%Latitude                   => Latitude
            if (.NOT. associated(Me%ExternalVar%Latitude))           &
               stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR06' 
            Me%ExternalVar%Longitude                  => Longitude
            if (.NOT. associated(Me%ExternalVar%Longitude))           &
               stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR07'                                      
            Me%ExternalVar%Ratios                     => Ratios
                if (.NOT. associated(Me%ExternalVar%Ratios))               &
               stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR08'
            Me%ExternalVar%GrossGrowthRate            => GrossGrowthRate 
                if (.NOT. associated(Me%ExternalVar%GrossGrowthRate ))      &
               stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR09'  
            if(present (GrossGrowthRateDiat))then
            Me%ExternalVar%GrossGrowthRateDiat        => GrossGrowthRateDiat 
                if (.NOT. associated(Me%ExternalVar%GrossGrowthRateDiat ))      &
               stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR10'       
            Me%ExternalVar%VolumenZ                   => VolumenZ
                if (.NOT. associated(Me%ExternalVar%VolumenZ))               &
               stop 'Subroutine ModifyCarbonateSystem - ModuleCarbonateSystem. ERR11'            
                    
                       
            !Associates array dimension (external) to array dimension variable of this module(Me%Array%) 
            Me%Array%ILB = ArraySize%ILB
            Me%Array%IUB = ArraySize%IUB
            
            !Store in individual variables info contained in Ratios array
            call Biogeochemical_ratios_and_parameters 
            
            !Checks along dimensions of open point 1Darray (grid array converted to 1D) the points with water
d1:         do index = Me%Array%ILB, Me%Array%IUB            
 i1:            if (present(OpenPoints)) then
   i2:              if (OpenPoints(index) == OpenPoint) then
                       CalcPoint = .true.
                    else
                       CalcPoint = .false.
                    endif i2
                else
                    CalcPoint = .true.
                endif i1               
               !If the array element contains water (CalcPoint = true), do computations in it 
 i3:           if (CalcPoint) then 
                call CS_computations(index)  
               endif i3
            end do d1
           
            !When calculations are done, nullify associated external var for the next time step calculation
            nullify(Me%ExternalVar%Salinity   )
            nullify(Me%ExternalVar%Temperature)
            nullify(Me%ExternalVar%Mass       )  
            nullify(Me%ExternalVar%Latitude   )  
            nullify(Me%ExternalVar%Longitude  )  
            nullify(Me%ExternalVar%Ratios     ) 
            nullify(Me%ExternalVar%GrossGrowthRate)
            if(present (GrossGrowthRateDiat))then
            nullify(Me%ExternalVar%GrossGrowthRate) 
            endif
                        
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_           
 !----------------------------------------------------------------------------
  end subroutine ModifyCarbonateSystem
 !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

                                     
                                     
 !>@author Marta López, Maretec
 !>@Brief: Stores in individual variables the ratios and parameters contained in 
 !> Me%ExternalVar%Ratios 1D array. Depending on waterquality configuration, stores
 !> Me%PelagicModuleConfiguration the WQ organism activated                                    
 !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
 subroutine Biogeochemical_ratios_and_parameters
 !-------------------------------------------------------------------------

   Me%ExternalRatio%Diatm_are_calculated    = Me%ExternalVar%Ratios(1) 
   Me%ExternalRatio%Cilia_are_calculated    = Me%ExternalVar%Ratios(2) 
   Me%ExternalParam%pelagic_module_dt       = Me%ExternalVar%Ratios(3) 
   Me%ExternalParam%pelagic_module_dt_day   = Me%ExternalVar%Ratios(4)    
   Me%ExternalRatio%NC_phyto                = Me%ExternalVar%Ratios(5)  
   Me%ExternalRatio%PC_phyto                = Me%ExternalVar%Ratios(6)  
   Me%ExternalRatio%NC_zoo                  = Me%ExternalVar%Ratios(7)  
   Me%ExternalRatio%PC_zoo                  = Me%ExternalVar%Ratios(8)   
   Me%ExternalParam%MinOxygen               = Me%ExternalVar%Ratios(9)
   Me%ExternalParam%TNitrification          = Me%ExternalVar%Ratios(10)
   Me%ExternalParam%NitrificationSatConst   = Me%ExternalVar%Ratios(11)
   Me%ExternalParam%KNitrificationRateK1    = Me%ExternalVar%Ratios(12)
   Me%ExternalParam%KNitrificationRateK2    = Me%ExternalVar%Ratios(13)
   Me%ExternalParam%KDenitrificationRate    = Me%ExternalVar%Ratios(14)
   Me%ExternalParam%TDenitrification        = Me%ExternalVar%Ratios(15)
   Me%ExternalParam%DenitrificationSatConst = Me%ExternalVar%Ratios(16)
   Me%ExternalParam%NSatConst               = Me%ExternalVar%Ratios(17)
   Me%ExternalParam%PhytoEndogRepConst      = Me%ExternalVar%Ratios(18)
   Me%ExternalParam%PhotorespFactor         = Me%ExternalVar%Ratios(19)
   Me%ExternalParam%ZK1                     = Me%ExternalVar%Ratios(20)
   Me%ExternalParam%ZK2                     = Me%ExternalVar%Ratios(21)
   Me%ExternalParam%ZK3                     = Me%ExternalVar%Ratios(22)
   Me%ExternalParam%ZK4                     = Me%ExternalVar%Ratios(23)
   Me%ExternalParam%Toptzoomin              = Me%ExternalVar%Ratios(24)
   Me%ExternalParam%Toptzoomax              = Me%ExternalVar%Ratios(25)
   Me%ExternalParam%Tzoomin                 = Me%ExternalVar%Ratios(26)
   Me%ExternalParam%Tzoomax                 = Me%ExternalVar%Ratios(27)
   Me%ExternalParam%ZooReferenceRespRate    = Me%ExternalVar%Ratios(28)
   
   
i1:  if(Me%ExternalRatio%Diatm_are_calculated == 1)then
  i2:   if(Me%ExternalRatio%Cilia_are_calculated == 0)then   
           Me%ExternalRatio%NC_diatm                = Me%ExternalVar%Ratios(29)  
           Me%ExternalRatio%PC_diatm                = Me%ExternalVar%Ratios(30)
           Me%ExternalParam%DiaNSatConst            = Me%ExternalVar%Ratios(31)
           
         elseif(Me%ExternalRatio%Cilia_are_calculated == 1)then
           Me%ExternalRatio%NC_diatm                = Me%ExternalVar%Ratios(29)  
           Me%ExternalRatio%PC_diatm                = Me%ExternalVar%Ratios(30) 
           Me%ExternalParam%DiaNSatConst            = Me%ExternalVar%Ratios(31)
           Me%ExternalRatio%NC_cilia                = Me%ExternalVar%Ratios(32)  
           Me%ExternalRatio%PC_cilia                = Me%ExternalVar%Ratios(33) 
         endif i2 
     elseif (Me%ExternalRatio%Cilia_are_calculated == 1)then
          Me%ExternalRatio%NC_cilia                = Me%ExternalVar%Ratios(29)  
          Me%ExternalRatio%PC_cilia                = Me%ExternalVar%Ratios(30) 
     endif i1
   
     
   !WQ configuration
       if (Me%ExternalRatio%Diatm_are_calculated == 0) then 
       
             if (Me%ExternalRatio%Cilia_are_calculated == 0) then 
                 
                 Me%PelagicModuleConfiguration = 1
                 
             elseif (Me%ExternalRatio%Cilia_are_calculated == 1) then  
                 
                 Me%PelagicModuleConfiguration = 3
             endif   
             
        elseif (Me%ExternalRatio%Diatm_are_calculated == 1) then 

             if (Me%ExternalRatio%Cilia_are_calculated == 0) then 
                 
                 Me%PelagicModuleConfiguration = 2
                 
             elseif (Me%ExternalRatio%Cilia_are_calculated == 1) then  
                 
                 Me%PelagicModuleConfiguration = 4
             endif                
       endif 
 !----------------------------------------------------------------------------
 end subroutine Biogeochemical_ratios_and_parameters
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
 
                                     
 !>@author Marta López Maretec
 !>@Brief: Calls the subroutines that do computations
 !>@param[in] index
 !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
                                     
  subroutine CS_computations(index)   
  
  !Arguments------------------------------------------------------------------  
   integer, intent(IN) :: index   !Element of the 1Darray
  !---------------------------------------------------------------------------- 
   
   SELECT CASE(Me%CSConfiguration)  !Depending on the choosen config on the .dat file
       
     CASE(1)                                           !Alk_param and DIC_no_calc
       
       call ComputeAlkalinity_param      (index)
       call ComputeDIC_no_calc           (index) 
       
     CASE(2)                                           !Alk_param and DIC_calc
       
        call ComputeAlkalinity_param     (index)      
        call ComputeDIC_calc             (index) 
        
     CASE(3)                                           !Alk_bio and DIC_no_cacl
       
         call ComputeAlkalinity_bio      (index)
         call ComputeDIC_no_calc         (index) 
         
     CASE(4)                                           !Alk_bio and DIC_calc
       
         call ComputeAlkalinity_bio_calc (index)
         call ComputeDIC_calc            (index)  
         
   END SELECT  
   
  ! Rest of carbonate system parameters calculation 
  ! i3: if () then             
  !      call mocsy (index)       
  !     endif i3
 !----------------------------------------------------------------------------
  end subroutine  CS_computations
 !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
 
  
  
 !>@author Marta López, Maretec
 !>@Brief:
 !> Calculates alkalinity values through algorithms proposed by Lee.et al 2006.
 !> Alkalinity units: umol/kg 
 !>@param[in] index 
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
  
  subroutine ComputeAlkalinity_param(index)
  
    !Arguments-------------------------------------------------------------
     integer, intent(IN) :: index       
    !Local-----------------------------------------------------------------     
     real                         :: ALK                !Alkalinity
     real                         :: Lat                !Latitude local
     real                         :: Long               !Longitude local
     real                         :: Temp               !Temperature  local       
     character(len=StringLength)  :: Ocean_Case         = null_str 
     !     'North_Atlantic'         30ºN-80ºN, 0ºC<SST<20ºC; 31<SSS<37 
     !     'North_Pacific'              >30ºN,     SST<20ºC; 31<SSS<35 
     !     'Subtropics'             30ºS-30ºN,     SST>20ºC; 31<SSS<38 
     !     'Southern_Ocean'         70ºS-30ºS,     SST<20ºC; 33<SSS<36
     !     'Eq_upwelling_Pacific'   20ºS-20ºN,75ºW-110ºW;10ºS-10ºN,110ºW-140ºW, ...             
     !     'Out_of_range_S'                        ... SST>18ºC; 31<SSS<36.5
     !     'Out_of_range_N'     
    !----------------------------------------------------------------------     
     Lat  = Me%ExternalVar%Latitude(index) 
     Long = Me%ExternalVar%Longitude(index)
     Temp = Me%ExternalVar%Temperature(index)     
    !----------------------------------------------------------------------  
     !Lee et al., 2006. Global relationships of TA  with salinity and temperature in surface waters:
     !  "For predicting At from SS and SST values in a given zone, the primary criterion for choosing 
     !  an appropiate AT algorithms is the SST. If the temperature of a water sample in the zone of 
     !  interest is out of the range of the applicable equation, the equation for an adjacent zone sharing
     !  a common boundary shoud be choosen" 
        
     !First select the ocean region
i1:     IF (((Lat .GE. 30.) .AND. (lat < 80.)) .and. ((Long .GE. -100.) .AND. (Long < 100.))) THEN
            Ocean_Case = 'North_Atlantic'        
        ELSEIF ((Lat .GE. 30.) .and. ((Long .GE. 100.) .AND. (Long < -100.))) THEN    
            Ocean_Case = 'North_Pacific'    
        ELSEIF ((Lat .GE. -30.) .AND. (lat < 30.)) THEN         
   i2:          IF (((Lat .GE. -20.) .AND. (lat < 20.)) .and. ((Long .GE. -110.) .AND. (Long < -75.))) THEN            
                    Ocean_Case = 'Eq_upwelling_Pacific'            
                ELSEIF (((Lat .GE. -10.) .AND. (lat < 10.)) .and. ((Long .GE. -140.) .AND. (Long < -110.))) THEN            
                    Ocean_Case = 'Eq_upwelling_Pacific'            
                ELSE                    
                    Ocean_Case = 'Subtropics'                
                ENDIF i2                 
        ELSEIF ((Lat .GE. -70.) .AND. (lat < -30.)) THEN    
            Ocean_Case = 'Southern_Ocean'            
        ELSEIF (Lat < -70.) THEN            
            Ocean_Case = 'Out_of_range_S'           
        ELSE              
            Ocean_Case = 'Out_of_range_N'            
        ENDIF i1
    
     !Second, evaluate if the temperature corresponds to the allowed region interval, if not use the equation 
     !from the adjacent zone     
   SELECT CASE(Ocean_Case)
         
     CASE ('North_Atlantic')
        !write(*,*)'Ocean_Case is North_Atlantic' 
        if (temp < 20.) then    
          ALK = 2305. + 53.97 * (Me%ExternalVar%Salinity(index) - 35.) + 2.74 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                   & - 1.16 *(Me%ExternalVar%Temperature(index) - 20.) - 0.040*((Me%ExternalVar%Temperature(index) - 20.)** 2.)
        else
        !write(*,*)'Ha cambiado a subtropics por la temperatura' 
         ALK = 2305. + 58.66 * (Me%ExternalVar%Salinity(index) - 35.) + 2.32 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                     & - 1.41 *(Me%ExternalVar%Temperature(index) - 20.) + 0.040*((Me%ExternalVar%Temperature(index) - 29.)** 2.) 
        endif    
          
        
     CASE ('North_Pacific')
      !write(*,*)'Ocean_Case is North_Pacific'
        if (temp < 20.) then 
         ALK = 2305. + 53.23 * (Me%ExternalVar%Salinity(index) - 35.) + 1.85 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                   & - 14.72 * (Me%ExternalVar%Temperature (index) - 20.)  &
                   & - 0.158 * ((Me%ExternalVar%Temperature(index) - 20.) ** 2.) &          
                   & + 0.062 * (Me%ExternalVar%Temperature (index) - 20.) * Me%ExternalVar%Longitude(index) 
        else
         ALK = 2305. + 58.66 * (Me%ExternalVar%Salinity(index) - 35.) + 2.32 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                     & - 1.41 *(Me%ExternalVar%Temperature(index) - 20.) + 0.040*((Me%ExternalVar%Temperature(index) - 29.)** 2.) 
        endif    
              
        
     CASE ('Subtropics')
      !write(*,*)'Ocean_Case is Suptropics'
i3:     if (temp > 20.) then         
          ALK = 2305. + 58.66 * (Me%ExternalVar%Salinity(index) - 35.) + 2.32 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                     & - 1.41 *(Me%ExternalVar%Temperature(index) - 20.) + 0.040*((Me%ExternalVar%Temperature(index) - 29.)** 2.) 
        else 
    
  i4:     if ((Long .GE. -100.) .AND. (Long < 100.)) then
            ALK = 2305. + 53.97 * (Me%ExternalVar%Salinity(index) - 35.) + 2.74 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                   & - 1.16 *(Me%ExternalVar%Temperature(index) - 20.) - 0.040*((Me%ExternalVar%Temperature(index) - 20.)** 2.)  
              
          else
            ALK = 2305. + 53.23 * (Me%ExternalVar%Salinity    (index) - 35.) + 1.85 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                      & - 14.72 * (Me%ExternalVar%Temperature (index) - 20.)        &
                      & - 0.158 * ((Me%ExternalVar%Temperature(index) - 20.) ** 2.) &          
                      & + 0.062 * (Me%ExternalVar%Temperature (index) - 20.) * Me%ExternalVar%Longitude(index)             
          endif i4
          
        endif i3         
          
                
     CASE ('Eq_upwelling_Pacific')
       !write(*,*)'Ocean_Case is Eq_upwelling_Pacific'
          ALK = 2294. + 64.88 * (Me%ExternalVar%Salinity(index) - 35.) + 0.39 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                    & -  4.52 * (Me%ExternalVar%Temperature (index) - 20.) &
                    & - 0.232 * ((Me%ExternalVar%Temperature(index) - 29.)** 2.)           
          
     CASE ('Southern_Ocean')
      !write(*,*)'Ocean_Case is Southern_Ocean'
        if (temp < 20.) then       
          ALK = 2305. + 52.48 * (Me%ExternalVar%Salinity(index) - 35.) + 2.85 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                     & - 0.49 *(Me%ExternalVar%Temperature(index) - 20.) + 0.086 *((Me%ExternalVar%Temperature(index) - 20.)** 2.)  
        else
          ALK = 2305. + 58.66 * (Me%ExternalVar%Salinity(index) - 35.) + 2.32 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                    & - 1.41 * (Me%ExternalVar%Temperature(index) - 20.) &
                    & + 0.040*((Me%ExternalVar%Temperature(index) - 29.)** 2.)            
        endif
        
        
     CASE('Out_of_range_S') 
          write(*,*)'The point with latitude', Lat
          write(*,*)'and longitude', Long
          write(*,*)'has not a representative algorithm for Alkalinity. It is going to be calculated through the algorithm of the '
          write(*,*)'Southern Ocean.It is recommended to use the biological alkalinity calculation option.Check it in the .datfile'
          ALK = 2305. + 52.48 * (Me%ExternalVar%Salinity(index) - 35.) + 2.85 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                     & - 0.49 *(Me%ExternalVar%Temperature(index) - 20.) + 0.086 *((Me%ExternalVar%Temperature(index) - 20.)** 2.)  
          
          
     CASE('Out_of_range_N') 
         write(*,*)'The point with latitude', Lat
         write(*,*)'and longitude', Long
         write(*,*)'has not a representative algorithm for Alkalinity. It is going to be calculated through the algorithm for the'
         write(*,*)'North Atlantic. It is recommended to use the biological alkalinity calculation option'          
           ALK = 2305. + 53.97 * (Me%ExternalVar%Salinity(index) - 35.) + 2.74 *((Me%ExternalVar%Salinity(index) - 35.) ** 2.) &
                     & - 1.16 *(Me%ExternalVar%Temperature(index) - 20.) - 0.040*((Me%ExternalVar%Temperature(index) - 20.)** 2.) 
           
     END SELECT
     
     ! Entender como es enviada a waterproperties, que no divida por el volumen! porque 
     ! aqui ya esta en umol/kg. Hay que multil por densidad??
     Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_p , Index) = ALK 
     !write(*,*) 'Alcalinidad (umol/kg) = ', ALK 
     !write(*,*)'latitud=', Lat
!----------------------------------------------------------------------------
  end subroutine ComputeAlkalinity_param
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  
 !>@author Marta López, Maretec
 !>@Brief:
 !> Calculates alkalinity values taking into account the changes due to 
 !> biological activity
 !> Units here: ueq (micro equivalents). Once it is sent to WaterProperties, units 
 !> are ueq/L.  
 !>@param[in] index, 
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
 subroutine ComputeAlkalinity_bio(index)
 
   !Arguments-----------------------------------------------------------------
     integer, intent(IN) :: index  
   !Local--------------------------------------------------------------------   
     real           :: Nitrif1       !Nitrification, first step 
     real           :: Nitrif2       !Nitrification, second step 
     real           :: Denit         !Denitrification  
     real           :: phyto_resp    !Phytoplankton organic matter aerobic mineralization
     real           :: zoo_resp      !Zooplankton organic matter aerobic mineralization
     real           :: diatm_resp    !Diatoms organic matter aerobic mineralization
     real           :: cil_resp      !Ciliates organic matter aerobic mineralization
     real           :: phyto_am_uptk !Phytoplankton ammonia uptake
     real           :: phyto_na_uptk !Phytoplankton nitrate uptake
     real           :: diatm_am_uptk !Diatoms ammonia uptake
     real           :: diatm_na_uptk !Diatoms nitrate uptake
     
     real           :: yN_pd         !Mean N/C stoichometry ratio used by phyto and diatoms
     real           :: yP_pd         !Mean N/C stoichometry ratio used by phyto and diatoms
     real           :: yN_p          !Phytoplankton N/C stoichometry ratio used in pelagic module
     real           :: yP_p          !Phytoplankton P/C stoichometry ratio used in pelagic module
     real           :: yN_d          !Diatoms N/C stoichometry ratio used in pelagic module
     real           :: yP_d          !Diatoms P/C stoichometry ratio used in pelagic module
     real           :: yN_z          !Zooplankton N/C stoichometry ratio used in pelagic module
     real           :: yP_z          !Zooplankton P/C stoichometry ratio used in pelagic module
     real           :: yN_c          !Ciliates N/C stoichometry ratio used in pelagic module
     real           :: yP_c          !Ciliates P/C stoichometry ratio used in pelagic module
  !----------------------------------------------------------------------------------------------    
     yN_p = Me%ExternalRatio%NC_phyto    
     yP_p = Me%ExternalRatio%PC_phyto
     yN_d = Me%ExternalRatio%NC_diatm
     yP_d = Me%ExternalRatio%PC_diatm
     yN_z = Me%ExternalRatio%NC_zoo
     yP_z = Me%ExternalRatio%PC_zoo  
     yN_c = Me%ExternalRatio%NC_cilia
     yP_c = Me%ExternalRatio%PC_cilia   
  !------------------------------------------------------------------------------------------------    
          
   call Compute_biogeoch_rates_nitr_denit (index, Nitrif1, Nitrif2, Denit)     
 
   select case(Me%PelagicModuleConfiguration)
        
   case(1) !Phyto and Zoo
       
        call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp) 
        call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk)       
               
        
        Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) &    
                                                           + (0.8 + yN_p - yP_p) * Denit                      & ! <- sources 
                                                           + (yN_p + yP_p) * phyto_na_uptk                    & 
                                                           + (yN_p - yP_p) * phyto_resp                       & 
                                                           + (yN_z - yP_z) * zoo_resp                         & 
                                                           - Nitrif1                                          &  ! <- sinks
                                                           - Nitrif2                                          &
                                                           - (yP_p + yN_p) * phyto_am_uptk
   case(2) !Phyto, Zoo and diatoms   
       
          yN_pd = (yN_p + yN_d) / 2. 
          yP_pd = (yP_p + yP_d) / 2.
          
          call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp, &
                                           diato_respiration = diatm_resp) 
          
          call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk, &
                                                      diatm_am_uptk, diatm_na_uptk )
          
          Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index)     &    
                                                              + (0.8 + yN_pd - yP_pd) * Denit                       &  ! <- sources
                                                              + (yN_p + yP_p) * phyto_na_uptk                       & 
                                                              + (yN_d + yP_d) * diatm_na_uptk                       &
                                                              + (yN_p - yP_p) * phyto_resp                          & 
                                                              + (yN_z - yP_z) * zoo_resp                            & 
                                                              + (yN_d - yP_d) * diatm_resp                          &                                                      
                                                              - Nitrif1                                             &  ! <- sinks
                                                              - Nitrif2                                             &
                                                              - (yP_p + yN_p) * phyto_am_uptk                       &   
                                                              - (yP_d + yN_d) * diatm_am_uptk)                          
    
   case(3) !Phyto, Zoo and ciliates
       
           call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp,& 
                                            cilia_respiration = cil_resp) 
          
           call Compute_biogeoch_rates_nitrogen_uptake(phyto_am_uptk, phyto_na_uptk)
                    
           Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index)  &    
                                                               + (0.8 + yN_p - yP_p) * Denit                      & ! <- sources 
                                                               + (yN_p + yP_p) * phyto_na_uptk                    & 
                                                               + (yN_p - yP_p) * phyto_resp                       & 
                                                               + (yN_z - yP_z) * zoo_resp                         & 
                                                               + (yN_c - yP_c) * cil_resp                         & 
                                                               - Nitrif1                                          &  ! <- sinks
                                                               - Nitrif2                                          &
                                                               - (yP_p + yN_p) * phyto_am_uptk                                                              
       
    case(4) !Phyto, Zoo, ciliates and diatoms
          
           yN_pd = (yN_p + yN_d) / 2. 
           yP_pd = (yP_p + yP_d) / 2.          
          
          call Compute_biogeoch_rates_resp(index,  phyto_resp,  zoo_resp,          &
                                                   diato_respiration = diatm_resp, & 
                                                   cilia_respiration = cil_resp  )
          
          call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk,  &
                                                             diatm_am_uptk, diatm_na_uptk )
          
          Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index)    &    
                                                              + (0.8 + yN_pd - yP_pd) * Denit                      &  ! <- sources
                                                              + (yN_p + yP_p) * phyto_na_uptk                      & 
                                                              + (yN_d + yP_d) * diatm_na_uptk                      &
                                                              + (yN_p - yP_p) * phyto_resp                         & 
                                                              + (yN_z - yP_z) * zoo_resp                           & 
                                                              + (yN_d - yP_d) * diatm_resp                         &  
                                                              + (yN_c - yP_c) * cil_resp                           & 
                                                              - Nitrif1                                            &  ! <- sinks
                                                              - Nitrif2                                            &
                                                              - (yP_p + yN_p) * phyto_am_uptk                      &   
                                                              - (yP_d + yN_d) * diatm_am_uptk)  
          
          
    end select
     ! To couple module LIFE!
       !Facer um get como o getWQparameters do WaterQuality,para crear um array com os mesmos valores e dimension. 
       !Acrescentar ese get na subroutina GetWQRatio do interface. No couple carbonatesystem do wproperties facer escolher
       !que get usar dependendo do modulo que tenha sido activo. Também facer com o dt (ver codigo do carbonatesystem no 
       !wproperties . Quando ficar pronto na subrt read data deste modulo mudar o stop cuando o life estar activo      
!----------------------------------------------------------------------------     
 end subroutine ComputeAlkalinity_bio
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 
 
 !>@author Marta López, Maretec
 !>@Brief:
 !> Calculates alkalinity values taking into account the changes due to 
 !> biological activity, including calcium carbonate processes
 !> Units here: ueq (micro equivalents). Once it is sent to WaterProperties, units 
 !> are ueq/L.  
 !> param[in] index
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
 subroutine ComputeAlkalinity_bio_calc(index)
 
   !Arguments-----------------------------------------------------------------
     integer, intent(IN) :: index  
   !Local--------------------------------------------------------------------   
     real           :: Nitrif1         !Nitrif1 
     real           :: Nitrif2         !Nitrif2 
     real           :: Denit           !Denitrification  
     real           :: phyto_resp      !Phytoplankton organic matter aerobic mineralization
     real           :: zoo_resp        !Zooplankton organic matter aerobic mineralization
     real           :: diatm_resp      !Diatoms organic matter aerobic mineralization
     real           :: cil_resp        !Ciliates organic matter aerobic mineralization
     real           :: phyto_am_uptk   !Phytoplankton ammonia uptake
     real           :: phyto_na_uptk   !Phytoplankton nitrate uptake
     real           :: diatm_am_uptk   !Diatoms ammonia uptake
     real           :: diatm_na_uptk   !Diatoms nitrate uptake
     !real           :: caco3disolution !CaCO3 disolution 
     !real           :: caco3precipitat !CaCO3 precipitation
     
     
     real           :: yN_pd         !Mean N/C stoichometry ratio used by phyto and diatoms
     real           :: yP_pd         !Mean N/C stoichometry ratio used by phyto and diatoms
     real           :: yN_p          !Phytoplankton N/C stoichometry ratio used in pelagic module
     real           :: yP_p          !Phytoplankton P/C stoichometry ratio used in pelagic module
     real           :: yN_d          !Diatoms N/C stoichometry ratio used in pelagic module
     real           :: yP_d          !Diatoms P/C stoichometry ratio used in pelagic module
     real           :: yN_z          !Zooplankton N/C stoichometry ratio used in pelagic module
     real           :: yP_z          !Zooplankton P/C stoichometry ratio used in pelagic module
   !----------------------------------------------------------------------------------------------    
     yN_p = Me%ExternalRatio%NC_phyto    
     yP_p = Me%ExternalRatio%PC_phyto
     yN_d = Me%ExternalRatio%NC_diatm
     yP_d = Me%ExternalRatio%PC_diatm
     yN_z = Me%ExternalRatio%NC_zoo
     yP_z = Me%ExternalRatio%PC_zoo    
   !------------------------------------------------------------------------------------------------    
          
   call Compute_biogeoch_rates_nitr_denit (index, Nitrif1, Nitrif2, Denit) 
   !call Compute_biogeoch_rates_caco3(index, caco3disolution, caco3precipitat)
   
   select case(Me%PelagicModuleConfiguration)
        
   case(1) !Phyto and Zoo
       
       call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp) 
       call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk)
     
      Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) & 
                                                              + (0.8 + yN_p - yP_p) * Denit                 & ! <- sources                                  
                                                              + (yN_p + yP_p) * phyto_na_uptk               & 
                                                              + (yN_p - yP_p) * phyto_resp                  & 
                                                              + (yN_z - yP_z) * zoo_resp                    & 
                                                              !+ 2. * CaCo3Disolution                       &
                                                              - Nitrif1                                     &  ! <- sinks
                                                              - Nitrif2                                     &
                                                              -(yP_p + yN_p) * phyto_am_uptk 
                                                              !- 2. * CaCo3Precipitation        
                                                           
   case(2) !Phyto, Zoo and diatoms   
       
          yN_pd = (yN_p + yN_d) / 2. 
          yP_pd = (yP_p + yP_d) / 2.
          
          call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp, &
                                           diato_respiration = diatm_resp) 
          
          call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk, &
                                                      diatm_am_uptk, diatm_na_uptk )
          
          Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index)    &    
                                                              + (0.8 + yN_pd - yP_pd) * Denit                        & ! <- sources 
                                                              + (yN_p + yP_p) * phyto_na_uptk                      & 
                                                              + (yN_d + yP_d) * diatm_na_uptk                      &
                                                              + (yN_p - yP_p) * phyto_resp                         & 
                                                              + (yN_z - yP_z) * zoo_resp                           & 
                                                              + (yN_d - yP_d) * diatm_resp                         &  
                                                             !+  2. * CaCo3Disolution                              &
                                                              - Nitrif1                                            &  ! <- sinks
                                                              - Nitrif2                                            &
                                                              - (yP_p + yN_p) * phyto_am_uptk                      !&    
                                                              - (yP_d + yN_d) * diatm_am_uptk)                     &
                                                             !- 2. * CaCo3Precipitation  
    
    case(3) !Phyto, Zoo and ciliates
                                                              
          call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp,& 
                                           cilia_respiration = cil_resp) 
          
          call Compute_biogeoch_rates_nitrogen_uptake(phyto_am_uptk, phyto_na_uptk)
                    
          Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) &    
                                                              + (0.8 + yN_p - yP_p) * Denit                     & ! <- sources 
                                                              + (yN_p + yP_p) * phyto_na_uptk                   & 
                                                              + (yN_p - yP_p) * phyto_resp                      & 
                                                              + (yN_z - yP_z) * zoo_resp                        & 
                                                              + (yN_c - yP_c) * cil_resp                        & 
                                                             !+  2. * CaCo3Disolution                           &
                                                              - Nitrif1                                         &  ! <- sinks
                                                              - Nitrif2                                         &           
                                                              - (yP_p + yN_p) * phyto_am_uptk                   &    
                                                             !- 2. * CaCo3Precipitation                           
       
    case(4) !Phyto, Zoo, ciliates and diatoms
          
           yN_pd = (yN_p + yN_d) / 2. 
           yP_pd = (yP_p + yP_d) / 2.
          
           call Compute_biogeoch_rates_resp(index, phyto_resp,  zoo_resp,          &
                                                   diato_respiration = diatm_resp, & 
                                                   cilia_respiration = cil_resp  )
          
           call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk,  &
                                                              diatm_am_uptk, diatm_na_uptk )
          
          Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index) = Me%ExternalVar%Mass(Me%PropIndex%ALK_cs_b, Index)    &    
                                                              + (0.8 + yN_pd - yP_pd) * Denit                      &  ! <- sources
                                                              + (yN_p + yP_p) * phyto_na_uptk                      & 
                                                              + (yN_d + yP_d) * diatm_na_uptk                      &
                                                              + (yN_p - yP_p) * phyto_resp                         & 
                                                              + (yN_z - yP_z) * zoo_resp                           & 
                                                              + (yN_d - yP_d) * diatm_resp                         &  
                                                              + (yN_c - yP_c) * cil_resp                           & 
                                                              !+  2. * CaCo3Disolution                              &
                                                              - Nitrif1                                            &  ! <- sinks
                                                              - Nitrif2                                            &
                                                              - (yP_p + yN_p) * phyto_am_uptk                      &   
                                                              - (yP_d + yN_d) * diatm_am_uptk)                     &
                                                             !- 2. * CaCo3Precipitation           
          
   end select
       
!----------------------------------------------------------------------------     
 end subroutine ComputeAlkalinity_bio_calc
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

                                               
 
 !>@author Marta López, Maretec
 !>@Brief:
 !> Calculates dissolved inorganic carbon (DIC) values taking into account the changes
 !> due to biological activity, including calcium carbonate processes
 !> Units here: umolC (micromolesC ). Once it is sent to WaterProperties, units 
 !> are umolC/L.  
 !> param[in] index                          
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
subroutine ComputeDIC_calc(index)
 !Arguments-------------------------------------------------------------------
    integer, intent(IN) :: index  
 !Local--------------------------------------------------------------------   
     real           :: phyto_resp      !Phytoplankton organic matter aerobic mineralization
     real           :: zoo_resp        !Zooplankton organic matter aerobic mineralization
     real           :: diatm_resp      !Diatoms organic matter aerobic mineralization
     real           :: cil_resp        !Ciliates organic matter aerobic mineralization
     real           :: phyto_am_uptk   !Phytoplankton ammonia uptake
     real           :: phyto_na_uptk   !Phytoplankton nitrate uptake
     real           :: diatm_am_uptk   !Diatoms ammonia uptake
     real           :: diatm_na_uptk   !Diatoms nitrate uptake   
     !real           :: caco3disolution    !CaCO3 disolution 
     !real           :: caco3precipitation !CaCO3 precipitation  
     real           :: Denit           !Denitrification           
     real           :: yN_pd         !Mean N/C stoichometry ratio used by phyto and diatoms    
     real           :: yN_p          !Phytoplankton N/C stoichometry ratio used in pelagic module
     real           :: yN_d          !Diatoms N/C stoichometry ratio used in pelagic module
   !----------------------------------------------------------------------------------------------    
    !yN_p = Me%ExternalRatio%NC_phyto  
    !yN_d = Me%ExternalRatio%NC_diatm
  !---------------------------------------------------------------------------------------------- 
     
   !call Compute_biogeoch_rates_nitr_denit (index, Nitrif1, Nitrif2, Denit) 
   !call Compute_biogeoch_rates_caco3(index, caco3disolution, caco3precipitation)
   
   select case(Me%PelagicModuleConfiguration)
        
    case(1) !Phyto and Zoo
       
       call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp) 
       call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk)   
          Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_c, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_c, Index) &    
                                                           !  + caco3disolution                  &           ! <- sources 
                                                           !  + (Denit / 0.8)                    &    
                                                             + phyto_resp                       & 
                                                             + zoo_resp                         & 
                                                            !- caco3precipit                    &            ! <- sinks   
                                                             - phyto_am_uptk                    &
                                                             - phyto_na_uptk
        
   case(2) !Phyto, Zoo and diatoms        

          yN_pd = (yN_p + yN_d) / 2. 
        
          call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp, &
                                                   diato_respiration = diatm_resp)           
          call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk , &
                                                              diatm_am_uptk, diatm_na_uptk )
          
          Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_c, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_c, Index)     & 
                                                             !+ caco3disolution                    &           ! <- sources 
                                                              + phyto_resp                         &  
                                                              + zoo_resp                           & 
                                                              + diatm_resp                         &                                                      
                                                             !- caco3precipit                      &           ! <- sinks                                                          
                                                              - phyto_am_uptk                      &     
                                                              - diatm_am_uptk                      &        
                                                              - phyto_na_uptk                      &   
                                                              - diatm_na_uptk                     
    
    case(3) !Phyto, Zoo and ciliates
          call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp,& 
                                           cilia_respiration = cil_resp)           
          call Compute_biogeoch_rates_nitrogen_uptake(phyto_am_uptk, phyto_na_uptk)
                    
          Me%ExternalVar%Mass(Me%PropIndex%%DIC_cs_c, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_c, Index) &    
                                                               !+ caco3disolution                     &          ! <- sources 
                                                               + phyto_resp                          &  
                                                               + zoo_resp                            & 
                                                               + cil_resp                            &                                                      
                                                              !- caco3precipit                       &           ! <- sinks                                                          
                                                               - phyto_am_uptk                       &    
                                                               - phyto_na_uptk                       &   
                                                                                                                       
       
    case(4) !Phyto, Zoo, ciliates and diatoms    
          
            yN_pd = (yN_p + yN_d) / 2.         
          
           call Compute_biogeoch_rates_resp(index,  phyto_resp,  zoo_resp,          &
                                                    diato_respiration = diatm_resp, & 
                                                    cilia_respiration = cil_resp  )
          
           call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk,  &
                                                               diatm_am_uptk, diatm_na_uptk )
           
           Me%ExternalVar%Mass(Me%PropIndex%%DIC_cs_c, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_c, Index)  &    
                                                              !+ caco3disolution                    &         ! <- sources                       
                                                              + phyto_resp                         & 
                                                              + zoo_resp                           & 
                                                              + diatm_resp                         &  
                                                              + cil_resp                           & 
                                                             !- caco3precipit                      &           ! <- sinks                                                          
                                                              - phyto_am_uptk                      &    
                                                              - diatm_am_uptk                      &   
                                                              - phyto_na_uptk                      &  
                                                              - diatm_na_uptk            
          
   end select                                                           

!----------------------------------------------------------------------------
end subroutine ComputeDIC_calc
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


 !>@author Marta López, Maretec
 !>@Brief:
 !> Calculates dissolved inorganic carbon (DIC) values not taking into account the changes
 !> due to calcium carbonate processes
 !> Units here: umolC (micromolesC ). Once it is sent to WaterProperties, units 
 !> are umolC/L.  
 !> param[in] index 
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 subroutine ComputeDIC_no_calc(index)
 
 !Arguments-------------------------------------------------------------------
     integer, intent(IN) :: index  
 !Local--------------------------------------------------------------------   
     real           :: phyto_resp      !Phytoplankton organic matter aerobic mineralization
     real           :: zoo_resp        !Zooplankton organic matter aerobic mineralization
     real           :: diatm_resp      !Diatoms organic matter aerobic mineralization
     real           :: cil_resp        !Ciliates organic matter aerobic mineralization
     real           :: phyto_am_uptk   !Phytoplankton ammonia uptake
     real           :: phyto_na_uptk   !Phytoplankton nitrate uptake
     real           :: diatm_am_uptk   !Diatoms ammonia uptake
     real           :: diatm_na_uptk   !Diatoms nitrate uptake
     real           :: yN_pd           !Mean N/C stoichometry ratio used by phyto and diatoms    
     real           :: yN_p            !Phytoplankton N/C stoichometry ratio used in pelagic module
     real           :: yN_d            !Diatoms N/C stoichometry ratio used in pelagic module     
  !---------------------------------------------------------------------------------------------- 
 
  select case(Me%PelagicModuleConfiguration)
        
   case(1) !Phyto and Zoo
       
        call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp) 
        call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk)   
        Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index) & 
                                                            + phyto_resp                       &  ! <- sources 
                                                            + zoo_resp                         & 
                                                            - phyto_am_uptk                    &  ! <- sinks   
                                                            - phyto_na_uptk
   
   case(2) !Phyto, Zoo and diatoms        

         yN_pd = (yN_p + yN_d) / 2. 
         
         call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp,  &
                                                              diato_respiration = diatm_resp)           
         call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk, &
                                                              diatm_am_uptk, diatm_na_uptk )
          
         Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index)     &
                                                             + phyto_resp                        &   ! <- sources 
                                                             + zoo_resp                          & 
                                                             + diatm_resp                        &                      
                                                             - phyto_am_uptk                     &   ! <- sinks 
                                                             - diatm_am_uptk                     &        
                                                             - phyto_na_uptk                     &   
                                                             - diatm_na_uptk                     
    
   case(3) !Phyto, Zoo and ciliates
       
         call Compute_biogeoch_rates_resp(index, phyto_resp, zoo_resp,& 
                                           cilia_respiration = cil_resp)           
         call Compute_biogeoch_rates_nitrogen_uptake(phyto_am_uptk, phyto_na_uptk)
                    
         Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index) &       
                                                             + phyto_resp                          &     ! <- sources 
                                                             + zoo_resp                            & 
                                                             + cil_resp                            &                                                                  
                                                             - phyto_am_uptk                       &     ! <- sinks  
                                                             - phyto_na_uptk                       &   
                                                                                                                       
       
   case(4) !Phyto, Zoo, ciliates and diatoms 
             
          yN_pd = (yN_p + yN_d) / 2.         
          
          call Compute_biogeoch_rates_resp(index, phyto_resp,  zoo_resp,          &
                                                  diato_respiration = diatm_resp, & 
                                                  cilia_respiration = cil_resp  )
          
          call Compute_biogeoch_rates_nitrogen_uptake(index, phyto_am_uptk, phyto_na_uptk,  &
                                                             diatm_am_uptk, diatm_na_uptk )
           
          Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index) = Me%ExternalVar%Mass(Me%PropIndex%DIC_cs_nc, Index)  & 
                                                               + phyto_resp                         &    ! <- sources 
                                                               + zoo_resp                           & 
                                                               + diatm_resp                         &  
                                                               + cil_resp                           &                                                                                                              
                                                               - phyto_am_uptk                      &    ! <- sinks      
                                                               - diatm_am_uptk                      &   
                                                               - phyto_na_uptk                      &  
                                                               - diatm_na_uptk            
          
   end select                                                         

!----------------------------------------------------------------------------
 end subroutine ComputeDIC_no_calc
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 !>@author Marta López, Maretec
 !>@Brief: Calculate several biogeochemical rates for each grid cell, in every
 !>time step
 !>@param[in] index
 !>@param[out] nitrification1, nitrification2, denitrification
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
 
 subroutine Compute_biogeoch_rates_nitr_denit(index, nitrification1, &
                                                     nitrification2, &
                                                     denitrification )
 ! Arguments -----------------------------------------------------------------
     integer, intent(IN )  :: index 
     real   , intent(OUT)  :: nitrification1
     real   , intent(OUT)  :: nitrification2
     real   , intent(OUT)  :: denitrification         
 ! Local ---------------------------------------------------------------------
     real           :: x5
     real           :: x1
     real           :: NitrificationRateK1
     real           :: NitrificationRateK2
     real           :: DenitrificationRate
     real           :: Pelagic_module_DTDay 
     real           :: Pelagic_module_DT   
     integer        :: O, AM, NA, NI
    
 !--------------------------------------------------------------------------       
     Pelagic_module_DTDay = Me%ExternalParam%pelagic_module_dt_day
        Pelagic_module_DT = Me%ExternalParam%pelagic_module_dt     
                       AM = Me%PropIndex%AM
                       NA = Me%PropIndex%NA
                       NI = Me%PropIndex%NI
                        O = Me%PropIndex%Oxygen
 !--------------------------------------------------------------------------                   

     
      x5 = MAX(Me%ExternalVar%Mass(O, index),Me%ExternalParam%MinOxygen)                                 &
           /(Me%ExternalParam%NitrificationSatConst + Me%ExternalVar%Mass(O, index))
      
      x1 = Me%ExternalParam%DenitrificationSatConst / (Me%ExternalParam%DenitrificationSatConst          &
                                        + MAX(Me%ExternalVar%Mass(O, index),Me%ExternalParam%MinOxygen))
      
            
      NitrificationRateK1 = Me%ExternalParam%KNitrificationRateK1 * Me%ExternalParam%TNitrification      &
                            **(Me%ExternalVar%Temperature(index) - 20.0) * x5
     
      NitrificationRateK2 = Me%ExternalParam%KNitrificationRateK2 * Me%ExternalParam%TNitrification      &
                            **(Me%ExternalVar%Temperature(index) - 20.0) * x5    
      
      DenitrificationRate = Me%ExternalParam%KDenitrificationRate * Me%ExternalParam%TDenitrification    &
                            ** (Me%ExternalVar%Temperature(index) - 20.) * x1               
      
       
      !  [mgN/l]      =        DtDay        *      1/Days          *          [mgN/l] 
      nitrification1  = NitrificationRateK1 * Pelagic_module_DTDay * Me%ExternalVar%Mass(AM,Index)
      nitrification2  = NitrificationRateK2 * Pelagic_module_DTDay * Me%ExternalVar%Mass(NI,Index)
      denitrification = DenitrificationRate * Pelagic_module_DTDay * Me%ExternalVar%Mass(NA,Index)
      
      ![mgN/dt]       =     [mgN/l]     *      [l/cell]                  /        [dt] 
      nitrification1  = nitrification1  * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT
      nitrification2  = nitrification2  * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT
      denitrification = denitrification * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT    
            
      ![mmolN/dt]     =   [mgN/dt]      /     [mgN/mmolN]          
      nitrification1  = nitrification1  / Me%AuxParam%N_AtomicMass
      nitrification2  = nitrification2  / Me%AuxParam%N_AtomicMass
      denitrification = denitrification / Me%AuxParam%N_AtomicMass
       
      ![umolN/dt]     = [mmolN/dt]      *  10^3
      nitrification1  = nitrification1  / 0.001
      nitrification2  = nitrification2  / 0.001
      denitrification = denitrification / 0.001
      
 !----------------------------------------------------------------------------
 end subroutine Compute_biogeoch_rates_nitr_denit
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
 

 !>@author Marta López, Maretec
 !>@Brief: Calculate several biogeochemical rates for each grid cell, in every
 !>time step
 !>@param[in] index
 !>@param[out] phyto, zoo, diatoms and ciliates respiration in umolC
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
 
  subroutine Compute_biogeoch_rates_resp (index, phyto_respiration &
                                               , zoopl_respiration &
                                               , diato_respiration & 
                                               , cilia_respiration )
  ! Arguments -----------------------------------------------------------------
     integer, intent(IN )            :: index 
     real   , intent(OUT)            :: phyto_respiration
     real   , intent(OUT)            :: zoopl_respiration
     real   , intent(OUT), optional  :: diato_respiration
     real   , intent(OUT), optional  :: cilia_respiration        
  ! Local ---------------------------------------------------------------------     
     real           :: PhyRespirationRate
     real           :: ZooRespirationRate
     real           :: DiaRespirationRate
     real           :: CilRespirationRate     
     real           :: Pelagic_module_DT
     real           :: Pelagic_module_DTDay          
     integer        :: Phyto, Zoo, O
     real           :: TZooLimitationFactor 
     real           :: PhotorespirationRate, PhytoEndogenousRepiration
     real           :: s1, s2, ya, yb, xa, xb     
  !-------------------------------------------------------------------------- 
     Pelagic_module_DTDay = Me%ExternalParam%pelagic_module_dt_day
        Pelagic_module_DT = Me%ExternalParam%pelagic_module_dt
                   Phyto  = Me%PropIndex%Phyto 
                     Zoo  = Me%PropIndex%Zoo
                        O = Me%PropIndex%Oxygen
  !---------------------------------------------------------------------------
     
   !Case 1: Phyto and Zoo 
     
    !ZOO 
        s1 = (1. / (Me%ExternalParam%TOptZooMin - Me%ExternalParam%TZooMin)) * &
               log((Me%ExternalParam%ZK2 * (1.0 - Me%ExternalParam%ZK1))       &
                 / (Me%ExternalParam%ZK1 * (1.0 - Me%ExternalParam%ZK2)))

        s2 = (1. / (Me%ExternalParam%TZooMax - Me%ExternalParam%TOptZooMax)) * &
               log((Me%ExternalParam%ZK3 * (1.0 - Me%ExternalParam%ZK4))       &
                 / (Me%ExternalParam%ZK4 * (1.0 - Me%ExternalParam%ZK3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Me%ExternalParam%TZooMin))
        yb = exp(s2 * (Me%ExternalParam%TZooMax - Me%ExternalVar%Temperature(index)))

        xa = (Me%ExternalParam%ZK1 * ya) / (1.0 + Me%ExternalParam%ZK1 * (ya - 1.0))
        xb = (Me%ExternalParam%ZK4 * yb) / (1.0 + Me%ExternalParam%ZK4 * (yb - 1.0))

        TZooLimitationFactor = xa * xb 
     
      !PHYTO      
        PhytoEndogenousRepiration = Me%ExternalParam%PhytoEndogRepConst * exp(0.069 * Me%ExternalVar%Temperature(index))
      
         if (MAX(Me%ExternalVar%Mass(O, index),Me%ExternalParam%MinOxygen).eq.Me%ExternalParam%MinOxygen)then
           PhytoEndogenousRepiration = -1.0 / null_real
         endif    
         PhotorespirationRate = Me%ExternalParam%PhotorespFactor * Me%ExternalVar%GrossGrowthRate(index)     
            
     !!!!
     ZooRespirationRate = Me%ExternalParam%ZooReferenceRespRate * TZooLimitationFactor
     PhyRespirationRate = PhytoEndogenousRepiration + PhotorespirationRate  
     
    ![mgC/l]            =       (DtDay)      *      (1/Days)         *               [mgC/l]   
     zoopl_respiration  = ZooRespirationRate * Pelagic_module_DTDay  *  Me%ExternalVar%Mass(Zoo  ,Index) 
     phyto_respiration  = PhyRespirationRate * Pelagic_module_DTDay  *  Me%ExternalVar%Mass(Phyto,Index)  
     
    ![mgC/dt]          =      [mgC/l]      *      [l/cell]                  /        [dt] 
     phyto_respiration = phyto_respiration * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT 
     zoopl_respiration = zoopl_respiration * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT 
    
    ![mmolC/dt]        =   [mgC/dt]        /     [mgC/mmolC]          
     phyto_respiration = phyto_respiration / Me%AuxParam%C_AtomicMass 
     zoopl_respiration = zoopl_respiration / Me%AuxParam%C_AtomicMass
     
    ! [umolC/dt]        = [mmolC/dt]        * 10^3
     phyto_respiration = phyto_respiration / 0.001
     zoopl_respiration = zoopl_respiration / 0.001 

  if(Me%PelagicModuleConfiguration .ne. 1) then
      
     select case (Me%PelagicModuleConfiguration)
     case(2)
          !DiaRespirationRate  = * Me%ExternalVar%GrossGrowthRateDiat(index) 
         
          !diato_respiration = DiaRespirationRate * Pelagic_module_DTDay  *  Me%ExternalVar%Mass(Diatoms ,Index)
          !diato_respiration = diato_respiration  * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT
          !diato_respiration = diato_respiration  / (Me%AuxParam%C_AtomicMass * 0.001)
     case(3)
          !CiliaRespirationRate =
          
          !cilia_respiration = CiliaRespirationRate * Pelagic_module_DTDay  *  Me%ExternalVar%Mass(Ciliates ,Index)
          !cilia_respiration = cilia_respiration * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT
          !cilia_respiration = cilia_respiration / (Me%AuxParam%C_AtomicMass * 0.001)
     case(4)
          !DiaRespirationRate   =  * Me%ExternalVar%GrossGrowthRateDiat(index) 
          !CiliaRespirationRate =
         
          !diato_respiration = DiaRespirationRate * Pelagic_module_DTDay  *  Me%ExternalVar%Mass(Diatoms,Index)
          !diato_respiration = diato_respiration  * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT
          !diato_respiration = diato_respiration / (Me%AuxParam%C_AtomicMass * 0.001)
         
          !cilia_respiration = CiliaRespirationRate * Pelagic_module_DTDay  *  Me%ExternalVar%Mass(Ciliates,Index)
          !cilia_respiration = cilia_respiration * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT
          !cilia_respiration = cilia_respiration / (Me%AuxParam%C_AtomicMass * 0.001)
         
          
     end select  
      
  endif       
  
!-----------------------------------------------------------------------------
 end subroutine Compute_biogeoch_rates_resp
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



 !>@author Marta López, Maretec
 !>@Brief: Calculate ammonia and nitrate phytoplankton and (if present diatoms)
 !> uptake due to primary production in each grid cell, for every time step   
 !>@param[in]  index
 !>@param[out] phyto_ammonia, phyto_ammonia, diato_ammonia, diato_nitrate 
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
 
 subroutine Compute_biogeoch_rates_nitrogen_uptake (index , phyto_ammonia &
                                                          , phyto_nitrate &
                                                          , diato_ammonia & 
                                                          , diato_nitrate )
 ! Arguments -----------------------------------------------------------------
     integer, intent(IN )            :: index 
     real   , intent(OUT)            :: phyto_ammonia
     real   , intent(OUT)            :: phyto_nitrate
     real   , intent(OUT), optional  :: diato_ammonia
     real   , intent(OUT), optional  :: diato_nitrate   
!-----------------------------------------------------------------------------
     
     real    :: AmmoniaPreferenceFactor
     real    :: DiaAmmoniaPreferenceFactor
     real    :: x1,   x2,  x3,  x4
     real    :: dx2, dx3, dx4
     real    :: Pelagic_module_DT
     real    :: Pelagic_module_DTDay
     real    :: NSatConst 
     real    :: DiaNSatConst 
     real    :: phyto_ammonia_uptk
     real    :: phyto_nitrate_uptk
     real    :: diato_ammonia_uptk
     real    :: diato_nitrate_uptk
     integer :: AM, NA, Phyto, Zoo
!-----------------------------------------------------------------------------  
     
     Pelagic_module_DTDay = Me%ExternalParam%pelagic_module_dt_day
        Pelagic_module_DT = Me%ExternalParam%pelagic_module_dt
                NSatConst = Me%ExternalParam%NSatConst
             DiaNSatConst = Me%ExternalParam%DiaNSatConst
                       AM = Me%PropIndex%AM
                       NA = Me%PropIndex%NA
                    Phyto = Me%PropIndex%Phyto
                      Zoo = Me%PropIndex%Zoo 
!-----------------------------------------------------------------------------  

     ! Phyto            
     x1 = Me%ExternalVar%Mass(AM,Index) * Me%ExternalVar%Mass(NA,Index)
     x2 = (NSatConst + Me%ExternalVar%Mass(AM, Index)) * (NSatConst + Me%ExternalVar%Mass(NA,Index))
     x3 = NSatConst * Me%ExternalVar%Mass(AM, Index)
     x4 = (Me%ExternalVar%Mass(AM ,Index) + Me%ExternalVar%Mass(NA,Index)) * (NSatConst + Me%ExternalVar%Mass(NA,Index))

       if ((x1 .EQ. 0.0) .AND. (x3 .EQ. 0.0)) then
           AmmoniaPreferenceFactor = 0.0                 
       else 
           AmmoniaPreferenceFactor = (x1 / x2) + (x3 / x4)
       end if         
   
     !    (DtDay)       =      (dimensionless)           *      (d-1)     
     phyto_ammonia_uptk = AmmoniaPreferenceFactor        * Me%ExternalVar%GrossGrowthRateDiat(index)                  
     phyto_nitrate_uptk = (1. - AmmoniaPreferenceFactor) * Me%ExternalVar%GrossGrowthRateDiat(index)
       

     ! [mgC/l]         =       (DtDay)      *      (1/Days)        *               [mgC/l]  
     phyto_ammonia     = phyto_ammonia_uptk * Pelagic_module_DTDay * Me%ExternalVar%Mass(Me%PropIndex%phyto,Index)
     phyto_nitrate     = phyto_nitrate_uptk * Pelagic_module_DTDay * Me%ExternalVar%Mass(Me%PropIndex%phyto,Index)   
           
     ! [mgC/dt]        =      [mgC/l]    *      (l/cell)                  /        (dt) 
     phyto_ammonia     =  phyto_ammonia  * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT 
     phyto_nitrate     =  phyto_nitrate  * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT 
    
     ! [mmolC/dt]      =   [mgC/dt]        /     [mgC/mmolC]          
     phyto_ammonia     = phyto_ammonia     / Me%AuxParam%C_AtomicMass 
     phyto_nitrate     = phyto_nitrate     / Me%AuxParam%C_AtomicMass
     
     ! [umolC/dt]      = [mmolC/dt]    * 10^3
     phyto_ammonia     = phyto_ammonia / 0.001
     phyto_nitrate     = phyto_nitrate / 0.001
       
       
     ! Diatoms  
       
     if(present(diato_ammonia))then
          
         x2 = (DiaNSatConst + Me%ExternalVar%Mass(AM,Index)) &
                        * (DiaNSatConst + Me%ExternalVar%Mass(NA,Index))
         x3 = DiaNSatConst * Me%ExternalVar%Mass(AM,Index)
         x4 = (Me%ExternalVar%Mass(AM,Index) + Me%ExternalVar%Mass(NA,Index))  &
                         * (DiaNSatConst + Me%ExternalVar%Mass(NA,Index))

        if ((x1 .EQ. 0.0) .AND. (x3 .EQ. 0.0)) then
                DiaAmmoniaPreferenceFactor = 0.0                 
        else 
                DiaAmmoniaPreferenceFactor = (x1 / x2) + (x3 / x4)
        end if     
       
     !    (DtDay)         =      (dimensionless)               *      (d-1)   
      diato_ammonia_uptk  = DiaAmmoniaPreferenceFactor         !* Me%ExternalVar%DiaGrossGrowRate (index)
      diato_nitrate_uptk  = (1. - DiaAmmoniaPreferenceFactor)  !* Me%ExternalVar%DiaGrossGrowRate (index)        
          
     ! [mgC/l]         =       (DtDay)      *      (1/Days)        *               [mgC/l]  
     diato_ammonia     = diato_ammonia_uptk * Pelagic_module_DTDay * Me%ExternalVar%Mass(Me%PropIndex%Diatoms,Index)
     diato_nitrate     = diato_nitrate_uptk * Pelagic_module_DTDay * Me%ExternalVar%Mass(Me%PropIndex%Diatoms,Index)   
           
     ! [mgC/dt]        =      [mgC/l]    *      (l/cell)                  /        (dt) 
     diato_ammonia     =  diato_ammonia  * Me%ExternalVar%VolumenZ(index)  / Pelagic_module_DT 
     phyto_nitrate     =  diato_nitrate  * Me%ExternalVar%VolumenZ(index) / Pelagic_module_DT 
    
     ! [mmolC/dt]      =   [mgC/dt]        /     [mgC/mmolC]          
     diato_ammonia     = diato_ammonia     / Me%AuxParam%C_AtomicMass 
     diato_nitrate     = diato_nitrate     / Me%AuxParam%C_AtomicMass
     
     ! [umolC/dt]      = [mmolC/dt]    * 10^3
     diato_ammonia     = diato_ammonia / 0.001
     diato_nitrate     = diato_nitrate / 0.001
          
     endif
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
  end subroutine Compute_biogeoch_rates_nitrogen_uptake
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                                                          
                                                          
                                                          
 !>@author Marta López, Maretec
 !>@Brief: Calculate calcium carbonate precipitation and dissolution
 !>in each grid cell, for every time step. 
 !> Assumptions: calcium carbonate only in calcite form, fixed  ...                                                           
 !>@param[in]  index
 !>@param[out] 
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
 
 subroutine Compute_biogeoch_rates_caco3 (index , caco3disolution &
                                                , caco3precipitation ) 
                                                          
 ! Arguments -----------------------------------------------------------------
     integer, intent(IN )            :: index 
     real   , intent(OUT)            :: caco3disolution
     real   , intent(OUT)            :: caco3precipitation 
!----------------------------------------------------------------------------- 
     real    :: Pelagic_module_DT
     real    :: Pelagic_module_DTDay 
     integer :: Phyto
!-----------------------------------------------------------------------------      
   Pelagic_module_DTDay = Me%ExternalParam%pelagic_module_dt_day
      Pelagic_module_DT = Me%ExternalParam%pelagic_module_dt 
                  Phyto = Me%PropIndex%Phyto                    
!-----------------------------------------------------------------------------  


 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
 end subroutine Compute_biogeoch_rates_caco3
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                                                         
                                                          
                                                          
                                                          


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
 !subroutine mocsy(index)
! This subroutine is adapted from mocsy package https://github.com/jamesorr/mocsy
! All the info about variables, compute options and others can be found on the link above
! integer, intent(IN) :: index
 
!Arguments-------------------------------------------------------------

     !integer, intent(IN) :: index 
!----------------------------------------------------------------------------
     
     !call (index)
     
!----------------------------------------------------------------------------

 !end subroutine mocsy
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 
 
 
 
 

 
 
 
 
 
 
 
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     
  subroutine KillCarbonateSystem(ObjCarbonateSystemID, STAT)

    !Arguments------------------------------------------------------------
    integer                             :: ObjCarbonateSystemID              
    integer, optional, intent(OUT)      :: STAT
    !External-------------------------------------------------------------
    integer                             :: ready_, nUsers           
    !Local----------------------------------------------------------------
    integer                             :: STAT_           
    !---------------------------------------------------------------------
    STAT_ = UNKNOWN_        
    call Ready(ObjCarbonateSystemID, ready_)  
    
if1 :   if (ready_ .NE. OFF_ERR_) then
            nUsers = DeassociateInstance(mCarbonateSystem_,  Me%InstanceID)  
            if (nUsers == 0) then
                call DeallocateInstance
                ObjCarbonateSystemID = 0
                STAT_ = SUCCESS_
            endif
        else 
            STAT_ = ready_
        end if if1   
        
     if (present(STAT)) STAT = STAT_
    !------------------------------------------------------------------------
    end subroutine KillCarbonateSystem
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
    
    
     
    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     
    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_CarbonateSystem), pointer          :: AuxObjCarbonateSystem
        type (T_CarbonateSystem), pointer          :: PreviousObjCarbonateSystem

        !Updates pointers
        if (Me%InstanceID == FirstObjCarbonateSystem%InstanceID) then
            FirstObjCarbonateSystem => FirstObjCarbonateSystem%Next
        else
            PreviousObjCarbonateSystem => FirstObjCarbonateSystem
            AuxObjCarbonateSystem      => FirstObjCarbonateSystem%Next
            do while (AuxObjCarbonateSystem%InstanceID /= Me%InstanceID)
                PreviousObjCarbonateSystem => AuxObjCarbonateSystem
                AuxObjCarbonateSystem      => AuxObjCarbonateSystem%Next
            enddo

            !Now update linked list
            PreviousObjCarbonateSystem%Next => AuxObjCarbonateSystem%Next

            !Deallocates instance
            deallocate (Me)
            nullify    (Me) 

        endif
     !----------------------------------------------------------------------       
    end subroutine DeallocateInstance
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

    subroutine Ready (ObjCarbonateSystemID, ready_) 
    
    !Arguments-------------------------------------------------------------
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
 !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
        
    
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