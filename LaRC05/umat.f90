 !DIR$ FREEFORM
    
include "plyTools.f90"
include "failureCriteria.f90"
include "materialResponse.f90"
    
!-------------------------------------------------------------------!
!          UMAT subroutine, defines the material behaviour          !
!-------------------------------------------------------------------!
SUBROUTINE UMAT( STRESS, STATEV, DDSDDE, SSE, SPD, SCD,             &
    &            RPL, DDSDDT, DRPLDE, DRPLDT,                       &
    &            STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,           &
    &            PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,   &
    &            PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT,       &
    &            DFGRD0, DFGRD1,                                    &
    &            NOEL, NPT, LAYER, KSPT, JSTEP, KINC                )
!   ------|     --------------------------------------------------- !
!   STRESS      stress tensor at the beginning TBU (to be updated)  !
!   STATEV      solution-dependent state variables TBU (as before)  ! 
!   DDSDDE      Jacobian: DDSDDE(I,J) = DS(I)/DE(J)                 !
!      SSE      specific elastic strain energy TBU                  !
!      SPD      specific plastic dissipation energy TBU             !
!      SCD      specific elastic “creep” dissipation energy TBU     !
!      RPL      volumetric heat generation by mechanical work    [T]!
!   DDSDDT      d(stress increments)/d(temperature)              [T]!
!   DRPLDE      d(RPL)/d(strain increments)                      [T]! 
!   DRPLDT      d(RPL)/d(temperature)                            [T]!
!    STRAN      total strains at the beginning of the increment     !
!   DSTRAN      array of strain increments                          !
!  TIME(1)      step time at the beginning of the current increment !
!  TIME(2)      total time                                          !
!    DTIME      time increment                                      ! 
!     TEMP      temperature at the start of the increment           !      
!    DTEMP      increment of temperature                            !
!   PREDEF      interpolated values of predefined field             !             
!    DPRED      array of increments of predefined field variables   !
!   CMNAME      user-defined material name, left justified          !
!      NDI      number of direct stress components at this point    !
!     NSHR      number of eng shear stress components at this point !
!    NTENS      size of the stress/strain array (NDI + NSHR)        !
!   NSTATV      number of sdv associated with this material         !    
!    PROPS      user-specified array of material                    !
!   NPROPS      number of PROPS                                     !
!   COORDS      coordinates of this point                           !
!     DROT      rotation increment matrix                           !
!   PNEWDT      ratio of suggested new time increment to DTIME      !
!   CELENT      characteristic element length                       !
!   DFDRD0      deformation grad at the beginning of the increment  !
!   DFGRD1      deformation gradient at the end of the increment    !
!     NOEL      element number                                      !
!      NPT      integration point number                            !
!    LAYER      layer number (for composite and layered solids)     !
!     KSPT      section point number within the current layer       !
!    JSTEP      step number                                         !
!     KINC      increment number                                    !
!  -------      --------------------------------------------------- !
!               Variables to be defined:                            !
!               - DDSDDE(NTENS,NTENS)                               !
!               - STRESS(NTENS)                                     !
!               - STATEV(NSTATV)                                    !
!               - SSE                                               !
!               - SPD                                               !
!               - SCD                                               !
!  -------      --------------------------------------------------- !
!                                     Miguel Matos, Silvestre Pinho !
!                                           miguelmatos@outlook.com !
!  -------      --------------------------------------------------- !
!                               created on 03.06.2019, Miguel Matos !
!                                tested on 03.07.2019, Miguel Matos !
!  -------      --------------------------------------------------- !

    use plyTools
    use failureCriteria
    use materialResponse

    !implicit none    

    include 'aba_param.inc'
    
    ! parameters
    integer, parameter      :: RWP = KIND(1.0D0)
    real(RWP), parameter    :: ZERO = 0.00D0, ONE = 1.00D0
    
    ! variables for table collections
    integer, parameter              ::  maxProps = 20,              &
        &                               maxIndependentVars = 1,     &
        &                               minNProps   = 18  
    
    ! input variables   -------------------------------------------------
    character(len=80), intent(in)   ::  CMNAME
    
    integer, intent(in)             ::  NDI, NSHR, NTENS, NSTATV, &
        &                               NPROPS, NOEL, NPT, LAYER, &
        &                               KSPT, KINC, JSTEP(4)
    
    real(RWP), intent(inout)        ::  STRESS(NTENS), STATEV(NTENS),   &
        &                               DDSDDE(NTENS,NTENS),            &
        &                               SSE, SPD, SCD,                  &
        &                               DDSDDT(NTENS), DRPLDE(NTENS),   &
        &                               RPL, DRPLDT,                    &
        &                               PNEWDT
    
    real(RWP), intent(in)           ::  STRAN(NTENS), DSTRAN(NTENS),    &
        &                               TIME(2), DTIME, TEMP, DTEMP,    &
        &                               PREDEF, DPRED, PROPS(NPROPS),   &
        &                               COORDS(3), DROT(3,3), CELENT,   &
        &                               DFGRD0(3,3), DFGRD1(3,3)
    
    ! -------------------------------------------------------------------
    integer, parameter      ::  SDV_POS_FI1 = 1,                        &
        &                       SDV_POS_FI2 = 2,                        &
        &                       SDV_POS_FI3 = 3,                        &
        &                       SDV_POS_FI4 = 4,                        &
        &                       SDV_POS_FIMAX = 5,                      &
        &                       SDV_POS_dMatrix = 6,                    &
        &                       SDV_POS_dSplit  = 7,                    &
        &                       SDV_POS_dFibre  = 8,                    &
        &                       SDV_POS_dKink   = 9
    
    
    ! variables for output
    integer, parameter      ::  OUTSIZE = 10
    integer                 ::  INTV(OUTSIZE)
    real(RWP)               ::  REALV(OUTSIZE)
    character*8             ::  CHARV(OUTSIZE)    
    character(len=128)      ::  OUT_TEXT
    
    ! local variables
    logical, parameter      ::  OUTDEBUG = .False.  ! Debug to .msg file 
    logical, parameter      ::  ENABLE_FI_DEN_MOD = .True.    
    
    ! local variables
    real(RWP)               ::  current2DStrain(3), elastic2DStress(3)
    real(RWP)               ::  plyIndexes(4), plyDamageVariables(4),   &   
        &                       plyDamage(3),                           &
        &                       dEnergy(2)
    ! ply 
    type(PlyDefinition)     :: currentPly
    
    ! variables for table collections
    integer                 ::  jError, numProps, numTCs
    real(RWP)               ::  rProps(maxProps),                       &
        &                       dPropDVar(maxProps*maxIndependentVars), &
        &                       rIndepVars(maxIndependentVars)
    character(len=90)       ::  collectionName
    
    
    !--------------------------------------------------------------------
    ! debug only
    if (OUTDEBUG) then
        INTV(1) = JSTEP(1)
        INTV(2) = KINC
        INTV(3) = NOEL
        INTV(4) = NPT
        INTV(5) = LAYER
        INTV(6) = KSPT
        write(OUT_TEXT,'(A)')   "* %I.%I UVARM EL %I PT %I LAY %I SPT %I"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
        REALV(1) = CELENT
        write(OUT_TEXT,'(A)')   "* Lc = %R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! error checks
    if (NTENS .NE. 3) then
        write(OUT_TEXT,'(A)')   "* This UMAT is intended for 2D elements"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! initializations
    collectionName  = "MATTC_"//trim(CMNAME)
    
    ! enable table collections
    call setTableCollection(collectionName, jError)
    if (jError .ne. 0) then
        write(OUT_TEXT,'(A)')   "* ERROR initializing table collections"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
        call XIT    ! failed initialization
    end if
       
    ! access tables
    call getPropertyTable("LARC05PROPERTIES", rIndepVars, ZERO, field, &
        &                 numProps, rProps, dPropDVar, 0, jError)
    
    if (jError .ne. 0) call XIT    ! failed initialization
    
    if (numProps .LT. minNProps) then
        INTV(1) = numProps
        INTV(2) = minNProps
        write(OUT_TEXT,'(A)')   "* ERROR: %I properties required; %I found"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
        call XIT    ! failed initialization
    end if
    
    ! get old ply Indexes
    plyIndexes(1) = STATEV(SDV_POS_FI1)
    plyIndexes(2) = STATEV(SDV_POS_FI2)
    plyIndexes(3) = STATEV(SDV_POS_FI3)
    plyIndexes(4) = STATEV(SDV_POS_FI4)
    
    ! get ply old damage variables
    plyDamageVariables(1) = STATEV(SDV_POS_dMatrix)
    plyDamageVariables(2) = STATEV(SDV_POS_dSplit)
    plyDamageVariables(3) = STATEV(SDV_POS_dFibre)
    plyDamageVariables(4) = STATEV(SDV_POS_dKink)
        
    !--------------------------------------------------------------------    
    
    ! initialize ply 
    call initializePly(currentPly, rProps, numProps)
    
    ! get actual strain
    current2DStrain = STRAN + DSTRAN
    
    if (OUTDEBUG) then
        REALV(1) = STRESS(1)
        REALV(2) = STRESS(2)
        REALV(3) = STRESS(3)
        
        write(OUT_TEXT,'(A)')   "* Si: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    if (OUTDEBUG) then
        REALV(1) = current2DStrain(1)
        REALV(2) = current2DStrain(2)
        REALV(3) = current2DStrain(3)
        
        write(OUT_TEXT,'(A)')   "* EE: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! get "elastic" stresses and failure criteria
    elastic2DStress = computeStresses2D(current2DStrain, currentPly)
    
    if (OUTDEBUG) then
        REALV(1) = elastic2DStress(1)
        REALV(2) = elastic2DStress(2)
        REALV(3) = elastic2DStress(3)
        
        write(OUT_TEXT,'(A)')   "* SE: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! compute elastic indexes
    plyIndexes   = plyCriteria2D( currentPly, elastic2DStress,      &
        &                         plyIndexes,                       &
                                  limitFIDen = ENABLE_FI_DEN_MOD    )
    
    if (OUTDEBUG) then
        REALV(1) = plyIndexes(1)
        REALV(2) = plyIndexes(2)
        REALV(3) = plyIndexes(3)
        REALV(4) = plyIndexes(4)
        
        write(OUT_TEXT,'(A)')   "* FI: %R,%R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! compute damage variables
    plyDamage = computeDamageVariables( current2DStrain,                &
        &                               elastic2DStress,                &
        &                               plyIndexes, plyDamageVariables, &
        &                               currentPly%Gvalues,             &
        &                               characteristicLength=CELENT     )
    ! store old stresses 
    elastic2DStress = STRESS
    
    ! compute real stress
    STRESS = computeStresses2D( current2DStrain, currentPly,            &
        &                       materialDamage = plyDamage,             &
        &                       constitutiveMatrix = DDSDDE )
    
    if (OUTDEBUG) then
        REALV(1) = STRESS(1)
        REALV(2) = STRESS(2)
        REALV(3) = STRESS(3)
        write(OUT_TEXT,'(A)')   "* So: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    if (OUTDEBUG) then
        REALV(1) = DDSDDE(1,1)
        REALV(2) = DDSDDE(1,2)
        REALV(3) = DDSDDE(1,3)
        write(OUT_TEXT,'(A)')   "* DDSDDE1: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
        REALV(1) = DDSDDE(2,1)
        REALV(2) = DDSDDE(2,2)
        REALV(3) = DDSDDE(2,3)
        write(OUT_TEXT,'(A)')   "* DDSDDE2: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
        REALV(1) = DDSDDE(3,1)
        REALV(2) = DDSDDE(3,2)
        REALV(3) = DDSDDE(3,3)
        write(OUT_TEXT,'(A)')   "* DDSDDE3: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    !-------------------------------------------------------------------- 
    ! compute energy (elastic2DStress has OLDSTRESS)
    dEnergy = computeEnergyVariation( STRAN, DSTRAN,                    &
        &                             elastic2DStress, STRESS,          &
        &                             NTENS                             )
    
    ! elastic energy
    SSE = SSE + dEnergy(1)
    ! plastic dissipation [total - elastic]
    SPD = SPD + dEnergy(2)
    
    !-------------------------------------------------------------------- 
    ! update SDVs
    STATEV(SDV_POS_dMatrix) = plyDamageVariables(1) 
    STATEV(SDV_POS_dSplit)  = plyDamageVariables(2) 
    STATEV(SDV_POS_dFibre)  = plyDamageVariables(3)
    STATEV(SDV_POS_dKink)   = plyDamageVariables(4) 
    
    ! get old ply Indexes
    STATEV(SDV_POS_FI1) = plyIndexes(1)
    STATEV(SDV_POS_FI2) = plyIndexes(2)
    STATEV(SDV_POS_FI3) = plyIndexes(3)
    STATEV(SDV_POS_FI4) = plyIndexes(4)
    
    if (SDV_POS_FIMAX .GT. 0) then
        STATEV(SDV_POS_FIMAX) = MAX( plyIndexes(1), plyIndexes(2), plyIndexes(3), plyIndexes(4) )
    end if
     
    ! complete    

END SUBROUTINE UMAT