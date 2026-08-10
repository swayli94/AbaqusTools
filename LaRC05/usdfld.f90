include "plyTools.f90"
include "failureCriteria.f90"
include "materialResponse.f90"
    
SUBROUTINE USDFLD( FIELD, STATEV, PNEWDT, DIRECT, T, CELENT,        &
    &              TIME, DTIME, CMNAME, ORNAME, NFIELD, NSTATV,     &
    &              NOEL, NPT, LAYER, KSPT, KSTEP, KINC, NDI, NSHR,  &
    &              COORD, JMAC, JMATYP, MATLAYO, LACCFLA            )
!   ------|     --------------------------------------------------- !
!    FIELD      field variables at the current material point   [u] !
!   STATEV      solution-dependent state variables                  !
!   PNEWDT      ratio suggested time increment/current increment    !
!   DIRECT      cos material dirs in global direcions (dirs)        !
!        T      cos material dirs relative to element dirs          !
!   CELENT      characteristic element length
!  TIME(1)      step time at the beginning of the current increment !
!  TIME(2)      total time                                          !
!    DTIME      time increment                                      !
!   CMNAME      user-defined material name, left justified          !
!   ORNAME      user-defined local orientation name, left just.     !
!   NFIELD      number of field variables at each point             !
!   NSTATV      number of solution-dependent state variables        !
!     NOEL      element number                                      !
!      NPT      integration point number                            !
!    LAYER      layer number (for composite and layered solids)     !
!     KSPT      section point number within the current layer       !
!    KSTEP      step number                                         !
!     KINC      increment number                                    !
!      NDI      number of direct stress components at this point    !
!     NSHR      number of shear sress components at this point      !
!    COORD      cordinates at this material point                   !
!   ------      --------------------------------------------------- !
!     JMAC      variable to be passed into GETVRM                   !
!   JMATYP      variable to be passed into GETVRM                   !
!  MATLAYO      variable to be passed into GETVRM                   !
!  LACCFLA      variable to be passed into GETVRM                   !
!  -------      --------------------------------------------------- !
!                                     Miguel Matos, Silvestre Pinho !
!                                           miguelmatos@outlook.com !
!  -------      --------------------------------------------------- !
!                               created on 12.06.2019, Miguel Matos !
!                                tested on 03.07.2019, Miguel Matos !
!  -------      --------------------------------------------------- !

    use plyTools
    use failureCriteria
    use materialResponse
    
    include 'aba_param.inc'
    
    ! input variables
    integer, parameter      ::  RWP = KIND(1.0D0)
    real(RWP), parameter    ::  ZERO = 0.000D0, ONE = 1.000D0
    
    real(RWP), parameter    ::  FIELD_MAX_VAL = 0.999D0
    
    ! variables for table collections
    integer, parameter      ::  maxProps = 20,                      &
        &                       maxIndependentVars = 1,             &
        &                       minNProps   = 18  

    ! input vars
    integer, intent(in)     ::  NFIELD, NOEL, NPT, LAYER, KSPT,     &
        &                       KSTEP, KINC, NDI, NSHR, NSTATV
    
    character(len=80), intent(in)   ::  CMNAME, ORNAME
    real(RWP), intent(inout)        ::  FIELD(NFIELD),              &
        &                               STATEV(NSTATV), PNEWDT
    real(RWP), intent(in)           ::  DIRECT(3,3), T(3,3),        &
        &                               CELENT, TIME(2), DTIME,     &
        &                               COORD(3)

    ! access variables
    real(RWP)               :: ARRAY(15)
    integer                 :: JARRAY(15), JRCD
    character*3             :: FLGRAY(15)
    
    ! variables for output
    integer,parameter       :: OUTSIZE=10
    integer                 :: INTV(OUTSIZE)
    real(RWP)               :: REALV(OUTSIZE)
    character*8             :: CHARV(OUTSIZE)    
    character(len=128)      :: OUT_TEXT
    
    ! variables for table collections
    integer                 ::  jError, numProps, numTCs
    real(RWP)               ::  rProps(maxProps),                       &
        &                       dPropDVar(maxProps*maxIndependentVars), &
        &                       rIndepVars(maxIndependentVars)
    character(len=90)       ::  collectionName
    
    ! sdv variable position
    integer, parameter      ::  SDV_POS_FI1     = 1,                &
        &                       SDV_POS_FI2     = 2,                &
        &                       SDV_POS_FI3     = 3,                &
        &                       SDV_POS_FI4     = 4,                &
        &                       SDV_POS_FIMAX   = 5,                &
        &                       SDV_POS_dMatrix = 6,                &
        &                       SDV_POS_dSplit  = 7,                &
        &                       SDV_POS_dFibre  = 8,                &
        &                       SDV_POS_dKink   = 9
    
    
    ! locals    ****    ****    ****    ****    ****    ****    ****
    ! ply 
    type(PlyDefinition)     ::  currentPly
    ! reals
    real(RWP)               ::  ply2DStrain(3), elastic2DStress(3), &
        &                       plyIndexes(4), plyDamageVariables(4),&
        &                       plyHealth(3)   
    ! bools
    logical, parameter      ::  ENABLE_FI_DEN_MOD = .True.  ! do not change
    logical, parameter      ::  OUTDEBUG = .False.          ! debug to msg file
        
    if (jError .ne. 0) then
        write(OUT_TEXT,'(A)')   "(USDFLD) Starting"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    
    ! --------------------------------------------------------------!
    ! input checks
    if ((NFIELD .NE. 3) .OR. (NSTATV .LT. 9) ) then
        write(OUT_TEXT,'(A)')   "(USDFLD) 3 field variables and at least 9 SDVs are required"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! collection name
    collectionName  = "MATTC_"//trim(CMNAME)
    
    ! enable table collections
    call setTableCollection(collectionName, jError)
    if (jError .ne. 0) then
        write(OUT_TEXT,'(A)')   "(USDFLD) ERROR initializing table collections"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
        call XIT    ! failed initialization
    end if
       
    ! access tables
    call getPropertyTable("LARC05PROPERTIES", rIndepVars, ZERO, field, &
        &                 numProps, rProps, dPropDVar, 0, jError)
    if (jError .ne. 0) then
        write(OUT_TEXT,'(A)')   "(USDFLD) LARC05PROPERTIES table not found"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! correct number of properties
    if (numProps .LT. minNProps) then
        INTV(1) = numProps
        INTV(2) = minNProps
        write(OUT_TEXT,'(A)')   "(USDFLD): %I properties required; %I found"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! initialize ply
    call initializePly(currentPly, rProps, numProps)
    
    ! get history
    plyIndexes(:) = ZERO
    plyIndexes(1) = STATEV(SDV_POS_FI1)
    plyIndexes(2) = STATEV(SDV_POS_FI2)
    plyIndexes(3) = STATEV(SDV_POS_FI3)
    plyIndexes(4) = STATEV(SDV_POS_FI4)
    
    ! get ply old damage variables
    plyDamageVariables(:) = ZERO
    plyDamageVariables(1) = STATEV(SDV_POS_dMatrix)
    plyDamageVariables(2) = STATEV(SDV_POS_dSplit)
    plyDamageVariables(3) = STATEV(SDV_POS_dFibre)
    plyDamageVariables(4) = STATEV(SDV_POS_dKink)
    
    ! get ply strains
    CALL GETVRM('E', ARRAY, JARRAY, FLGRAY,             &
        &       JRCD, JMAC, JMATYP, MATLAYO, LACCFLA    )
    
    ! copy stresses to local variable
    ply2DStrain(:) = ZERO
    ply2DStrain(1) = dble(ARRAY(1))
    ply2DStrain(2) = dble(ARRAY(2))
    ply2DStrain(3) = dble(ARRAY(4))
    
    ! debug output
    if (OUTDEBUG) then
        REALV(1) = ply2DStrain(1)
        REALV(2) = ply2DStrain(2)
        REALV(3) = ply2DStrain(3)
        write(OUT_TEXT,'(A)')   "(USDFLD) E : %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! elastic 2D stress
    elastic2DStress = computeStresses2D(ply2DStrain, currentPly)
    ! debug output
    if (OUTDEBUG) then
        REALV(1) = elastic2DStress(1)
        REALV(2) = elastic2DStress(2)
        REALV(3) = elastic2DStress(3)
        write(OUT_TEXT,'(A)')   "(USDFLD) SE: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! compute elastic indexes
    plyIndexes   = plyCriteria2D( currentPly, elastic2DStress,      &
        &                         plyIndexes,                       &
                                  limitFIDen = ENABLE_FI_DEN_MOD    )
    ! debug
    if (OUTDEBUG) then
        REALV(1) = plyIndexes(1)
        REALV(2) = plyIndexes(2)
        REALV(3) = plyIndexes(3)
        REALV(4) = plyIndexes(4)
        write(OUT_TEXT,'(A)')   "(USDFLD) FI: %R,%R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    ! compute damage variables
    plyHealth = computeDamageVariables( ply2DStrain,                    &
        &                               elastic2DStress,                &
        &                               plyIndexes, plyDamageVariables, &
        &                               currentPly%Gvalues,             &
        &                               characteristicLength=CELENT     )
    
    ! update field
    FIELD(1) = MIN(ONE - plyHealth(1),FIELD_MAX_VAL)
    FIELD(2) = MIN(ONE - plyHealth(2),FIELD_MAX_VAL)
    FIELD(3) = MIN(ONE - plyHealth(3),FIELD_MAX_VAL)
    
    ! debug
    if (OUTDEBUG) then
        REALV(1) = FIELD(1)
        REALV(2) = FIELD(2)
        REALV(3) = FIELD(3)
        write(OUT_TEXT,'(A)')   "(USDFLD) fv: %R,%R,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    
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
    ! dealloc
    
    ! done
    
END SUBROUTINE USDFLD
