include "plyTools.f90"
include "failureCriteria.f90"
include "materialResponse.f90"

SUBROUTINE USDFLD( FIELD, STATEV, PNEWDT, DIRECT, T, CELENT,            &
    &              TIME, DTIME, CMNAME, ORNAME, NFIELD, NSTATV,         &
    &              NOEL, NPT, LAYER, KSPT, KSTEP, KINC, NDI, NSHR,      &
    &              COORD, JMAC, JMATYP, MATLAYO, LACCFLA                )
!   ------|     --------------------------------------------------- ----!
!    FIELD      field variables at the current material point   [u]     !
!   STATEV      solution-dependent state variables                      !
!   PNEWDT      ratio suggested time increment/current increment        !
!   DIRECT      cos material dirs in global direcions (dirs)            !
!        T      cos material dirs relative to element dirs              !
!   CELENT      characteristic element length                           !
!  TIME(1)      step time at the beginning of the current increment     !
!  TIME(2)      total time                                              !   
!    DTIME      time increment                                          !
!   CMNAME      user-defined material name, left justified              !
!   ORNAME      user-defined local orientation name, left just.         !
!   NFIELD      number of field variables at each point                 !
!   NSTATV      number of solution-dependent state variables            !
!     NOEL      element number                                          !
!      NPT      integration point number                                !
!    LAYER      layer number (for composite and layered solids)         !
!     KSPT      section point number within the current layer           !
!    KSTEP      step number                                             !
!     KINC      increment number                                        !
!      NDI      number of direct stress components at this point        !
!     NSHR      number of shear sress components at this point          !
!    COORD      cordinates at this material point                       !
!   ------      --------------------------------------------------------!
!     JMAC      variable to be passed into GETVRM                       !
!   JMATYP      variable to be passed into GETVRM                       !
!  MATLAYO      variable to be passed into GETVRM                       !
!  LACCFLA      variable to be passed into GETVRM                       !
!  -------      --------------------------------------------------------!
!                                     Miguel Matos, Silvestre Pinho     !
!                                           miguelmatos@outlook.com     !
!                                     Miguel Matos                      !
!                                           miguel.matos@airbus.com     !
!  -------      --------------------------------------------------------!
!  LaRC05       damage model for 2D/3D elements                         !
!  -------      --------------------------------------------------------!
!                               created on 19.11.2019, Miguel Matos     !
!  -------      --------------------------------------------------------!

    use plyTools
    use failureCriteria
    use materialResponse
    
    !implicit none
    include 'aba_param.inc'
    
    ! input variables
    integer, parameter      ::  RWP = KIND(1.0D0)
    real(RWP), parameter    ::  ZERO = 0.000D0, ONE = 1.000D0
    
    real(RWP), parameter    ::  MAX_FIELD_VAL = 0.9990D0, & !0.9990D0
        &                       MIN_FIELD_VAL = ZERO
    
    ! sdv variable position
    integer, parameter      ::  SDV_NUMBER_2D    = 8,               &
        &                       SDV_NUMBER_3D    = 10,              &
        &                       FV_NUMBER_2D     = 3,               &
        &                       FV_NUMBER_3D     = 4,               &
        &                       MIN_PROP_NUMBER_2D = 18,            &
        &                       MIN_PROP_NUMBER_3D = 25,            &
        &                       F_VAL_TEN = ZERO,                   &
        &                       F_VAL_COMP = -ONE

    ! bools
    logical, parameter      ::  OUTDEBUG = .False.      ! debug to msg file
    logical, parameter      ::  SIMPLEDEBUG = .False.   ! simple debug to the msg file
    logical, parameter      ::  ROUND_RESULT = .False.  ! round the final result
    logical, parameter      ::  READ_OLD_FIELD = .True. ! read old field values
    
    ! variables for table collections
    integer, parameter      ::  maxProps = 28,                      &
        &                       maxIndependentVars = 1

    ! input vars
    integer, intent(in)     ::  NFIELD, NOEL, NPT, LAYER, KSPT,     &
        &                       KSTEP,  KINC, NDI, NSHR,  NSTATV
    character(len=80), intent(in)       ::  CMNAME, ORNAME
    real(RWP), intent(inout)            ::  FIELD(NFIELD),          &   
        &                                   STATEV(NSTATV), PNEWDT
    real(RWP), intent(in)   ::  DIRECT(3,3), T(3,3),                &
        &                       CELENT, TIME(2), DTIME,             &
        &                       COORD(3)

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
    
    ! 3D ELEMENT FLAG
    logical                 ::  isElement2D = .False.
    
    ! locals    ****    ****    ****    ****    ****    ****    ****
    ! ply 
    type(PlyDefinition)  ::  currentPly
    ! cycles
    integer             ::  i, NSTR=-1
    integer             ::  nFieldTenComp = -1, nFIs = 4, nDMG = 3, &
        &                   nFIMax = -1, NFIELDDMG = 3
    ! reals
    real(RWP)           ::  strainVector(6),                            &
        &                   failureIndexes(5), damageVariables(5),      &
        &                   healthVariables(6),                         &
        &                   fieldValue,                                 &
        &                   oldFieldValues(NFIELD), fieldValues(NFIELD)

    ! initializations
    oldFieldValues(:) = ZERO
    NFIELDDMG = NFIELD

    ! --------------------------------------------------------------!
    ! debug print
    if (OUTDEBUG) then
        INTV(1) = KSTEP
        INTV(2) = KINC
        INTV(3) = NOEL
        INTV(4) = NPT
        INTV(5) = LAYER
        INTV(6) = KSPT
        write(OUT_TEXT,'(A)')   "(USDFLD) inc %I.%I elm %I SP %I.%I.%I"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
        REALV(1) = CELENT   
        write(OUT_TEXT,'(A)')   "(USDFLD) Le: %R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
    end if
    
    ! element type id
    NSTR = NDI + NSHR
    if (NSTR .EQ. 3) then
        isElement2D = .True.
    else if (NSTR .EQ. 6) then
        isElement2D = .False.
    else
        INTV(1) = NSTR
        write(OUT_TEXT,'(A)')   "(USDFLD) Expecting 3 or 6 stress components. Got %I"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,CHARV)
    end if
           
    ! number checks only for the first increment
    if ((KSTEP * KINC) .LT. 2) then
        if (OUTDEBUG) then
            INTV(1) = 3
            if (isElement2D) INTV(1) = 2
            write(OUT_TEXT,'(A)') "%I-D element identified"
        end if
        ! 2D
        if (isElement2D) then
            if (NFIELD .LT. FV_NUMBER_2D) then
            INTV(1) = NFIELD
            INTV(2) = FV_NUMBER_2D
            write(OUT_TEXT,'(A)')   &
            &   "(USDFLD) Incorrect number of field variables. Got %I, expecting %I"
            call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,CHARV)
            end if
            if (NSTATV .LT. SDV_NUMBER_2D) then
            INTV(1) = NSTATV
            INTV(2) = SDV_NUMBER_2D
            write(OUT_TEXT,'(A)')   &
            &   "(USDFLD) Incorrect number of SDVs. Got %I, expecting %I"
            call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,CHARV)
            end if
        else
            if (NFIELD .LT. FV_NUMBER_3D) then
            INTV(1) = NFIELD
            INTV(2) = FV_NUMBER_3D
            write(OUT_TEXT,'(A)')   &
            &   "(USDFLD) Incorrect number of field variables. Got %I, expecting %I"
            call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,CHARV)
            end if
            if (NSTATV .LT. SDV_NUMBER_3D) then
            INTV(1) = NSTATV
            INTV(2) = SDV_NUMBER_3D
            write(OUT_TEXT,'(A)')   &
            &   "(USDFLD) Incorrect number of SDVs. Got %I, expecting %I"
            call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,CHARV)
            end if
        end if
    end if
    ! end of input checks

    ! set numbering of entries
    if (isElement2D) then
        nFIs = 4
        nDMG = 3
        if (NSTATV .GT. SDV_NUMBER_2D) then
            nFIMax = 2*nFIs + 1
        end if
        if (NFIELD .GT. FV_NUMBER_2D) then
            NFIELDDMG = NFIELD-1 
            nFieldTenComp = FV_NUMBER_2D + 1
        end if
    else
        nFIs = 5
        nDMG = 4
        if (NSTATV .GT. SDV_NUMBER_3D) then
            nFIMax = 2*nFIs + 1
        end if
        if (NFIELD .GT. FV_NUMBER_3D) then
            NFIELDDMG = NFIELD-1 
            nFieldTenComp = FV_NUMBER_3D + 1
        end if
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
    if (numProps .LT. MIN_PROP_NUMBER_2D) then
        INTV(1) = numProps
        INTV(2) = MIN_PROP_NUMBER_2D
        write(OUT_TEXT,'(A)')   "(USDFLD): %I properties required; %I found"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
    end if    

    
    ! initialize ply
    call initializePly(currentPly, rProps, numProps)
    
    ! tension/compression
    if (nFieldTenComp .GT. 0) then
        ! get last stress
        ! this is needed in advance as it is used to compute elastic S 
        CALL GETVRM('S', ARRAY, JARRAY, FLGRAY,             &
        &       JRCD, JMAC, JMATYP, MATLAYO, LACCFLA)
        
        !NZSTRESS is defined in failureCriteria.f90
        if (ARRAY(1) .LT. -NZSTRESS) then
            fieldValue = F_VAL_COMP
            currentPly%TENCOMP = .False.
        else
           fieldValue = F_VAL_TEN
        end if
        ! debug
        if (OUTDEBUG) then
            REALV(1) = fieldValue
            write(OUT_TEXT,'(A)')   "(USDFLD) Tension/compression: %R"
            call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
        end if
        ! write 
        FIELD(nFieldTenComp) = fieldValue
    end if
    
    ! get ply strains
    CALL GETVRM('E', ARRAY, JARRAY, FLGRAY,             &
        &       JRCD, JMAC, JMATYP, MATLAYO, LACCFLA    )
    strainVector(:) = ZERO
    do i=1,6
        strainVector(i) = DBLE(ARRAY(i))
    end do
    
    ! debug
    if (OUTDEBUG) then  
        ! Fields
        OUT_TEXT = "(USDFLD) E : %R,%R,%R"
        REALV(1:3) = strainVector(1:3)
        write(OUT_TEXT,'(A,A)') trim(OUT_TEXT)
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
        REALV(1:3) = strainVector(4:6)
        write(OUT_TEXT,'(A,A)') trim(OUT_TEXT)
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
    end if
    
    ! get history   ->  simplified version as of 03.12.2019
    ! failure indexes
    failureIndexes(:) = ZERO
    ! mode damage variables
    damageVariables(:) = ZERO
    do i=1,nFIs
        failureIndexes(i)  = STATEV(i)
        damageVariables(i) = STATEV(nFIs+i)
    end do
      
    healthVariables(:) = ONE
      
    ! update all damage variables
    call updateDamageVariables( currentPly,                             &
        &                       strainVector, 6,                        &
        &                       failureIndexes, 5,                      &
        &                       damageVariables,                        &
        &                       healthVariables, 6,                     &
        &                       celent,                                 & 
        &                       isElement2D,                            &
        &                       .False.,                                &
        &                       -ONE                                    )
    ! 3D failure angles are disabled for damage propagation, as change in
    ! angle during damage propagation would not respect the softening law

    ! retrieve the old damage/field variables
    if (READ_OLD_FIELD) then
        CALL GETVRM('FV', ARRAY, JARRAY, FLGRAY,             &
            &       JRCD, JMAC, JMATYP, MATLAYO, LACCFLA    )
        do i=1,NFIELD
            oldFieldValues(i) = ARRAY(i)
        end do
    end if

    ! debug for initial field values
    if (OUTDEBUG) then  
        OUT_TEXT = ""
        do i=1,NFIELD
            OUT_TEXT = trim(OUT_TEXT)// "%R,"
            REALV(i) = oldFieldValues(i)
        end do
        write(OUT_TEXT,'(A,A)')   "(USDFLD) Fd0: ",trim(OUT_TEXT)
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
    end if
    
    ! update field   
    if (isElement2D) then
        fieldValues(1) = ONE- healthVariables(1)
        fieldValues(2) = ONE- healthVariables(2)
        fieldValues(3) = ONE- healthVariables(4)
    else
        fieldValues(1) = ONE- healthVariables(1)
        fieldValues(2) = ONE- healthVariables(2)
        fieldValues(3) = ONE- healthVariables(3)
        fieldValues(4) = ONE- MIN(healthVariables(4), healthVariables(5), healthVariables(6))
    end if
    
    ! limit field
    do i=1,NFIELDDMG
        fieldValues(i) = MIN(fieldValues(i),MAX_FIELD_VAL)
        fieldValues(i) = MAX(fieldValues(i),MIN_FIELD_VAL)
    end do

    ! viscosity - field values go from 0 to 1 ( = damage)
    if (READ_OLD_FIELD .and. (currentPly%viscosity .GT. ZERO)) then
        fieldValues = viscousRegularization(fieldValues,          &
            &                               oldFieldValues,       &
            &                               NFIELD,               &
            &                               DTIME,                &
            &                               currentPly%viscosity  )
    end if

    ! write back
    do i=1,NFIELDDMG
        ! ensure damage does not decrease
        fieldValue = MAX(oldFieldValues(i),fieldValues(i))
        if (ROUND_RESULT) then
            fieldValue = NINT(fieldValue*1.00D+3) * 1.000D-3
        end if
        FIELD(i) = fieldValue
    end do
        
    ! Update solution dependent variables
    ! 1. failure indexes and damage variables
    STATEV(:) = ZERO
    do i=1,nFIs
        STATEV(i) = failureIndexes(i)
        STATEV(i+nFIs) = damageVariables(i)
    end do
    if (nFIMax .GT. 0) then
        STATEV(nFIMax) = MAXVAL( failureIndexes )
    end if
    
    ! debug
    if (OUTDEBUG .OR. SIMPLEDEBUG) then  
        ! Fields
        OUT_TEXT = ""
        do i=1,NFIELD
            OUT_TEXT = trim(OUT_TEXT)// "%R,"
            REALV(i) = FIELD(i)
        end do
        write(OUT_TEXT,'(A,A)')   "(USDFLD) Fd: ",trim(OUT_TEXT)
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
    end if
    
    if (OUTDEBUG) then      
        ! Failure indexes
        OUT_TEXT = ""
        do i=1,nFIs
            OUT_TEXT = trim(OUT_TEXT)// "%R,"
            REALV(i) = failureIndexes(i)
        end do
        write(OUT_TEXT,'(A,A)')   "(USDFLD) FI: ",trim(OUT_TEXT)
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
        ! Damage variables
        OUT_TEXT = ""
        do i=1,5
            OUT_TEXT = trim(OUT_TEXT)// "%R,"
            REALV(i) = damageVariables(i)
        end do
        write(OUT_TEXT,'(A,A)')   "(USDFLD) dI: ",trim(OUT_TEXT)
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
        write(OUT_TEXT,'(A,A)')   "(USDFLD) -------------------------- "
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,CHARV)
    end if
   
    ! done
    
END SUBROUTINE USDFLD
