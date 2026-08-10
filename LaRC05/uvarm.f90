include "plyTools.f90"
include "failureCriteria.f90"
    
SUBROUTINE UVARM( UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,NUVARM,    &
    &             NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,    &
    &             JMAC,JMATYP,MATLAYO,LACCFLA)
!   ------|     --------------------------------------------------- !
!     UVAR      field variables at the current material point   [u] !
!   DIRECT      cos material dirs in global direcions (dirs)        !
!        T      cos material dirs relative to element dirs          !
!  TIME(1)      step time at the beginning of the current increment !
!  TIME(2)      total time                                          !
!    DTIME      time increment                                      !
!   CMNAME      user-defined material name, left justified          !
!   ORNAME      user-defined local orientation name, left just.     !
!   NUVARM      number of user-defined output variables             !
!     NOEL      element number                                      !
!      NPT      integration point number                            !
!    LAYER      layer number (for composite and layered solids)     !
!     KSPT      section point number within the current layer       !
!    KSTEP      step number                                         !
!     KINC      increment number                                    !
!      NDI      number of direct stress components at this point    !
!     NSHR      number of shear sress components at this point      !
!    COORD      cordinates at this material point
!   ------      --------------------------------------------------- !
!     JMAC      variable to be passed into GETVRM                   !
!    JMTYP      variable to be passed into GETVRM                   !
!  MATLAYO      variable to be passed into GETVRM                   !
!  LACCFLA      variable to be passed into GETVRM                   !
!  -------      --------------------------------------------------- !
!   log
!       20.05.2019: Added history dependency
!       21.05.2019: Tested in Abq (run and value checks)
!       25.05.2019: Added enhanced UVARM
!       10.06.2019: Added support for 3D 
!       13.06.2019: Modified UVARM7
!       27.06.2019: Compatibility between output of 2D and 3D FIs
!          03.2020: Updated ply definition
!  -------      --------------------------------------------------- !
!                                     Miguel Matos, Silvestre Pinho !
!                                           miguelmatos@outlook.com !
!  -------      --------------------------------------------------- !
!  -------      --------------------------------------------------- !

    use plyTools
    use failureCriteria
    
    include 'aba_param.inc'

    ! parameters
    integer, parameter      :: RWP = KIND(1.0000D0)
    real(RWP), parameter    :: ZERO = 0.00D0, ONE = 1.00D0
    ! variables for table collections
    integer, parameter              :: maxProps = 18,             &
        &                              maxIndependentVars = 1

    ! input variables   -------------------------------------------------
    character(len=80)       :: CMNAME, ORNAME
    
    real(RWP)               :: UVAR(NUVARM), DIRECT(3,3), T(3,3),       &
                            &  TIME(2), COORD(3)

    real(RWP)               :: DTIME

    integer                 :: NUVARM, NOEL, NPT, LAYER, KSPT, KSTEP,   &
                            &  KINC, NDI, NSHR
    ! -------------------------------------------------------------------
        
    ! local variables
    logical, parameter      :: OUTDEBUG = .False. 
    real(RWP), parameter    :: LIMIT_UVARM5 = ONE

    ! access variables
    real(RWP)               :: ARRAY(15)
    integer                 :: JARRAY(15), JRCD
    character(len=3)        :: FLGRAY(15)
    
    ! variables for output
    integer,parameter       :: OUTSIZE=10
    integer                 :: INTV(OUTSIZE)
    real(RWP)               :: REALV(OUTSIZE)
    character(len=8)        :: CHARV(OUTSIZE)    
    character(len=128)      :: OUT_TEXT

    ! variables for table collections
    integer                 :: jError, numProps, numTCs
    real(RWP)               :: rProps(maxProps), dPropDVar(maxProps*maxIndependentVars), rIndepVars(maxIndependentVars)
    character(len=90)       :: collectionName
    
    ! ply properties, stresses and indexes
    type(PlyDefinition)     :: currentPly
    real(RWP)               :: plyStresses(6), plyIndexes(5)  
    logical                 :: isElement3D, computeMIIndex
    
    ! locals
    integer                 :: i, maxFailureMode = 1, numIndexes
    real(RWP)               :: maxFailureIndex = ZERO
    
    ! definition checks
    if (NUVARM .LT. 4) then
        write(OUT_TEXT,'(A)')   "** This UVARM requires at least 4 variables"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
        ! this should terminate the analysis
    end if
    
    ! element dim
    isElement3D = (NDI+NSHR == 6)
    computeMIIndex = .False.
    
    ! debug only
    if (OUTDEBUG) then
        INTV(1) = KSTEP
        INTV(2) = KINC
        INTV(3) = NOEL
        INTV(4) = NPT
        INTV(5) = LAYER
        INTV(6) = KSPT
        write(OUT_TEXT,'(A)')   "(UVARM) %I.%I UVARM EL %I PT %I LAY %I SPT %I"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
        if (isElement3D) then
            write(OUT_TEXT,'(A)')   "(UVARM) 3D Element identified"
        else
            write(OUT_TEXT,'(A)')   "(UVARM) 2D Element identified"
        end if
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
        ! Uvarm V
        write(OUT_TEXT,'(A)')   "(UVARM) UVARM V2"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
        
    end if

    ! initializations
    collectionName  = "MATTC_"//trim(CMNAME)
    
    ! enable table collections
    call setTableCollection(collectionName, jError)
    if (jError .ne. 0) then
        write(OUT_TEXT,'(A)')   "(UVARM) * ERROR initializing table collections"
        call STDB_ABQERR(-3,OUT_TEXT,INTV,REALV,"none")
        call XIT    ! failed initialization
    end if
    
    ! access tables
    call getPropertyTable("LARC05PROPERTIES", rIndepVars, ZERO, field, &
        &                 numProps, rProps, dPropDVar, 0, jError)
    
    if (jError .ne. 0) call XIT    ! failed initialization
      
    ! initialize ply properties from read property list
    call initializePly(currentPly, rProps, numProps)
    
    
    !   :DELETE:
    if (OUTDEBUG) then
        if (NOEL*NPT*LAYER*KSPT .EQ. 657) then
            call plyPrint(currentPly)
        end if
    end if
    
    
    ! check properties for NCF
    isElement3D = (isElement3D .AND. (currentPly%Zt .GT. ZERO) .AND. NUVARM .GE. 6)
    ! debug only
    if (isElement3D) then
        numIndexes = 5
        computeMIIndex = .True.
        if (OUTDEBUG) then
            write(OUT_TEXT,'(A)')   "(UVARM): Computing 5 failure indexes in 3D"
            call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
            REALV(1) = currentPly%Zt    ! :DELETE:
            REALV(2) = currentPly%ILSS
            write(OUT_TEXT,'(A)')   "(UVARM): Zt %R, ILSS %R"
            call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
        end if
    else
        numIndexes = 4
    end if
    
    ! get ply stresses
    CALL GETVRM('S', ARRAY, JARRAY, FLGRAY,             &
        &       JRCD, JMAC, JMATYP, MATLAYO, LACCFLA    )
    
    ! copy stresses to local variable
    plyStresses(:) = ZERO
    do i=1,6
        plyStresses(i) = dble(ARRAY(i))
    end do
    
    ! evaluate ply indexes
    plyIndexes(:) = ZERO
    do i=1,4
        plyIndexes(i) = dble(UVAR(i))
    end do    
    if (computeMIIndex) then
        if (NUVARM .GT. 5) then
            plyIndexes(5) = dble(UVAR(5))
        end if
    end if
    
    ! new version
    plyIndexes = completeCriteria( currentPly, plyStresses, plyIndexes, &
        &                          6, 5,                                &
        &                          isElement3D, matrixInterfaceIndex=computeMIIndex)
    
    
    ! debug output
    if (OUTDEBUG) then
        REALV(1) = plyIndexes(1)
        REALV(2) = plyIndexes(2)
        REALV(3) = plyIndexes(3)
        REALV(4) = plyIndexes(4)
        REALV(5) = plyIndexes(5)
        write(OUT_TEXT,'(A)')   "* fI: %R, %R, %R, %R ,%R"
        call STDB_ABQERR(-1,OUT_TEXT,INTV,REALV,"none")
    end if
    
    
    UVAR(:) = ZERO
    
    ! copy back to UVAR
    do i=1,4
        UVAR(i) = plyIndexes(i)
    end do
    
    ! UVARM 5 -> matrix interface failure mode
    if (NUVARM .GT. 4) then
       if (computeMIIndex) then
           UVAR(5) = plyIndexes(5)
        end if
    end if
        
    ! UVARM 6 -> MAX INDEX
    if (NUVARM .GT. 5) then
        maxFailureIndex = maxval(plyIndexes)
        UVAR(6) = maxFailureIndex
        if (LIMIT_UVARM5 .GT. ZERO) then
            UVAR(6) = min(LIMIT_UVARM5, UVAR(6))
        end if
    end if
    
    ! UVARM 7
    if (NUVARM .GT. 6) then
    ! get max and mode        
        maxFailureIndex = ZERO
        maxFailureMode  = 0
        do i=1,numIndexes
            if (plyIndexes(i) .GT. maxFailureIndex) then
                maxFailureMode = i
                maxFailureIndex = plyIndexes(i)
            end if
        end do
        
        ! limit failure index
        maxFailureIndex = min(ONE,maxFailureIndex)
        
        UVAR(7) = 1000.000D0 +  10.0000D0 * nint(10.0d0 * maxFailureIndex) + &
                    &   float(maxFailureMode)  ! +1000 added on 13.06.2019
    end if
  
    
END SUBROUTINE UVARM