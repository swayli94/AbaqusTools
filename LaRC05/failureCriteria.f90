!DIR$ FREEFORM
!  failureCriteria.f90 
!   
!  FUNCTIONS:
!   completeCriteria    -   computes all (5) failure indexes for a ply
!                           deals with 2D and 3D S or elements
!   plyCriteria2D       -   simplified 2D failure indexes for a ply
!   rotateStress        -   rotates stress vectors
!   matrixTractions     -   computes the matrix tractions for a crack angle
!
!****************************************************************************
!                   developed by Miguel A.S. Matos and Silvestre T. Pinho
!                                       Imperial College London, Aeronautics
!                                                       created 17.04.2019
!****************************************************************************
!   log:
!       04.06.2019  Modified plyCriteria to export failure angles [reverted]
!       05.06.2019  Reverted
!       06.06.2019  Parameter changes for performance
!       10.06.2019  Support for 3D and NCF criteria    
!       02.07.2019  plyCriteria2D increased support for large stress values
!       03.07.2019  Final comment and DEBUG mode disable
!       -- --   --
!       11.11.2019  completeCriteria modified to support damage control
!       20.11.2019  removed deprecated plyCriteria
!       22.11.2019  added NZSTRESS threshold for positive/negative criteria
!                   modification limited to completeCriteria
!
!   last error free build on:
!       03.07.2019
!       22.11.2019
!****************************************************************************

module failureCriteria

    use plyTools
    
    implicit none

    ! global definitions
    integer, parameter, private     ::  RWP = KIND(1.0D0)
    integer, parameter, private     ::  n23 = 6, n13 = 5
    real(RWP), parameter, private   ::  ZERO = 0.000000D0, ONE  = 1.000000D0,   &
        &                               HALF = 0.500000D0, NZNEG = -1.000D-12,  &
        &                               TEN  = 10.00000D0              
    ! this one is not private
    real(RWP), parameter            ::  NZSTRESS = 1.000D-6
    real(RWP), parameter            ::  NEG_STRESS_THS = -1.000D-3
    ! math constants
    real(RWP), parameter,private    ::  PI = ACOS(-1.000000D0),                 &
        &                               DEGTORAD = PI/180.000000D0

    ! model parameters
    real(RWP), parameter, private   ::  KINK_ANGLE_SPACING = 30.000D0

    contains
  
    !_rotateStress______________________________________________________
    ! rotates the stress vector around the specified axis
    !___________________________________________________________________
    ! consider using rotsig
    function rotateStress(inputStress,angleRad,axis)
    !-------------------------------------------------------------------!
    !   inputStress   stress vector (6)
    !   angleRad      angle to rotate in radians
    !   axis          for axis (1,2,3) for (x,y,z), respectively
    !-------------------------------------------------------------------!
    !   edited on 26.11.2019 to include parametric nS23, nS31
    !-------------------------------------------------------------------!

        implicit none
    
        real(RWP),intent(in)    ::  inputStress(6), angleRad
        integer, intent(in)     ::  axis
        real(RWP)               ::  rotateStress(6) 

        real(RWP)               ::  c1, s1, c2, s2, cs

        !i = (/2, 3, 1, 5, 6, 4/) rotate around x
        !i = (/3, 1, 2, 6, 4, 5/) rotate around y
        !i = (/1, 2, 3, 4, 5, 6/) rotate around z
        
        ! initialization
        rotateStress = inputStress
        
        ! if the angle is zero just leave
        if (angleRad .eq. 0.00d0) then
            return 
        end if

        
        ! sin and cos
        c1 = cos(angleRad)
        c2 = c1*c1
        s1 = sin(angleRad)
        s2 = s1*s1
        cs = c1*s1


        ! rotate
        select case (axis)
        case (1)
            rotateStress(2) =  c2*inputStress(2) + s2*inputStress(3) + 2*cs*inputStress(n23)
            rotateStress(3) =  s2*inputStress(2) + c2*inputStress(3) - 2*cs*inputStress(n23)
            rotateStress(1) =  inputStress(1)
            rotateStress(n23) = -cs*inputStress(2) + cs*inputStress(3) + (c2-s2)*inputStress(n23)
            rotateStress(n13) =  cs*inputStress(n13) - s1*inputStress(4)
            rotateStress(4) =  s1*inputStress(n13) + c1*inputStress(4)
        case (2)
            rotateStress(3) =  c2*inputStress(3) + s2*inputStress(1) + 2*cs*inputStress(n13)
            rotateStress(1) =  s2*inputStress(3) + c2*inputStress(1) - 2*cs*inputStress(n13)
            rotateStress(2) =  inputStress(2)
            rotateStress(n13) = -cs*inputStress(3) + cs*inputStress(1) + (c2-s2)*inputStress(n13)
            rotateStress(4) =  cs*inputStress(4) - s1*inputStress(n23)
            rotateStress(n23) =  s1*inputStress(4) + c1*inputStress(n23)
        case default
            rotateStress(1) =  c2*inputStress(1) + s2*inputStress(2) + 2*cs*inputStress(4)
            rotateStress(2) =  s2*inputStress(1) + c2*inputStress(2) - 2*cs*inputStress(4)
            rotateStress(3) =  inputStress(3)
            rotateStress(4) = -cs*inputStress(1) + cs*inputStress(2) + (c2-s2)*inputStress(4)
            rotateStress(n23) =  cs*inputStress(n23) - s1*inputStress(n13)
            rotateStress(n13) =  s1*inputStress(n23) + c1*inputStress(n13)
        end select

    end function
    
    
    !_matrixTractions___________________________________________________
    ! retrieves the traction components in the fracture plane
    !___________________________________________________________________
    function matrixTractions(plyStresses, angleRad)
    !-------------------------------------------------------------------!
    !   plyStresses     stress vector (6)
    !   angleRad        angle to rotate in radians
    !-------------------------------------------------------------------!
    !   modified on 26.11.2019 for nS23 and nS13
    !-------------------------------------------------------------------!
    
        implicit none

        ! arguments
        real(RWP), intent(in)   ::  plyStresses(6), angleRad
        real(RWP)               ::  matrixTractions(3)
        
        ! locals
        real(RWP)               ::  cosv, sinv, sin2v
        
        ! initialize
        cosv  = cos(angleRad)
        sinv  = sin(angleRad)
        sin2v = sin(2*angleRad)
        
        matrixTractions(1)  = plyStresses(2) * (cosv*cosv) &
            &               + plyStresses(3) * (sinv*sinv) &
            &               + plyStresses(n23) * sin2v
        matrixTractions(2)  = HALF * (plyStresses(3) - plyStresses(2)) * sin2v &
            &               + plyStresses(n23) * COS(2*angleRad)
        matrixTractions(3)  = plyStresses(4) * cosv + plyStresses(n13) * sinv
        
    end function
    
    
    !_plyCriteria2D__________________________________________________________
    ! evaluates all criteria related to a single ply, with 2D stress tensors
    !________________________________________________________________________
    function plyCriteria2D( plyProperties, ply2DStresses, oldPlyIndexes,    &
        &                   limitFIDen                                      )
    !-----------------------------------------------------------------------!
    !   plyProperties   PlyDefinition with ply material properties
    !   ply2DStresses   ply stresses in local reference frame
    !   oldPlyIndexes   ply strains in local reference frame
    !   limitFIDen      denominator limitation for damage propagation
    !-----------------------------------------------------------------------!
   
        implicit none

        ! arguments
        type(PlyDefinition), intent(in)     :: plyProperties
        real(RWP), intent(in)               :: ply2DStresses(3)
        real(RWP), intent(in)               :: oldPlyIndexes(4)
        logical, intent(in), optional       :: limitFIDen
        ! output
        real(RWP)                           :: plyCriteria2D(4)
        
        ! local settings
        real(RWP)       ::  matrixAngles(2), kinkAngles(2)

        ! local variables
        real(RWP)       ::  angle, cosv, sinv,                              &
            &               aux1, aux2, aux3, currentMax, currentIndex,     &
            &               phi,                                            &
            &               SN, SST, SSL,                                   &
            &               s1mod, s2mod
        real(RWP)       ::  auxStresses(6)
        logical         ::  fiLimitSn = .False.
        integer         ::  i
        
        ! initializations
        if (PRESENT(limitFIDen)) then
            fiLimitSn = limitFIDen
            if (DEBUG_MODE) print *, '(plyCriteria2D): sigmaN denominator control active'
        end if
        
        ! more settings
        matrixAngles    = (/ ZERO, plyProperties%a0/)
        kinkAngles      = (/ ZERO, 90.0000D0/)

        ! debug
        if (DEBUG_MODE) print *, '(plyCriteria2D)'
        plyCriteria2D = oldPlyIndexes
        
        ! ***************************************************************
        ! Matrix cracking       <- [/] ->
        if (DEBUG_MODE) print *, '(plyCriteria2D): Evaluating matrix cracking'
        currentMax = ZERO
        ! iterate failure angles
        do i=1, size(matrixAngles)
            angle   = matrixAngles(i)*DEGTORAD
            if (DEBUG_MODE) print *, '(plyCriteria2D): Evaluating matrix cracking: angle',angle
            cosv    = cos(angle)
            sinv    = sin(angle)
            ! potential fracture plane
            SN      = ply2DStresses(2) * (cosv*cosv) 
            SST     = - HALF * (ply2DStresses(2)) * SIN(2*angle)
            SSL     = ply2DStresses(3) * cosv
            !if (DEBUG_MODE) print '(E11.4, E11.4, E11.4)', SN, SST, SSL  
            ! calculate failure index
            currentIndex = plyEvaluateCriteria(plyProperties, SST, SSL, SN, limitSigmaN = fiLimitSn)
            if (DEBUG_MODE) print *, '(plyCriteria2D):                           : index',currentIndex
            currentMax   = max(currentMax,currentIndex)
        end do
        
        plyCriteria2D(1) = max( currentMax, oldPlyIndexes(1) )

        ! ***************************************************************
        ! Fibre tension         <- [-/-] ->
        if (DEBUG_MODE) print *, '(plyCriteria2D): Evaluating fibre tension xxx'
        currentIndex  = max(ZERO,ply2DStresses(1)) / plyProperties%Xt
        
        plyCriteria2D(3) = max( currentIndex, oldPlyIndexes(3) )

        ! ***************************************************************
        ! Matrix splitting and Fibre Kinking       ->[/]<-
        if (ply2DStresses(1) .lt. NEG_STRESS_THS) then
            if (DEBUG_MODE) print *, '(plyCriteria2D): Evaluating splitting and kinking'
            currentMax = ZERO
            ! find psi (kink band angle)
            aux1  = (plyProperties%G12 - plyProperties%Xc) * plyProperties%phiC
            if (DEBUG_MODE) print *, '(plyCriteria2D): Number of kink angles', size(kinkAngles)
            do i=1, size(kinkAngles)
                angle = kinkAngles(i)*DEGTORAD
                if (DEBUG_MODE) print *, '(plyCriteria2D): Trying kink band: angle', angle
                auxStresses(:) = ZERO
                auxStresses(1) = ply2DStresses(1)
                auxStresses(2) = ply2DStresses(2)
                auxStresses(4) = ply2DStresses(3)
                    
                auxStresses = rotateStress(auxStresses,angle,1) ! rotate around x ! Fixed 16.05.2019
                aux2  = ABS(auxStresses(4)) + aux1

                s1mod = MAX(-plyProperties%Xc,auxStresses(1))
                s2mod = MAX(-plyProperties%Yc, MIN(auxStresses(2), plyProperties%Yt))
                aux3  = plyProperties%G12 + s1mod - s2mod
                if (DEBUG_MODE) print *, '(plyCriteria): AUX3N = ', aux3
                phi   = SIGN(ONE,auxStresses(4)) * aux2/aux3    ! this sign will not return zero
                if (DEBUG_MODE) print *, '(plyCriteria2D): phi = ', phi
                if (DEBUG_MODE) print *, '(plyCriteria2D): Initial stresses'
                if (DEBUG_MODE) print '       (E10.2, E10.2,E10.2,E10.2,E10.2,E10.2)', ply2DStresses
                if (DEBUG_MODE) print *, '(plyCriteria2D): Rotated stresses'
                if (DEBUG_MODE) print '       (F10.3, F10.3, F10.3, F10.3, F10.3, F10.3)', auxStresses
                auxStresses = rotateStress(auxStresses,phi,3) ! rotate around z ! Fixed 16.05.2019
                if (DEBUG_MODE) print *, '(plyCriteria2D): Rotated stresses'
                if (DEBUG_MODE) print '       (F10.3, F10.3, F10.3, F10.3, F10.3, F10.3)', auxStresses
                    
                currentIndex = plyEvaluateCriteria(plyProperties,       &
                            &                      auxStresses(n23),    &
                            &                      auxStresses(4),      &
                            &                      auxStresses(2),      &
                            &                      limitSigmaN = fiLimitSn )
                    
                if (DEBUG_MODE) print *, '(plyCriteria2D):                 : index', currentIndex
                    
                currentMax = max(currentMax,currentIndex)
            end do
            if (ply2DStresses(1) .lt. -0.500d0*plyProperties%Xc) then
                ! fibre kinking
                if (DEBUG_MODE) print *, '(plyCriteria2D): Updating fibre kinking'
                plyCriteria2D(4) = max( currentMax, oldPlyIndexes(4) )
            else
                if (DEBUG_MODE) print *, '(plyCriteria2D): Updating matrix splitting'
                plyCriteria2D(2) = max( currentMax, oldPlyIndexes(2) )
            end if
        end if

    end function plyCriteria2D


    !_completeCriteria_______________________________________________________
    ! evaluates all criteria related to a single ply
    !________________________________________________________________________
    function completeCriteria( plyProperties, plyStresses, previousIndexes, &
        &                      NSTRS, NINDEXES,                             &
        &                      threeDimensions, matrixInterfaceIndex,       &
        &                      limitFIDen)
    !-----------------------------------------------------------------------!
    !   plyProperties           PlyDefinition w/ ply material properties
    !   plyStresses             ply stresses in local reference frame
    !   previousIndexes         history maximum indexes
    !   NSTRS                   number of stress components
    !   NINDEXES                number of indexes
    !   threeDimensions         flag for 3D kink and matrix angles
    !   matrixInterfaceIndex    flag for matrix interface index (FI5)
    !   limitFIDen              limit denominator - damage propagation
    !-----------------------------------------------------------------------!
    
        implicit none

        ! arguments
        integer, intent(in)             ::  NSTRS, NINDEXES
        type(PlyDefinition), intent(in) ::  plyProperties
        real(RWP), intent(in)           ::  plyStresses(NSTRS)
        real(RWP), intent(in)           ::  previousIndexes(NINDEXES)
        logical, intent(in)             ::  threeDimensions
        logical, intent(in), optional   ::  matrixInterfaceIndex, limitFIDen
        
        ! return
        real(RWP)                       ::  completeCriteria(NINDEXES)
        
        ! locals
        real(RWP), allocatable  ::  matrixAngles(:), kinkAngles(:)
        integer                 ::  alocStatus, ikink, ii=1
        
        ! more locals
        real(RWP)       ::  angle, cosv, sinv,                              &
            &               aux1, aux2, aux3, currentMax, currentIndex,     &
            &               phi,                                            &
            &               SN, SST, SSL,                                   &
            &               s1mod, s2mod
        real(RWP)       ::  auxStresses(6), stressVector(6)
        logical         ::  toComputeMI = .False., fiLimitSn = .False.
        ! --------------------------------------------------------------


        ! initialize
        if (DEBUG_MODE) print 19122,'completeCriteria','Initializing'
        
        auxStresses(:) = ZERO
        stressVector(:) = ZERO
        ! fix the stress vector if it is 2D
        if (NSTRS .EQ. 6) then
            stressVector(:) = plyStresses(:)
        else
            stressVector(1) = plyStresses(1)
            stressVector(2) = plyStresses(2)
            stressVector(4) = plyStresses(3)
        end if
        ! matrix interface FI flag
        if (PRESENT(matrixInterfaceIndex)) then
            toComputeMI = (matrixInterfaceIndex .AND. (NSTRS .EQ. 6))
        end if
        ! FI for damage propagation flag
        if (PRESENT(limitFIDen)) then
            fiLimitSn = limitFIDen
        end if
        
        ! initialize 
        do ii=1,NINDEXES
            completeCriteria(ii) = previousIndexes(ii)
        end do
        
        ! allocate       
        if (threeDimensions) then
            allocate(matrixAngles(17),stat=alocStatus)  
            ! explicit initialization is required
            !   or Abaqus will abort the analysis with no message
            matrixAngles(:)  = ZERO
            matrixAngles(2)  = 5.00D0
            matrixAngles(3)  = 10.00D0
            matrixAngles(4)  = 15.00D0
            matrixAngles(5)  = 30.00D0
            matrixAngles(6)  = 45.00D0
            matrixAngles(7)  = plyProperties%a0
            matrixAngles(8)  = 60.00D0 
            matrixAngles(9)  = 75.00D0
            matrixAngles(10) = 90.00D0
            matrixAngles(11) = 105.00D0
            matrixAngles(12) = 120.00D0
            matrixAngles(13) = 180.00D0 - plyProperties%a0
            matrixAngles(14) = 135.00D0
            matrixAngles(15) = 150.00D0
            matrixAngles(16) = 165.00D0
            matrixAngles(17) = 170.00D0
            
            ikink = int(180.0d0/KINK_ANGLE_SPACING)
            allocate(kinkAngles(ikink),stat=alocStatus)
            do ii=1,ikink
                kinkAngles(ii) = ZERO + float(ii-1) * KINK_ANGLE_SPACING
            end do

        else
            allocate(matrixAngles(2),stat=alocStatus) 
            allocate(kinkAngles(2),  stat=alocStatus)
            ! explicit initialization is required
            matrixAngles(:) = ZERO
            matrixAngles(2) = plyProperties%a0
            kinkAngles(:)   = ZERO
            kinkAngles(2)   = 90.0000D0
        end if

        
        
        ! ***************************************************************
        ! Matrix cracking       <- [/] ->
        if (DEBUG_MODE) print 19122,'completeCriteria','Evaluating matrix cracking'
        currentMax = ZERO
        ! iterate failure angles
        do ii=1, size(matrixAngles)
            ! angle
            angle   = matrixAngles(ii)*DEGTORAD
            if (DEBUG_MODE) print 19123,'completeCriteria','Evaluating matrix cracking angle', angle
            cosv    = cos(angle)
            sinv    = sin(angle)
            
            ! potential fracture plane
            SN      = stressVector(2) * (cosv*cosv) &
                    + stressVector(3) * (sinv*sinv) &
                    + stressVector(n23) * SIN(2*angle)
            SST     = HALF * (stressVector(3) - stressVector(2)) * SIN(2*angle) &
                    + stressVector(n23) * COS(2*angle)
            SSL     = stressVector(4) * cosv + stressVector(n13) * sinv
          
            ! calculate failure index
            currentIndex = plyEvaluateCriteria(plyProperties, SST, SSL, SN, limitSigmaN = fiLimitSn)
            if (DEBUG_MODE) print 19123,'completeCriteria','       -> index', currentIndex
            currentMax   = max(currentMax,currentIndex)
        end do
        
        ! update max
        completeCriteria(1) = max(currentMax, completeCriteria(1))

        ! ***************************************************************
        ! Fibre tension         <- [-/-] ->
        if (DEBUG_MODE) print 19122,'completeCriteria','Evaluating fibre tension'
        currentIndex  = max(ZERO,stressVector(1)) / plyProperties%Xt
        completeCriteria(3) = max(currentIndex, completeCriteria(3))
        
        ! ***************************************************************
        ! Matrix splitting and Fibre Kinking       ->[/]<-
        if (stressVector(1) .LT. NEG_STRESS_THS) then ! if in compression
            ! modified to tollerate near zero stresses 22.11.2019
            if (DEBUG_MODE) print 19122,'completeCriteria','Evaluating split and kinking'
            currentMax = ZERO
            ! find psi (kink band angle)
            aux1  = (plyProperties%G12 - plyProperties%Xc) * plyProperties%phiC

            do ii=1, size(kinkAngles)
                angle = kinkAngles(ii)*DEGTORAD
                if (DEBUG_MODE) print 19123,'completeCriteria','Evaluating kink band angle', angle
                auxStresses = rotateStress(stressVector,angle,1) ! rotate around x ! Fixed 16.05.2019
                aux2  = ABS(auxStresses(4)) + aux1
                s1mod = MAX(-plyProperties%Xc,auxStresses(1))
                s2mod = MAX(-plyProperties%Yc, MIN(auxStresses(2), plyProperties%Yt))
                aux3  = plyProperties%G12 + s1mod - s2mod
                phi   = SIGN(ONE,auxStresses(4)) * aux2/aux3    ! this sign will not return zero
                auxStresses = rotateStress(auxStresses,phi,3) ! rotate around z ! Fixed 16.05.2019
                currentIndex = plyEvaluateCriteria(plyProperties,       &
                            &                      auxStresses(n23),    & 
                            &                      auxStresses(4),      &
                            &                      auxStresses(2),      &
                            &                      limitSigmaN = fiLimitSn )
                    
                if (DEBUG_MODE) print 19123,'completeCriteria','       -> index', currentIndex
                currentMax = max(currentMax,currentIndex)
            end do
            ! place the index according to S11
            if (stressVector(1) .lt. -0.500d0*plyProperties%Xc) then
                ! fibre kinking
                if (DEBUG_MODE) print 19122,'completeCriteria','Updating fibre kinking'
                completeCriteria(4) = max( completeCriteria(4), currentMax)
            else
                ! matrix splitting 
                if (DEBUG_MODE) print 19122,'completeCriteria','Updating matrix splitting'
                completeCriteria(2) = max( completeCriteria(2), currentMax)
            end if
        end if
        
        ! ***************************************************************
        ! NFC criteria       ->[/]<-
        if ((NINDEXES .GE. 5) .AND. (toComputeMI)) then
            if (DEBUG_MODE) print 19122,'completeCriteria','Evaluating matrix interface'
            if (stressVector(3) .GT. NZSTRESS) then !NZSTRESS 22.11.2019
                aux1 = stressVector(5) / plyProperties%ILSS
                aux2 = stressVector(6) / plyProperties%ILSS
                aux3 = stressVector(3) / plyProperties%Zt
                currentIndex = (aux1*aux1) + (aux2*aux2) + (aux3*aux3)
                currentIndex = SQRT(currentIndex)
            else
                currentIndex = ZERO
            end if
            completeCriteria(5) = max( completeCriteria(5), currentIndex )
        end if
    
    ! deallocate
    deallocate(matrixAngles, stat=alocStatus)
    deallocate(kinkAngles,   stat=alocStatus)
        
19122 format("(",A,"): ", (A))
19123 format("(",A,"): ", (A),': ', 12(E11.3))
        
    end function completeCriteria

    !_ncfCriteria_______________________________________________________
    !   function to calculate the 5th failure index
    !___________________________________________________________________
    function ncfCriteria( plyProperties, stressVector)

        ! parameters
        type(PlyDefinition), intent(in) ::  plyProperties
        real(RWP), intent(in)           ::  stressVector(6)

        ! output
        real(RWP) ::    ncfCriteria
        ! locals
        real(RWP) ::    aux1, aux2, aux3

        if (stressVector(3) .GT. NZSTRESS) then
            aux1 = stressVector(5) / plyProperties%ILSS
            aux2 = stressVector(6) / plyProperties%ILSS
            aux3 = stressVector(3) / plyProperties%Zt
            ncfCriteria = (aux1*aux1) + (aux2*aux2) + (aux3*aux3)
            ncfCriteria = SQRT(ncfCriteria)
        else
            ncfCriteria = ZERO
        end if

    end function ncfCriteria


    !_maxVector_________________________________________________________
    !   returns the element-wise max of 2 vectors, with the size of the 
    !       first vector
    !___________________________________________________________________
    function max_vector( vector1, vector2 )

        real(RWP), intent(in)   ::  vector1(:), vector2(:)
        real(RWP)               ::  max_vector(size(vector1))
        
        integer :: i 
        
        max_vector(:) = vector1(:)
        
        do i=1,min(size(vector1),size(vector2))
            max_vector(i) = max(vector1(i),vector2(i))
        end do
                   
    end function max_vector
    
    
    !_zeros_____________________________________________________________
    !   returns a zero array with n entries
    !___________________________________________________________________
    function zeros( n )
    
        integer, intent(in) ::  n
        real(RWP)           ::  zeros(n)
        
        zeros(:) = ZERO
        
    end function zeros
        
    
    
    !_zzeros_like_vector________________________________________________
    !   returns a zero array with the same length as the input vector
    !___________________________________________________________________
    function zeros_like_vector( vector )
    
        real(RWP)           ::  vector(:)
        real(RWP)           ::  zeros_like_vector(size(vector))
        
        zeros_like_vector(:) = ZERO
        
    end function zeros_like_vector
     
    
    !_round_significant
    !
    !___________________________________________________________________
    function round_significant( x, sdigits, zeroThreshold)
    
        real(RWP), intent(in)           ::  x
        integer, intent(in)             ::  sdigits
        real(RWP), intent(in), optional ::  zeroThreshold
	
        real(RWP)   ::  round_significant
        real(RWP)   ::  th = NZSTRESS, dd, s = ONE
        
        round_significant = ZERO
        
        if (present(zeroThreshold)) then
            th = zeroThreshold
        end if
        
        if (ABS(x) .LT. th) then
            round_significant = ZERO
        else 
            if (x .LT. ZERO) then
                s = -s
            end if
            dd = TEN**(sdigits-CEILING(LOG10(s*x)))
            round_significant = s*float(NINT(x*dd))/dd
        end if
  
    end function round_significant
    
    !_round_significant_vector
    !
    !___________________________________________________________________
    function round_significant_vector( vector, sdigits, zeroThreshold)
        real(RWP), intent(in)           ::  vector(:)
        integer, intent(in)             ::  sdigits
        real(RWP), intent(in), optional ::  zeroThreshold
    
        real(RWP)   ::  round_significant_vector(size(vector))
        integer     ::  i
        
        round_significant_vector(:) = ZERO
        
        do i=1,size(vector)
            if (present(zeroThreshold)) then
                round_significant_vector(i) = &
                    & round_significant(vector(i),sdigits,zeroThreshold)
            else 
                round_significant_vector(i) = &
                    & round_significant(vector(i),sdigits)
            end if
        end do
        
    end function round_significant_vector
    
    
    
    !_cumulativeFailureIndexes
    !
    !___________________________________________________________________
    function cumulativeFailureIndexes( currentFIs,                      &
        &                              oldFIs,                          &
        &                              NINDEXES                         )
        
        
    integer, intent(in)     ::  NINDEXES
    real(RWP), intent(in)   ::  currentFIs(NINDEXES), oldFIs(NINDEXES)
    
    real(RWP)   :: cumulativeFailureIndexes(NINDEXES)
    
    real(RWP)   ::  aux2, aux4
    integer     ::  iActive
    
    ! initialize
    cumulativeFailureIndexes(1:NINDEXES) = currentFIs(1:NINDEXES)
    
    ! get both FIs
    aux2 = MAX( currentFIs(2), oldFIs(2) )
    aux4 = MAX( currentFIs(4), oldFIs(4) )
    
    ! if both active
    if ((aux2 .GT. ONE) .AND. (aux4 .GT. ONE)) then
        if (currentFIs(2) .GT. ONE) then
            iActive = 2
        else
            iActive = 4
        end if
        ! already active/activated now
        if (oldFIs(iActive) .GT. ONE) then
            aux2 = oldFIs(2) + oldFIs(4) - ONE
            cumulativeFailureIndexes(iActive) = oldFIs(iActive)
        else
            aux2 = oldFIs(2) + oldFIs(4)
            cumulativeFailureIndexes(iActive) = ONE
        end if
        ! increment
        aux4 = MAX(ZERO, currentFIs(iActive) - aux2)
        ! increment current index
        cumulativeFailureIndexes(iActive) = cumulativeFailureIndexes(iActive) + aux4
    end if     
    
    end function cumulativeFailureIndexes
    

    !_limitValue
    ! returns the variable limited by minValue and maxValue 
    !___________________________________________________________________
    function limitVariable(variable, minValue, maxValue)

    real(RWP), intent(in)           ::  variable
    real(RWP), intent(in), optional ::  minValue, maxValue

    real(RWP)   ::  limitVariable

    limitVariable = variable

    if (present(minValue)) then
        limitVariable = MAX(minValue, variable)
    end if
    if (present(maxValue)) then
        limitVariable = MIN(maxValue, variable)
    end if

    end function limitVariable


end module failureCriteria

