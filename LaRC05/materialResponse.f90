!DIR$ FREEFORM
!  materialResponse.f90 
!   
!  FUNCTIONS:
!   constitutiveMatrix2D    constitutive matrix
!   macauley                Mac Cauley brackets with a judge entry
!   computeStresses2D       compute 2D stress vector
!   computeDamageVariables  compute damage (returns health) variables
!   computeEnergyVariation  compute change in energy during an increment
!
!****************************************************************************
!****************************************************************************
!                   developed by Miguel A.S. Matos and Silvestre T. Pinho
!                                       Imperial College London, Aeronautics
!                                                       created 03.06.2019
!****************************************************************************
!                   continued at Airbus by Miguel Matos (mmatos)
!****************************************************************************
!   log:
!       12.06.2019  Fixes for lambda
!       02.07.2019  Fixes for better mixed mode representation for Gmatrix 
!       03.07.2019  Final comments, debug disabled
!       06.11.2019  Modifications to support material behaviour in 3D
!       21.11.2019  Finalized updateDamageVariables
!       22.11.2019  Added computeStressPlane: compute stresses in 2D 
!                   but able to handle S/E vectors of 3 or 6
!   last error free build on:
!       03.07.2019
!       22.11.2019
!****************************************************************************
    
module materialResponse

use plyTools
use failureCriteria

implicit none

integer, parameter, private     ::  N_DAMAGE_VARIABLES = 3
integer, parameter, private     ::  N_DAMAGE_VARIABLES_3D = 6
integer, parameter, private     ::  RWP = KIND(1.0D0)
!
real(RWP), parameter, private   ::  ONE = 1.0000D0, ZERO = 0.0000D0,    &
    &                               HALF = 0.50000D0,                   &
    &                               NEAR_ZERO_POS = 1.00D-6,            &
    &                               NEAR_ZERO_NEG = -1.00D-6
real(RWP), parameter, private   ::  MAX_LAMBDA_F = 1.200

real(RWP), parameter, private   ::  DAMAGE_LIMIT = 0.9999D0
real(RWP), parameter, private   ::  MIN_MATERIAL_HEALTH = 0.00001D0

real(RWP), parameter, private   ::  MIN_DMG_COMBINED = 0.0100D00

real(RWP), parameter, private   ::  REG_VISCOSITY = 0.0002500D0

contains
    
    !_constitutiveMatrix2D_____________________________________________________
    ! returns an array with non-zero entries of the constitutive matrix D:
    !   [D11, D12, D21, D22, D33]
    !________________________________________________________________________
    function constitutiveMatrix2D( plyProperties, materialDamage )
    !-----------------------------------------------------------------------!
    !   plyProperties   PlyDefinition with ply material properties
    !   materialDamage  material health variables (1 is undamaged)
    !-----------------------------------------------------------------------!
    
        implicit none

        ! parameters
        !integer, parameter      :: RWP  = KIND(1.0D0)
        !real(RWP), parameter    :: ONE  = 1.000000D0
        
        ! arguments        
        type(PlyDefinition), intent(in) :: plyProperties
        real(RWP), intent(in)           :: materialDamage(N_DAMAGE_VARIABLES)
        real(RWP)                       :: constitutiveMatrix2D(5)
        
        ! locals
        real(RWP)                       :: k0, E11
        
        E11 = getE11(plyProperties)
        
        ! constants
        k0  = plyProperties%nu12 * plyProperties%nu21
        k0 = k0 * materialDamage(1) * materialDamage(2)
        k0 = ONE - k0
        
        ! matrix entries
        constitutiveMatrix2D(1) = E11 * materialDamage(1) / k0
        constitutiveMatrix2D(2) = plyProperties%nu21 * E11 *  &
            &                       materialDamage(1) * materialDamage(2) / k0
        constitutiveMatrix2D(3) = plyProperties%nu12 * plyProperties%E22 * &
            &                       materialDamage(1) * materialDamage(2) / k0
        constitutiveMatrix2D(4) = plyProperties%E22 * materialDamage(2) / k0
        constitutiveMatrix2D(5) = plyProperties%G12 *  materialDamage(3)
        
    end function
    
    
    
    !_constitutiveMatrix3D___________________________________________________
    ! returns the constitutive matrix D
    !      1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12
    !   [D11, D12, D21, D22, D44, D13, D31, D23, D32, D33, D55, D66]
    !       note the order of the outputs (compatible with the 2D version)
    !________________________________________________________________________
    function constitutiveMatrix3D( plyProperties, materialHealth)
    !-----------------------------------------------------------------------!
    !   plyProperties   PlyDefinition with ply material properties
    !   materialHealth  material health variables (1 is undamaged) = (1-d)
    !
    !   20/19/2019      modified to include all damage variables
    !   20/19/2019      verified by hand
    !-----------------------------------------------------------------------!
    
        implicit none
        
        ! parameters
        ! none
        
        ! arguments
        type(PlyDefinition), intent(in) ::  plyProperties
        real(RWP), intent(in)           ::  materialHealth(N_DAMAGE_VARIABLES_3D)
        ! output
        real(RWP)                       ::  constitutiveMatrix3D(12)
        
        ! locals
        real(RWP)                       ::  k0 = ZERO
        real(RWP)                       ::  aux1,aux2,aux3,aux4, faux
        real(RWP)                       ::  nu31, nu32
        real(RWP)                       ::  E11
        
                    
        ! additional material properties
        E11  = getE11(plyProperties)
        nu31 = plyProperties%nu13 * plyProperties%E33 / plyProperties%E11
        nu32 = plyProperties%nu23 * plyProperties%E33 / plyProperties%E22
                
        ! initializations
        constitutiveMatrix3D(:) = ZERO
        
        ! constants
        aux1 = plyProperties%nu12 * plyProperties%nu21
        aux1 = aux1 * materialHealth(1) * materialHealth(2)
        aux2 = plyProperties%nu23 * nu32
        aux2 = aux2 * materialHealth(2) * materialHealth(3)
        aux3 = nu31 * plyProperties%nu13
        aux3 = aux3 * materialHealth(3) * materialHealth(1)
        aux4 = plyProperties%nu21 * nu32 * plyProperties%nu13
        aux4 = aux4 * materialHealth(1)*materialHealth(2)*materialHealth(3)
        ! k0
        k0   = aux1 + aux2 + aux3 + aux4
        k0   = ONE - k0
        
        ! diagonal components
        faux = ONE - aux2 !C11
        constitutiveMatrix3D(1) = E11 * materialHealth(1) * faux / k0
        faux = ONE - aux3 !C22
        constitutiveMatrix3D(4) = plyProperties%E22 * materialHealth(2) * faux/ k0
        faux = ONE - aux1 !C33
        constitutiveMatrix3D(10) = plyProperties%E33 * materialHealth(3) * faux / k0
        
        ! cross components
        !C12
        faux = materialHealth(3) * nu31 * materialHealth(2) * plyProperties%nu23
        faux = faux + materialHealth(2) * plyProperties%nu21
        constitutiveMatrix3D(2) = materialHealth(1) * E11 * faux / k0
        !C13
        faux = materialHealth(2) * plyProperties%nu21 * materialHealth(3) * nu32
        faux = faux + materialHealth(3) * nu31
        constitutiveMatrix3D(6) = materialHealth(1) * E11 * faux / k0
        
        !C21
        faux = materialHealth(3) * nu32 * materialHealth(1) * plyProperties%nu13
        faux = faux + materialHealth(1) * plyProperties%nu12
        constitutiveMatrix3D(3) = materialHealth(2) * plyProperties%E22 * faux / k0
        !C23
        faux = materialHealth(1) * plyProperties%nu12 * materialHealth(3) * nu31
        faux = faux + materialHealth(3) * nu32
        constitutiveMatrix3D(8) = plyProperties%E22 * materialHealth(2) * faux / k0
        
        !C31
        faux = materialHealth(1) * plyProperties%nu12 * materialHealth(2) * plyProperties%nu23
        faux = faux + materialHealth(1) * plyProperties%nu13
        constitutiveMatrix3D(7) = plyProperties%E33 * materialHealth(3) * faux / k0
        !C32
        faux = materialHealth(2) * plyProperties%nu21 * materialHealth(1) * plyProperties%nu13
        faux = faux + materialHealth(2) * plyProperties%nu23
        constitutiveMatrix3D(9) = plyProperties%E33 * materialHealth(3) * faux / k0
                
        ! shear components
        constitutiveMatrix3D(5)  = plyProperties%G12 *  materialHealth(4)
        constitutiveMatrix3D(11) = plyProperties%G13 *  materialHealth(5)
        constitutiveMatrix3D(12) = plyProperties%G23 *  materialHealth(6)
    
    end function
    
    
    !_macauley_brackets______________________________________________________
    ! Macauley brackets <value>+
    !   if judge is present value * <sign(judge)>
    !________________________________________________________________________
    function macauley( value, judge)
        
        implicit none
        
        real(RWP), intent(in)           ::  value
        real(RWP), intent(in), optional ::  judge
        
        real(RWP)                       ::  macauley
        
        macauley = value
        if ( PRESENT(judge) ) then
            if (judge .LT. ZERO) then
                macauley = ZERO
            end if
        else
            if (value .LT. ZERO) then
                macauley = ZERO
            end if
        end if
        
    end function
        
        
    !_macauley_brackets______________________________________________________
    ! Macauley brackets <value>+
    !   if judge is present value * <judge is True>
    !________________________________________________________________________
    function mabrackets( invalue, judge)
        
        implicit none
        
        real(RWP), intent(in)           ::  invalue
        logical,   intent(in), optional ::  judge
        
        
        real(RWP)                       ::  mabrackets
        
        ! initialize
        mabrackets = invalue
        ! look at the judge
        if ( PRESENT(judge) ) then
            if (.NOT. judge) then
                mabrackets = ZERO
            end if
        else
        ! new expression
            mabrackets = HALF * ( invalue + ABS(invalue) )
        end if
        
    end function
        
        
    !_convertVector2Dto3D______________________________________________________
    ! convets a vector from 2D to 3D
    !   [S11, S22, S12] -> [S11, S22, S33, S12, S23, S31]
    !________________________________________________________________________
    function convertVector2Dto3D( vector2D )
    
        implicit none
        
        real(RWP), intent(in)   ::  vector2D(3)
        
        real(RWP)               ::  convertVector2Dto3D(6)
        
        convertVector2Dto3D(:) = ZERO
        convertVector2Dto3D(1) = vector2D(1)
        convertVector2Dto3D(2) = vector2D(2)
        convertVector2Dto3D(4) = vector2D(3)
        
    end function
        
    !_convertVector3Dto2D______________________________________________________
    ! convets a vector from 3D to 2D
    !   [S11, S22, S33, S12, S23, S31] -> [S11, S22, S12]
    !________________________________________________________________________
    function convertVector3Dto2D( vector3D )
    
        implicit none
        
        real(RWP), intent(in)   ::  vector3D(6)
        
        real(RWP)               ::  convertVector3Dto2D(3)
    
        convertVector3Dto2D(:) = ZERO
        convertVector3Dto2D(1) = vector3D(1)
        convertVector3Dto2D(2) = vector3D(2)
        convertVector3Dto2D(3) = vector3D(4)
    
    end function
        
    !_computeStresses2D______________________________________________________
    ! returns the planar stress vector
    !   [S11, S22, S12]
    !________________________________________________________________________
    function computeStresses2D( strainVector, plyProperties, materialDamage,&
        &                       constitutiveMatrix )
    !-----------------------------------------------------------------------!
    !   plyProperties   PlyDefinition with ply material properties
    !   materialDamage  Material damage (1 is undamaged)
    !-----------------------------------------------------------------------!
    
        implicit none
        
        ! arguments        
        type(PlyDefinition), intent(in) ::  plyProperties
        real(RWP), intent(in)           ::  strainVector(3)
        real(RWP), intent(in), optional ::  materialDamage(N_DAMAGE_VARIABLES)
        real(RWP), intent(out),optional ::  constitutiveMatrix(3,3)
        real(RWP)                       ::  computeStresses2D(3)
        
        ! locals
        real(RWP)                       ::  dMatrixEntries(5),  &
            &                               currentMaterialDamage(N_DAMAGE_VARIABLES)
        
        if (PRESENT(materialDamage)) then
            currentMaterialDamage = materialDamage
        else
            currentMaterialDamage(:) = ONE
        end if
        
        ! get matrix
        dMatrixEntries = constitutiveMatrix2D( plyProperties, currentMaterialDamage )
        ! multiply
        computeStresses2D(:) = ZERO
        computeStresses2D(1) = dMatrixEntries(1) * strainVector(1) + dMatrixEntries(2) * strainVector(2)    ! S1
        computeStresses2D(2) = dMatrixEntries(3) * strainVector(1) + dMatrixEntries(4) * strainVector(2)    ! S2
        computeStresses2D(3) = dMatrixEntries(5) * strainVector(3) ! S3
        
        ! return the matrix if requested
        if ( PRESENT(constitutiveMatrix) ) then
            constitutiveMatrix(:,:) = ZERO
            constitutiveMatrix(1,1) = dMatrixEntries(1)
            constitutiveMatrix(1,2) = dMatrixEntries(2)
            constitutiveMatrix(2,1) = dMatrixEntries(3)
            constitutiveMatrix(2,2) = dMatrixEntries(4)
            constitutiveMatrix(3,3) = dMatrixEntries(5)
        end if           
        
    end function
    
    
    !_computeStressPlane______________________________________________
    ! returns the full stress tensor considering plane stress 
    !   [S11, S22, S12] or [S11,S22,S33,S12,S23,S31] (S33,S23,S31=0)
    !       depending on the specified NSTR
    !___________________________________________________________________
    function computeStressPlane( plyProperties,                         &
        &                        strainVector, NSTR,                    &
        &                        materialHealth,                        &
        &                        constitutiveMatrix                     )
    !-------------------------------------------------------------------!
    !   plyProperties   PlyDefinition with ply material properties
    !   materialHealth  Material damage (1 is undamaged)
    !-------------------------------------------------------------------!
    
        implicit none
        
        ! arguments        
        type(PlyDefinition), intent(in) ::  plyProperties
        integer, intent(in)             ::  NSTR
        real(RWP), intent(in)           ::  strainVector(NSTR)
        real(RWP), intent(in), optional ::  materialHealth(N_DAMAGE_VARIABLES)
        real(RWP), intent(out),optional ::  constitutiveMatrix(NSTR,NSTR)
        ! output
        real(RWP)                       ::  computeStressPlane(NSTR)
        
        ! locals
        real(RWP)   ::  dMatrixEntries(5),  &
            &           currentMaterialHealth(N_DAMAGE_VARIABLES)
        integer     ::  nS12 = 3
        
        ! fix nS12 for 3D
        if (NSTR .EQ. 6) then
            nS12 = 4
        end if
        
        if (DEBUG_MODE) print *, '(computeStressPlane)'
        
        ! initialize optoinals        
        if (PRESENT(materialHealth)) then
            currentMaterialHealth = materialHealth
        else
            currentMaterialHealth(:) = ONE
        end if
        
        ! get matrix
        dMatrixEntries = constitutiveMatrix2D( plyProperties, currentMaterialHealth )

        ! multiply
        computeStressPlane(:) = ZERO
        computeStressPlane(1) = dMatrixEntries(1) * strainVector(1) + &
            &                   dMatrixEntries(2) * strainVector(2)    ! S1
        computeStressPlane(2) = dMatrixEntries(3) * strainVector(1) + &
            &                   dMatrixEntries(4) * strainVector(2)    ! S2
        computeStressPlane(nS12) = dMatrixEntries(5) * strainVector(nS12) ! S3
        
        ! return the matrix if requested
        if ( PRESENT(constitutiveMatrix) ) then
            constitutiveMatrix(:,:) = ZERO
            constitutiveMatrix(1,1) = dMatrixEntries(1)
            constitutiveMatrix(1,2) = dMatrixEntries(2)
            constitutiveMatrix(2,1) = dMatrixEntries(3)
            constitutiveMatrix(2,2) = dMatrixEntries(4)
            constitutiveMatrix(nS12,nS12) = dMatrixEntries(nS12)
        end if           
        
    end function
    
    
    !_computeStresses______________________________________________________
    ! returns the complete stress vector
    !   [S11, S22, S33, S12, S23, S31]
    !________________________________________________________________________
    function computeStresses3D( strainVector, plyProperties, materialHealth,&
        &                       constitutiveMatrix )
    !-----------------------------------------------------------------------!
    !   strainVector        3D (6x1) strain vector
    !   plyProperties       PlyDefinition with ply material properties
    !   materialHealth      Material health (1 is undamaged) = 1 - damage
    !   constitutiveMatrix  if present, will output the constitutive matrix
    !-----------------------------------------------------------------------!
    
        implicit none
        
        ! arguments
        type(PlyDefinition), intent(in)     ::  plyProperties
        real(RWP), intent(in)               ::  strainVector(6)
        real(RWP), intent(in), optional     ::  materialHealth(N_DAMAGE_VARIABLES_3D)
        real(RWP), intent(out), optional    ::  constitutiveMatrix(6,6)
        
        ! output
        real(RWP)               ::  computeStresses3D(6)
        
        ! arguments
        real(RWP)               ::  dMatrixEntries(12),     &
            &                       currentMaterialHealth(N_DAMAGE_VARIABLES_3D)
            
        ! get material damage
        if (PRESENT(materialHealth)) then
            currentMaterialHealth = materialHealth
        else
            currentMaterialHealth(:) = ONE
        end if 
        
        ! retrieve constitutive matrix
        dMatrixEntries = constitutiveMatrix3D( plyProperties, currentMaterialHealth)
        !    1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12
        ! [D11, D12, D21, D22, D44, D13, D31, D23, D32, D33, D55, D66]
        
        computeStresses3D(:) = ZERO
        computeStresses3D(1) = dMatrixEntries(1) * strainVector(1)
        computeStresses3D(1) = computeStresses3D(1) + &
            &                  dMatrixEntries(2) * strainVector(2)
        computeStresses3D(1) = computeStresses3D(1) + &
            &                  dMatrixEntries(6) * strainVector(3)
        !
        computeStresses3D(2) = dMatrixEntries(4) * strainVector(2)
        computeStresses3D(2) = computeStresses3D(2) + &
            &                  dMatrixEntries(3) * strainVector(1)
        computeStresses3D(2) = computeStresses3D(2) + &
            &                  dMatrixEntries(8) * strainVector(3)
        !
        computeStresses3D(3) = dMatrixEntries(10) * strainVector(3)
        computeStresses3D(3) = computeStresses3D(3) + &
            &                  dMatrixEntries(7) * strainVector(1)
        computeStresses3D(3) = computeStresses3D(3) + &
            &                  dMatrixEntries(9) * strainVector(2)
                 
        ! shear components
        computeStresses3D(4) = dMatrixEntries(5) * strainVector(4)
        computeStresses3D(5) = dMatrixEntries(11) * strainVector(5)
        computeStresses3D(6) = dMatrixEntries(12) * strainVector(6)
        
        
        ! return the matrix if requested
        if ( PRESENT(constitutiveMatrix) ) then
            constitutiveMatrix(:,:) = ZERO
            constitutiveMatrix(1,1) = dMatrixEntries(1)
            constitutiveMatrix(1,2) = dMatrixEntries(2)
            constitutiveMatrix(1,3) = dMatrixEntries(6)
            constitutiveMatrix(2,1) = dMatrixEntries(3)
            constitutiveMatrix(2,2) = dMatrixEntries(4)
            constitutiveMatrix(2,3) = dMatrixEntries(8)
            constitutiveMatrix(3,1) = dMatrixEntries(7)
            constitutiveMatrix(3,2) = dMatrixEntries(9)
            constitutiveMatrix(3,3) = dMatrixEntries(10)
            constitutiveMatrix(4,4) = dMatrixEntries(5)
            constitutiveMatrix(5,5) = dMatrixEntries(11)
            constitutiveMatrix(6,6) = dMatrixEntries(12)
        end if
        
    end function computeStresses3D
    
            
    !_computeDamageVariables_________________________________________________
    ! computes all damage variables di, and returns the HEALTH variables
    !   hi = (1-di) for i=1,2,s
    !________________________________________________________________________
    function computeDamageVariables( strainVector, stressVector,            &
        &                            failureIndexes, damageVariables,       &
        &                            failureEnergy,                         &
        &                            characteristicLength,                  &
        &                            errorFlag                              )
    !-----------------------------------------------------------------------!
    !   strainVector    ply 2D strain vector
    !   stressVector    ply 2D elastic stress vector
    !   failureIndexes  all 4 failure indexes
    !   damageVariables current damage variables dm,ds,df,dk
    !   failureEnergy   critical energy release rates, and mixed mode exp
    !   characteristicLength    element characteristic length
    !   errorFlag       error flag. not currently in use: inc. size control
    !-----------------------------------------------------------------------!
    
    implicit none
    
    ! arguments
    real(RWP), intent(in)           ::  strainVector(3), stressVector(3),   &
        &                               failureIndexes(4),                  &
        &                               failureEnergy(5),                   &
        &                               characteristicLength
    real(RWP), intent(inout)        ::  damageVariables(4)
    ! optional
    integer, optional, intent(inout)::  errorFlag
    ! out
    real(RWP)                       ::  computeDamageVariables(3)
    
    ! local variables
    real(RWP)   ::  matrixElasticEnergy, fibreElasticEnergy,        & ! per Volume
        &           aux1, aux2

    real(RWP)   ::  currentLambda, failureLambda, currentDelta,         &
        &           currentDamage(4), currentFailureEnergy,             &
        &           onsetEnergy
    integer     ::  imode
      
    ! get elastic energy that can be disspipated by the mode
    fibreElasticEnergy  = HALF * DOT_PRODUCT(strainVector,stressVector)
    matrixElasticEnergy = HALF * macauley(strainVector(2)) * macauley(stressVector(2)) + &
        &                 HALF * strainVector(3) * stressVector(3)  
    
    ! initialize 
    computeDamageVariables(:) = ONE
    currentDamage             = damageVariables
    
    if (DEBUG_MODE) print *, '(computeDamageVariables)'
    
    ! elastic energy at onset
    do imode=1,4
        if (DEBUG_MODE) print *, '(computeDamageVariables): Mode',imode
        if (failureIndexes(imode) .GE. ONE) then
            currentLambda = failureIndexes(imode) * failureIndexes(imode)
            
            ! get onset energy
            if (imode .LT. 3) then  ! matrix
                onsetEnergy = matrixElasticEnergy / currentLambda
                if (DEBUG_MODE) print '(" ",A," :",F12.4)', '(computeDamageVariables): Umat/V  ', matrixElasticEnergy
            else
                if (DEBUG_MODE) print '(" ",A," :",F12.4)', '(computeDamageVariables): Ufib/V  ', fibreElasticEnergy
                onsetEnergy = fibreElasticEnergy / currentLambda
            end if            

            if (DEBUG_MODE) print '(" ",A," :",F12.4)', '(computeDamageVariables): Onset/V ', onsetEnergy
            !
            ! damage driving parameters
            currentFailureEnergy = failureEnergy(imode)
            if (imode .EQ. 1) then
                if (matrixElasticEnergy .GT. NEAR_ZERO_POS) then
                    ! Correction (power and denominator) as of 02.07.2019
                    aux1 = macauley(strainVector(2)*stressVector(2)) / (2.000D0*matrixElasticEnergy)
                    aux2 = strainVector(3)*stressVector(3) / (2.000D0*matrixElasticEnergy)
                    aux1 = aux1*aux1
                    aux2 = aux2*aux2
                    aux1 = (aux1 / failureEnergy(1))
                    aux2 = (aux2 / failureEnergy(2))
                    aux1 = aux1 ** failureEnergy(5)
                    aux2 = aux2 ** failureEnergy(5)
                    currentFailureEnergy = (aux1 + aux2) ** (-ONE/failureEnergy(5))
                end if
            end if
            if (DEBUG_MODE) print '(" ",A," :",F12.4)', '(computeDamageVariables): Ufail/V ', currentFailureEnergy
            ! compute current delta
            if (imode .LT. 3 .AND. (onsetEnergy .LT. NEAR_ZERO_POS) ) then
                failureLambda = ZERO
                currentDelta  = ZERO
            else
                failureLambda = currentFailureEnergy / onsetEnergy
                failureLambda = failureLambda / characteristicLength
                failureLambda = MAX( failureLambda, MAX_LAMBDA_F)
                failureLambda = failureLambda*failureLambda
                ! calculate damage variables
                currentDelta = (SQRT(failureLambda)-failureIndexes(imode)) / (SQRT(failureLambda)-ONE)
                currentDelta = ONE - currentDelta/failureIndexes(imode)
            end if

            if (DEBUG_MODE) print '(" ",A," :",F12.4)','(computeDamageVariables): lambda_f', failureLambda
            if (DEBUG_MODE) print '(" ",A," :",F12.4)', '(computeDamageVariables): lambda  ', currentLambda
            !
            if (DEBUG_MODE) print '(" ",A," :",F12.4)', '(computeDamageVariables): delta   ', currentDelta
            !
            ! damage variable
            currentDamage(imode) = max(currentDelta, damageVariables(imode))
            if (DEBUG_MODE) print '(" ",A," :",F12.4)', '(computeDamageVariables): damageV ', currentDamage(imode)
        end if
        ! limit damage variables
        if (DAMAGE_LIMIT .GT. ZERO) then
            currentDamage(imode) = min(DAMAGE_LIMIT, currentDamage(imode))
        end if
        ! d12 (affected by all)
        computeDamageVariables(3) = computeDamageVariables(3) * (ONE - currentDamage(imode))
        ! d11 and d22
        if (imode .LT. 3) then 
            ! matrix damage
            computeDamageVariables(2) = computeDamageVariables(2) * &
                &   (ONE - macauley(currentDamage(imode),judge=stressVector(2)))
        else
            ! fibre damage
            computeDamageVariables(1) = computeDamageVariables(1) * &
                &   (ONE - currentDamage(imode))
            computeDamageVariables(2) = computeDamageVariables(2) * &
                &   (ONE - currentDamage(imode))
        end if
        damageVariables(imode) = currentDamage(imode)
    end do
    
    if (MIN_MATERIAL_HEALTH .GT. ZERO) then 
        computeDamageVariables(1) = MAX(MIN_MATERIAL_HEALTH, computeDamageVariables(1) )
        computeDamageVariables(2) = MAX(MIN_MATERIAL_HEALTH, computeDamageVariables(2) )
        computeDamageVariables(3) = MAX(MIN_MATERIAL_HEALTH, computeDamageVariables(3) )
    end if

    if (present(errorFlag) .AND. DEBUG_MODE) print '(" ",A)', '(computeDamageVariables): Completed'

    end function
    
     
    !_computeCurrentEnergy___________________________________________________
    ! computes current energy
    !   NOT IN USE/FOR TESTING ONLY, PLEASE DELETE
    !   [elastic;total]        
    !________________________________________________________________________
    function computeCurrentEnergy( currentStrain, currentStress, elasticMatrix)
    !   currentStrain   2D strain vector (3)
    !   currentStress   2D stress vector
    !   elasticMatrix   2D constitutive matrix
    
    implicit none
    
    ! arguments
    real(RWP), intent(in)       ::  currentStrain(3),                       &
        &                           currentStress(3),                       &
        &                           elasticMatrix(3,3)
    
    real(RWP)                   ::  computeCurrentEnergy(2)
    
    ! locals
    real(RWP)                   ::  elasticStress(3)
    
    elasticStress = matmul(elasticMatrix,currentStrain)
    computeCurrentEnergy(1) = HALF * DOT_PRODUCT(elasticStress,currentStrain)
    computeCurrentEnergy(2) = HALF * DOT_PRODUCT(currentStress,currentStrain)
    
    end function    
    
    
    !_computeEnergyVariation________________________________________________
    ! computes the energy change in the current step
    !   elastic energy                          SSE = SSE + dEnergy(1)
    !   plastic dissipation [total - elastic]   SPD = SPD + dEnergy(2)
    !_______________________________________________________________________
    function computeEnergyVariation( initialStrain, dStrain,                &
        &                            initialStress, finalStress,            &
        &                            NSTR                                   )
    ! initialStrain     strain vector at the increment start (NSTR)
    ! dStrain           change in strain during the increment (NSTR)
    ! initialStress     stress vector at the start of the increment
    ! finalStress       stress vector at the end of the increment
    ! NSTR              number of stress/strain components (3 in 2D)
    
    implicit none
    
    integer, intent(in)     ::  NSTR
    real(RWP), intent(in)   ::  initialStrain(NSTR), dStrain(NSTR),         &
        &                       initialStress(NSTR), finalStress(NSTR)
    
    real(RWP)               ::  computeEnergyVariation(2)
    
    ! locals
    real(RWP)               ::  dE, dEElastic
    
    ! total energy = sigma2*deltaEpsilon + 1/2*deltaSigma*deltaEpsilon
    dE = HALF * DOT_PRODUCT(initialStress+finalStress,dStrain)
    ! available elastic energy
    dEElastic = HALF * DOT_PRODUCT(finalStress,initialStrain+dStrain) - &
        &       HALF * DOT_PRODUCT(initialStress,initialStrain)
    
    computeEnergyVariation(1) = dEElastic
    computeEnergyVariation(2) = dE - dEElastic

    end function
        

        
    !_updateDamageVariables_____________________________________________
    ! updates the damage variables for the current increment                
    !   history dependent variables:
    !       FI  (failure indexes 1,2,3,4,5 - no max!)
    !       DMG (damage variables in each of the failure modes)
    !   updated variables
    !       FI
    !       DMG
    !       H (health variables in each of the directions: 1,2,S)
    !___________________________________________________________________
    subroutine updateDamageVariables( plyProperties,                    &
        &                             strainVector, NSTR,               &   
        &                             failureIndexes, NINDEXES,         &
        &                             damageVariables,                  &  
        &                             healthVariables, NHEALTHVARS,     &
        &                             characteristicLength,             &
        &                             twoDimensionalElm,                &
        &                             threeDimensionalFI,               &
        &                             timeIncrement,                    &
        &                             force2DResponse                   )
                
    implicit none
    
    ! arguments
    integer, intent(in)                 ::  NSTR, NINDEXES, NHEALTHVARS
    type(PlyDefinition), intent(in)     ::  plyProperties
    real(RWP), intent(in)               ::  strainVector(NSTR)
    real(RWP), intent(inout)            ::  failureIndexes(NINDEXES),   &
        &                                   damageVariables(NINDEXES)
    real(RWP), intent(out)              ::  healthVariables(NHEALTHVARS)
    real(RWP), intent(in)               ::  characteristicLength
    logical, intent(in)                 ::  twoDimensionalElm,          &
        &                                   threeDimensionalFI
    real(RWP), intent(in), optional     ::  timeIncrement
    logical, intent(in), optional       ::  force2DResponse
    
    ! locals
    real(RWP)                           ::  elasticStress(NSTR),        &
        &                                   oldFailureIndexes(NINDEXES)
    real(RWP)                           ::  elastic2DStress(3)
    logical                             ::  computeMI = .False.,        &
        &                                   force2D   = .False.
    
    real(RWP)                           ::  dt = -ONE
    integer                             ::  ii
    
    ! debug
    if (DEBUG_MODE) print 19130, 'updateDamageVariables'
        
    ! initialize flags
    computeMI = .NOT. twoDimensionalElm

    ! initialize variables
    healthVariables(:) = ONE
    oldFailureIndexes(1:NINDEXES) = failureIndexes(1:NINDEXES)
       
    ! local initializations
    if (PRESENT(timeIncrement)) then
        dt = timeIncrement
    end if
    if (PRESENT(force2DResponse)) then
        force2D = force2DResponse
    end if
        
    ! Step 1 -> compute elastic stress
    if (twoDimensionalElm) then
        elasticStress = computeStressPlane( plyProperties,              &
            &                               strainVector, NSTR          )
    else
        elasticStress = computeStresses3D( strainVector, plyProperties)
    end if
    
    if (DEBUG_MODE) print 19132, 'updateDamageVariables', 'EE', strainVector
    if (DEBUG_MODE) print 19132, 'updateDamageVariables', 'SE', elasticStress
    
    ! Step 2 -> compute elastic FIs
    failureIndexes(:) = ZERO ! reset these
    if (twoDimensionalElm .or. force2D) then
        elastic2DStress(1) = elasticStress(1)
        elastic2DStress(2) = elasticStress(2)
        if (NSTR .EQ. 6) then
            elastic2DStress(3) = elasticStress(4)
        else
            elastic2DStress(3) = elasticStress(3)
        end if
        failureIndexes(1:4) = plyCriteria2D( plyProperties, elastic2DStress,   &
            &                           failureIndexes(1:4),                   &
            &                           limitFIDen = .True.               )
        
    else
        failureIndexes = completeCriteria( plyProperties, elasticStress,    &
            &                              failureIndexes,                  &
            &                              NSTR, NINDEXES,  &
            &                              threeDimensionalFI, computeMI,   &
            &                              limitFIDen = .True.              )
    end if    
    if (DEBUG_MODE) print 19132, 'updateDamageVariables', 'FI', failureIndexes
        
    ! Step 2.1 -> kink and split
    failureIndexes = cumulativeFailureIndexes( failureIndexes,          &
        &   oldFailureIndexes, NINDEXES)
        

    if (force2D .and. NINDEXES .GT. 4) then
        failureIndexes(5) = ZERO
    end if

    ! Step 3 -> compute damage variables
    healthVariables = fHealthVariables(plyProperties,                   &
        &                        strainVector, elasticStress, NSTR,     &
        &                        failureIndexes, damageVariables,       & 
        &                        NINDEXES, NHEALTHVARS,                 &
        &                        characteristicLength,                  &
        &                        (twoDimensionalElm .or. force2D),      &
        &                        dt                                     )

    if (force2D .and. NINDEXES .GT. 4) then
        failureIndexes(5) = ncfCriteria( plyProperties, elasticStress )
    end if

    if (DEBUG_MODE) print 19132, 'updateDamageVariables', 'h ', healthVariables

    ! update failure indexes
    do ii=1,NINDEXES
        failureIndexes(ii) = MAX( failureIndexes(ii), oldFailureIndexes(ii) )
    end do

19130 format("(",A,") ",A)
19132 format("(",A,") ",A,": ", 10(F12.4))

    ! end
    end subroutine updateDamageVariables
        
        
    !_computeHealthVariables____________________________________________
    ! computes all damage variables di, and returns the HEALTH variables 
    ! considers a 3D tensor
    !   hi = (1-di) for i=1,2,s
    !___________________________________________________________________
    function computeHealthVariables(plyProperties,                      &
        &                           strainVector, stressVector, NSTR,   &
        &                           failureIndexes, damageVariables,    & 
        &                           oldFailureIndexes,                  &
        &                           NINDEXES, NVARIABLES,               &
        &                           characteristicLength,               &
        &                           twoDimensional,                     &
        &                           timeIncrement                       )
    !-------------------------------------------------------------------!
    !   plyProperties           plyDefinition object
    !   strainVector            ply 3D strain vector
    !   stressVector            ply 3D elastic stress vector
    !   failureIndexes          instantaneous elastic FIs
    !   damageVariables         current damage variables dm,ds,df,dk
    !   oldFailureIdx           history dependent old FIs
    !   characteristicLength    element characteristic length   
    !-------------------------------------------------------------------!
    
    implicit none
    
    ! arguments
    type(PlyDefinition), intent(in) ::  plyProperties   
    ! in
    integer, intent(in)             ::  NSTR, NINDEXES, NVARIABLES
    real(RWP), intent(in)           ::  strainVector(NSTR),             &
        &                               stressVector(NSTR),             &
        &                               oldFailureIndexes(NINDEXES),    &
        &                               characteristicLength
    logical, intent(in)             ::  twoDimensional 
    ! inout
    real(RWP), intent(inout)        ::  damageVariables(NINDEXES)
    real(RWP), intent(inout)        ::  failureIndexes(NINDEXES)
    real(RWP), intent(in)           ::  timeIncrement
    
    ! out
    real(RWP)   ::  computeHealthVariables(NVARIABLES)
    
    ! local variables
    real(RWP)   ::  matrixElasticEnergy, fibreElasticEnergy,            &
        &           matrixElasticEnergy3 = ZERO
    real(RWP)   ::  aux1, aux2, dMatrix = ZERO, dFibre = ZERO
    real(RWP)   ::  newFailureIndexes(NINDEXES)

    real(RWP)   ::  currentLambda, currentDelta,                        &
        &           currentDamage(NINDEXES),                            &
        &           currentFailureEnergy,                               &
        &           onsetEnergy,                                        &
        &           ultimateIndex = ZERO, currentFI = ZERO
    integer     ::  imode
    integer     ::  nS12 = 4

    ! flag
    logical     ::  includeMI = .False.
    logical     ::  kinkAndSplit = .False.
    logical     ::  matrixCompression = .False.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    
    ! Initializations
    newFailureIndexes(:) = ZERO
    newFailureIndexes(1:NINDEXES) = failureIndexes(1:NINDEXES)
    computeHealthVariables(:) = ONE
    currentDamage(:) = damageVariables(:)
    
    ! Matrix compression
    matrixCompression = (stressVector(2) .LT. -NZSTRESS)
    
    ! Debug
    if (DEBUG_MODE) print 19130, 'computeHealthVariables', 'v1'
    if (DEBUG_MODE .AND. twoDimensional .AND. matrixCompression) &
        &   print 19130, 'computeHealthVariables', '2D mtx compression'
    if (DEBUG_MODE) print 19131, 'computeHealthVariables', 'Vector size', NSTR
    if (DEBUG_MODE) print 19131, 'computeHealthVariables', 'Size', NVARIABLES
    
    ! Position of in plane shear
    if (NSTR .EQ. 6) then
        nS12 = 4
    else
        nS12 = 3
    end if
    
    ! Include matrix interface
    if ((NVARIABLES .GT. 3) .AND. (NINDEXES .GT. 4) .AND. (NSTR .EQ. 6)) then
        includeMI = .True.
    end if  
      
    
    ! Quick skip if no failure
    if (MAXVAL(failureIndexes) .GT. ONE) then
    
    ! get elastic energy that can be disspipated by the mode
    fibreElasticEnergy    = HALF * DOT_PRODUCT(strainVector,stressVector)
    matrixElasticEnergy   = HALF * mabrackets(strainVector(2)) * mabrackets(stressVector(2)) +  &
        &                   HALF * strainVector(nS12) * stressVector(nS12)
    ! the following two are initialized as ZERO
    if (includeMI) then
    matrixElasticEnergy3 = HALF * mabrackets(strainVector(3)) * mabrackets(stressVector(3)) +  &
        &                  HALF * strainVector(5) * stressVector(5) +                      &
        &                  HALF * strainVector(6) * stressVector(6)
    end if
    ! mod for matrix compressive elastic energy
    if (twoDimensional .AND. matrixCompression) then
    matrixElasticEnergy = HALF * mabrackets(-strainVector(2)) * mabrackets(-stressVector(2)) +  &
        &                 HALF * strainVector(nS12) * stressVector(nS12)
    end if

    ! Energy debug
    if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'Ufib/V  ', fibreElasticEnergy
    if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'Umtx/V  ', matrixElasticEnergy
    if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'Umtx2/V ', matrixElasticEnergy3

    ! Combined mode kink and split    
    aux1 = MAX(failureIndexes(2), oldFailureIndexes(2))
    aux2 = MAX(failureIndexes(4), oldFailureIndexes(4))
    if ((aux1 .GT. ONE) .AND. (aux2 .GT. ONE)) then
        kinkAndSplit = .True.
        print *,'(computeHealthVariables) k&S  :', kinkAndSplit, 'FI2:',aux1,', FI4', aux2
    end if    

    ! cycle through all failure modes
    do imode=1,NINDEXES
        ! get failure index
        currentFI = failureIndexes(imode)
        ! continue only if FI > 1
        if (DEBUG_MODE) print 19131, 'computeHealthVariables', 'Mode', imode
        if (currentFI .LT. ONE) then
            if (DEBUG_MODE) print 19130, 'computeHealthVariables', ' -> skipped'
            cycle
        else
            ! Fix failure indexes
            ! this is the correction implemented to account for 
            !   mode transition = between 2 and 4, which have the same 
            !   equation, meaning the modes
            if ((imode==2 .OR. imode==4) .AND. (kinkAndSplit)) then
                ! if it is now activated, theres a +1 discontinuity
                aux1 = oldFailureIndexes(2) + oldFailureIndexes(4) - ONE
                if (oldFailureIndexes(imode) .LT. ONE) then
                    aux1 = aux1 + ONE
                    currentFI = ONE
                else
                    currentFI = oldFailureIndexes(imode)
                end if
                aux2 = MAX(ZERO,failureIndexes(imode) - aux1)
                print 19130, 'computeHealthVariables', '>>mixed mode is active<<'
                if (DEBUG_MODE) then
                    print 19130, 'computeHealthVariables', '>>mixed mode is active<<'
                    print 19132, 'computeHealthVariables', 'FI_old  ', currentFI
                    print 19132, 'computeHealthVariables', 'FI_comb ', aux1
                    print 19132, 'computeHealthVariables', 'FI_inc  ', aux2
                end if
                ! account for mixed mode
                currentFI = currentFI + aux2
                newFailureIndexes(imode) = currentFI
            end if        
            ! get lambda value
            currentLambda = currentFI * currentFI
            ! mode elastic energy and failure energy
            if (imode == 1) then
                onsetEnergy = matrixElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G1cMat
                if (matrixElasticEnergy .GT. NEAR_ZERO_POS) then
                    ! new version
                    aux1 = HALF * mabrackets(strainVector(2) * stressVector(2))
                    aux2 = HALF * strainVector(nS12) * stressVector(nS12)
                    currentFailureEnergy = combinedMatrixG( plyProperties, &
                        &   aux1, aux2, matrixElasticEnergy,               &
                        &   matrixCompression )
                end if
            else if (imode == 2) then
                onsetEnergy = matrixElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G2cMat
            else if (imode == 3) then
                onsetEnergy = fibreElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G1cFibT
            else if (imode == 4) then
                onsetEnergy = fibreElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G1cFibK
            else if (imode == 5) then
                onsetEnergy = matrixElasticEnergy3 / currentLambda
                currentFailureEnergy = plyProperties%G1cMat
                if (matrixElasticEnergy3 .GT. NEAR_ZERO_POS) then
                    aux1 = HALF * mabrackets(strainVector(3)) * mabrackets(stressVector(3))
                    aux2 = HALF * strainVector(5) * stressVector(5) + &
                        &  HALF * strainVector(6) * stressVector(6) 
                    currentFailureEnergy = combinedMatrixG( plyProperties, &
                        &   aux1, aux2, matrixElasticEnergy3, .False.)
                end if
            else
                onsetEnergy = ZERO
            end if
  
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'Onset/V ', onsetEnergy
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'Uf/V    ', currentFailureEnergy
            
            ! damage driving parameters
            if (currentFailureEnergy .LT. NEAR_ZERO_POS) then
                if (DEBUG_MODE) print 19130, 'computeHealthVariables', ' -> no failure energy'
                currentDelta = ZERO
                cycle
            end if
            
            if (onsetEnergy .LT. NEAR_ZERO_POS) then 
                if (DEBUG_MODE) print 19130, 'computeHealthVariables', ' -> near zero onset'
                currentDelta  = ZERO
                cycle
            end if

            ultimateIndex = (currentFailureEnergy / onsetEnergy) / characteristicLength
            ultimateIndex = MAX( ultimateIndex, MAX_LAMBDA_F ) ! :TODO consider sqrt(1.2)
            currentDelta  = (ultimateIndex - currentFI) / (ultimateIndex - ONE)
            currentDelta  = ONE - (currentDelta/currentFI)

            
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'Le      ', characteristicLength
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'fi_f    ', ultimateIndex
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'fi      ', failureIndexes(imode)
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'fi*     ', currentFI
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'delta   ', currentDelta
        
            ! damage variable
            currentDamage(imode) = MAX( currentDelta, damageVariables(imode) )
            if (DEBUG_MODE) print 19132, 'computeHealthVariables', 'dmg     ', currentDamage(imode)
            
        end if
        
        ! limit damage 
        if (DAMAGE_LIMIT .GT. ZERO) then
            currentDamage(imode) = MIN( currentDamage(imode), DAMAGE_LIMIT )
        end if
        
    end do

    ! end of skip if no failure index is CURRENTLY > 1
    ! new modification to save unnecessary calculations 
    !   11.12.2019
    else
        if (DEBUG_MODE) print 19130, 'computeHealthVariables', &
            &   '>>>>>>>>>>> skipping all >>>>>>>>>>>'
        ! damage does not evolve 
        ! currentDamage(:) is initialized as damageVariables(:)
    end if


    ! insert regularization here !
    ! >                        < !
    ! __________________________ !
    if ((timeIncrement .GT. ZERO) .AND. (REG_VISCOSITY .GT. ZERO)) then
        if (DEBUG_MODE) print 19130, 'computeHealthVariables', &
            &   'Applying viscous regularization'
        currentDamage = viscousRegularization(currentDamage,            &
            &                                 damageVariables, NINDEXES,&
            &                                 timeIncrement,            &
            &                                 REG_VISCOSITY             )
        !------------- new as of 20.12.2019
    end if

    ! copy back
    damageVariables(:) = currentDamage(:)

    ! update the new failure indexes
    failureIndexes = max_vector(oldFailureIndexes,newFailureIndexes)
    
    if (kinkAndSplit) then
        print *,'(computeHealthVariables) k&S  :', kinkAndSplit
        print 19132, 'computeHealthVariables', 'FI_old', oldFailureIndexes
        print 19132, 'computeHealthVariables', 'FI_new', newFailureIndexes
        print 19132, 'computeHealthVariables', 'FI    ', failureIndexes
        print 19132, 'computeHealthVariables', 'S     ', stressVector
        print 19132, 'computeHealthVariables', 'E     ', strainVector
        print 19132, 'computeHealthVariables', 'dmg   ', damageVariables
    end if
    
    if (DEBUG_MODE) print 19133, 'computeHealthVariables', 'dmg     ',damageVariables
    
    
    ! damage variables
    dFibre  = ONE - (ONE - currentDamage(3)) * (ONE - currentDamage(4))
    dMatrix = MAX(currentDamage(1),currentDamage(2))
    ! 3D
    
    ! convert damage into material health (components)
    ! d1: fibre
    computeHealthVariables(1) = (ONE - dFibre)
    ! d2: matrix direction
    if (twoDimensional) then
        computeHealthVariables(2) = (ONE - dFibre) * (ONE - dMatrix)
    else
        computeHealthVariables(2) = (ONE - dFibre) * (ONE - &
            &   mabrackets(dMatrix,judge=(stressVector(2) .GT. -NZSTRESS)))
    end if
    ! dshear_12
    computeHealthVariables(nS12) = (ONE - dFibre) * (ONE - dMatrix)
    
    ! 3D components
    if ((.NOT. twoDimensional) .AND. includeMI) then
        dMatrix = MAX(dMatrix,currentDamage(5))
        !d3
        computeHealthVariables(3) = (ONE - dFibre) * (ONE -             &
            &   mabrackets(dMatrix,judge=(stressVector(3) .GT. -NZSTRESS)))
        !d23
        computeHealthVariables(5) = (ONE - dFibre) * (ONE - dMatrix)
        !d31
        computeHealthVariables(6) = computeHealthVariables(5)
    end if
        
   
    ! min health
    if (MIN_MATERIAL_HEALTH .GE. NEAR_ZERO_NEG) then 
        if (DEBUG_MODE) print 19133, 'computeHealthVariables', 'h min:  ', MIN_MATERIAL_HEALTH
        do imode=1,NVARIABLES
            computeHealthVariables(imode) = MAX(MIN_MATERIAL_HEALTH, computeHealthVariables(imode) )
        end do
    end if
    
19130 format("(",A,") ",A)
19131 format("(",A,") ",A,": ", 10(I0))     
19132 format("(",A,") ",A,": ", 10(F13.5))
19133 format("(",A,") ",A,": ", 10(E13.4))

    end function
    
    !_computeHealthVariables____________________________________________
    ! computes all damage variables di, and returns the HEALTH variables 
    ! considers a 3D tensor
    !   hi = (1-di) for i=1,2,s
    !___________________________________________________________________
    function fHealthVariables(plyProperties,                      &
        &                     strainVector, stressVector, NSTR,   &
        &                     failureIndexes, damageVariables,    & 
        &                     NINDEXES, NVARIABLES,               &
        &                     characteristicLength,               &
        &                     twoDimensional,                     &
        &                     timeIncrement                       )
    !-------------------------------------------------------------------!
    !   plyProperties           plyDefinition object
    !   strainVector            ply 3D strain vector
    !   stressVector            ply 3D elastic stress vector
    !   failureIndexes          instantaneous elastic FIs
    !   damageVariables         current damage variables dm,ds,df,dk
    !   oldFailureIdx           history dependent old FIs
    !   characteristicLength    element characteristic length   
    !-------------------------------------------------------------------!
    
    implicit none
    
    ! arguments
    type(PlyDefinition), intent(in) ::  plyProperties   
    ! in
    integer, intent(in)             ::  NSTR, NINDEXES, NVARIABLES
    real(RWP), intent(in)           ::  strainVector(NSTR),             &
        &                               stressVector(NSTR),             &
        &                               characteristicLength
    logical, intent(in)             ::  twoDimensional 
    ! inout
    real(RWP), intent(inout)        ::  damageVariables(NINDEXES)
    real(RWP), intent(inout)        ::  failureIndexes(NINDEXES)
    real(RWP), intent(in)           ::  timeIncrement
    
    ! out
    real(RWP)   ::  fHealthVariables(NVARIABLES)
    
    ! local variables
    real(RWP)   ::  matrixElasticEnergy, fibreElasticEnergy,            &
        &           matrixElasticEnergy3 = ZERO
    real(RWP)   ::  aux1, aux2, dMatrix = ZERO, dFibre = ZERO

    real(RWP)   ::  currentLambda, currentDelta,                        &
        &           currentDamage(NINDEXES),                            &
        &           currentFailureEnergy,                               &
        &           onsetEnergy, failureLambda,                         &
        &           failureEnergy(5)
    
    integer     ::  imode
    integer     ::  nS12 = 4

    ! flag
    logical     ::  includeMI = .False.
    logical     ::  matrixCompression = .False.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    
    ! Initializations
    fHealthVariables(:) = ONE
    currentDamage(:)    = damageVariables(:)
    
    ! Matrix compression
    matrixCompression = (stressVector(2) .LT. -NZSTRESS)
    dMatrix = ZERO
    dFibre = ZERO
       
    ! Position of in plane shear
    if (NSTR .EQ. 6) then
        nS12 = 4
    else
        nS12 = 3
    end if
    
    ! Include matrix interface
    if ((NVARIABLES .GT. 3) .AND. (NINDEXES .GT. 4) .AND. (NSTR .EQ. 6)) then
        includeMI = .True.
    end if  
      
    ! get elastic energy that can be disspipated by the mode
    fibreElasticEnergy    = HALF * DOT_PRODUCT(strainVector,stressVector)
    matrixElasticEnergy   = HALF * mabrackets(strainVector(2)) * mabrackets(stressVector(2)) +  &
        &                   HALF * strainVector(nS12) * stressVector(nS12)
    ! the following two are initialized as ZERO
    if (includeMI) then
    matrixElasticEnergy3 = HALF * mabrackets(strainVector(3)) * mabrackets(stressVector(3)) +  &
        &                  HALF * strainVector(5) * stressVector(5) +                      &
        &                  HALF * strainVector(6) * stressVector(6)
    end if
    ! mod for matrix compressive elastic energy
    if (twoDimensional .AND. matrixCompression) then
    matrixElasticEnergy = HALF * mabrackets(-strainVector(2)) * mabrackets(-stressVector(2)) +  &
        &                 HALF * strainVector(nS12) * stressVector(nS12)
    end if

    ! failure energies
    failureEnergy(1) = plyProperties%G1cMat
    failureEnergy(2) = plyProperties%G2cMat
    failureEnergy(3) = plyProperties%G1cFibT
    failureEnergy(4) = plyProperties%G1cFibK
    failureEnergy(5) = plyProperties%G1cMat
    
    ! cycle through all failure modes
    ! elastic energy at onset
    do imode=1,NINDEXES
        if (DEBUG_MODE) print '(A,I0)', '(fHealthVariables): Mode ',imode
        if (failureIndexes(imode) .GE. ONE) then
            currentLambda = failureIndexes(imode) * failureIndexes(imode)
            
            ! get onset energy
            if (imode .LT. 3) then  ! matrix
                onsetEnergy = matrixElasticEnergy / currentLambda
                if (DEBUG_MODE) print '(A," :",F12.4)', '(fHealthVariables): Umat/V  ', matrixElasticEnergy
            else if (imode .LT. 5) then ! fibre
                if (DEBUG_MODE) print '(A," :",F12.4)', '(fHealthVariables): Ufib/V  ', fibreElasticEnergy
                onsetEnergy = fibreElasticEnergy / currentLambda
            else ! matrix direction3 
                if (DEBUG_MODE) print '(A," :",F12.4)', '(fHealthVariables): Umt3/V  ', matrixElasticEnergy3
                onsetEnergy = matrixElasticEnergy3 / currentLambda
            end if            

            if (DEBUG_MODE) print '(A," :",F12.4)', '(fHealthVariables): Onset/V ', onsetEnergy
            !
            ! damage driving parameters
            currentFailureEnergy = failureEnergy(imode)
            if (imode .EQ. 1) then
                if (matrixElasticEnergy .GT. NEAR_ZERO_POS) then
                    ! new version
                    aux1 = HALF * mabrackets(strainVector(2) * stressVector(2))
                    aux2 = HALF * strainVector(nS12) * stressVector(nS12)
                    currentFailureEnergy = combinedMatrixG( plyProperties, &
                        &   aux1, aux2, matrixElasticEnergy,               &
                        &   matrixCompression )
                end if
            else if (imode .EQ. 5) then
                if (matrixElasticEnergy3 .GT. NEAR_ZERO_POS) then
                    ! new version
                    aux1 = HALF * mabrackets(strainVector(3)) * mabrackets(stressVector(3))
                    aux2 = HALF * strainVector(5) * stressVector(5)
                    aux2 = aux2 + HALF * strainVector(6) * stressVector(6) 
                    currentFailureEnergy = combinedMatrixG( plyProperties, &
                        &   aux1, aux2, matrixElasticEnergy3,              &
                        &   .False. )
                end if
            end if
            
            if (DEBUG_MODE) print '(A," :",F12.4)', '(fHealthVariables): Ufail/V ', currentFailureEnergy
            ! compute current delta
            if (onsetEnergy .LT. NEAR_ZERO_POS) then
                failureLambda = ZERO
                currentDelta  = ZERO
            else
                failureLambda = currentFailureEnergy / onsetEnergy
                failureLambda = failureLambda / characteristicLength
                failureLambda = MAX( failureLambda, MAX_LAMBDA_F)
                failureLambda = failureLambda*failureLambda
                ! calculate damage variables
                currentDelta = (SQRT(failureLambda)-failureIndexes(imode)) / (SQRT(failureLambda)-ONE)
                currentDelta = ONE - currentDelta/failureIndexes(imode)
            end if

            ! damage variable
            currentDamage(imode) = max(currentDelta, damageVariables(imode))
        
            if (DEBUG_MODE) then
                print '(A," :",E16.8)','(fHealthVariables): lambda_f', failureLambda
                print '(A," :",E16.8)','(fHealthVariables): Le      ', characteristicLength
                print '(A," :",E16.8)','(fHealthVariables): lambda  ', currentLambda
                print '(A," :",E16.8)','(fHealthVariables): delta   ', currentDelta
                print '(A," :",E16.8)','(fHealthVariables): damageV ', currentDamage(imode)
            end if

        end if
        ! limit damage variables
        if (DAMAGE_LIMIT .GT. ZERO) then
            currentDamage(imode) = min(DAMAGE_LIMIT, currentDamage(imode))
        end if
    end do
    
    ! viscosity
    if ((timeIncrement .GT. ZERO) .AND. (plyProperties%viscosity .GT. ZERO)) then
        if (DEBUG_MODE) then
            print '(A," :",E16.8)','(fHealthVariables): regularization', plyProperties%viscosity
        end if
        !viscousRegularization(newDamage, oldViscousDamage, NDMG,   &
        !&                          timeIncrement, viscosity             )
        currentDamage = viscousRegularization(currentDamage,            &
            &                                 damageVariables,          &
            &                                 NINDEXES,                 &
            &                                 timeIncrement,            &
            &                                 plyProperties%viscosity   )
    end if

    damageVariables(1:NINDEXES) = currentDamage(1:NINDEXES)
    
    fHealthVariables(:) = ONE
    !dMatrix = (ONE - currentDamage(1))
    !dMatrix = dMatrix * (ONE - currentDamage(2))
    dMatrix = ONE - MAX(currentDamage(1),currentDamage(2))
    dFibre  = (ONE - currentDamage(3))
    dFibre  = dFibre * (ONE - currentDamage(4)) ! bug until 08/01
    fHealthVariables(1) = fHealthVariables(1) * dFibre
    
    fHealthVariables(nS12) = fHealthVariables(nS12) * dFibre
    fHealthVariables(nS12) = fHealthVariables(nS12) * dMatrix
    
    ! todo: fix dMatrix for compression in 3D
    
    fHealthVariables(2) = fHealthVariables(2) * dFibre
    fHealthVariables(2) = fHealthVariables(2) * dMatrix
    
    if (MIN_MATERIAL_HEALTH .GT. ZERO) then 
        fHealthVariables(1) = MAX(MIN_MATERIAL_HEALTH, fHealthVariables(1) )
        fHealthVariables(2) = MAX(MIN_MATERIAL_HEALTH, fHealthVariables(2) )
        fHealthVariables(nS12) = MAX(MIN_MATERIAL_HEALTH, fHealthVariables(nS12) )
    end if
    
    ! 3D components
    if ((.NOT. twoDimensional) .AND. includeMI) then
        dMatrix = MIN(dMatrix,ONE-currentDamage(5)) ! inverse
        fHealthVariables(5) = fHealthVariables(5)*dFibre
        fHealthVariables(5) = fHealthVariables(5)*dMatrix
        fHealthVariables(6) = fHealthVariables(6)*dFibre
        fHealthVariables(6) = fHealthVariables(6)*dMatrix
        if (stressVector(3) .LE. NZSTRESS) then
            dMatrix = ONE   ! inverse
        end if

        fHealthVariables(3) = fHealthVariables(3)*dFibre
        fHealthVariables(3) = fHealthVariables(3)*dMatrix
        
        if (MIN_MATERIAL_HEALTH .GT. ZERO) then 
            fHealthVariables(3) = MAX(MIN_MATERIAL_HEALTH, fHealthVariables(3) )
            fHealthVariables(5) = MAX(MIN_MATERIAL_HEALTH, fHealthVariables(5) )
            fHealthVariables(6) = MAX(MIN_MATERIAL_HEALTH, fHealthVariables(6) )
        end if
    end if
    
    end function
    
    
    !_abortWithError____________________________________________________
    ! error handling subroutine
    !___________________________________________________________________
    subroutine abortWithError(errorMessage, errorNumber, condition)
    
        character(len=*), intent(in)    ::  errorMessage
        integer, intent(in), optional   ::  errorNumber
        logical, intent(in), optional   ::  condition
        
        logical                         ::  checks = .True.
        integer                         ::  errorCode = -1
    
        if (present(errorNumber))   errorCode = errorNumber
        if (present(condition))     checks = condition
        
        if (checks) then
            print *,'ERROR! ', trim(errorMessage)
            call exit(errorNumber)
        end if
        
    end subroutine
        
        
        
    function viscousRegularization(newDamage, oldViscousDamage, NDMG,   &
        &                          timeIncrement, viscosity             )
        
        integer, intent(in)     ::  NDMG
        real(RWP), intent(in)   ::  newDamage(NDMG),                    &
            &                       oldViscousDamage(NDMG)
        real(RWP), intent(in)   ::  timeIncrement, viscosity
        
        ! returns
        real(RWP)             ::    viscousRegularization(NDMG)
        
        ! locals
        real(RWP)   ::  aux1, aux2
        
        aux1 = timeIncrement / (viscosity + timeIncrement)
        aux2 = viscosity / (viscosity + timeIncrement)
                      
        ! visc reg
        viscousRegularization(:) = ZERO
        ! vector operation
        viscousRegularization = aux1 * newDamage + aux2 * oldViscousDamage
                
    end function viscousRegularization
     
     
    !_computeHealthVariables____________________________________________
    ! computes all damage variables di, and returns the HEALTH variables 
    ! considers a 3D tensor
    !   hi = (1-di) for i=1,2,s
    !___________________________________________________________________
    function computeHealth(plyProperties,                        &
        &                  strainVector, stressVector, NSTR,     &   
        &                  failureIndexes, damageVariables,      & 
        &                  oldFailureIndexes,                    &
        &                  NINDEXES, NVARIABLES,                 &
        &                  characteristicLength,                 &
        &                  twoDimensional,                       &
        &                  timeIncrement                         )
    !-------------------------------------------------------------------!
    !   plyProperties           plyDefinition object
    !   strainVector            ply 3D strain vector
    !   stressVector            ply 3D elastic stress vector
    !   failureIndexes          instantaneous elastic FIs
    !   damageVariables         current damage variables dm,ds,df,dk
    !   oldFailureIdx           history dependent old FIs
    !   characteristicLength    element characteristic length   
    !-------------------------------------------------------------------!
    
    implicit none
    
    ! arguments
    type(PlyDefinition), intent(in) ::  plyProperties   
    ! in
    integer, intent(in)             ::  NSTR, NINDEXES, NVARIABLES
    real(RWP), intent(in)           ::  strainVector(NSTR),             &
        &                               stressVector(NSTR),             &
        &                               oldFailureIndexes(NINDEXES),    &
        &                               characteristicLength
    logical, intent(in)             ::  twoDimensional 
    ! inout
    real(RWP), intent(inout)        ::  damageVariables(NINDEXES)
    real(RWP), intent(inout)        ::  failureIndexes(NINDEXES)
    real(RWP), intent(in)           ::  timeIncrement
    
    ! out
    real(RWP)   ::  computeHealth(NVARIABLES)
    
    ! local variables
    real(RWP)   ::  matrixElasticEnergy, fibreElasticEnergy,            &
        &           matrixElasticEnergy3 = ZERO
    real(RWP)   ::  aux1, aux2, dMatrix = ZERO, dFibre = ZERO
    real(RWP)   ::  newFailureIndexes(NINDEXES)

    real(RWP)   ::  currentLambda, currentDelta,                        &
        &           currentDamage(NINDEXES),                            &
        &           currentFailureEnergy,                               &
        &           onsetEnergy,                                        &
        &           ultimateIndex = ZERO, currentFI = ZERO
    integer     ::  imode
    integer     ::  nS12 = 4

    ! flag
    logical     ::  includeMI = .False.
    logical     ::  kinkAndSplit = .False.
    logical     ::  matrixCompression = .False.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    
    ! Initializations
    newFailureIndexes(:) = ZERO
    newFailureIndexes(1:NINDEXES) = failureIndexes(1:NINDEXES)
    computeHealth(:) = ONE
    currentDamage(:) = damageVariables(:)
    
    ! Matrix compression
    matrixCompression = (stressVector(2) .LT. -NZSTRESS)
    matrixCompression = .False.
    
    ! Debug
    if (DEBUG_MODE) print 19130, 'computeHealth', 'v1'
    if (DEBUG_MODE .AND. twoDimensional .AND. matrixCompression) &
        &   print 19130, 'computeHealth', '2D mtx compression'
    if (DEBUG_MODE) print 19131, 'computeHealth', 'Vector size', NSTR
    if (DEBUG_MODE) print 19131, 'computeHealth', 'Size', NVARIABLES
    
    ! Position of in plane shear
    if (NSTR .EQ. 6) then
        nS12 = 4
    else
        nS12 = 3
    end if
    
    ! Include matrix interface
    if ((NVARIABLES .GT. 3) .AND. (NINDEXES .GT. 4) .AND. (NSTR .EQ. 6)) then
        includeMI = .True.
    end if    

    ! Elastic energy that can be disspipated by the mode
    fibreElasticEnergy    = HALF * DOT_PRODUCT(strainVector,stressVector)
    matrixElasticEnergy   = HALF * mabrackets(strainVector(2) * stressVector(2))
    matrixElasticEnergy   = matrixElasticEnergy + &
        &                   HALF * strainVector(nS12) * stressVector(nS12)
        
    ! Energy in the 3rd direction
    if (includeMI) then
    matrixElasticEnergy3 = HALF * mabrackets(strainVector(3)) * mabrackets(stressVector(3)) +  &
        &                  HALF * strainVector(5) * stressVector(5) +                      &
        &                  HALF * strainVector(6) * stressVector(6)
    end if

    ! Energy debug
    if (DEBUG_MODE) print 19132, 'computeHealth', 'Ufib/V  ', fibreElasticEnergy
    if (DEBUG_MODE) print 19132, 'computeHealth', 'Umtx/V  ', matrixElasticEnergy
    if (DEBUG_MODE) print 19132, 'computeHealth', 'Umtx2/V ', matrixElasticEnergy3

    ! Combined mode kink and split    
    aux1 = MAX(failureIndexes(2), oldFailureIndexes(2))
    aux2 = MAX(failureIndexes(4), oldFailureIndexes(4))
    if ((aux1 .GT. ONE) .AND. (aux2 .GT. ONE)) then
        kinkAndSplit = .True.
        print *,'(computeHealth) k&S  :', kinkAndSplit, 'FI2:',aux1,', FI4', aux2
    end if  
    kinkAndSplit = .False.  

    ! cycle through all failure modes
    do imode=1,NINDEXES
        ! get failure index
        currentFI = failureIndexes(imode)
        ! continue only if FI > 1
        if (DEBUG_MODE) print 19131, 'computeHealth', 'Mode', imode
        if (currentFI .GT. ONE) then
            ! fix failure indexes
            ! this is the correction implemented to account for 
            !   mode transition = between 2 and 4, which have the same 
            !   equation, meaning the modes
            if ((imode==2 .OR. imode==4) .AND. (kinkAndSplit)) then
                ! if it is now activated, theres a +1 discontinuity
                if (oldFailureIndexes(imode) .LT. ONE) then
                    aux1 = oldFailureIndexes(2) + oldFailureIndexes(4)
                    currentFI = ONE
                else
                    aux1 = oldFailureIndexes(2) + oldFailureIndexes(4) - ONE
                    currentFI = oldFailureIndexes(imode)
                end if
                aux2 = MAX(ZERO,failureIndexes(imode) - aux1)
                print 19130, 'computeHealth', '>> mixed mode is active <<'
                if (DEBUG_MODE) then
                    print 19130, 'computeHealth', '>>mixed mode is active<<'
                    print 19132, 'computeHealth', 'FI_old  ', currentFI
                    print 19132, 'computeHealth', 'FI_comb ', aux1
                    print 19132, 'computeHealth', 'FI_inc  ', aux2
                end if
                ! account for mixed mode
                currentFI = currentFI + aux2
                newFailureIndexes(imode) = currentFI
            end if        
            ! get lambda value
            currentLambda = currentFI * currentFI
            ! get mode energy and onset
            select case (imode)
            case (1)
                onsetEnergy = matrixElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G1cMat
                if (matrixElasticEnergy .GT. NEAR_ZERO_POS) then
                    ! new version
                    aux1 = HALF * mabrackets(strainVector(2) * stressVector(2))
                    aux2 = HALF * strainVector(nS12) * stressVector(nS12)
                    currentFailureEnergy = combinedMatrixG( plyProperties, &
                        &   aux1, aux2, matrixElasticEnergy,               &
                        &   matrixCompression )
                end if
            case (2)
                onsetEnergy = matrixElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G2cMat
            case (3)
                onsetEnergy = fibreElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G1cFibT
            case (4)
                onsetEnergy = fibreElasticEnergy / currentLambda
                currentFailureEnergy = plyProperties%G1cFibK
            case (5) 
                onsetEnergy = matrixElasticEnergy3 / currentLambda
                currentFailureEnergy = plyProperties%G1cMat
                if (matrixElasticEnergy3 .GT. NEAR_ZERO_POS) then
                    aux1 = HALF * mabrackets(strainVector(3)) * mabrackets(stressVector(3))
                    aux2 = HALF * strainVector(5) * stressVector(5) + &
                        &  HALF * strainVector(6) * stressVector(6) 
                    currentFailureEnergy = combinedMatrixG( plyProperties, &
                        &   aux1, aux2, matrixElasticEnergy3, .False.)
                end if
            case default
                onsetEnergy = ZERO
                currentFailureEnergy = ZERO
            end select

            if (DEBUG_MODE) print 19132, 'computeHealth', 'Onset/V ', onsetEnergy
            if (DEBUG_MODE) print 19132, 'computeHealth', 'Uf/V    ', currentFailureEnergy
            
            ! checks
            if ((currentFailureEnergy .LT. NEAR_ZERO_POS) .or. &
                & (onsetEnergy .LT. NEAR_ZERO_POS)) then
                if (DEBUG_MODE) print 19130, 'computeHealth', ' -> zero failure or onset'
                currentDelta = ZERO
            else
                ultimateIndex = (currentFailureEnergy / onsetEnergy)
                ultimateIndex = ultimateIndex / characteristicLength
                ultimateIndex = MAX( ultimateIndex, SQRT(MAX_LAMBDA_F) )
                currentDelta  = (ultimateIndex - currentFI) / (ultimateIndex - ONE)
                currentDelta  = ONE - (currentDelta/currentFI)
            end if
            
            ! damage variable
            currentDamage(imode) = MAX( currentDelta, damageVariables(imode) )
            
            ! limit damage 
            if (DAMAGE_LIMIT .GT. ZERO) then
                currentDamage(imode) = MIN( currentDamage(imode), DAMAGE_LIMIT )
            end if
            
            ! debug
            if (DEBUG_MODE) then
                print 19132, 'computeHealth', 'Le      ', characteristicLength
                print 19132, 'computeHealth', 'fi_f    ', ultimateIndex
                print 19132, 'computeHealth', 'fi      ', failureIndexes(imode)
                print 19132, 'computeHealth', 'fi*     ', currentFI
                print 19132, 'computeHealth', 'delta   ', currentDelta
                print 19132, 'computeHealth', 'dmg     ', currentDamage(imode)
            end if
        end if
    end do

    ! insert regularization here !
    ! >                        < !
    ! __________________________ !
    if ((timeIncrement .GT. ZERO) .AND. (REG_VISCOSITY .GT. ZERO)) then
        if (DEBUG_MODE) print 19130, 'computeHealth', &
            &   'Applying viscous regularization'
        currentDamage = viscousRegularization(currentDamage,            &
            &                                 damageVariables, NINDEXES,&
            &                                 timeIncrement,            &
            &                                 REG_VISCOSITY             )
        !------------- new as of 20.12.2019
    end if

    ! copy back
    damageVariables(:) = currentDamage(:)

    ! update the new failure indexes
    failureIndexes = max_vector(oldFailureIndexes,newFailureIndexes)

    ! health variables
    computeHealth(:) = ONE
    dMatrix = (ONE - currentDamage(1))
    dMatrix = dMatrix * (ONE - currentDamage(2))
    dFibre  = (ONE - currentDamage(3))
    dFibre  = dFibre * (ONE - currentDamage(4))
    ! direction 1
    computeHealth(1) = computeHealth(1) * dFibre
    ! direction 2
    computeHealth(2) = computeHealth(2) * dFibre
    computeHealth(2) = computeHealth(2) * dMatrix
    ! shear 12
    computeHealth(3) = computeHealth(3) * dFibre
    computeHealth(3) = computeHealth(3) * dMatrix
    
    ! min health
    if (MIN_MATERIAL_HEALTH .GT. NEAR_ZERO_NEG) then 
        do imode=1,NVARIABLES
            computeHealth(imode) = &
                &   MAX(MIN_MATERIAL_HEALTH, computeHealth(imode))
        end do
    end if
    
19130 format("(",A,") ",A)
19131 format("(",A,") ",A,": ", 10(I0))     
19132 format("(",A,") ",A,": ", 10(F13.5))

    end function
        
end module    
