!DIR$ FREEFORM
!  plyTools.f90 
!
!  OBJECTS:
!  PlyDefinition - Defines a ply class and its properties
! 
!  SUBROUTINES:
!  initializePly    - initializes the object from PROPS(NPROPS)
!  plyPrint         - prints the ply definition and all of its properties
!   
!  FUNCTIONS:
!  plyCriticalAlignmnent    - computes phi_c
!  plyEvaluateCriteria      - evaluates the common criterion
!
!****************************************************************************
!****************************************************************************
!                     developed by Miguel A.S. Matos and Silvestre T. Pinho
!                                      Imperial College London, Aeronautics
!                                                        created 17.04.2019
!****************************************************************************
!****************************************************************************
!   log:
!       04.06.2019  Added failure energies
!       10.06.2019  Redefined failure energies
!       03.07.2019  Final checks
!       20.11.2019  Added 3D
!       03.12.2019  Completed 3D full support
!       04.12.2019  Added DEBUG_MODE as a public variable
!   last error free build on:
!       03.07.2019
!****************************************************************************
    
    
module plyTools

    implicit none
    private
    
    public  ::  PlyDefinition,              &
        &       initializePly,              &
        &       plyCriticalAlignment,       &
        &       plyEvaluateCriteria,        &
        &       plyPrint,                   &
        &       plyFromFile,                &
        &       plyAdditional3DProperties,  &
        &       combinedMatrixG,            &
        &       getE11          
    
    logical, parameter, public  ::  DEBUG_MODE = .False.
    
    !types (for abaqus)
    integer,   parameter    ::  RWP = KIND(1.0D0)
    real(RWP), parameter    ::  ZERO = 0.0000D0, ONE = 1.0000D0,        &
        &                       TWO  = 2.0000D0,                        &
        &                       SQ2 = SQRT(2.000D0), HALF = 0.5000D0,   &
        &                       ZERONEG = -1.0000D-9, ZEROPOS = 1.0000D-9
    real(RWP), parameter    ::  PI = ACOS(-1.0D0)
    real(RWP), parameter    ::  DEGTORAD = PI/180.0D0
    
    integer, parameter      ::  PLY_TYPE_UD = 1,                &
        &                       PLY_TYPE_THICK_EMBEDDED = 2,    &
        &                       PLY_TYPE_THIN_EMBEDDED  = 3,    &
        &                       PLY_TYPE_THIN_OUTER     = 4,    &
        &                       PLY_TYPE_EMBEDDED_AUTO  = 5
    
    real(RWP), parameter    ::  DEFAULT_ALPHA_0 = 53.00D0,      &
        &                       DEFAULT_ALPHA_G = 1.21
    
    
    ! local definitions (relates properties with input indexes)
    integer,parameter ::    IPROP_E11     =  1, &   ! fiber direction modulus
        &                   IPROP_E22     =  2, &   ! transverse modulus
        &                   IPROP_nu12    =  3, &   ! Poisson's ratio
        &                   IPROP_G12     =  4, &   ! In plane shear modulus
        &                   IPROP_Xt      =  5, &   ! longitudinal tensile strength
        &                   IPROP_Xc      =  6, &   ! longitudinal compressive strength
        &                   IPROP_Yt      =  7, &   ! transverse tensile strength
        &                   IPROP_Yc      =  8, &   ! transverse compressive strength
        &                   IPROP_Sl      =  9, &   ! longitudinal shear strength
        &                   IPROP_a0      = 10, &   ! fracture angle plane
        &                   IPROP_nL      = 11, &   ! longitudinal shear friction coefficient
        &                   IPROP_nT      = 12, &   ! transverse shear friction coefficient
        &                   IPROP_ILSS    = 13, &   ! Interlaminar shear strength
        &                   IPROP_Zt      = 14, &   ! Tensile strength in the 3rd direction
        &                   IPROP_G1cMat  = 15, &   ! Matrix mode I critical energy release rate
        &                   IPROP_G2cMat  = 16, &   ! Matrix mode II critical energy release rate
        &                   IPROP_G1cFibT = 17, &   ! Fibre tensile critical eneregy release rate 
        &                   IPROP_G1cFibK = 18, &   ! Fibre kinking critical energy release rate
        &                   IPROP_GExpMM  = 19, &   ! matrix mixed-mode power law exponent
        &                   IPROP_E11comp = 20, &   !
        &                   IPROP_E33     = 21, &   !
        &                   IPROP_nu13    = 22, &   !
        &                   IPROP_nu23    = 23, &   !
        &                   IPROP_G13     = 24, &   !
        &                   IPROP_G23     = 25, &   !
        &                   IPROP_STAB    = 26, &   !
        &                   IPROP_PTYP    = 27, &   !
        &                   IPROP_TH      = 28, &   ! Ply 
        &                   IPROP_G1c     = 29, &   !
        &                   IPROP_G2c     = 30      !
        
    ! Order:
    !    1,    2,     3,   4,  5,  6,   7,   8,  9, 10, 11, 12,   13, 14,   
    !  E11,  E22,  nu12, G12, Xt, Xc,  Yt,  Yc, Sl, a0, nL, nT, ILSS, Zt, 
    !   15,   16,    17,     18,   19,    20, 21,  22,  23,  
    ! G1cM, G2cM, GfibT,  GfibK, GExp, PTYPE, TH, G1c, G2c
    
    integer, parameter  ::  DEFAULT_PTYP = PLY_TYPE_UD    ! default ply type (1 means UD)
    
    ! PlyDefinition object
    type PlyDefinition
        real(RWP)   ::  E11, E22,               &   ! elastic moduli
            &           nu12, nu21,             &   ! Poisson's ratio
            &           G12,                    &   ! shear moduli
            &           Xt, Xc, Yt, Yc,         &   ! strengths
            &           Sl, a0, nL, nT, St,     &   ! additional properties
            &           YtIS, SlIS, StIS,       &   ! in-situ effects
            &           phiC,                   &   ! 
            &           ILSS, Zt,               &   !
            &           Gvalues(5),             &   ! Failure energies
            &           G1cMat,G2cMat,          &   
            &           G1cFibT,G1cFibK,        &   
            &           GalphaM,                &
            &           E11comp,                &   ! Compression modulus
            &           E33,                    &   ! Through the thick
            &           nu23, nu13,             &   ! Poisson's
            &           G13, G23,               &   ! more shear compts
            &           viscosity                   ! viscous regularization
        
        integer     ::  PTYP                        ! ply type
        logical     ::  TENCOMP ! tension compression flag
        
    end type PlyDefinition
    
    ! content
    contains

    !_initializePly_____________________________________________________
    ! initialization of this class based on Abaqus PROPS vector
    !___________________________________________________________________
    subroutine initializePly(this, PROPS, NPROPS)
    
        ! input variables
        type(PlyDefinition), intent(inout)  ::  this
        integer, intent(in)                 ::  NPROPS
        real(RWP), intent(in)               ::  PROPS(NPROPS)
        
        ! locals
        real(RWP)   :: aux1, aux2, a0, AA, G1c, G2c, TH
                
        ! property definition
        this%E11  = PROPS(IPROP_E11)    !1
        this%E22  = PROPS(IPROP_E22)    !2
        this%nu12 = PROPS(IPROP_nu12)   !3 
        this%G12  = PROPS(IPROP_G12)    !4
        this%Xt   = PROPS(IPROP_Xt)     !5
        this%Xc   = PROPS(IPROP_Xc)     !6
        this%Yt   = PROPS(IPROP_Yt)     !7
        this%Yc   = PROPS(IPROP_Yc)     !8
        this%Sl   = PROPS(IPROP_Sl)     !9
        this%nL   = PROPS(IPROP_nL)     !10
        this%nT   = PROPS(IPROP_nT)     !10
        this%a0   = PROPS(IPROP_a0)     !11
        
        ! fix properties
        if (this%a0 .LT. ZERONEG) then
            this%a0 = DEFAULT_ALPHA_0
        end if
        a0 = this%a0 * DEGTORAD    
        if (this%nL .LT. ZERONEG) then
            aux1 = this%Sl * COS(2*a0)
            aux2 = this%Yc * COS(a0) * COS(a0)
            this%nL = - aux1/aux2
        end if
        if (this%nT .LT. ZERONEG) then
            this%nT = -ONE/TAN(2*a0)
        end if
        
        ! 3D NCF failure indexes
        if (NPROPS >= IPROP_Zt) then
            this%ILSS = PROPS(IPROP_ILSS)
            this%Zt   = PROPS(IPROP_Zt)
            if (this%Zt .LT. ZERO) then
                this%Zt = HALF * this%Yt
            end if
        else
            this%ILSS = -ONE
            this%Zt   = -ONE
        end if
            
        if (NPROPS >= IPROP_nu23) then
            this%nu23 = PROPS(IPROP_nu23)
        else
            this%nu23 = ZERO
        end if
        
        ! Failure energies
        this%Gvalues(:) = ZERO
        this%G1cMat = ZERO
        this%G2cMat = ZERO
        this%G1cFibT = ZERO
        this%G1cFibK = ZERO
        this%GalphaM = ZERO
        if (NPROPS >= IPROP_G1cFibK) then
            this%Gvalues(1) = PROPS(IPROP_G1cMat)
            this%Gvalues(2) = PROPS(IPROP_G2cMat)
            this%Gvalues(3) = PROPS(IPROP_G1cFibT)
            this%Gvalues(4) = PROPS(IPROP_G1cFibK)
            if (NPROPS >= IPROP_GExpMM) then
                this%Gvalues(5) = PROPS(IPROP_GExpMM)
                if (this%Gvalues(5) < ZERO) then
                    this%Gvalues(5) = DEFAULT_ALPHA_G
                end if
            else 
                this%Gvalues(5) = DEFAULT_ALPHA_G
            end if
            ! populate new dedicated variables
            this%G1cMat = PROPS(IPROP_G1cMat)
            this%G2cMat = PROPS(IPROP_G2cMat)
            this%G1cFibT = PROPS(IPROP_G1cFibT)
            this%G1cFibK = PROPS(IPROP_G1cFibK)
            this%GalphaM = PROPS(IPROP_GExpMM)
            if (NPROPS >= IPROP_GExpMM) then
                this%GalphaM = PROPS(IPROP_GExpMM)
                if (this%GalphaM < ZERO) then
                    this%GalphaM = DEFAULT_ALPHA_G
                end if
            else 
                this%GalphaM = DEFAULT_ALPHA_G
            end if
        end if
        
        ! Defaults - to be fixed later on
        this%E11comp = -ONE
        this%E33     = -ONE
        this%nu13    = -ONE
        this%nu23    = -ONE
        this%G13     = -ONE
        this%G23     = -ONE
        this%viscosity = -ONE
        
        ! Compression modulus
        if (NPROPS >= IPROP_E11comp) then
            this%E11comp = PROPS(IPROP_E11comp)
        end if
        if (this%E11comp .LT. ZERONEG) then
            this%E11comp = this%E11
        end if
        
        ! 3D properties
        ! E33
        if (NPROPS >= IPROP_E33) then
            this%E33 = PROPS(IPROP_E33)
        end if
        if (this%E33 .LT. ZERONEG) then
            this%E33 = this%E22
        end if
        ! nu13
        if (NPROPS >= IPROP_nu13) then
            this%nu13 = PROPS(IPROP_nu13)
        end if
        if (this%nu13 .LT. ZERONEG) then
            this%nu13 = this%nu12
        end if
        ! nu23        
        if (NPROPS >= IPROP_nu23) then
            this%nu23 = PROPS(IPROP_nu23)
        end if
        if (this%nu23 .LT. ZERONEG) then
            this%nu23 = HALF
        end if        
        ! G13
        if (NPROPS >= IPROP_G13) then
            this%G13 = PROPS(IPROP_G13)
        end if
        if (this%G13 .LT. ZERONEG) then
            this%G13 = this%G12
        end if
        ! G23
        if (NPROPS >= IPROP_G23) then
            this%G23 = PROPS(IPROP_G23)
        end if
        if (this%G23 .LT. ZERONEG) then
            this%G23 = this%E33/(TWO*(ONE+this%nu23))
        end if
        
        ! Viscous stabilization
        if (NPROPS >= IPROP_STAB) then
            this%viscosity = PROPS(IPROP_STAB)
        end if
        
        ! Ply type       
        if (NPROPS >= IPROP_PTYP) then
            this%PTYP = int(PROPS(IPROP_PTYP))
        else
            this%PTYP = DEFAULT_PTYP
        end if
                
        
        ! Shear strength        
        this%St = this%Yc * COS(a0) * (sin(a0) + (cos(a0)/tan(2*a0)) )
        
        ! poisson
        this%nu21 = this%nu12 * this%E22 / getE11(this)
        
        ! in situ effects
        select case (this%PTYP)
        case (1)    ! UD
            this%StIS = this%St
            this%YtIS = this%Yt
            this%SlIS = this%Sl
        case (2)
            this%YtIS = 1.1200D0 * SQ2 * this%Yt
            this%SlIS = SQ2 * this%Sl
            this%StIS = SQ2 * this%St
        case (3)
            G1c  = PROPS(IPROP_G1c)
            G2c  = PROPS(IPROP_G2c)
            TH   = PROPS(IPROP_TH)
            AA   = 2.00D0 * (ONE/this%E22 - this%nu21/this%E11)
            this%YtIS = SQRT( 8 * G1c / (PI * TH * AA) )
            this%SlIS = SQRT( 8 * G2c * this%G12 / (PI * TH) )
            this%StIS = this%SlIS
        case (4)
            G1c  = PROPS(IPROP_G1c)
            G2c  = PROPS(IPROP_G2c)
            TH   = PROPS(IPROP_TH)
            AA   = 2.00D0 * (ONE/this%E22 - this%nu21/this%E11)
            this%YtIS = SQRT( 4 * G1c / (PI * TH * AA) )
            this%SlIS = SQRT( 4 * G2c * this%G12 / (PI * TH) )
            this%StIS = this%SlIS * SQ2
        case default
            G1c  = PROPS(IPROP_G1c)
            G2c  = PROPS(IPROP_G2c)
            TH   = PROPS(IPROP_TH)
            AA   = 2.00D0 * (ONE/this%E22 - this%nu21/this%E11)
            aux1 = 1.1200D0 * SQ2 * this%Yt
            aux2 = SQRT( 8 * G1c / (PI * TH * AA) )
            this%YtIS = MAX( aux1, aux2 )
            aux1 = SQ2 * this%Sl
            aux2 = SQRT( 8 * G2c * this%G12 / (PI * TH) )
            this%SlIS = MAX( aux1, aux2 )
            aux1 = SQ2 * this%St
            this%StIS = MAX( aux1, aux2 )

        end select 
        
        this%phiC = plyCriticalAlignment(this)
            
        this%TENCOMP = .True.
            
    end subroutine

    
    !_plyCriticalAlignment______________________________________________
    ! calculates the ply material Phi_c
    !___________________________________________________________________
    function plyCriticalAlignment(this)
        ! input 
        type(PlyDefinition), intent(in)  ::  this
        ! output
        real(RWP)   :: plyCriticalAlignment
        ! locals
        real(RWP)   :: aux1,aux2
        ! calculate
        aux1 = this%Sl/this%Xc
        aux2 = 2.000D0*(aux1 + this%nL)
        aux1 = 1.00D0 - SQRT(1.00d0-2.00d0*aux2*aux1)
        plyCriticalAlignment = ATAN(aux1/aux2)  
    end function
    
    !_combinedMatrixG___________________________________________________
    ! returns Gmatrix for a combined uniaxial and shear loading
    !   see Ref: 
    !___________________________________________________________________
    function combinedMatrixG(this, uniaxialEnergy, shearEnergy,         &
        &                    totalEnergy, matrixCompression)
        ! input 
        type(PlyDefinition), intent(in)  ::  this
        real(RWP), intent(in)   ::  uniaxialEnergy, shearEnergy, totalEnergy
        logical, intent(in)     ::  matrixCompression
        
        ! output
        real(RWP)   ::  combinedMatrixG
        ! locals
        real(RWP)   ::  aux1, aux2, Guni
        
        Guni = this%G1cMat
        if (matrixCompression) then
            Guni = this%G2cMat / COS(this%a0*DEGTORAD)
        end if
        
        ! G 
        aux1 = uniaxialEnergy / totalEnergy
        aux2 = shearEnergy / totalEnergy
        ! fail safe
        if ((aux1 .LT. ZEROPOS) .AND. (aux2 .LT. ZEROPOS)) then
            combinedMatrixG = this%G1cMat
        else 
            ! issue raised on 06.02.2020
            ! the reference included a square which may ?? overpredict the 
            !   mixed mode quantities
            !aux1 = aux1*aux1    ! This will overpredict the mixed mode !
            !aux2 = aux2*aux2    ! This will overpredict the mixed mode !
            aux1 = (aux1 / Guni)
            aux2 = (aux2 / this%G2cMat)
            aux1 = aux1 ** this%GalphaM
            aux2 = aux2 ** this%GalphaM
            ! return
            combinedMatrixG = (aux1 + aux2) ** (-ONE/this%GalphaM)
        end if
    end function
    
    !_plyAdditional3DProperties_________________________________________
    ! returns all additional 3D required properties
    ! lamina assumptions
    !   [E33, v13, v31, v23, v32, G23, G13]
    !___________________________________________________________________
    function plyAdditional3DProperties(this)
        ! input
        type(PlyDefinition), intent(in) :: this
        ! output
        real(RWP)   ::  plyAdditional3DProperties(7)
        ! initialization
        plyAdditional3DProperties(:) = ZERO
        ! E33
        plyAdditional3DProperties(1) = this%E33
        ! v13
        plyAdditional3DProperties(2) = this%nu13
        ! v31 <- ok
        plyAdditional3DProperties(3) = this%nu13 * this%E33 / this%E11
        ! v23
        plyAdditional3DProperties(4) = this%nu23
        ! v32 <- ok
        plyAdditional3DProperties(5) = this%nu23 * this%E33 / this%E22
        ! G23
        plyAdditional3DProperties(6) = this%G23
        ! G13
        plyAdditional3DProperties(7) = this%G13
    end function        
    
    !_plyEvaluateCriteria____________________________________________________
    ! evaluates the LaRC05 failure criteria, used by modes indexes 1,2 and 4
    !________________________________________________________________________
    function plyEvaluateCriteria(this,SST,SSL,SN,limitSigmaN)
    
        implicit none
        
        ! input
        type(PlyDefinition), intent(in) ::  this
        real(RWP), intent(in)           ::  SST, SSL, SN
        logical, intent(in), optional   ::  limitSigmaN
        
        ! output
        real(RWP)   ::  plyEvaluateCriteria
        
        ! locals
        real(RWP)   ::  aux1, aux2, aux3, SNDEN
        
        ! added 18.06.2019
        ! modified 13.02.2020
        !   so that, the beneficial effect of normal stress
        !   is limited by its value when pure compressive failure
        !   happens Yc * (cos(ALPHA0)^2)
        SNDEN = SN
        if ( PRESENT(limitSigmaN) ) then
            if (limitSigmaN) then
                aux1 = COS(this%a0*DEGTORAD)
                aux1 = aux1*aux1
                SNDEN = MAX(-this%Yc*aux1,MIN(ZERO,SN))
            end if
        end if
        
        !            
        aux3 = max(SN/(this%YtIS),ZERO)
        aux2 = SSL / (this%SlIS - this%nL * SNDEN)
        aux1 = SST / (this%StIS - this%nT * SNDEN)
        
        aux1 = aux1*aux1
        aux2 = aux2*aux2
        aux3 = aux3*aux3
        !
        plyEvaluateCriteria = sqrt(aux1 + aux2 + aux3)
    end function
    
    
    !_getE11____________________________________________________________
    ! retrieves E11 as a function S11 (or a POS/NEG flag)
    !___________________________________________________________________
    function getE11(this)
    
        implicit none
        
        ! input
        type(PlyDefinition), intent(in) ::  this
        ! output
        real(RWP)                       ::  getE11

        getE11 = this%E11
        if (.not. this%TENCOMP) then
            getE11 = this%E11comp
        end if
    
    end function getE11
    
    !_plyPrint_______________________________________________________________
    ! prints the ply defintion
    !________________________________________________________________________
    subroutine plyPrint(this)
        ! local variables
        type(PlyDefinition), intent(in)     ::  this
1313    format ( A, F12.2, A)
1314    format ( A, E12.4, A)
1315    format ( A, I12, A)
        ! print
        print *, 'Ply: '
        print 1313, ' 1. E11+ ',this%E11
        print 1313, ' 2. E22  ',this%E22
        print 1313, ' 3. v12  ',this%nu12
        print 1313, ' 4. G12  ',this%G12
        print 1313, ' 5. Xt   ',this%Xt
        print 1313, ' 6. Xc   ',this%Xc
        print 1313, ' 7. Yt   ',this%Yt
        print 1313, ' 7. YtIS ',this%YtIS
        print 1313, ' 8. Yc   ',this%Yc
        print 1313, ' 9. Sl   ',this%Sl
        print 1313, ' 9. SlIS ',this%SlIS
        print 1313, '  . St   ',this%St
        print 1313, '  . StIS ',this%StIS
        print 1313, '10. a0   ',this%a0,' degrees'
        print 1314, '11. nL   ',this%nL
        print 1314, '12. nT   ',this%nT
        print 1314, '  . pC   ',this%phiC,' rad'
        
        print 1314, '13. ILSS ',this%ILSS
        print 1314, '14. Zt   ',this%Zt
        
        print 1314, '15. G1mat',this%G1cMat
        print 1314, '16. G2mat',this%G2cMat
        print 1314, '17. G1ft ',this%G1cFibT
        print 1314, '18. G1fki',this%G1cFibK
        print 1314, '19. GExp ',this%GalphaM
        
        print 1313, '20. E11- ',this%E11comp
        print 1313, '21. E33  ',this%E33
        print 1313, '22. v13  ',this%nu13
        print 1313, '23. v23  ',this%nu23
        print 1313, '24. G13  ',this%G13
        print 1313, '25. G23  ',this%G23

        print 1314, '26. mu   ',this%viscosity
        
        print 1315, '27. PTYP ',this%PTYP


        print *
    end subroutine
    
    
    !_plyFromFile____________________________________________________________
    ! retrieves the ply properties from a file
    ! the file should be in the following format:
    !   | First line is unread
    !   | NPROPS (number of properties)
    !   | PROPS  (comma separated properties)
    !________________________________________________________________________
    subroutine plyFromFile( this, filename, readstat)
    
        implicit none
        
        type(PlyDefinition), intent(inout)  ::  this
        character(len=*), intent(in)        ::  filename
        integer, intent(out), optional      ::  readstat
        
        real(RWP), allocatable              ::  PROPS(:)
        integer                             ::  NPROPS=1, allocstat=0
        logical                             ::  existsFile
        
        ! initialize variables
        if (present(readstat) ) readstat = 0
        
        !  check if file exists
        inquire(file=filename, exist=existsFile)
        if (.NOT. existsFile) then  
            print *, 'ERROR: file ', trim(filename), ' not found'
            if (present(readstat))  readstat = 1
            stop
        end if
        
        ! open file
        open (unit=10102, file=filename, action='read')
        ! Skip first line
        read (10102,*)
        ! Read NPROPS
        read (10102,*) NPROPS
        allocate(PROPS(NPROPS), stat=allocstat)
        if (allocstat .GT. 0) then
            print *, 'ERROR: failed allocation'
            if (present(readstat))  readstat = 1
            stop
        end if
        read (10102,*) PROPS
        !print *, PROPS
        close(10102)

        ! Initialize ply object based on PROPS(NPROPS)
        call initializePly(this,PROPS,NPROPS)
        
        ! deallocate and return
        deallocate(PROPS, stat=allocstat)
    
    end subroutine
    
end module plyTools

    
