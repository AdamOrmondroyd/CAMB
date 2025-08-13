    module DarkEnergyFK
    use precision
    use DarkEnergyInterface
    use classes
    implicit none

    private

    type, extends(TDarkEnergyModel) :: TDarkEnergyFK
        integer :: n_w = 1                           ! Number of w-values
        real(dl), allocatable :: w_knots(:)         ! [w0, w1, ..., w_{n-1}]
        real(dl), allocatable :: a_knots(:)          ! [a1, a2, ..., a_{n-2}]
        real(dl) :: a_min = 1.0e-3_dl                ! Minimum scale factor
        ! Internal full knot arrays (computed from inputs)
        real(dl), allocatable :: full_a(:), full_w(:)
        real(dl), allocatable :: knot_log_density(:)
        logical :: initialized = .false.
        real(dl) :: c_Gamma_ppf = 0.4_dl  ! PPF anisotropy parameter
    contains
        ! Public DE methods we must override
        ! NOTE: copying TAxionEffectiveFluid
        procedure :: ReadParams => TDarkEnergyFK_ReadParams
        procedure, nopass :: PythonClass => TDarkEnergyFK_PythonClass
        procedure, nopass :: SelfPointer => TDarkEnergyFK_SelfPointer
        procedure :: Init => TDarkEnergyFK_Init
        ! TODO:
        procedure :: PerturbedStressEnergy => TDarkEnergyFK_PerturbedStressEnergy
        ! NOTE: diff_rhopi_Add_Term is inherited from TDarkEnergyEqnOfState
        ! end todo
        procedure :: w_de => TDarkEnergyFK_w_de
        procedure :: grho_de => TDarkEnergyFK_grho_de
        procedure :: Effective_w_wa => TDarkEnergyFK_Effective_w_wa
        procedure :: PrintFeedback => TDarkEnergyFK_PrintFeedback
        procedure :: SetFlexKnots => TDarkEnergyFK_SetFlexKnots
        ! Flexknot-specific methods
        procedure, private :: BuildFullKnots
        procedure, private :: FindSegment
        procedure, private :: ComputeDensityCache
        procedure, private :: IntegrateSegment
    end type TDarkEnergyFK

    public TDarkEnergyFK
    contains

    subroutine TDarkEnergyFK_ReadParams(this, Ini)
        use IniObjects
        class(TDarkEnergyFK) :: this
        class(TIniFile), intent(in) :: Ini
        call this%TDarkEnergyModel%ReadParams(Ini)
        ! TODO: Implementation for reading from .ini files can be added here
        error stop 'FlexKnot INI parameters not yet implemented - use Python interface'

    end subroutine TDarkEnergyFK_ReadParams

    function TDarkEnergyFK_PythonClass()
        character(LEN=:), allocatable :: TDarkEnergyFK_PythonClass
        TDarkEnergyFK_PythonClass = 'DarkEnergyFK'
    end function TDarkEnergyFK_PythonClass

    subroutine TDarkEnergyFK_SelfPointer(cptr, P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type(TDarkEnergyFK), pointer :: PType
        class(TPythonInterfacedClass), pointer :: P
        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TDarkEnergyFK_SelfPointer

    subroutine TDarkEnergyFK_Init(this, State)
        use classes
        use config  ! TODO: why?
        class(TDarkEnergyFK), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        call this%TDarkEnergyModel%Init(State)

        if (this%is_cosmological_constant) then
            this%num_perturb_equations = 0
        else
            this%num_perturb_equations = 1  ! PPF-style
      end if
    end subroutine TDarkEnergyFK_Init

    subroutine TDarkEnergyFK_SetFlexKnots(this, n_w, w_knots, a_knots)
        class(TDarkEnergyFK), intent(inout) :: this
        integer, intent(in) :: n_w
        real(dl), intent(in) :: w_knots(n_w)
        real(dl), intent(in) :: a_knots(*)  ! Assumed size to handle n_w=1 case
        integer :: i

        ! Input validation
        if (n_w < 1) then
            error stop 'FlexKnot: n_w must be >= 1'
        end if

        this%n_w = n_w

        ! Allocate arrays safely
        if (allocated(this%w_knots)) deallocate(this%w_knots)
        if (allocated(this%a_knots)) deallocate(this%a_knots)

        allocate(this%w_knots(n_w))
        this%w_knots = w_knots

        if (n_w > 2) then
            allocate(this%a_knots(n_w-2))
            do i = 1, n_w-2
                this%a_knots(i) = a_knots(i)
            end do

            ! Validate: a_knots should be decreasing and in (0,1)
            if (n_w > 3) then
                do i = 1, n_w-3
                    if (this%a_knots(i) <= this%a_knots(i+1)) then
                        error stop 'FlexKnot: a_knots must be in decreasing order'
                    end if
                end do
            end if

            if (this%a_knots(1) >= 1.0_dl .or. this%a_knots(n_w-2) <= this%a_min) then
                error stop 'FlexKnot: a_knots must be in range (a_min, 1.0)'
            end if
        end if

        call this%BuildFullKnots()
        call this%ComputeDensityCache()
        this%is_cosmological_constant = (n_w == 1 .and. abs(w_knots(1) + 1.0_dl) < 1e-6_dl)
        this%initialized = .true.

    end subroutine TDarkEnergyFK_SetFlexKnots

    subroutine BuildFullKnots(this)
        class(TDarkEnergyFK), intent(inout) :: this
        integer :: i

        ! Allocate full knot arrays
        if (allocated(this%full_a)) deallocate(this%full_a)
        if (allocated(this%full_w)) deallocate(this%full_w)
        if (allocated(this%knot_log_density)) deallocate(this%knot_log_density)

        allocate(this%full_a(this%n_w), this%full_w(this%n_w))
        allocate(this%knot_log_density(this%n_w))

        ! Build full knot arrays
        this%full_a(1) = 1.0_dl                    ! Today
        this%full_w(1) = this%w_knots(1)          ! w0

        if (this%n_w > 2) then
            ! Interior knots
            do i = 1, this%n_w-2
                this%full_a(i+1) = this%a_knots(i)
                this%full_w(i+1) = this%w_knots(i+1)
            end do
        end if

        this%full_a(this%n_w) = this%a_min         ! Early time
        this%full_w(this%n_w) = this%w_knots(this%n_w)  ! w_{n-1}

    end subroutine BuildFullKnots

    subroutine FindSegment(this, a, idx, a1, a2, w1, w2)
        class(TDarkEnergyFK), intent(in) :: this
        real(dl), intent(in) :: a
        integer, intent(out) :: idx
        real(dl), intent(out) :: a1, a2, w1, w2
        integer :: i

        ! Find segment (full_a is in decreasing order)
        idx = this%n_w  ! Default to last segment
        do i = 1, this%n_w - 1
            if (a >= this%full_a(i+1)) then
                idx = i
                exit
            end if
        end do

        if (idx == this%n_w) then
            ! Extrapolate with constant w at early times
            a1 = this%full_a(this%n_w)
            a2 = this%full_a(this%n_w)
            w1 = this%full_w(this%n_w)
            w2 = this%full_w(this%n_w)
        else
            a1 = this%full_a(idx)
            a2 = this%full_a(idx+1)
            w1 = this%full_w(idx)
            w2 = this%full_w(idx+1)
        end if

    end subroutine FindSegment

    function TDarkEnergyFK_w_de(this, a) result(w)
        class(TDarkEnergyFK) :: this
        real(dl), intent(in) :: a
        real(dl) :: w
        integer :: idx
        real(dl) :: a1, a2, w1, w2

        if (.not. this%initialized) then
            w = -1.0_dl
            return
        end if

        call this%FindSegment(a, idx, a1, a2, w1, w2)

        if (abs(a2 - a1) < 1e-15_dl) then
            w = w1  ! Constant segment
        else
            w = w1 + (w2 - w1) * (a - a1) / (a2 - a1)  ! Linear interpolation
        end if

    end function TDarkEnergyFK_w_de

    subroutine TDarkEnergyFK_Effective_w_wa(this, w, wa)
        class(TDarkEnergyFK), intent(inout) :: this
        real(dl), intent(out) :: w, wa
        integer :: idx
        real(dl) :: a1, a2, w1, w2, slope

        if (.not. this%initialized) then
            w = -1.0_dl
            wa = 0.0_dl
            return
        end if

        ! Find segment containing a=1.0 (today)
        call this%FindSegment(1.0_dl, idx, a1, a2, w1, w2)

        if (abs(a2 - a1) < 1e-15_dl) then
            w = w1
            wa = 0.0_dl
        else
            ! Convert w(a) = w1 + slope*(a-a1) to w + wa*(1-a) form
            slope = (w2 - w1) / (a2 - a1)
            w = w1 + slope * (1.0_dl - a1)  ! Value at a=1
            wa = -slope                     ! dw/d(1-a) = -dw/da
        end if

    end subroutine TDarkEnergyFK_Effective_w_wa

    subroutine IntegrateSegment(this, a_start, a_end, w_start, w_end, integral)
        class(TDarkEnergyFK), intent(in) :: this
        real(dl), intent(in) :: a_start, a_end, w_start, w_end
        real(dl), intent(out) :: integral
        real(dl) :: slope, log_ratio

        if (abs(a_end - a_start) < 1e-15_dl) then
            integral = 0.0_dl
            return
        end if

        slope = (w_end - w_start) / (a_end - a_start)
        log_ratio = log(a_end / a_start)

        ! Analytical integral of -3 * ∫[(1+w(a))/a] da
        ! where w(a) = w_start + slope*(a - a_start)
        integral = -3.0_dl * ((1.0_dl + w_start) * log_ratio + &
                             slope * (a_end - a_start - a_start * log_ratio))

    end subroutine IntegrateSegment

    subroutine ComputeDensityCache(this)
        class(TDarkEnergyFK), intent(inout) :: this
        integer :: i
        real(dl) :: integral

        ! Compute log(a^4 * rho_de / rho_de(a=1)) at each knot
        ! Following the pattern from TDarkEnergyEqnOfState_SetwTable
        
        ! The reference formula is: ln(a^4) - 3*int[1 to a] (1+w) d ln a
        ! For our piecewise linear w(a), integrate segment by segment
        
        this%knot_log_density(1) = 0.0_dl  ! At a=1, this is 0 by definition
        
        do i = 2, this%n_w
            ! Integrate this segment from full_a(i-1) to full_a(i)
            call this%IntegrateSegment(this%full_a(i-1), this%full_a(i), &
                                      this%full_w(i-1), this%full_w(i), integral)
            
            ! Accumulate: log(a^4 rho_de(a_i)) relative to log(a^4 rho_de(a=1))
            this%knot_log_density(i) = this%knot_log_density(i-1) + integral + &
                                       4.0_dl * log(this%full_a(i) / this%full_a(i-1))
        end do

    end subroutine ComputeDensityCache

    function TDarkEnergyFK_grho_de(this, a) result(grho_de)  !relative density (8 pi G a^4 rho_de /grhov)
        class(TDarkEnergyFK) :: this
        real(dl), intent(in) :: a
        real(dl) :: grho_de
        integer :: idx
        real(dl) :: a1, a2, w1, w2, log_rho, extra_integral, w_at_a

        if (a == 0.0_dl) then
            grho_de = 0.0_dl
            return
        end if

        if (a >= 1.0_dl) then
            grho_de = 1.0_dl
            return
        end if

        if (.not. this%initialized) then
            grho_de = 1.0_dl
            return
        end if

        call this%FindSegment(a, idx, a1, a2, w1, w2)

        if (idx == this%n_w) then
            ! Constant extrapolation beyond last knot
            ! ln(a^4 rho) = ln(a_knot^4 rho_knot) + [4 - 3(1+w)] * ln(a/a_knot)
            log_rho = this%knot_log_density(this%n_w) + &
                     (4.0_dl - 3.0_dl * (1.0_dl + this%full_w(this%n_w))) * log(a / this%full_a(this%n_w))
        else
            ! Linear interpolation within segment
            w_at_a = w1 + (w2-w1)*(a-a1)/(a2-a1)
            call this%IntegrateSegment(a1, a, w1, w_at_a, extra_integral)
            log_rho = this%knot_log_density(idx) + extra_integral + &
                     4.0_dl * log(a / a1)
        end if

        grho_de = exp(log_rho)

    end function TDarkEnergyFK_grho_de

    subroutine TDarkEnergyFK_PrintFeedback(this, FeedbackLevel)
        class(TDarkEnergyFK) :: this
        integer, intent(in) :: FeedbackLevel
        integer :: i

        if (FeedbackLevel > 0 .and. this%initialized) then
            write(*,'("FlexKnot Dark Energy:")')
            write(*,'("  w-values: ",20F8.4)') this%w_knots
            if (this%n_w > 2) then
                write(*,'("  a-knots:  ",20F8.4)') this%a_knots
            end if
            write(*,'("  Full knots:")')
            do i = 1, this%n_w
                write(*,'("    (",F8.5,", ",F8.4,")")') this%full_a(i), this%full_w(i)
            end do
        end if

    end subroutine TDarkEnergyFK_PrintFeedback

    subroutine TDarkEnergyFK_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyFK), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: Gamma, S_Gamma, ckH, Gammadot, Fa, sigma
    real(dl) :: vT, grhoT, k2

    k2=k**2
    !ppf
    grhoT = grho - grhov_t
    vT = dgq / (grhoT + gpres_noDE)
    Gamma = ay(w_ix)

    !sigma for ppf
    sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / kf1 - &
        k * Gamma
    sigma = sigma / adotoa

    S_Gamma = grhov_t * (1 + w) * (vT + sigma) * k / adotoa / 2._dl / k2
    ckH = this%c_Gamma_ppf * k / adotoa

    if (ckH * ckH > 1000) then
        ! Was ckH^2 > 30 originally, but this is better behaved (closer to fluid)
        ! for some extreme models (thanks Yanhui Yang, Simeon Bird 2024)
        Gamma = 0
        Gammadot = 0.d0
    else
        Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
        Gammadot = Gammadot * adotoa
    endif
    ayprime(w_ix) = Gammadot !Set this here, and don't use PerturbationEvolve

    Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / kf1
    dgqe = S_Gamma - Gammadot / adotoa - Gamma
    dgqe = -dgqe / Fa * 2._dl * k * adotoa + vT * grhov_t * (1 + w)
    dgrhoe = -2 * k2 * kf1 * Gamma - 3 / k * adotoa * dgqe

    end subroutine TDarkEnergyFK_PerturbedStressEnergy

    end module DarkEnergyFK
