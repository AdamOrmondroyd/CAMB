from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool


class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int)]

    def validate_params(self):
        return True


class DarkEnergyEqnOfState(DarkEnergyModel):
    """
    Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or call set_w_a_table to set another tabulated w(a). If tabulated w(a) is used, w and wa are set
    to approximate values at z=0.

    See :meth:`.model.CAMBparams.set_initial_power_function` for a convenience constructor function to
    set a general interpolated P(k) model from a python function.

    """
    _fortran_class_module_ = 'DarkEnergyInterface'
    _fortran_class_name_ = 'TDarkEnergyEqnOfState'

    _fields_ = [
        ("w", c_double, "w(0)"),
        ("wa", c_double, "-dw/da(0)"),
        ("cs2", c_double, "fluid rest-frame sound speed squared"),
        ("use_tabulated_w", c_bool, "using an interpolated tabulated w(a) rather than w, wa above"),
        ("__no_perturbations", c_bool, "turn off perturbations (unphysical, so hidden in Python)")
    ]

    _methods_ = [('SetWTable', [numpy_1d, numpy_1d, POINTER(c_int)])]

    def set_params(self, w=-1.0, wa=0, cs2=1.0):
        """
         Set the parameters so that P(a)/rho(a) = w(a) = w + (1-a)*wa

        :param w: w(0)
        :param wa: -dw/da(0)
        :param cs2: fluid rest-frame sound speed squared
        """
        self.w = w
        self.wa = wa
        self.cs2 = cs2
        self.validate_params()

    def validate_params(self):
        if not self.use_tabulated_w and self.wa + self.w > 0:
            raise CAMBError('dark energy model has w + wa > 0, giving w>0 at high redshift')

    def set_w_a_table(self, a, w):
        """
        Set w(a) from numerical values (used as cubic spline). Note this is quite slow.

        :param a: array of scale factors
        :param w: array of w(a)
        :return: self
        """
        if len(a) != len(w):
            raise ValueError('Dark energy w(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):
            raise ValueError('Dark energy w(a) arrays must end at a=1')
        if np.any(a <= 0):
            raise ValueError('Dark energy w(a) table cannot be set for a<=0')

        a = np.ascontiguousarray(a, dtype=np.float64)
        w = np.ascontiguousarray(w, dtype=np.float64)

        self.f_SetWTable(a, w, byref(c_int(len(a))))
        return self

    def __getstate__(self):
        if self.use_tabulated_w:
            raise TypeError("Cannot save class with splines")
        return super().__getstate__()


@fortran_class
class DarkEnergyFluid(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization using the constant sound-speed single fluid model
    (as for single-field quintessence).

    """

    _fortran_class_module_ = 'DarkEnergyFluid'
    _fortran_class_name_ = 'TDarkEnergyFluid'

    def validate_params(self):
        super().validate_params()
        if not self.use_tabulated_w:
            if self.wa and (self.w < -1 - 1e-6 or 1 + self.w + self.wa < - 1e-6):
                raise CAMBError('fluid dark energy model does not support w crossing -1')

    def set_w_a_table(self, a, w):
        # check w array has elements that do not cross -1
        if np.sign(1 + np.max(w)) - np.sign(1 + np.min(w)) == 2:
            raise CAMBError('fluid dark energy model does not support w crossing -1')
        super().set_w_a_table(a, w)


@fortran_class
class DarkEnergyPPF(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization in the PPF perturbation approximation
    (`arXiv:0808.3125 <https://arxiv.org/abs/0808.3125>`_)
    Use inherited methods to set parameters or interpolation table.

    Note PPF is not a physical model and just designed to allow crossing -1 in an ad hoc smooth way. For models
    with w>-1 but far from cosmological constant, it can give quite different answers to the fluid model with c_s^2=1.

    """
    # cannot declare c_Gamma_ppf directly here as have not defined all fields in DarkEnergyEqnOfState (TCubicSpline)
    _fortran_class_module_ = 'DarkEnergyPPF'
    _fortran_class_name_ = 'TDarkEnergyPPF'


@fortran_class
class AxionEffectiveFluid(DarkEnergyModel):
    """
    Example implementation of a specific (early) dark energy fluid model
    (`arXiv:1806.10608 <https://arxiv.org/abs/1806.10608>`_).
    Not well tested, but should serve to demonstrate how to make your own custom classes.
    """
    _fields_ = [
        ("w_n", c_double, "effective equation of state parameter"),
        ("fde_zc", c_double, "energy density fraction at z=zc"),
        ("zc", c_double, "decay transition redshift (not same as peak of energy density fraction)"),
        ("theta_i", c_double, "initial condition field value")]

    _fortran_class_name_ = 'TAxionEffectiveFluid'
    _fortran_class_module_ = 'DarkEnergyFluid'

    def set_params(self, w_n, fde_zc, zc, theta_i=None):
        self.w_n = w_n
        self.fde_zc = fde_zc
        self.zc = zc
        if theta_i is not None:
            self.theta_i = theta_i


# base class for scalar field quintessence models
class Quintessence(DarkEnergyModel):
    r"""
    Abstract base class for single scalar field quintessence models.

    For each model the field value and derivative are stored and splined at sampled scale factor values.

    To implement a new model, need to define a new derived class in Fortran,
    defining Vofphi and setting up initial conditions and interpolation tables (see TEarlyQuintessence as example).

    """
    _fields_ = [
        ("DebugLevel", c_int),
        ("astart", c_double),
        ("integrate_tol", c_double),
        ("sampled_a", AllocatableArrayDouble),
        ("phi_a", AllocatableArrayDouble),
        ("phidot_a", AllocatableArrayDouble),
        ("__npoints_linear", c_int),
        ("__npoints_log", c_int),
        ("__dloga", c_double),
        ("__da", c_double),
        ("__log_astart", c_double),
        ("__max_a_log", c_double),
        ("__ddphi_a", AllocatableArrayDouble),
        ("__ddphidot_a", AllocatableArrayDouble),
        ("__state", f_pointer)
    ]
    _fortran_class_module_ = 'Quintessence'

    def __getstate__(self):
        raise TypeError("Cannot save class with splines")


@fortran_class
class EarlyQuintessence(Quintessence):
    r"""
    Example early quintessence (axion-like, as `arXiv:1908.06995 <https://arxiv.org/abs/1908.06995>`_) with potential

     V(\phi) = m^2f^2 (1 - cos(\phi/f))^n + \Lambda_{cosmological constant}

    """

    _fields_ = [
        ("n", c_double, "power index for potential"),
        ("f", c_double, r"f/Mpl (sqrt(8\piG)f); only used for initial search value when use_zc is True"),
        ("m", c_double, "mass parameter in reduced Planck mass units; "
                        "only used for initial search value when use_zc is True"),
        ("theta_i", c_double, "phi/f initial field value"),
        ("frac_lambda0", c_double, "fraction of dark energy in cosmological constant today (approximated as 1)"),
        ("use_zc", c_bool, "solve for f, m to get specific critical redshift zc and fde_zc"),
        ("zc", c_double, "redshift of peak fractional early dark energy density"),
        ("fde_zc", c_double, "fraction of early dark energy density to total at peak"),
        ("npoints", c_int, "number of points for background integration spacing"),
        ("min_steps_per_osc", c_int, "minimum number of steps per background oscillation scale"),
        ("fde", AllocatableArrayDouble, "after initialized, the calculated background early dark energy "
                                        "fractions at sampled_a"),
        ("__ddfde", AllocatableArrayDouble)
    ]
    _fortran_class_name_ = 'TEarlyQuintessence'

    def set_params(self, n, f=0.05, m=5e-54, theta_i=0.0, use_zc=True, zc=None, fde_zc=None):
        self.n = n
        self.f = f
        self.m = m
        self.theta_i = theta_i
        self.use_zc = use_zc
        if use_zc:
            if zc is None or fde_zc is None:
                raise ValueError("must set zc and fde_zc if using 'use_zc'")
            self.zc = zc
            self.fde_zc = fde_zc


@fortran_class
class DarkEnergyFK(DarkEnergyModel):
    """
    FlexKnot dark energy model using linear splines between nodes.

    This model parameterizes w(a) as piecewise linear segments (flexknots) between
    specified knot points. The parameterization uses:
        w_knots = [w₀, w₁, ..., w_{n-1}]  (n values)
        a_knots = [a₁, a₂, ..., a_{n-2}]   (n-2 values)

    Creating knots at: (1.0, w₀), (a₁, w₁), ..., (a_{n-2}, w_{n-2}), (a_min, w_{n-1})
    where a_min ≈ 10⁻⁸ represents the early universe limit.

    This avoids cubic spline overhead while providing exact linear segments that
    work perfectly with the Effective_w_wa method since flexknots are exactly
    linear, making w + wa(1-a) not an approximation but the exact form.
    """

    _fortran_class_module_ = 'DarkEnergyFK'
    _fortran_class_name_ = 'TDarkEnergyFK'

    _fields_ = [
        ("__n_w", c_int, "number of w-values"),
        ("__a_min", c_double, "minimum scale factor")
    ]

    _methods_ = [('SetFlexKnots', [POINTER(c_int), numpy_1d, numpy_1d])]

    def set_flexknots(self, w_knots, a_knots=None):
        """
        Set flexknot parameterization for w(a).

        Parameters:
        -----------
        w_knots : array_like
            Array of w-values [w₀, w₁, ..., w_{n-1}] where:
            - w₀ is w(a=1) today
            - w_{n-1} is w(a→0) in early universe

        a_knots : array_like, optional
            Array of interior knot positions [a₁, a₂, ..., a_{n-2}].
            Must be in decreasing order and in range (a_min, 1.0).
            If None and len(w_knots) > 1, creates logarithmically spaced knots.

        Returns:
        --------
        self : DarkEnergyFlexKnot
            Returns self for method chaining

        Examples:
        ---------
        # Cosmological constant
        de = DarkEnergyFlexKnot()
        de.set_flexknots([-1.0])

        # Simple w0-wa model equivalent
        de.set_flexknots([-1.0, -0.5], [0.5])  # Linear from w=-1 at a=1 to w=-0.5 at a=0

        # More complex flexknot model
        w_vals = [-0.9, -1.1, -0.8, -0.6]
        a_knots = [0.7, 0.3]  # Creates knots at (1.0,-0.9), (0.7,-1.1), (0.3,-0.8), (a_min,-0.6)
        de.set_flexknots(w_vals, a_knots)
        """
        w_knots = np.ascontiguousarray(w_knots, dtype=np.float64)
        n_w = len(w_knots)

        if n_w == 1:
            # Cosmological constant case
            a_knots_array = np.array([0.0], dtype=np.float64)  # Dummy value, not used
        else:
            if a_knots is None:
                # Create logarithmically spaced interior knots
                if n_w == 2:
                    a_knots_array = np.array([0.0], dtype=np.float64)  # Dummy, no interior knots needed
                else:
                    # Create n_w-2 logarithmically spaced knots between a_min and 1
                    a_min = 1e-8
                    a_knots_array = np.logspace(np.log10(a_min), 0, n_w-1)[1:-1][::-1]  # Reverse for decreasing order
                    a_knots_array = np.ascontiguousarray(a_knots_array, dtype=np.float64)
            else:
                a_knots = np.asarray(a_knots, dtype=np.float64)
                if len(a_knots) != n_w - 2:
                    raise ValueError(f"Expected {n_w-2} a_knots for {n_w} w_knots, got {len(a_knots)}")

                # Validate knots are in decreasing order
                if n_w > 2 and not np.all(a_knots[:-1] > a_knots[1:]):
                    raise ValueError("a_knots must be in decreasing order")

                # Validate range
                a_min = 1e-8
                if n_w > 2 and (np.any(a_knots >= 1.0) or np.any(a_knots <= a_min)):
                    raise ValueError(f"a_knots must be in range ({a_min}, 1.0)")

                a_knots_array = np.ascontiguousarray(a_knots, dtype=np.float64)

        # Call Fortran subroutine
        self.f_SetFlexKnots(byref(c_int(n_w)), w_knots, a_knots_array)

        return self

    def set_params(
        self, w0,
        a1=None, w1=None,
        a2=None, w2=None,
        a3=None, w3=None,
        a4=None, w4=None,
        a5=None, w5=None,
        a6=None, w6=None,
        a7=None, w7=None,
        a8=None, w8=None,
        a9=None, w9=None,
        a10=None, w10=None,
        a11=None, w11=None,
        a12=None, w12=None,
        a13=None, w13=None,
        a14=None, w14=None,
        a15=None, w15=None,
        a16=None, w16=None,
        a17=None, w17=None,
        a18=None, w18=None,
        a19=None, w19=None,
        a20=None, w20=None,
        a21=None, w21=None,
        a22=None, w22=None,
        a23=None, w23=None,
        wn=None,
    ):
        """
        Convenience method that calls set_flexknots.
        Provided for consistency with other dark energy classes.
        """
        args = locals()
        w_knots = np.array([
            args[f'w{i}'] for i in range(24) if args[f'w{i}'] is not None
        ], dtype=np.float64)
        if wn is not None:
            w_knots = np.concatenate([w_knots, [wn]])

        a_knots = np.array([
            args[f'a{i}'] for i in range(1, 24) if args[f'a{i}'] is not None
        ], dtype=np.float64)
        return self.set_flexknots(w_knots, a_knots)

    def __getstate__(self):
        # FlexKnot can be pickled since it doesn't use splines
        return super().__getstate__()


# short names for models that support w/wa
F2003Class._class_names.update({'fluid': DarkEnergyFluid, 'ppf': DarkEnergyPPF, 'flexknot': DarkEnergyFK})
