! constants.f90
! author: sunder

!-----------------------------------------------------------------------
! Number of physical constants and conversion factors frequently used
! in computations
!-----------------------------------------------------------------------

module constants
    implicit none
    public

    ! ---------------------------------------------------------------------------------------------------------
    ! Mathematical Constants
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: m_e                       = 2.71828182845904523536028747135d0      ! e
    double precision, parameter :: m_log2e                   = 1.44269504088896340735992468100d0      ! log_2 (e)
    double precision, parameter :: m_log10e                  = 0.43429448190325182765112891892d0      ! log_10 (e)
    double precision, parameter :: m_sqrt2                   = 1.41421356237309504880168872421d0      ! sqrt(2)
    double precision, parameter :: m_sqrt1_2                 = 0.70710678118654752440084436210d0      ! sqrt(1/2)
    double precision, parameter :: m_sqrt3                   = 1.73205080756887729352744634151d0      ! sqrt(3)
    double precision, parameter :: m_pi                      = 3.14159265358979323846264338328d0      ! pi
    double precision, parameter :: m_pi_2                    = 1.57079632679489661923132169164d0      ! pi/2
    double precision, parameter :: m_pi_4                    = 0.78539816339744830961566084582d0      ! pi/4
    double precision, parameter :: m_sqrtpi                  = 1.77245385090551602729816748334d0      ! sqrt(pi)
    double precision, parameter :: m_1_pi                    = 0.31830988618379067153776752675d0      ! 1/pi
    double precision, parameter :: m_2_pi                    = 0.63661977236758134307553505349d0      ! 2/pi
    double precision, parameter :: m_ln10                    = 2.30258509299404568401799145468d0      ! ln(10)
    double precision, parameter :: m_ln2                     = 0.69314718055994530941723212146d0      ! ln(2)
    double precision, parameter :: m_lnpi                    = 1.14472988584940017414342735135d0      ! ln(pi)
    double precision, parameter :: m_euler                   = 0.57721566490153286060651209008d0      ! Euler constant

    ! ---------------------------------------------------------------------------------------------------------
    ! Fundamental Physical Constants (All units are in MKSA)
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: speed_of_light            = 2.99792458d8                           ! Speed of light in vacuum
    double precision, parameter :: vacuum_permeability       = 1.25663706144d-6                       ! Permeability of free space
    double precision, parameter :: vacuum_permittivity       = 8.854187817d-12                        ! Permittivity of free space
    double precision, parameter :: plancks_constant_h        = 6.62606896d-34                         ! Planck’s constant, h
    double precision, parameter :: plancks_constant_hbar     = 1.05457162825d-34                      ! Planck’s constant divided by 2\pi, h_bar
    double precision, parameter :: num_avogadro              = 6.02214199d23                          ! Avogadro’s number, 1 mol
    double precision, parameter :: faraday                   = 9.64853429775d4                        ! The molar charge of 1 Faraday
    double precision, parameter :: boltzman                  = 1.3806504d-23                          ! Boltzmann constant, k
    double precision, parameter :: molar_gas                 = 8.314472d0                             ! Molar gas constant
    double precision, parameter :: standard_gas_volume       = 2.2710981d-2                           ! Standard gas volume
    double precision, parameter :: stefan_boltzman_constant  = 5.67040047374d-8                       ! Stefan-Boltzmann radiation constant
    double precision, parameter :: gauss                     = 1.0d-4                                 ! Magnetic field of 1 Gauss

    ! ---------------------------------------------------------------------------------------------------------
    ! Astronomy and Astrophysics
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: astronomical_unit         = 1.49597870691d11                       ! The length of 1 astronomical unit (mean earth-sun distance)
    double precision, parameter :: gravitational_constant    = 6.673d-11                              ! Universal gravitational constant
    double precision, parameter :: light_year                = 9.46053620707d15                       ! The distance of 1 light-year
    double precision, parameter :: parsec                    = 3.08567758135d16                       ! The distance of 1 parsec
    double precision, parameter :: grav_accel                = 9.80665d0                              ! Standard gravitational acceleration on Earth, g
    double precision, parameter :: solar_mass                = 1.98892d30                             ! Mass of the Sun

    ! ---------------------------------------------------------------------------------------------------------
    ! Atomic and Nuclear Physics
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: electron_charge           = 1.602176487d-19                        ! Charge of the electron
    double precision, parameter :: electron_volt             = 1.602176487d-19                        ! Energy of 1 electron volt
    double precision, parameter :: unified_atomic_mass       = 1.660538782d-27                        ! Unified atomic mass
    double precision, parameter :: mass_electron             = 9.10938188d-31                         ! Mass of the electron
    double precision, parameter :: mass_muon                 = 1.88353109d-28                         ! Mass of the muon
    double precision, parameter :: mass_proton               = 1.67262158d-27                         ! Mass of the proton
    double precision, parameter :: mass_neutron              = 1.67492716d-27                         ! Mass of the neutron
    double precision, parameter :: num_fine_structure        = 7.297352533d-3                         ! Electromagnetic fine structure constant
    double precision, parameter :: rydberg                   = 2.17987196968d-18                      ! Rydberg constant
    double precision, parameter :: bohr_radius               = 5.291772083d-11                        ! Bohr radius
    double precision, parameter :: angstrom                  = 1.0d-10                                ! Length of 1 angstrom
    double precision, parameter :: barn                      = 1.0d-28                                ! Area of 1 barn
    double precision, parameter :: bohr_magneton             = 9.27400899d-24                         ! Bohr Magneton
    double precision, parameter :: nuclear_magneton          = 5.05078317d-27                         ! Nuclear Magneton
    double precision, parameter :: electron_magnetic_moment  = 9.28476362d-24                         ! Absolute value of the magnetic moment of the electron
    double precision, parameter :: proton_magnetic_moment    = 1.410606633d-26                        ! Magnetic moment of the proton
    double precision, parameter :: thomson_cross_section     = 6.65245893699d-29                      ! Thomson cross section
    double precision, parameter :: debye                     = 3.33564095198d-30                      ! Electric dipole moment of 1 Debye

    ! ---------------------------------------------------------------------------------------------------------
    ! Measurement of Time
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: minute                    = 60.d0                                   ! Number of seconds in 1 minute
    double precision, parameter :: hour                      = 3.6d3                                  ! Number of seconds in 1 hour
    double precision, parameter :: day                       = 8.64d4                                 ! Number of seconds in 1 day
    double precision, parameter :: week                      = 6.048d5                                ! Number of seconds in 1 week

    ! ---------------------------------------------------------------------------------------------------------
    ! Imperial units
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: inch                      = 2.54d-2                                ! Length of 1 inch
    double precision, parameter :: foot                      = 3.048d-1                               ! Length of 1 foot
    double precision, parameter :: yard                      = 9.144d-1                               ! Length of 1 yard
    double precision, parameter :: mile                      = 1.609344d3                             ! Length of 1 mile
    double precision, parameter :: mil                       = 2.54d-5                                ! Length of 1 mil (1/1000 of inch)

    ! ---------------------------------------------------------------------------------------------------------
    ! Speed and Nautical Units
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: kilometers_per_hour       = 2.77777777778d-1                       ! Speed of 1 kilometer per hour
    double precision, parameter :: miles_per_hour            = 4.4704d-1                              ! Speed of 1 mile per hour
    double precision, parameter :: nautical_mile             = 1.852d3                                ! One Nautical mile
    double precision, parameter :: fathom                    = 1.8288d0                               ! Length of 1 fathom
    double precision, parameter :: knot                      = 5.14444444444d-1                       ! Speed of 1 knot

    ! ---------------------------------------------------------------------------------------------------------
    ! Printers Units
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: point                     = 3.52777777778d-4                       ! Length of 1 printer’s point (1/72 inch)
    double precision, parameter :: texpoint                  = 3.51459803515d-4                       ! length of 1 TeX point (1/72.27 inch)

    ! ---------------------------------------------------------------------------------------------------------
    ! Volume, Area and Length
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: micron                    = 1.0d-6                                 ! Length of 1 micron
    double precision, parameter :: hectare                   = 1.0d4                                  ! Area of 1 hectare
    double precision, parameter :: acre                      = 4.04685642241d3                        ! Area of 1 acre
    double precision, parameter :: liter                     = 1.0d-3                                 ! Volume of 1 liter
    double precision, parameter :: us_gallon                 = 3.78541178402d-3                       ! Volume of 1 US gallon
    double precision, parameter :: canadian_gallon           = 4.54609d-3                             ! Volume of 1 Canadian gallon
    double precision, parameter :: uk_gallon                 = 4.546092d-3                            ! Volume of 1 UK gallon
    double precision, parameter :: quart                     = 9.46352946004d-4                       ! Volume of 1 quart
    double precision, parameter :: pint                      = 4.73176473002d-4                       ! Volume of 1 pint

    ! ---------------------------------------------------------------------------------------------------------
    ! Mass and Weight
    ! ---------------------------------------------------------------------------------------------------------

    double precision, parameter :: pound_mass                = 4.5359237d-1                           ! Mass of 1 pound
    double precision, parameter :: ounce_mass                = 2.8349523125d-2                        ! Mass of 1 ounce
    double precision, parameter :: ton                       = 9.0718474d2                            ! Mass of 1 ton
    double precision, parameter :: metric_ton                = 1.0d3                                  ! Mass of 1 metric ton
    double precision, parameter :: uk_ton                    = 1.0160469088d3                         ! Mass of 1 UK ton
    double precision, parameter :: troy_ounce                = 3.1103475d-2                           ! Mass of 1 troy ounce
    double precision, parameter :: carat                     = 2.0d-4                                 ! Mass of 1 carat
    double precision, parameter :: gram_force                = 9.80665d-3                             ! Weight force of 1 gram
    double precision, parameter :: pound_force               = 4.44822161526d0                        ! Weight force of 1 pound
    double precision, parameter :: kilopound_force           = 4.44822161526d3                        ! Weight force of 1 kilo-pound
    double precision, parameter :: poundal                   = 1.38255d-1                             ! Weight force of 1 poundal

    ! ---------------------------------------------------------------------------------------------------------
    ! Thermal energy and power
    ! ---------------------------------------------------------------------------------------------------------





end module constants

