module inputs
!Inputs to the program.

implicit none

	!m_1(2) = Mass of nucleon 1(2) in MeV
	real(8), parameter, public :: m_1 = 938.926d0
	real(8), parameter, public :: m_2 = 938.926d0

	!q = |q| Momentum variable.
	!n_theta = Number of Gauss-Legendre nodes to evaluate V (the potential) at i.e. the "b" in Ax = b.
	!T = Isospin. The value is 0 or 1
	real(8), parameter, public :: q = 375.2851d0
	integer, parameter, public :: n_theta = 40
	integer, parameter, public :: T = 1

	!Pcm = bvec{P}-->P z-hat. By definition it's one half the CM momentum.
	!k_fermi = The Fermi momentum. Units are [fm^(-1)]
	!Q_pauli = 0 for medium and 1 for free space.
	real(8), parameter, public :: Pcm = q
	real(8), parameter, public :: k_fermi = 1.4d0
	real(8), parameter, public :: Q_pauli = 0.0d0

	!n_V_avg = Number of integrating nodes for the phi integrated potential on [0,2pi].
	!n_0_2q,n_2q_inf = Number of integrating nodes for the integral over momentum (k variable). The interval is [0,2q] and [2q,inf] respectivly.
	!	See Eq. 15 of T matrix rev. report. **MAKE SURE n_0_2q IS EVEN**
	!n_theta_k = Number of integrating nodes for the integral over momentum (theta_k) space on [0,pi].
	integer, parameter, public :: n_V_avg = 20
	integer, parameter, public :: n_0_2q = 30
	integer, parameter, public :: n_2q_inf = 100
	integer, parameter, public :: n_theta_k = 80

	!Cutoff_inf = Value for infinity i.e. the [2q,infinity] integral.
	real(8), parameter, public :: cutoff_inf = 2000.0d0

end module
