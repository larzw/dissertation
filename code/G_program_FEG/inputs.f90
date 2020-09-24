module inputs
!Inputs to the program.

implicit none

	!m_1(2) = Mass of nucleon 1(2) in MeV
	real(8), parameter, public :: m_1 = 938.926d0
	real(8), parameter, public :: m_2 = 938.926d0

	!q = |q|.
	!T = Isospin. The value is 0 or 1
	real(8), parameter, public :: q = 306.42d0
	integer, parameter, public :: T = 0

	!Pcm = One half the CM momentum.
	!k_fermi_1 and k_fermi_2 = Fermi momentum in fm^(-1) for nucleon 1 and 2.
	!Q_pauli = 0 for medium and 1 for free space.
	real(8), parameter, public :: Pcm = q
	real(8), parameter, public :: k_fermi_1 = 1.4d0
	real(8), parameter, public :: k_fermi_2 = 1.4d0
	real(8), parameter, public :: Q_pauli = 1.0d0

	!n_V_avg = Nodes/weights for phi integrated potential over [0,2pi].
	!n_0_2q,n_2q_inf = Nodes/weights for [0,inf] integral. The interval is split into [0,2q] and [2q,inf].
	!**MAKE SURE n_0_2q IS EVEN**
	!n_theta_k = Nodes/weights for [0,pi] integral
	integer, parameter, public :: n_V_avg = 20
	integer, parameter, public :: n_0_2q = 2
	integer, parameter, public :: n_2q_inf = 2
	integer, parameter, public :: n_theta_k = 2

	!cutoff0_inf = Value for infinity i.e. the [2q,infinity] integral. Used for ^0T isospin = 0.
	!cutoff_inf = Value for infinity i.e. the [2q,infinity] integral. Used on all the others.
	real(8), parameter, public :: cutoff0_inf = 20000.0d0
	real(8), parameter, public :: cutoff_inf = 2000.0d0

end module
