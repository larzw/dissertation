module ab
!Setup A and b i.e. we are solving Ax=b

use gauss_legendre
use inputs
use bonn
implicit none

	!n_k = Number of nodes and weights for [0,inf] integral
	integer, parameter, public :: n_k = n_0_2q + n_2q_inf

	!ss = Length of xi vector/function.
	integer, parameter, public :: ss = n_theta_k*(1 + n_k)

	!m_avg = Used in the denomenator of the scatering equation.
	!i.e. E_q - E_k + i*epsilon ---> sqrt( m_avg^2 + q^2 ) - (m_avg^2 + k^2) + i*epsilon
	real(8), parameter, public :: m_avg = ( m_1 + m_2 )/2.0d0

contains

real(8) function V_pip(qp,thetap,qpp,thetapp,lambda1p,lambda2p,lambda1pp,lambda2pp,LAMBDA,T)
!Computes the phi integrated potential (FEG method) over the interval [0,2pi].

implicit none

	!T = Isospin
	!(qp,thetap), (qpp,thetapp) = vec{q},vec{q'}.
	!lambda1p,lambda2p,lambda1pp,lambda2pp = < | >,< |sigma^(1)*sigma^(2)| >
	!LAMBDA = 0,1,-1 value that appears in FEG method.
	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,lambda1p,lambda2p,lambda1pp,lambda2pp,LAMBDA

	!s = Summing variable in the integral approximation
	!c1 = Counter
	!x,w = Nodes and weights for quadrature
	complex(8) :: s
	integer :: c1
	real(8), dimension(n_V_avg) :: x,w
		
	call gauss_legendre_quadrature(n_V_avg,0.0d0,2.0d0*pi,x,w)

	s = (0.0d0,0.0d0)
	do c1 = 1,n_V_avg

		!phi = 0 in the potential
		s = s + w(c1)*cdexp( i*LAMBDA*x(c1) )*V(qp,thetap,0.0d0,qpp,thetapp,x(c1),lambda1p,lambda2p,lambda1pp,lambda2pp,T)

	enddo

	!The imaginary part is zero.
	V_pip = ( 1.0d0/(2.0d0*pi) )*real(s)

end function

real(8) function V_1(qp,thetap,q,T)
!<++|V|++>

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V_1 = V(qp,thetap,0.0d0,q,0.0d0,0.0d0,0.5d0,0.5d0,0.5d0,0.5d0,T)

end function

real(8) function VLAMBDA_1(qp,thetap,qpp,thetapp,LAMBDA,T)
!<++|V^Lambda|++>

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA_1 = V_pip(qp,thetap,qpp,thetapp,0.5d0,0.5d0,0.5d0,0.5d0,LAMBDA,T)

end function

real(8) function V_2(qp,thetap,q,T)
!<++|V|-->

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V_2 = V(qp,thetap,0.0d0,q,0.0d0,0.0d0,0.5d0,0.5d0,-0.5d0,-0.5d0,T)

end function

real(8) function VLAMBDA_2(qp,thetap,qpp,thetapp,LAMBDA,T)
!<++|V^Lambda|-->

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA_2 = V_pip(qp,thetap,qpp,thetapp,0.5d0,0.5d0,-0.5d0,-0.5d0,LAMBDA,T)

end function

real(8) function V_3(qp,thetap,q,T)
!<+-|V|+->

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V_3 = V(qp,thetap,0.0d0,q,0.0d0,0.0d0,0.5d0,-0.5d0,0.5d0,-0.5d0,T)

end function

real(8) function VLAMBDA_3(qp,thetap,qpp,thetapp,LAMBDA,T)
!<+-|V^Lambda|+->

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA_3 = V_pip(qp,thetap,qpp,thetapp,0.5d0,-0.5d0,0.5d0,-0.5d0,LAMBDA,T)

end function

real(8) function V_4(qp,thetap,q,T)
!<+-|V|-+>

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V_4 = V(qp,thetap,0.0d0,q,0.0d0,0.0d0,0.5d0,-0.5d0,-0.5d0,0.5d0,T)

end function

real(8) function VLAMBDA_4(qp,thetap,qpp,thetapp,LAMBDA,T)
!<+-|V^Lambda|-+>

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA_4 = V_pip(qp,thetap,qpp,thetapp,0.5d0,-0.5d0,-0.5d0,0.5d0,LAMBDA,T)

end function

real(8) function V_5(qp,thetap,q,T)
!<++|V|+->

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V_5 = V(qp,thetap,0.0d0,q,0.0d0,0.0d0,0.5d0,0.5d0,0.5d0,-0.5d0,T)

end function

real(8) function VLAMBDA_5(qp,thetap,qpp,thetapp,LAMBDA,T)
!<++|V^Lambda|+->

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA_5 = V_pip(qp,thetap,qpp,thetapp,0.5d0,0.5d0,0.5d0,-0.5d0,LAMBDA,T)

end function

real(8) function V_6(qp,thetap,q,T)
!<+-|V|++>

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V_6 = V(qp,thetap,0.0d0,q,0.0d0,0.0d0,0.5d0,-0.5d0,0.5d0,0.5d0,T)

end function

real(8) function VLAMBDA_6(qp,thetap,qpp,thetapp,LAMBDA,T)
!<+-|V^Lambda|++>

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA_6 = V_pip(qp,thetap,qpp,thetapp,0.5d0,-0.5d0,0.5d0,0.5d0,LAMBDA,T)

end function

real(8) function V0(qp,thetap,q,T)
!^{0}V

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V0 = V_1(qp,thetap,q,T) - V_2(qp,thetap,q,T)

end function

real(8) function VLAMBDA0(qp,thetap,qpp,thetapp,LAMBDA,T)
!^{0}V^Lambda

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA0 = VLAMBDA_1(qp,thetap,qpp,thetapp,LAMBDA,T) - VLAMBDA_2(qp,thetap,qpp,thetapp,LAMBDA,T)

end function

real(8) function V12(qp,thetap,q,T)
!^{12}V

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,q

	V12 = V_1(qp,thetap,q,T) + V_2(qp,thetap,q,T)

end function

real(8) function VLAMBDA12(qp,thetap,qpp,thetapp,LAMBDA,T)
!^{12}V^Lambda

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA12 = VLAMBDA_1(qp,thetap,qpp,thetapp,LAMBDA,T) + VLAMBDA_2(qp,thetap,qpp,thetapp,LAMBDA,T)

end function

real(8) function VLAMBDA1(qp,thetap,qpp,thetapp,LAMBDA,T)
!^{1}V^Lambda
!Note: This is used purly for convience. I'm using Rup's linear combinations.

implicit none

	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA

	VLAMBDA1 = VLAMBDA_3(qp,thetap,qpp,thetapp,LAMBDA,T) - VLAMBDA_4(qp,thetap,qpp,thetapp,LAMBDA,T)

end function

subroutine principle_value_quadrature(q,xx,ww,cutoff_inf)
!Generates quadradure rule to preform the priciple value integral i.e.,
!pv int[f(x)/(x-q),{x,0,inf}] with a singularity at x = q.
!
!To preform the principle value integral,
!pv int[f(x)/(x-q),{x,0,inf}] = int[f(x)/(x-q),{x,0,2q}] + int[f(x)/(x-q),{x,2q,inf}]
!Then employ Gauss-Legendre quadrature normally with the exception that
!YOU MAKE SURE TO CHOOSE AN EVEN NUMBER OF PONTS FOR THE [0,2q] INTEGRAL.	
	
implicit none

	!q = |q|
	!cutoff_inf = Cutoff for [2q,inf] integral.
	!xx,ww = Nodes/Weights	
	real(8), intent(in) :: q,cutoff_inf
	real(8), dimension(n_k), intent(out) :: xx,ww
		
	!x_1,w_1 = Nodes/weights for [0,2q] integral.
	!x_2,w_2 = Nodes/weights for [2q,inf] integral.
	real(8), dimension(n_0_2q) :: x_1,w_1
	real(8), dimension(n_2q_inf) :: x_2,w_2
		
	call gauss_legendre_quadrature(n_0_2q,0.0d0,2.0d0*q,x_1,w_1)
	call gauss_legendre_quadrature(n_2q_inf,2.0d0*q,cutoff_inf,x_2,w_2)
		
	xx(1:n_0_2q) = x_1
	xx(n_0_2q + 1:n_k) = x_2
	ww(1:n_0_2q) = w_1
	ww(n_0_2q + 1:n_k) = w_2

end subroutine

subroutine xi(q,xi_k,cutoff_inf)
!Values of vec{q}' to evaluate at.

implicit none

	!q = |q|
	!cutoff_inf = Cutoff for [2q,inf] integral.
	!xi_k = Values of vec{q}' to evaluate at.
	real(8), intent(in) :: q,cutoff_inf
	real(8), dimension(ss,2), intent(out) :: xi_k

	!c1,c2 = Counters
	!p,j = Variables in Eq. 18 of T matrix report.
	!x_theta_k, w_theta_k,x_k,w_k = Nodes and weights for the theta_k, and k integral.
	integer :: c1,c2,p,j
	real(8), dimension(n_theta_k) :: x_theta_k,w_theta_k
	real(8), dimension(n_k) :: x_k,w_k

	call gauss_legendre_quadrature(n_theta_k,0.0d0,pi,x_theta_k,w_theta_k)
	call principle_value_quadrature(q,x_k,w_k,cutoff_inf)

	!c_j(q) vector see Eq. 19 of T matrix report.
	do c1 = 1,n_theta_k

		xi_k(c1,1) = q
		xi_k(c1,2) = x_theta_k(c1)

	enddo

	!b_l vector see Eq. 19 of T matrix report.
	do c2 = 1,n_k*n_theta_k

		p = floor(  ( dble(c2-1)/dble(n_theta_k) ) + 1.0d0  ) 	
		j = c2 - n_theta_k*(p - 1)

		xi_k(c2 + n_theta_k,1) = x_k(p)
		xi_k(c2 + n_theta_k,2) = x_theta_k(j)

	enddo

end subroutine

real(8) function pauli(qpp,thetapp)
!The Pauli operator i.e. Q(q",theta",P,kF1,kF2)
!Notes: Pcm, k_fermi_1, k_fermi_2, and Q_pauli are globaly defined variables found in input.f90

implicit none

	!qpp,thetapp = vec{q}".
	real(8), intent(in) :: qpp,thetapp

	!P_plus_qpp and P_minus_qpp = |P + q"| and |P - q"|
	real(8) :: P_plus_qpp,P_minus_qpp

	P_plus_qpp = dsqrt(  (qpp**2.0d0) + (Pcm**2.0d0) + 2.0d0*Pcm*qpp*dcos(thetapp)  )
	P_minus_qpp = dsqrt(  (qpp**2.0d0) + (Pcm**2.0d0) - 2.0d0*Pcm*qpp*dcos(thetapp)  )

	!I have to convert kF[fm^(-1)] to MeV by multiplying by hbar*c.
	if( (P_plus_qpp > (k_fermi_1*197.327d0)) .and. (P_minus_qpp > (k_fermi_2*197.327d0)) ) then

		pauli = 1.0d0

	else

		pauli = Q_pauli

	endif

end function

complex(8) function alpha(k,j,q,LAMBDA,T,Vn,cutoff_inf)
!Computes alpha

implicit none

	interface

		real(8) function Vn(qp,thetap,qpp,thetapp,LAMBDA,T)
			
				integer, intent(in) :: T
				real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA
			
		end function

	end interface

	!k,j = Row,Column of matrix.
	!T = Isospin.
	!q =|q|.
	!LAMBDA = 0,1,-1 value that appears in FEG method.
	!cutoff_inf = Cutoff for [2q,inf] integral.
	integer, intent(in) :: k,j,T
	real(8), intent(in) :: q,LAMBDA,cutoff_inf

	!E_q = I.e. E_q = sqrt( q^2 + m_avg^2 ).
	!x_theta_k, w_theta_k = Nodes and weights for the theta_k integral.
	!xi_k = Values to evaluate vec{q}' at.
	real(8) :: E_q
	real(8), dimension(n_theta_k) :: x_theta_k,w_theta_k
	real(8), dimension(ss,2):: xi_k

	call xi(q,xi_k,cutoff_inf)
	call gauss_legendre_quadrature(n_theta_k,0.0d0,pi,x_theta_k,w_theta_k)

	!Use "m_avg" i.e. not m or m_1 or m_2.
	E_q = dsqrt( q**2.0d0 + m_avg**2.0d0 )

	alpha = -i*(pi**2.0d0)*q*E_q*w_theta_k(j)*dsin( x_theta_k(j) )*Vn(xi_k(k,1),xi_k(k,2),q,x_theta_k(j),LAMBDA,T)&
			&*pauli(q,x_theta_k(j))

end function

real(8) function beta(k,L,q,LAMBDA,T,Vn,cutoff_inf)
!Computes beta

implicit none

	interface

		real(8) function Vn(qp,thetap,qpp,thetapp,LAMBDA,T)

			integer, intent(in) :: T
			real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA
			
		end function

	end interface

	!k,L = Row,Column of matrix.
	!T = Isospin.
	!q =|q|.
	!LAMBDA = 0,1,-1 value that appears in FEG method.
	!cutoff_inf = Cutoff for [2q,inf] integral.
	integer, intent(in) :: k,L,T
	real(8), intent(in) :: q,LAMBDA,cutoff_inf

	!p,j = Variables in Eq. 18 of T matrix report.
	!E_q,E_yp = I.e. E_q = sqrt(m_avg^2 + q^2).
	!temp1 = Temporary variable.
	!x_theta_k, w_theta_k,x_k,w_k = Nodes and weights for the theta_k, and k integral.
	integer :: p,j
	real(8) :: E_q,E_yp,temp1
	real(8), dimension(n_theta_k) :: x_theta_k,w_theta_k
	real(8), dimension(n_k) :: x_k,w_k
	real(8), dimension(ss,2):: xi_k

	call xi(q,xi_k,cutoff_inf)

	call gauss_legendre_quadrature(n_theta_k,0.0d0,pi,x_theta_k,w_theta_k)
	call principle_value_quadrature(q,x_k,w_k,cutoff_inf)

	p = floor(  ( dble(L-1)/dble(n_theta_k) ) + 1.0d0  ) 	
	j = L - n_theta_k*(p - 1)

	!Use "m_avg" i.e. not m or m_1 or m_2.
	E_q = dsqrt( q**2.0d0 + m_avg**2.0d0 )
	E_yp = dsqrt(  ( x_k(p) )**2.0d0 + m_avg**2.0d0  )

	temp1 = ( x_k(p)**2.0d0 )*dsin( x_theta_k(j) )
	beta = pi*w_k(p)*w_theta_k(j)*( temp1/(E_q - E_yp) )*Vn(xi_k(k,1),xi_k(k,2),x_k(p),x_theta_k(j),LAMBDA,T)&
			&*pauli(x_k(p),x_theta_k(j))

end function

subroutine build_A(A,q,LAMBDA,T,Vn,cutoff_inf)
!Computes matrix A i.e. Ax=b

implicit none

	interface

		real(8) function Vn(qp,thetap,qpp,thetapp,LAMBDA,T)

			integer, intent(in) :: T
			real(8), intent(in) :: qp,thetap,qpp,thetapp,LAMBDA
			
		end function

	end interface

	!T = Isospin.
	!q = |q|.
	!LAMBDA = 0,1,-1 value that appears in FEG method.
	!cutoff_inf = Cutoff for [2q,inf] integral.
	!A = Matrix representing the kernal.
	integer, intent(in) :: T
	real(8), intent(in) :: q,LAMBDA,cutoff_inf 
	complex, dimension(ss,ss), intent(out) :: A

	!c1,c2,c3,c4 = Counters
	integer :: c1,c2,c3,c4

	!Put column on the outside of loop. Since this is more efficient in fortran.

!$omp parallel
!$omp do

	!Alpha block
	do c2 = 1,n_theta_k

		do c1 = 1,ss

			!They are automatically typecast
			A(c1,c2) = alpha(c1,c2,q,LAMBDA,T,Vn,cutoff_inf)

		enddo

	enddo

!$omp end do
!$omp do

	!Beta block
	do c4 = 1,(n_k*n_theta_k)

		do c3 = 1,ss

			!They are automatically typecast
			A(c3,n_theta_k + c4) = beta(c3,c4,q,LAMBDA,T,Vn,cutoff_inf)

		enddo

	enddo

!$omp end do
!$omp end parallel	

end subroutine

subroutine build_V(V,q,T,Vn,cutoff_inf)
!Builds the potential vector. i.e. the b in Ax=b

implicit none

	interface

		real(8) function Vn(qp,thetap,q,T)

			integer, intent(in) :: T
			real(8), intent(in) :: qp,thetap,q
			
		end function

	end interface

	!q = |q|
	!cutoff_inf = Cutoff for [2q,inf] integral.
	!T = Isospin
	!V = Potential vector
	real(8), intent(in) :: q,cutoff_inf
	integer, intent(in) :: T
	complex, dimension(ss), intent(out) :: V

	!c1 = Counter
	!xi_k = Where to evaluate vec{q}' at.
	integer :: c1
	real(8), dimension(ss,2) :: xi_k

	call xi(q,xi_k,cutoff_inf)

	do c1 = 1,ss

		V(c1) = Vn(xi_k(c1,1),xi_k(c1,2),q,T)

	enddo

end subroutine

end module
