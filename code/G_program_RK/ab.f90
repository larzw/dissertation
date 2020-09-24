module ab
!Setup A and b i.e. we are solving Ax=b

use gauss_legendre
use inputs
use bonn
implicit none

	!n_k = Number of y_p/W_p nodes. See Eq.15 of T matrix rev. report.
	integer, parameter, public :: n_k = n_0_2q + n_2q_inf

	!ss = Length of xi vector.
	integer, parameter, public :: ss = n_theta_k*(1 + n_k)

	!m_avg = Used in the denomenator of the scatering equation.
	!	i.e. E_q - E_k + i*epsilon ---> sqrt( m_avg^2 + q^2 ) - (m_avg^2 + k^2) + i*epsilon
	real(8), parameter, public :: m_avg = ( m_1 + m_2 )/2.0d0

contains

real(8) function V_pip(q,theta,qp,thetap,lambda1p,lambda2p,lambda1,lambda2,T)
!Computes the phi integrated potential over the interval [0,2*pi].
!See R.A. Rice and Y.E. Kim, Few-Body Systems 14,127-148 (1993)

implicit none

	!T = Isospin
	!(q,theta,phi),(qp,thetap,phip) = vec{q},vec{q'} momentum variables.
	!lambda1p,lambda2p,lambda1,lambda2 = < | >,< |sigma^(1)*sigma^(2)| >
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap,lambda1p,lambda2p,lambda1,lambda2

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
		s = s + w(c1)*V(q,theta,0.0d0,qp,thetap,x(c1),lambda1p,lambda2p,lambda1,lambda2,T)

	enddo

	!The imaginary part is zero.
	V_pip = ( 1.0d0/(2.0d0*pi) )*real(s)

end function

real(8) function V_1(q,theta,qp,thetap,T)
!<++|V|++>

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V_1 = V_pip(q,theta,qp,thetap,0.5d0,0.5d0,0.5d0,0.5d0,T)

end function

real(8) function V_2(q,theta,qp,thetap,T)
!<++|V|-->

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V_2 = V_pip(q,theta,qp,thetap,0.5d0,0.5d0,-0.5d0,-0.5d0,T)

end function

real(8) function V_3(q,theta,qp,thetap,T)
!<+-|V|+->

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V_3 = V_pip(q,theta,qp,thetap,0.5d0,-0.5d0,0.5d0,-0.5d0,T)

end function

real(8) function V_4(q,theta,qp,thetap,T)
!<+-|V|-+>

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V_4 = V_pip(q,theta,qp,thetap,0.5d0,-0.5d0,-0.5d0,0.5d0,T)

end function

real(8) function V_5(q,theta,qp,thetap,T)
!<++|V|+->

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V_5 = V_pip(q,theta,qp,thetap,0.5d0,0.5d0,0.5d0,-0.5d0,T)

end function

real(8) function V_6(q,theta,qp,thetap,T)
!<+-|V|++>

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V_6 = V_pip(q,theta,qp,thetap,0.5d0,-0.5d0,0.5d0,0.5d0,T)

end function

real(8) function V0(q,theta,qp,thetap,T)
!^{0}V

implicit none

	!See V_phi_int()
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V0 = V_1(q,theta,qp,thetap,T) - V_2(q,theta,qp,thetap,T)

end function

real(8) function V1(q,theta,qp,thetap,T)
!^{1}V

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V1 = V_3(q,theta,qp,thetap,T) + V_4(q,theta,qp,thetap,T)

end function

real(8) function V12(q,theta,qp,thetap,T)
!^{12}V

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V12 = V_1(q,theta,qp,thetap,T) + V_2(q,theta,qp,thetap,T)

end function

real(8) function V34(q,theta,qp,thetap,T)
!^{34}V

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V34 = V_3(q,theta,qp,thetap,T) - V_4(q,theta,qp,thetap,T)

end function

real(8) function V55(q,theta,qp,thetap,T)
!^{55}V

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V55 = 2.0d0*V_5(q,theta,qp,thetap,T)

end function

real(8) function V66(q,theta,qp,thetap,T)
!^{66}V

implicit none

	!See V_pip
	integer, intent(in) :: T
	real(8), intent(in) :: q,theta,qp,thetap

	V66 = 2.0d0*V_6(q,theta,qp,thetap,T)

end function

subroutine principle_value_quadrature(q,xx,ww)
!Generates quadradure rule to preform the priciple value integral i.e.:
!	pv int[f(x)/(x-q),{x,0,inf}] where f(x) has a singularity at x=q.
!To preform the principle value integral:
!	pv int[f(x)/(x-q),{x,0,inf}] = int[f(x)/(x-q),{x,0,2q}] + int[f(x)/(x-q),{x,2q,inf}]
!	Then employ Gauss-Legendre quadrature normally with the exception that
!	YOU MAKE SURE TO CHOOSE AN EVEN NUMBER OF PONTS FOR THE [0,2q] INTEGRAL.	
	
implicit none

	!q = Momentum variable
	!xx,ww = Nodes/Weights	
	real(8), intent(in) :: q
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

subroutine xi(q,xi_k)
!Produces the xi_k(q) vector found in Eq. 19 of T matrix report.

implicit none

	!q = Momentum variable
	!xi_k = Xi vector appearing in Eq. 19 of T matirx report.
	real(8), intent(in) :: q
	real(8), dimension(ss,2), intent(out) :: xi_k

	!c1,c2 = Counters
	!p,j = Variables in Eq. 18 of T matrix report.
	!x_theta_k, w_theta_k,x_k,w_k = Nodes and weights for the theta_k, and k integral.
	integer :: c1,c2,p,j
	real(8), dimension(n_theta_k) :: x_theta_k,w_theta_k
	real(8), dimension(n_k) :: x_k,w_k

	call gauss_legendre_quadrature(n_theta_k,0.0d0,pi,x_theta_k,w_theta_k)
	call principle_value_quadrature(q,x_k,w_k)

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
!The Pauli operator i.e. Q(q",theta",P,kF)
!
!Notes:
!	Pcm, k_fermi, and Q_pauli are globaly defined variables found in input.f90

implicit none

	!qpp,thetapp = Momentum variables. I sometimes refere to vec{k} as vec{q}". 
	real(8), intent(in) :: qpp,thetapp

	!Temporary variables
	real(8) :: temp1,temp2

	temp1 = dsqrt(  (qpp**2.0d0) + (Pcm**2.0d0) + 2.0d0*Pcm*qpp*dcos(thetapp)  )
	temp2 = dsqrt(  (qpp**2.0d0) + (Pcm**2.0d0) - 2.0d0*Pcm*qpp*dcos(thetapp)  )

	!I have to convert kF[fm^(-1)] to MeV by multiplying by hbar*c.
	if( (temp1 > (k_fermi*197.327d0)) .and. (temp2 > (k_fermi*197.327d0)) ) then

		pauli = 1.0d0

	else

		pauli = Q_pauli

	endif

end function

complex(8) function alpha(k,j,q,T,Vn)
!Computes alpha

implicit none

	interface

		real(8) function Vn(q,theta,qp,thetap,T)
			
			integer, intent(in) :: T
			real(8), intent(in) :: q,theta,qp,thetap
			
		end function

	end interface

	!k,j = Row,Column of matrix.
	!T = Isospin.
	!q = Momentum variable.
	integer, intent(in) :: k,j,T
	real(8), intent(in) :: q

	!E_q = I.e. E_q = sqrt( q^2 + m_avg^2 ).
	!x_theta_k, w_theta_k = Nodes and weights for the theta_k integral.
	!	See Eq. 3 of T matrix report.
	!xi_k = Values to evaluate vec{q'} at. See Eq. 20 of T matrix rev. report.
	real(8) :: E_q
	real(8), dimension(n_theta_k) :: x_theta_k,w_theta_k
	real(8), dimension(ss,2):: xi_k

	call xi(q,xi_k)

	call gauss_legendre_quadrature(n_theta_k,0.0d0,pi,x_theta_k,w_theta_k)

	!Use "m_avg" i.e. not m or m_1 or m_2.
	E_q = dsqrt( q**2.0d0 + m_avg**2.0d0 )

	alpha = -i*(pi**2.0d0)*q*E_q*w_theta_k(j)*dsin( x_theta_k(j) )*Vn(q,x_theta_k(j),xi_k(k,1),xi_k(k,2),T)*pauli(q,x_theta_k(j))

end function

real(8) function beta(k,L,q,T,Vn)
!Computes beta

implicit none

	interface

		real(8) function Vn(q,theta,qp,thetap,T)

			integer, intent(in) :: T
			real(8), intent(in) :: q,theta,qp,thetap
			
		end function

	end interface

	!k,L = Row,Column of matrix.
	!T = Isospin.
	!q = Momentum variable.
	integer, intent(in) :: k,L,T
	real(8), intent(in) :: q

	!p,j = Variables in Eq. 18 of T matrix report.
	!E_q,E_yp = I.e. E_q = sqrt(m_avg^2 + q^2).
	!temp1 = Temporary variable.
	!x_theta_k, w_theta_k,x_k,w_k = Nodes and weights for the theta_k, and k integral.
	!	See Eq. 3 of T matrix report.
	integer :: p,j
	real(8) :: E_q,E_yp,temp1
	real(8), dimension(n_theta_k) :: x_theta_k,w_theta_k
	real(8), dimension(n_k) :: x_k,w_k
	real(8), dimension(ss,2):: xi_k

	call xi(q,xi_k)

	call gauss_legendre_quadrature(n_theta_k,0.0d0,pi,x_theta_k,w_theta_k)
	call principle_value_quadrature(q,x_k,w_k)

	p = floor(  ( dble(L-1)/dble(n_theta_k) ) + 1.0d0  ) 	
	j = L - n_theta_k*(p - 1)

	!Use "m_avg" i.e. not m or m_1 or m_2.
	E_q = dsqrt( q**2.0d0 + m_avg**2.0d0 )
	E_yp = dsqrt(  ( x_k(p) )**2.0d0 + m_avg**2.0d0  )

	temp1 = ( x_k(p)**2.0d0 )*dsin( x_theta_k(j) )
	beta = pi*w_k(p)*w_theta_k(j)*( temp1/(E_q - E_yp) )*Vn(x_k(p),x_theta_k(j),xi_k(k,1),xi_k(k,2),T)*pauli(x_k(p),x_theta_k(j))

end function

subroutine build_A(A,q,T,Vn)
!Computes matrix A i.e. Ax=b

implicit none

	interface

		real(8) function Vn(q,theta,qp,thetap,T)

			integer, intent(in) :: T
			real(8), intent(in) :: q,theta,qp,thetap
			
		end function

	end interface

	!T = Isospin.
	!q = Momentum variable.
	!A = Matrix in Eq. 23 of T matrix report.
	integer, intent(in) :: T
	real(8), intent(in) :: q
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
			A(c1,c2) = alpha(c1,c2,q,T,Vn)

		enddo

	enddo

!$omp end do
!$omp do

	!Beta block
	do c4 = 1,(n_k*n_theta_k)

		do c3 = 1,ss

			!They are automatically typecast
			A(c3,n_theta_k + c4) = beta(c3,c4,q,T,Vn)

		enddo

	enddo

!$omp end do
!$omp end parallel	

end subroutine

subroutine build_V(V,q,T,Vn)
!Builds the potential vector. i.e. the b in Ax=b

implicit none

	interface

		real(8) function Vn(q,theta,qp,thetap,T)

			integer, intent(in) :: T
			real(8), intent(in) :: q,theta,qp,thetap
			
		end function

	end interface

	!q= Momentum variable
	!T = Isospin
	!V = Potential vector
	real(8), intent(in) :: q
	integer, intent(in) :: T
	complex, dimension(ss,n_theta), intent(out) :: V

	!c1 = Counter
	!xi_k = See Eq. 19 of T matrix report.
	!x,w = Nodes/Weights
	integer :: c1,c2
	real(8), dimension(ss,2) :: xi_k
	real(8), dimension(n_theta) :: x,w

	call xi(q,xi_k)
	call gauss_legendre_quadrature(n_theta,0.0d0,pi,x,w)

	!Put column on the outside of loop. Since this is more efficient in fortran.
	do c2 = 1,n_theta

		do c1 = 1,ss

			!They are automatically typecast
			V(c1,c2) = Vn(q,x(c2),xi_k(c1,1),xi_k(c1,2),T)

		enddo

	enddo

end subroutine

end module
