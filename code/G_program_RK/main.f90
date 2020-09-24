program main
!Calculate the effective interaction.

use gauss_legendre
use bonn
use ab
use x
implicit none

	integer :: c1,k1,c2,k2,c3,k3,c4,k4,c5,k5,c6,k6
	real(8), dimension(ss,2) :: xi_k
	real(8), dimension(n_theta) :: xx,ww

	complex, dimension(ss,n_theta) :: B0
	complex, dimension(ss,n_theta) :: B1
	complex, dimension(ss,n_theta) :: B66,B55,B34,B12
	call T0(B0,q,T)
	call T1(B1,q,T)
	call Tnn(B66,B55,B34,B12,q,T)

	call xi(q,xi_k)
	call gauss_legendre_quadrature(n_theta,0.0d0,pi,xx,ww)

	open(unit = 10, file = 'calc/g0_1.d')
	open(unit = 20, file = 'calc/g1_1.d')
	open(unit = 30, file = 'calc/g66_1.d')
	open(unit = 40, file = 'calc/g55_1.d')
	open(unit = 50, file = 'calc/g34_1.d')
	open(unit = 60, file = 'calc/g12_1.d')

	!^0G
	do k1 = 1,n_theta

		do c1 = 1,ss

			write(10,*) real(B0(c1,k1)),aimag(B0(c1,k1)),xi_k(c1,1),xi_k(c1,2),xx(k1)

		enddo

	enddo

	!^1G
	do k2 = 1,n_theta

		do c2 = 1,ss

			write(20,*) real(B1(c2,k2)),aimag(B1(c2,k2)),xi_k(c2,1),xi_k(c2,2),xx(k2)

		enddo

	enddo

	!^66G
	do k3 = 1,n_theta

		do c3 = 1,ss

			write(30,*) real(B66(c3,k3)),aimag(B66(c3,k3)),xi_k(c3,1),xi_k(c3,2),xx(k3)

		enddo

	enddo

	!^55G
	do k4 = 1,n_theta

		do c4 = 1,ss

			write(40,*) real(B55(c4,k4)),aimag(B55(c4,k4)),xi_k(c4,1),xi_k(c4,2),xx(k4)

		enddo

	enddo

	!^34G
	do k5 = 1,n_theta

		do c5 = 1,ss

			write(50,*) real(B34(c5,k5)),aimag(B34(c5,k5)),xi_k(c5,1),xi_k(c5,2),xx(k5)

		enddo

	enddo

	!^12G
	do k6 = 1,n_theta

		do c6 = 1,ss

			write(60,*) real(B12(c6,k6)),aimag(B12(c6,k6)),xi_k(c6,1),xi_k(c6,2),xx(k6)

		enddo

	enddo

	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
	close(60)

end program
