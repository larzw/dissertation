program main
!Calculate the effective interaction.

use gauss_legendre
use bonn
use ab
use x
implicit none

	integer :: c1
	real(8), dimension(ss,2) :: xi0_k,xi_k

	complex, dimension(ss) :: B0
	complex, dimension(ss) :: B12,B6
	complex, dimension(ss) :: B3,B4,B5

	call T0(B0,q,T,cutoff0_inf)
	call T_12_6(B12,B6,q,T)
	call T_3_4_5(B3,B4,B5,q,T)

	open(unit = 10, file = 'calc/g0_0.d')
	open(unit = 20, file = 'calc/g12_0.d')
	open(unit = 30, file = 'calc/g6_0.d')
	open(unit = 40, file = 'calc/g3_0.d')
	open(unit = 50, file = 'calc/g4_0.d')
	open(unit = 60, file = 'calc/g5_0.d')

! 	open(unit = 10, file = 'calc/g0_1.d')
! 	open(unit = 20, file = 'calc/g12_1.d')
! 	open(unit = 30, file = 'calc/g6_1.d')
! 	open(unit = 40, file = 'calc/g3_1.d')
! 	open(unit = 50, file = 'calc/g4_1.d')
! 	open(unit = 60, file = 'calc/g5_1.d')


	call xi(q,xi0_k,cutoff0_inf)
	call xi(q,xi_k,cutoff_inf)

	!^0g
	do c1 = 1,ss

		write(10,*) real( B0(c1) ),aimag( B0(c1) ),xi0_k(c1,1),xi0_k(c1,2)

	enddo

	!^12g,g_6
	do c1 = 1,ss

		write(20,*) real( B12(c1) ),aimag( B12(c1) ),xi_k(c1,1),xi_k(c1,2)
		write(30,*) real( B6(c1) ),aimag( B6(c1) ),xi_k(c1,1),xi_k(c1,2)

	enddo

	!g_3,g_4,g_5
	do c1 = 1,ss

		write(40,*) real( B3(c1) ),aimag( B3(c1) ),xi_k(c1,1),xi_k(c1,2)
		write(50,*) real( B4(c1) ),aimag( B4(c1) ),xi_k(c1,1),xi_k(c1,2)
		write(60,*) real( B5(c1) ),aimag( B5(c1) ),xi_k(c1,1),xi_k(c1,2)

	enddo


	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
	close(60)

end program