module x
!Solve for x i.e. we are solving x = ( A^(-1) )*b

use inputs
use ab
implicit none

contains

subroutine T0(B0,q,T)
!^{0}T

implicit none

    !Lapack parameters.
    integer, parameter :: n = ss
    integer, parameter :: nrhs = n_theta
    integer, parameter :: lda = n
    integer, parameter :: ldb = n
    integer, parameter :: lwork = n
    integer :: info
    real, dimension(lwork) :: work
    integer, dimension(n) :: ipiv

	!q = Momentum variable.
	!T = Isospin
	!B0 = Right hand side of Ax=b. Then it gets replaced by the solution vector ^{0}T.
	real(8), intent(in) :: q
	integer, intent(in) :: T
	complex, dimension(ldb,nrhs), intent(out) :: B0

	!c1 = Counter
	!A0 = The matrix in the equation Ax=b.
	integer :: c1
	complex, dimension(:,:), allocatable :: A0

	allocate( A0(lda,n) )

	!Build the 0^A matrix
	call build_A(A0,q,T,V0)

	!Compute 1 - ^{0}A i.e. "^{0}B" in T matrix report.
	A0 = -A0
	do c1 = 1,n

		A0(c1,c1) = 1.0 + A0(c1,c1)

	enddo

	!Build right hand side vector in the equation Ax=b.
	call build_V(B0,q,T,V0)

	!Call to lapack to solve Ax=b using LU-factorization. The solution is written into B0.
	call cgesv(n,nrhs,A0,lda,ipiv,B0,ldb,info)

	deallocate(A0)

end subroutine

subroutine T1(B1,q,T)
!^{1}T

implicit none

	!Lapack parameters.
    integer, parameter :: n = ss
    integer, parameter :: nrhs = n_theta
    integer, parameter :: lda = n
    integer, parameter :: ldb = n
    integer, parameter :: lwork = n
    integer :: info
    real, dimension(lwork) :: work
    integer, dimension(n) :: ipiv

	!q = Momentum variable.
	!T = Isospin
	!B1 = Right hand side of Ax=b. Then it gets replaced by the solution vector ^{1}T.
	real(8), intent(in) :: q
	integer, intent(in) :: T
	complex, dimension(ldb,nrhs), intent(out) :: B1

	!c1 = Counter
	!A0 = The matrix in the equation Ax=b.
	integer :: c1
	complex, dimension(:,:), allocatable :: A1

	allocate( A1(lda,n) )

	!Build the 1^A matrix
	call build_A(A1,q,T,V1)

	!Compute 1 - ^{1}A i.e. "^{1}B" in T matrix report.
	A1 = -A1
	do c1 = 1,n

		A1(c1,c1) = 1.0 + A1(c1,c1)

	enddo

	!Build right hand side vector in the equation Ax=b.
	call build_V(B1,q,T,V1)

	!Call to lapack to solve Ax=b using LU-factorization. The solution is written into B1.
	call cgesv(n,nrhs,A1,lda,ipiv,B1,ldb,info)

	deallocate(A1)

end subroutine

subroutine Tnn(B66,B55,B34,B12,q,T)
!^{66}T,^{55}T,^{34}T,^{12}T

implicit none

  	!Lapack parameters. 
    integer, parameter :: n = 2*ss
    integer, parameter :: nrhs = 2*n_theta
    integer, parameter :: lda = n
    integer, parameter :: ldb = n
    integer, parameter :: lwork = n
    integer :: info
    real, dimension(lwork) :: work
    integer, dimension(n) :: ipiv

	!q = Momentum variable.
	!T = Isospin
	!B66...B12 = Solution vector ^{66}T...^{12}T.
	real(8), intent(in) :: q
	integer, intent(in) :: T
	complex, dimension(ss,n_theta), intent(out) :: B66,B55,B34,B12

	!c1,c2 = Counters
	!Ann = The A matrix in the equation Ax=b.
	!Vnn = The right hand side (b matrix) in the equation Ax=b.
	!A11...A22 = Matrix entries in Ann.
	!VV11...VV22 = Matrix entries in Vnn.
	integer :: c1,c2
	complex, dimension(:,:), allocatable :: Ann
	complex, dimension(:,:), allocatable :: Vnn
	complex, dimension(:,:), allocatable :: A11
	complex, dimension(:,:), allocatable :: A12
	complex, dimension(:,:), allocatable :: A21
	complex, dimension(:,:), allocatable :: A22
	complex, dimension(:,:), allocatable :: VV11
    complex, dimension(:,:), allocatable :: VV12
    complex, dimension(:,:), allocatable :: VV21
    complex, dimension(:,:), allocatable :: VV22

	allocate( Ann(lda,n) )
	allocate( Vnn(ldb,nrhs) )

	!Build A11 and put it in Ann
	allocate( A11(ss,ss) )
	call build_A(A11,q,T,V12)
	A11 = -A11
	do c1 = 1,ss

		A11(c1,c1) = 1.0 + A11(c1,c1)

	enddo
	Ann(1:ss,1:ss) = A11
	deallocate(A11)

	!Build A22 and put it in Ann
    allocate( A22(ss,ss) )
    call build_A(A22,q,T,V34)
    A22 = -A22
    do c2 = 1,ss

        A22(c2,c2) = 1.0 + A22(c2,c2)

    enddo
	Ann( (ss + 1):2*ss, (ss + 1):2*ss ) = A22
    deallocate(A22)

	!Build A12 and put it in Ann
    allocate( A12(ss,ss) )
    call build_A(A12,q,T,V55)
    Ann( 1:ss, (ss + 1):2*ss ) = -A12
    deallocate(A12)

	!Build A21 and put it in Ann
    allocate( A21(ss,ss) )
    call build_A(A21,q,T,V66)
    Ann( (ss + 1):2*ss, 1:ss ) = -A21
    deallocate(A21)

	!Build VV11 and put it in Vnn
	allocate( VV11(ss,n_theta) )
	call build_V(VV11,q,T,V12)
	Vnn(1:ss,1:n_theta) = VV11
	deallocate(VV11)

	!Build VV12 and put it in Vnn
	allocate( VV12(ss,n_theta) )
   	call build_V(VV12,q,T,V55)
   	Vnn( 1:ss,(n_theta + 1):2*n_theta ) = VV12
   	deallocate(VV12)

	!Build VV21 and put it in Vnn
	allocate( VV21(ss,n_theta) )
   	call build_V(VV21,q,T,V66)
   	Vnn( (ss + 1):2*ss,1:n_theta ) = VV21
   	deallocate(VV21)

	!Build VV22 and put it in Vnn
   	allocate( VV22(ss,n_theta) )
   	call build_V(VV22,q,T,V34)
   	Vnn( (ss + 1):2*ss,(n_theta + 1):2*n_theta ) = VV22
   	deallocate(VV22)

	!Call to lapack to solve Ax=b using LU-factorization. The solution is written into Vnn which is then written into B12...B66.
	call cgesv(n,nrhs,Ann,lda,ipiv,Vnn,ldb,info)
	B12 = Vnn(1:ss, 1:n_theta)
	B55 = Vnn(1:ss, (n_theta + 1 ):2*n_theta )
	B66 = Vnn( (ss + 1):2*ss,1:n_theta )
	B34 = Vnn( (ss + 1):2*ss, (n_theta + 1):2*n_theta )

	deallocate(Ann)
	deallocate(Vnn)

end subroutine

end module
