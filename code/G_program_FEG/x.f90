module x
!Solve for x i.e. we are solving x = ( A^(-1) )*b

use inputs
use ab
implicit none

contains

subroutine T0(B0,q,T,cutoff_inf)
!^{0}T

implicit none

    !Lapack parameters.
    integer, parameter :: n = ss
    integer, parameter :: nrhs = 1
    integer, parameter :: lda = n
    integer, parameter :: ldb = n
    integer, parameter :: lwork = n
    integer :: info
    real, dimension(lwork) :: work
    integer, dimension(n) :: ipiv

	!q = |q|.
	!cutoff_inf = Cutoff for [2q,inf] integral.
	!T = Isospin
	!B0 = Right hand side of Ax=b. Then it gets replaced by the solution vector ^{0}T.
	real(8), intent(in) :: q,cutoff_inf
	integer, intent(in) :: T
	complex, dimension(ldb,nrhs), intent(out) :: B0

	!c1 = Counter
	!A0 = The matrix in the equation Ax=b.
	integer :: c1
	complex, dimension(:,:), allocatable :: A0

	allocate( A0(lda,n) )

	!Build the 0^A matrix
	call build_A(A0,q,0.0d0,T,VLAMBDA0,cutoff_inf)

	!Compute 1 - ^{0}A i.e. "^{0}B" in T matrix report.
	A0 = -A0
	do c1 = 1,n

		A0(c1,c1) = 1.0 + A0(c1,c1)

	enddo

	!Build right hand side vector in the equation Ax=b.
	call build_V(B0,q,T,V0,cutoff_inf)

	!Call to lapack to solve Ax=b using LU-factorization. The solution is written into B0.
	call cgesv(n,nrhs,A0,lda,ipiv,B0,ldb,info)

	deallocate(A0)

end subroutine

subroutine T_12_6(B12,B6,q,T)
!^{12}T,T_6

implicit none

  	!Lapack parameters. 
    integer, parameter :: n = 2*ss
    integer, parameter :: nrhs = 1
    integer, parameter :: lda = n
    integer, parameter :: ldb = n
    integer, parameter :: lwork = n
    integer :: info
    real, dimension(lwork) :: work
    integer, dimension(n) :: ipiv

	!q = |q|
	!T = Isospin
	!B12,B6 = Solution vector ^{12}T,T_6.
	real(8), intent(in) :: q
	integer, intent(in) :: T
	complex, dimension(ss), intent(out) :: B12,B6

	!c1,c2 = Counters
	!Ann = The A matrix in the equation Ax=b.
	!Vnn = The right hand side (b matrix) in the equation Ax=b.
	!A11...A22 = Matrix entries in Ann.
	!VV1...VV2 = Vector entries in Vnn.
	integer :: c1,c2
	complex, dimension(:,:), allocatable :: Ann
	complex, dimension(:), allocatable :: Vnn
	complex, dimension(:,:), allocatable :: A11
	complex, dimension(:,:), allocatable :: A12
	complex, dimension(:,:), allocatable :: A21
	complex, dimension(:,:), allocatable :: A22
	complex, dimension(:), allocatable :: VV1
    complex, dimension(:), allocatable :: VV2

	allocate( Ann(lda,n) )
	allocate( Vnn(ldb) )

	!Build A11 and put it in Ann
	allocate( A11(ss,ss) )
	call build_A(A11,q,0.0d0,T,VLAMBDA12,cutoff_inf)
	A11 = -A11
	do c1 = 1,ss

		A11(c1,c1) = 1.0 + A11(c1,c1)

	enddo
	Ann(1:ss,1:ss) = A11
	deallocate(A11)

	!Build A12 and put it in Ann
    allocate( A12(ss,ss) )
	call build_A(A12,q,0.0d0,T,VLAMBDA_5,cutoff_inf)
    Ann( 1:ss, (ss + 1):2*ss ) = -4.0d0*A12
    deallocate(A12)

	!Build A21 and put it in Ann
    allocate( A21(ss,ss) )
	call build_A(A21,q,0.0d0,T,VLAMBDA_6,cutoff_inf)
    Ann( (ss + 1):2*ss, 1:ss ) = -A21
    deallocate(A21)

	!Build A22 and put it in Ann
    allocate( A22(ss,ss) )
	call build_A(A22,q,0.0d0,T,VLAMBDA1,cutoff_inf) 
    A22 = -A22
    do c2 = 1,ss

        A22(c2,c2) = 1.0 + A22(c2,c2)

    enddo
	Ann( (ss + 1):2*ss, (ss + 1):2*ss ) = A22
    deallocate(A22)

	!Build VV1 and put it in Vnn
	allocate( VV1(ss) )
	call build_V(VV1,q,T,V12,cutoff_inf)
	Vnn(1:ss) = VV1
	deallocate(VV1)

	!Build VV2 and put it in Vnn
	allocate( VV2(ss) )
	call build_V(VV2,q,T,V_6,cutoff_inf)
   	Vnn( (ss + 1):2*ss) = VV2
   	deallocate(VV2)

	!Call to lapack to solve Ax=b using LU-factorization. The solution is written into Vnn which is then written into B12,B6.
	call cgesv(n,nrhs,Ann,lda,ipiv,Vnn,ldb,info)
	B12 = Vnn(1:ss)
	B6 = Vnn( (ss + 1):2*ss )

	deallocate(Ann)
	deallocate(Vnn)

end subroutine

subroutine T_3_4_5(B3,B4,B5,q,T)
!T_3,T_4,T_5

implicit none

  	!Lapack parameters. 
    integer, parameter :: n = 3*ss
    integer, parameter :: nrhs = 1
    integer, parameter :: lda = n
    integer, parameter :: ldb = n
    integer, parameter :: lwork = n
    integer :: info
    real, dimension(lwork) :: work
    integer, dimension(n) :: ipiv

	!q = |q|
	!T = Isospin
	!B3,B4,B5 = Solution vector T_3,T_4,T_5.
	real(8), intent(in) :: q
	integer, intent(in) :: T
	complex, dimension(ss), intent(out) :: B3,B4,B5

	!c1,c2 = Counters
	!Ann = The A matrix in the equation Ax=b.
	!Vnn = The right hand side (b matrix) in the equation Ax=b.
	!A11...A33 = Matrix entries in Ann.
	!VV1,VV2,VV3 = Vector entries in Vnn.
	integer :: c1,c2
	complex, dimension(:,:), allocatable :: Ann
	complex, dimension(:), allocatable :: Vnn
	complex, dimension(:,:), allocatable :: A11
	complex, dimension(:,:), allocatable :: A12
	complex, dimension(:,:), allocatable :: A13
	complex, dimension(:,:), allocatable :: A21
	complex, dimension(:,:), allocatable :: A22
	complex, dimension(:,:), allocatable :: A23
	complex, dimension(:,:), allocatable :: A31
	complex, dimension(:,:), allocatable :: A32
	complex, dimension(:,:), allocatable :: A33
	complex, dimension(:), allocatable :: VV1
    complex, dimension(:), allocatable :: VV2
	complex, dimension(:), allocatable :: VV3

	allocate( Ann(lda,n) )
	allocate( Vnn(ldb) )

	!Build A11 and put it in Ann
	allocate( A11(ss,ss) )
	call build_A(A11,q,1.0d0,T,VLAMBDA_3,cutoff_inf)
	A11 = -A11
	do c1 = 1,ss

		A11(c1,c1) = 1.0 + A11(c1,c1)

	enddo
	Ann(1:ss,1:ss) = A11
	deallocate(A11)

	!Build A12 and put it in Ann
	allocate( A12(ss,ss) )
	call build_A(A12,q,1.0d0,T,VLAMBDA_4,cutoff_inf)
	A12 = -A12
	Ann(1:ss,(ss+1):2*ss) = A12
	deallocate(A12)

	!Build A13 and put it in Ann
	allocate( A13(ss,ss) )
	call build_A(A13,q,1.0d0,T,VLAMBDA_6,cutoff_inf)
	A13 = -2.0d0*A13
	Ann(1:ss,(2*ss+1):3*ss) = A13
	deallocate(A13)

	!Build A21 and put it in Ann
	allocate( A21(ss,ss) )
	call build_A(A21,q,-1.0d0,T,VLAMBDA_4,cutoff_inf)
	A21 = -A21
	Ann((ss+1):2*ss,1:ss) = A21
	deallocate(A21)

	!Build A22 and put it in Ann
	allocate( A22(ss,ss) )
	call build_A(A22,q,-1.0d0,T,VLAMBDA_3,cutoff_inf)
	A22 = -A22
	do c1 = 1,ss

		A22(c1,c1) = 1.0 + A22(c1,c1)

	enddo
	Ann((ss+1):2*ss,(ss+1):2*ss) = A22
	deallocate(A22)

	!Build A23 and put it in Ann
	allocate( A23(ss,ss) )
	call build_A(A23,q,-1.0d0,T,VLAMBDA_6,cutoff_inf)
	A23 = 2.0d0*A23
	Ann((ss+1):2*ss,(2*ss+1):3*ss) = A23
	deallocate(A23)

	!Build A31 and put it in Ann
	allocate( A31(ss,ss) )
	call build_A(A31,q,1.0d0,T,VLAMBDA_5,cutoff_inf)
	A31 = -A31
	Ann((2*ss+1):3*ss,1:ss) = A31
	deallocate(A31)

	!Build A32 and put it in Ann
	allocate( A32(ss,ss) )
	call build_A(A32,q,-1.0d0,T,VLAMBDA_5,cutoff_inf)
	Ann((2*ss+1):3*ss,(ss+1):2*ss) = A32
	deallocate(A32)

	!Build A33 and put it in Ann
	allocate( A33(ss,ss) )
	call build_A(A33,q,1.0d0,T,VLAMBDA12,cutoff_inf)
	A33 = -A33
	do c1 = 1,ss

		A33(c1,c1) = 1.0 + A33(c1,c1)

	enddo
	Ann((2*ss+1):3*ss,(2*ss+1):3*ss) = A33
	deallocate(A33)

	!Build VV1 and put it in Vnn
	allocate( VV1(ss) )
	call build_V(VV1,q,T,V_3,cutoff_inf)
	Vnn(1:ss) = VV1
	deallocate(VV1)

	!Build VV2 and put it in Vnn
	allocate( VV2(ss) )
	call build_V(VV2,q,T,V_4,cutoff_inf)
	Vnn((ss+1):2*ss) = VV2
	deallocate(VV2)

	!Build VV3 and put it in Vnn
	allocate( VV3(ss) )
	call build_V(VV3,q,T,V_5,cutoff_inf)
	Vnn((2*ss+1):3*ss) = VV3
	deallocate(VV3)

	!Call to lapack to solve Ax=b using LU-factorization. The solution is written into Vnn which is then written into B3,B4,B5.
	call cgesv(n,nrhs,Ann,lda,ipiv,Vnn,ldb,info)
	B3 = Vnn(1:ss)
	B4 = Vnn((ss+1):2*ss)
	B5 = Vnn((2*ss+1):3*ss)

	deallocate(Ann)
	deallocate(Vnn)

end subroutine

end module