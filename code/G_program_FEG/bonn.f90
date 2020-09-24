module bonn
!Computes the Bonn potential from L. White and F. Sammarruca, "Solution of the Bethe-Goldstone Equation without Partial Wave Decomposition".
!
!References:
!   Ref.[1] = R. Machleidt et al., The Bonn meson-exchange model for the nucleon-nucleon interation.
!   Ref.[2] = R. Machleidt, Adv. Nucl. Phys. 19, 1989, "The Meson Theory of Nuclear Forces and Nuclear Structure".
!   Ref.[3] = R. Machleidt, "Computational Nuclear Physics 2: Nuclear Reactions", Editors: K. Langanke et al.
!	Ref.[4] = L. White and F. Sammarruca, "Solution of the Bethe-Goldstone Equation without Partial Wave Decomposition".
!   In the comments below alpha = pion,eta,sigma,delta,omega,rho.

use inputs
implicit none

    !Stores the Bonn potential parameters. The parameters came from Ref.[2].
 
   !___________Bonn B_________________

    !m = Free space nucleon mass in MeV
	real(8), parameter, public :: m = 938.926d0

    !Mass of pi,eta,sigma,delta,omega,rho mesons in MeV.
	real(8), parameter, private :: m_pi = 138.03d0
	real(8), parameter, private :: m_eta = 548.8d0
	real(8), parameter, private :: m_sigma = 550.0d0
	real(8), parameter, private :: m_delta = 983.0d0
	real(8), parameter, private :: m_omega = 782.6d0
	real(8), parameter, private :: m_rho = 769.0d0

    !Coupling constant of pi,eta,sigma,delta,omega,rho in MeV.
    !The numbers below (x) are:  x = ( (g_alpha)**2 )/(4*pi) and g_alpha appears in Ref.[1] p.75-76.
	real(8), parameter, private :: g_pi = 14.6d0
	real(8), parameter, private :: g_eta = 5.0d0
	real(8), parameter, private :: g_sigma = 8.0769d0
	real(8), parameter, private :: g_delta = 3.1155d0
	real(8), parameter, private :: g_omega = 20.0d0
	real(8), parameter, private :: g_rho = 0.95d0

    !Value for f_rho/g_rho and f_omega/g_omega.
	real(8), parameter, private :: coupling_ratio_rho = 6.1d0
	real(8), parameter, private :: coupling_ratio_omega = 0.0d0  

    !n_alpha appering in the form factor.
	real(8), parameter, private :: n_pi = 1.0d0
	real(8), parameter, private :: n_eta = 1.0d0
	real(8), parameter, private :: n_sigma = 1.0d0
	real(8), parameter, private :: n_delta = 1.0d0
	real(8), parameter, private :: n_omega = 1.0d0
	real(8), parameter, private :: n_rho = 1.0d0

    !lambda_alpha appering in the form factor in MeV.
	real(8), parameter, private :: lambda_pi = 1200.0d0
	real(8), parameter, private :: lambda_eta = 1500.0d0
	real(8), parameter, private :: lambda_sigma = 2000.0d0
	real(8), parameter, private :: lambda_delta = 1500.0d0
	real(8), parameter, private :: lambda_omega = 1500.0d0
	real(8), parameter, private :: lambda_rho = 1300.0d0

    !________________End of Bonn B_____________________


    !pi = Value for the constant pi
    real(8), parameter, public :: pi = 3.1415926535897932385d0

    !i = sqrt(-1)
	complex(8), parameter, public :: i = (0.0d0,1.0d0)


contains

subroutine helicity_elements(thetap,phip,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_innerprod,lambda_sigma_element)
!Computes <lambda_1' lambda_2' | lambda_1 lambda_2> and <lambda_1' lambda_2' | vec{sigma}^(1)*vec{sigma}^(2) |lambda_1 lambda_2>.

implicit none

    !(thetap,phip), (theta,phi) = Angular variables for vec{q},vec{q'}.
    !lambda1p,lambda2p,lambda1,lambda2 = Matrix elements i.e. < | >, < |sigma^(1)*sigma^(2)| >.
    !Note, values for lambda1p,lambda2p,lambda1,lambda2 are  +(-) = 0.5(-0.5) i.e. spin up and down.
	real(8), intent(in) :: thetap,phip,theta,phi
	real(8), intent(in) :: lambda1p,lambda2p,lambda1,lambda2

    !lambda_innerprod = Value for <lambda_1' lambda_2'|lambda_1 lambda_2>.
    !lambda_sigma_element = Value for <lambda_1' lambda_2'|vec{sigma}^(1)*vec{sigma}^(2)|lambda_1 lambda_2>.
	complex(8), intent(out) :: lambda_innerprod,lambda_sigma_element

    !lambda = Array contaning [lambda_1',lambda_2',lambda_1,lambda_2]
	real(8), dimension(4) :: lambda

    lambda = (/lambda1p,lambda2p,lambda1,lambda2/)

	if ( compare(lambda,(/0.5d0,0.5d0,0.5d0,0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,-0.5d0,-0.5d0,-0.5d0/)) == 1) then  
    !<++|++>

		lambda_innerprod = 0.5d0*( 1.0d0 + dcos(thetap)*dcos(theta) + dsin(thetap)*dsin(theta)*dcos(phip - phi) )
		lambda_sigma_element = lambda_innerprod - 2.0d0

	elseif ( compare(lambda,(/0.5d0,0.5d0,0.5d0,-0.5d0/)) == 1 .or. compare(lambda,(/0.5d0,0.5d0,-0.5d0,0.5d0/)) == 1&
    &.or. compare(lambda,(/-0.5d0,-0.5d0,0.5d0,-0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,-0.5d0,-0.5d0,0.5d0/)) == 1 ) then
    !<++|+->

		lambda_innerprod = 0.5d0*( dcos(thetap)*dsin(theta) - dsin(thetap)*dcos(theta)*dcos(phip - phi) ) -&
							&0.5d0*dsin(thetap)*i*dsin(phip - phi)
 
		lambda_sigma_element = lambda_innerprod

		if ( compare(lambda,(/0.5d0,0.5d0,-0.5d0,0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,-0.5d0,-0.5d0,0.5d0/)) == 1 ) then

			lambda_innerprod = dcmplx( -real(lambda_innerprod),aimag(lambda_innerprod) )
			lambda_sigma_element = dcmplx( -real(lambda_sigma_element),aimag(lambda_sigma_element) )

		endif

	elseif ( compare(lambda,(/0.5d0,0.5d0,-0.5d0,-0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,-0.5d0,0.5d0,0.5d0/)) == 1 ) then
    !<++|-->

		lambda_innerprod = 0.5d0*( -1.0d0 + dcos(thetap)*dcos(theta) + dsin(thetap)*dsin(theta)*dcos(phip - phi) )
		lambda_sigma_element = lambda_innerprod + 2.0d0

	elseif ( compare(lambda,(/0.5d0,-0.5d0,0.5d0,0.5d0/)) == 1 .or. compare(lambda,(/0.5d0,-0.5d0,-0.5d0,-0.5d0/)) == 1&
    &.or. compare(lambda,(/-0.5d0,0.5d0,-0.5d0,-0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,0.5d0,0.5d0,0.5d0/)) == 1 ) then
    !<+-|++>

		lambda_innerprod = 0.5d0*( dsin(thetap)*dcos(theta) - dsin(theta)*dcos(thetap)*dcos(phip - phi) ) -&
							&0.5d0*dsin(theta)*i*dsin(phip - phi)

		lambda_sigma_element = lambda_innerprod

		if ( compare(lambda,(/-0.5d0,0.5d0,-0.5d0,-0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,0.5d0,0.5d0,0.5d0/)) == 1 ) then

			lambda_innerprod = dcmplx( -real(lambda_innerprod),aimag(lambda_innerprod) )
			lambda_sigma_element = dcmplx( -real(lambda_sigma_element),aimag(lambda_sigma_element) )

		endif

	elseif ( compare(lambda,(/0.5d0,-0.5d0,0.5d0,-0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,0.5d0,-0.5d0,0.5d0/)) == 1 ) then
    !<+-|+->

		lambda_innerprod = 0.5d0*(  dsin(thetap)*dsin(theta) + ( 1.0d0 + dcos(thetap)*dcos(theta) )*dcos(phip - phi)  ) +&
							&0.5d0*i*( dcos(thetap) + dcos(theta) )*dsin(phip - phi)

		lambda_sigma_element = lambda_innerprod

		if ( compare(lambda,(/-0.5d0,0.5d0,-0.5d0,0.5d0/)) == 1 ) then

			lambda_innerprod = dcmplx( real(lambda_innerprod),-aimag(lambda_innerprod) )
			lambda_sigma_element = dcmplx( real(lambda_sigma_element),-aimag(lambda_sigma_element) )

		endif

	elseif ( compare(lambda,(/0.5d0,-0.5d0,-0.5d0,0.5d0/)) == 1 .or. compare(lambda,(/-0.5d0,0.5d0,0.5d0,-0.5d0/)) == 1 ) then
    !<+-|-+>

		lambda_innerprod = 0.5d0*(  -dsin(thetap)*dsin(theta) + ( 1.0d0 - dcos(thetap)*dcos(theta) )*dcos(phip - phi)  ) +&
							&0.5d0*i*( dcos(thetap) - dcos(theta) )*dsin(phip - phi)

		lambda_sigma_element = lambda_innerprod

		if ( compare(lambda,(/-0.5d0,0.5d0,0.5d0,-0.5d0/)) == 1 ) then

			lambda_innerprod = dcmplx( real(lambda_innerprod),-aimag(lambda_innerprod) )
			lambda_sigma_element = dcmplx( real(lambda_sigma_element),-aimag(lambda_sigma_element) )

		endif

	endif

end subroutine

integer function compare(A,B)
!Check to see if array A = array B.
!If array A = array B then it outputs 1. Otherwise it outputs 0.

implicit none

    !A,B = array to compare
	real(8), dimension(4), intent(in) :: A,B

    !c1,x = Counters
	integer :: c1,x

	x = 0
	do c1 = 1,4

		if ( A(c1) == B(c1) ) then

			x = x + 1

		endif

	enddo

	if ( x == 4 ) then

		compare = 1
		
	else

		compare = 0

	endif

end function

complex(8) function V_pv(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,&
					&lambda1,lambda2,lambda_alpha,n_alpha,g_alpha,m_alpha,m_1,m_2,m)
!V_pv, used for pi and eta mesons.
!Notes:
!   (f_{ps})**2/( (m_{ps}**2)(2pi)**3 ) --> x/( 8(m*pi)**2 ) Where x is the coupling constant listed as g_alpha. See page 347 Ref.[2].

implicit none
	
    !(qp,thetap,phip), (q,theta,phi) = vec{q},vec{q'}.
    !lambda1p,lambda2p,lambda1,lambda2 = Value for < | >, < |sigma^(1)*sigma^(2)| >.
    !Note, values for lambda1p,lambda2p,lambda1,lambda2 are 0.5d0 = + and -0.5d0 = - i.e spin up and down.
    !lambda_alpha,n_alpha,g_alpha,m_alpha,m_1,m_2,m = See parameters at begining of code for explanation.
	real(8), intent(in) :: qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_alpha,n_alpha,g_alpha,m_alpha,m_1,m_2,m

    !qp_q_squared = |vec{q'} - vec{q}|**2. It will be used in the form factor.
    !F = Form factor p.74 Ref.[1].
    !E1,E2,E1p,E2p = Relativistic energies (see Ref.[4]).
    !W1,W2,W1p,W2p = E1+M1,etc... (see Ref.[4])
    !lambda_innerprod = <lambda_1' lambda_2' | lambda_1 lambda_2>.
    !lambda_sigma_element = <lambda_1' lambda_2' | vec{sigma}^(1)*vec{sigma}^(2) |lambda_1 lambda_2>.
    !temp1,temp2,..., and beta,L1p,L1,L2p,L2,ans = Temporary variables
	real(8) :: qp_q_squared,F,E1,E2,E1p,E2p,W1,W2,W1p,W2p
	complex(8) :: lambda_innerprod,lambda_sigma_element
	real(8) :: temp1,temp2,temp3,temp4,temp5,beta,L1p,L1,L2p,L2

	qp_q_squared = qp**2.0d0 + q**2.0d0 - 2.0d0*qp*q*( dsin(theta)*dsin(thetap)*dcos(phip-phi) + dcos(thetap)*dcos(theta) )

	F = (   ( (lambda_alpha)**2.0d0 - (m_alpha)**2.0d0 )/(  (lambda_alpha)**2.0d0 + qp_q_squared  )   )**n_alpha
	E1 = dsqrt( m_1**2.0d0 + q**2.0d0 )
	E2 = dsqrt( m_2**2.0d0 + q**2.0d0 )
	E1p = dsqrt( m_1**2.0d0 + qp**2.0d0 )
	E2p = dsqrt( m_2**2.0d0 + qp**2.0d0 )	
	W1 = E1 + m_1
	W2 = E2 + m_2
	W1p = E1p + m_1
	W2p = E2p + m_2
	
	call helicity_elements(thetap,phip,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_innerprod,lambda_sigma_element)

    !Temporary variables
	temp1 = g_alpha/( 8.0d0*((m*pi)**2.0d0) )
	temp2 = dsqrt( (W1p*W2p*W1*W2)/(E1*E1p*E2*E2p) )
	temp3 = (F**2.0d0)/( qp_q_squared + m_alpha**2.0d0 )
	beta = temp1*temp2*temp3
	L1p = (lambda1p*qp)/W1p
	L1 = (lambda1*q)/W1
	L2p = (lambda2p*qp)/W2p
	L2 = (lambda2*q)/W2
	temp4 = E1p - E1
	temp5 = E2p - E2

	V_pv = beta*m_1*m_2&
	&*(  4.0d0*( L1p - L1 )*( L2p - L2 ) + ( (temp4*temp5)/(m_1*m_2) )*( L1p + L1 )*( L2p + L2 )&
	&+ ( (2.0d0*temp4)/m_1 )*( L1p + L1 )*( L2p - L2 ) + ( (2.0d0*temp5)/m_2 )*( L1p - L1 )*( L2p + L2 )  )*lambda_innerprod

end function

complex(8) function V_s(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,&
					&lambda1,lambda2,lambda_alpha,n_alpha,g_alpha,m_alpha,m_1,m_2)
!V_s, used for sigma and delta mesons.
!Notes:
!   (g_{ps})**2/( (2pi)**3 ) --> x/( 2*pi**2 ) Where x is the coupling constant listed as g_alpha.
!
implicit none
	
    !(qp,thetap,phip),(q,theta,phi) = vec{q},vec{q'}.
    !lambda1p,lambda2p,lambda1,lambda2 = < | >, < |sigma^(1)*sigma^(2)| >.
    !Note, values for lambda1p,lambda2p,lambda1,lambda2 are 0.5d0 = + and -0.5d0 = - i.e spin up and down.
    !lambda_alpha,n_alpha,g_alpha,m_alpha,m_1,m_2 = See parameters at begining of code for explanation.
	real(8), intent(in) :: qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_alpha,n_alpha,g_alpha,m_alpha,m_1,m_2

    !qp_q_squared = |vec{q'} - vec{q}|**2. It will be used in the form factor.
    !F = Form factor p.74 Ref.[1].
    !E1,E2,E1p,E2p = Relativistic energies (see Ref.[4]).
    !W1,W2,W1p,W2p = E1+M1,etc... (see Ref.[4])
    !lambda_innerprod = <lambda_1' lambda_2' | lambda_1 lambda_2>
    !lambda_sigma_element = <lambda_1' lambda_2' | vec{sigma}^(1)*vec{sigma}^(2) |lambda_1 lambda_2>
    !temp1,temp2,..., and beta,L1p,L1,L2p,L2,ans = Temporary variables.
	real(8) :: qp_q_squared,F,E1,E2,E1p,E2p,W1,W2,W1p,W2p
	complex(8) :: lambda_innerprod,lambda_sigma_element
	real(8) :: temp1,temp2,temp3,beta,L1p,L1,L2p,L2

	qp_q_squared = qp**2.0d0 + q**2.0d0 - 2.0d0*qp*q*( dsin(theta)*dsin(thetap)*dcos(phip-phi) + dcos(thetap)*dcos(theta) )

	F = (   ( (lambda_alpha)**2.0d0 - (m_alpha)**2.0d0 )/(  (lambda_alpha)**2.0d0 + qp_q_squared  )   )**n_alpha
	E1 = dsqrt( m_1**2.0d0 + q**2.0d0 )
	E2 = dsqrt( m_2**2.0d0 + q**2.0d0 )
	E1p = dsqrt( m_1**2.0d0 + qp**2.0d0 )
	E2p = dsqrt( m_2**2.0d0 + qp**2.0d0 )	
	W1 = E1 + m_1
	W2 = E2 + m_2
	W1p = E1p + m_1
	W2p = E2p + m_2
	
	call helicity_elements(thetap,phip,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_innerprod,lambda_sigma_element)

    !Temporary variables
	temp1 = -g_alpha/( 2.0d0*(pi**2.0d0) )
	temp2 = dsqrt( (W1p*W2p*W1*W2)/(E1*E1p*E2*E2p) )
	temp3 = (F**2.0d0)/( qp_q_squared + m_alpha**2.0d0 )
	beta = temp1*temp2*temp3
	L1p = (lambda1p*qp)/W1p
	L1 = (lambda1*q)/W1
	L2p = (lambda2p*qp)/W2p
	L2 = (lambda2*q)/W2

	V_s = 0.25d0*beta*( 1.0d0 - 4.0d0*L1p*L1 )*( 1.0d0 - 4.0d0*L2p*L2 )*lambda_innerprod

end function

complex(8) function V_v(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,&
					&lambda1,lambda2,lambda_alpha,n_alpha,g_alpha&
	                &,coupling_ratio_alpha,m_alpha,m_1,m_2,m)
!V_v, used for omega and rho mesons.
!Notes:
!   (g_{ps})**2/( (2pi)**3 ) --> x/( 2*pi**2 ) Where x is the coupling constant listed as g_alpha.
!	For vector boson exchange, the potential is the sum of three terms (V_v = V_vv + V_tt + V_vt).

implicit none
	
    !(qp,thetap,phip), (q,theta,phi) = vec{q},vec{q'}.
    !lambda1p,lambda2p,lambda1,lambda2 = < | >,< |sigma^(1)*sigma^(2)| >.
    !Note, values for lambda1p,lambda2p,lambda1,lambda2 are 0.5d0 = + and -0.5d0 = - i.e spin up and down.
    !lambda_alpha,n_alpha,g_alpha,coupling_ratio_alpha,m_alpha,m_1,m_2,m = See parameters at the begining of code for explanation.
	real(8), intent(in) :: qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_alpha,n_alpha,g_alpha
	real(8), intent(in) :: coupling_ratio_alpha,m_alpha,m_1,m_2,m

    !qp_dot_q = vec{q'} dot vec{q}.
    !qp_q_squared = |vec{q'} - vec{q}|**2. It will be used in the form factor.
    !F = Form factor p.74 Ref.[1].
    !E1,E2,E1p,E2p = Relativistic energies (see Ref.[4]).
    !W1,W2,W1p,W2p = E1+M1,etc... (see Ref.[4])
    !lambda_innerprod = <lambda_1' lambda_2' | lambda_1 lambda_2>.
    !lambda_sigma_element = <lambda_1' lambda_2' | vec{sigma}^(1)*vec{sigma}^(2) |lambda_1 lambda_2>.
    !temp1,temp2,L1p,L1,L2p,L1 and temp1vv(vt,tt),temp2vv(vt,tt),..., and betavv(vt,tt),V_vv(vt,tt) = Temporary variables.
	real(8) :: qp_dot_q,qp_q_squared,F,E1,E2,E1p,E2p,W1,W2,W1p,W2p,temp1,temp2,L1p,L1,L2p,L2
	complex(8) :: lambda_innerprod,lambda_sigma_element
	real(8) :: temp1vv,betavv
	complex(8):: V_vv
	real(8) :: temp1vt,temp2vt,temp3vt,temp4vt,temp5vt,temp6vt,betavt
	complex(8) :: V_vt
	real(8) :: temp1tt,temp2tt,temp3tt,temp4tt,temp5tt,temp6tt,temp7tt,betatt
	complex(8) :: V_tt

	qp_dot_q = qp*q*( dsin(thetap)*dsin(theta)*dcos(phip-phi) + dcos(thetap)*dcos(theta) )
	qp_q_squared = qp**2.0d0 + q**2.0d0 - 2.0d0*qp_dot_q

	F = (   ( (lambda_alpha)**2.0d0 - (m_alpha)**2.0d0 )/(  (lambda_alpha)**2.0d0 + qp_q_squared  )   )**n_alpha
	E1 = dsqrt( m_1**2.0d0 + q**2.0d0 )
	E2 = dsqrt( m_2**2.0d0 + q**2.0d0 )
	E1p = dsqrt( m_1**2.0d0 + qp**2.0d0 )
	E2p = dsqrt( m_2**2.0d0 + qp**2.0d0 )	
	W1 = E1 + m_1
	W2 = E2 + m_2
	W1p = E1p + m_1
	W2p = E2p + m_2
	
	call helicity_elements(thetap,phip,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_innerprod,lambda_sigma_element)

    !Temporary variables
	temp1 = dsqrt( (W1p*W2p*W1*W2)/(E1*E1p*E2*E2p) )
	temp2 = (F**2.0d0)/( qp_q_squared + m_alpha**2.0d0 )
	L1p = (lambda1p*qp)/W1p
	L1 = (lambda1*q)/W1
	L2p = (lambda2p*qp)/W2p
	L2 = (lambda2*q)/W2

!V_vv

    !Temporary variables
	temp1vv = g_alpha/( 2.0d0*(pi**2.0d0) )
	betavv = 0.25d0*temp1*temp2*temp1vv 		

	V_vv = betavv*(  ( 1.0d0 + 4.0d0*L1p*L1 )*( 1.0d0 + 4.0d0*L2p*L2 )*lambda_innerprod&
	&- 4.0d0*( L1p + L1 )*( L2p + L2 )*lambda_sigma_element  )

!V_vt

    !Temporary variables
	temp1vt = g_alpha/( 2.0d0*(pi**2.0d0) )
	betavt = temp1*temp2*temp1vt
	temp2vt = ( W1p + W2p + W1 + W2 )/m
	temp3vt = ( E1p + E2p + E1 + E2 - 2.0d0*(m_1 + m_2) )/(2.0d0*m)
	temp4vt = (m_1 + m_2)/m
	temp5vt = ( E1p - E1 )/m
	temp6vt = (E2p - E2)/m

	V_vt = (coupling_ratio_alpha/2.0d0)*betavt*(   (  temp2vt*( 8.0d0*L1p*L1*L2p*L2 ) - temp3vt  )*lambda_innerprod&
	&- (  (2.0d0*temp4vt)*(L1p + L1)*(L2p + L2) + temp5vt*(L1p - L1)*(L2p + L2)&
	&+ temp6vt*(L1p + L1)*(L2p - L2)  )*lambda_sigma_element   )

!V_tt

    !Temporary variables
	temp1tt = g_alpha/( 2.0d0*(pi**2.0d0) )
	betatt = temp1*temp2*temp1tt
	temp2tt = (m_1*m_2)/(m**2.0d0)
    temp3tt = ( E1p + E2p + E1 + E2 )/( 2.0d0*(m**2.0d0) )
	temp4tt = (E1p + E1)*(E2p + E2) - (E1p - E1)**2.0d0 - (E2p - E2)**2.0d0 + qp**2.0d0 + q**2.0d0 + 2.0d0*qp_dot_q	
	temp5tt = ( m_1*(E2p - E2) )/(m**2.0d0)
	temp6tt = ( m_2*(E1p - E1) )/(m**2.0d0)
	temp7tt = ( (E1p - E1)*(E2p - E2) )/(m**2)

	V_tt = ( (coupling_ratio_alpha/2.0d0)**2.0d0 )*betatt*(    (   temp2tt*(1.0d0 + 4.0d0*L1p*L1)*(1.0d0 + 4.0d0*L2p*L2)&
	&- temp3tt*( m_1*(1.0d0 + 4.0d0*L1p*L1)*(1.0d0 - 4.0d0*L2p*L2) + m_2*(1.0d0 - 4.0d0*L1p*L1)*(1.0d0 + 4.0d0*L2p*L2) )&
	&+ (1.0d0/(m**2.0d0))*(1.0d0 - 4.0d0*L1p*L1)*(1.0d0 - 4.0d0*L2p*L2)*( 2.0d0*m_1*m_2 + 0.25d0*temp4tt )   )*lambda_innerprod&
	&- (   4.0d0*temp2tt*(L1p + L1)*(L2p + L2) + 2.0d0*temp5tt*(L1p + L1)*(L2p - L2)&
	&+ 2.0d0*temp6tt*(L1p - L1)*(L2p + L2) + temp7tt*(L1p - L1)*(L2p - L2)   )*lambda_sigma_element    )

!V_v

	V_v = V_vv + V_vt + V_tt

end function

complex(8) function V(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,T)
!V_Bonn

implicit none

    !T = Total isospin. T can be 0 or 1.
    !(qp,thetap,phip), (q,theta,phi) = vec{q},vec{q'}.
    !lambda1p,lambda2p,lambda1,lambda2 = < | >,< |sigma^(1)*sigma^(2)| >.
    !Note, values for lambda1p,lambda2p,lambda1,lambda2 are 0.5d0 = + and -0.5d0 = - i.e spin up and down.
	integer, intent(in) :: T
	real(8), intent(in) :: qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2

    !V_pi,V_eta,V_sigma,V_delta,V_omega,V_rho = Different meson contributions to the potential.
    !x = Isospin coefficient. i.e. x = -3 for T = 0 and x = 1 for T = 1.
	complex(8) :: V_pi,V_eta,V_sigma,V_delta,V_omega,V_rho
	real(8) :: x

    !Pi (isovector)
	V_pi = V_pv(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_pi,n_pi,g_pi,m_pi,m_1,m_2,m)

    !Eta (isoscalar)
	V_eta = V_pv(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_eta,n_eta,g_eta,m_eta,m_1,m_2,m)

    !Sigma (isoscalar)
	V_sigma = V_s(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_sigma,n_sigma,g_sigma,m_sigma,m_1,m_2)

    !Delta (isovector)
	V_delta = V_s(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_delta,n_delta,g_delta,m_delta,m_1,m_2)

    !Omega (isoscalar)
	V_omega = V_v(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_omega,n_omega,g_omega&
	&,coupling_ratio_omega,m_omega,m_1,m_2,m)

    !Rho (isovector)
	V_rho = V_v(qp,thetap,phip,q,theta,phi,lambda1p,lambda2p,lambda1,lambda2,lambda_rho,n_rho,g_rho&
	&,coupling_ratio_rho,m_rho,m_1,m_2,m)

    !Isospin coefficient
	if( T == 1 ) then

		x = 1.0d0

	else

		x = -3.0d0

	endif

    !V_Bonn
	V = V_eta + V_sigma + V_omega + x*( V_pi + V_delta + V_rho )

end function

end module
