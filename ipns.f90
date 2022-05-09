!============================================================================!
        program ipns 
!  
!  Purpose:   Incompressible parabolized Navier-Stokes Solver
!  
!  Notes:     This is an incompressible parabolized solver for the
!             parabolic cylinder problem.  The formulation follows
!             that of R.T. Davis in JFM.
!
!             For some reason, as you march downstream, there is
!             a tendency for the Newton iteration to fail.  I 
!             have implemented quite a bit of logic to delay
!             failure, but it never can quite make it to 
!             infinity.  Be careful about using too many points
!             in the streamwise direction.  I have pretty good
!             luck with 32 points in \xi, but using more points
!             causes the Newton iteration to quite earlier.
!
!             Note that fourth-order differencing is used in the
!             wall normal direction and second-order midpnt 
!             integration is used in the steamwise direction.
!
!   Author:   Scott Collis
!
!   Revised:  3-18-96
!   
!============================================================================!
        use constants
        use diff
        implicit none

!.... flow data

        real, allocatable :: q(:,:,:)
        real, allocatable :: g(:), h(:), rhs(:,:), lhs(:,:,:,:)
        real, allocatable :: g1(:), g2(:), h1(:), h2(:)
        real, allocatable :: dg(:), dh(:)
        real, allocatable :: gbest(:), hbest(:)
        real :: norm, normold

!.... flow parameters

        real :: R               ! Reynolds number
        real :: Pr=1.0          ! Prandtl number

!.... mesh

        real, allocatable :: N(:), eta(:), N1(:), N2(:)
        real, allocatable :: s(:), xi(:), jac(:)
        real, external    :: s1, s2
        
        real    :: dN, ds, A, xi_max
        integer :: neta, nxi
        
        real, external :: calc_xi
        
        real :: time = 0.0, delt = 5.0, alpha = 1.0, scale = 1.0

!.... local vars

        integer :: i, j, k, iter, niter
        real    :: b1, b2, b3, b4, b5
        
        real    :: gold1=0.0, gold2=0.0, h1old1=0.0, h1old2=0.0
        integer :: count=0
        real    :: sm, xim, s1m, gtmp=0.0, isign=-1.0
        
!.... control flags 

        logical :: midpnt=.true.
        logical :: stretch=.true.

!.... Lapack LU solve

        integer, allocatable :: ipiv(:)
        integer :: info
        
!.... physical mesh and B-spline

        real, allocatable :: xy(:,:,:), xknot(:), yknot(:), bs(:,:,:), v(:,:,:)
        integer :: nx, ny, nz, korder, nxi2
        real :: tmp, xl, yl, etal, xil, sl, Nl
        real :: m1l, m2l, n1l, n2l, jl, fl, f1l, f2l, gl, g1l, g2l
        real :: hl, h1l, h2l
        real, external :: BS2VL, BS2DR
        
!       real, parameter :: eps1 = 1.0e-10       ! for Newton
!       real, parameter :: eps2 = 1.0e-07       ! for wall vorticity

        real, parameter :: eps1 = 1.0e-8        ! for Newton
        real, parameter :: eps2 = 1.0e-3        ! for wall vorticity
        
        real :: gamma = 1.4, gamma1 = 0.4, cv = 716.5
        real :: Ma, psil, vorl, ul, vl, tl, rhol, pl
!=============================================================================!
!       I n p u t  a n d   A l l o c a t e
!=============================================================================!
        write(*,"('Enter R ==> ',$)")
        read(*,*) R

        write(*,"('Enter Nxi, Neta ==> ',$)")
        read(*,*) nxi, neta

        write(*,"('Enter Niter ==> ',$)")
        read(*,*) niter

        delt = zero
!       write(*,"('Enter delt (0 for pure Newton) ==> ',$)")
!       read(*,*) delt
                
        allocate( s(nxi), xi(nxi), jac(neta) )
        allocate( q(neta,nxi,2) )
        allocate( N(neta), eta(neta), g(neta), h(neta), g1(neta), g2(neta), &
                  h1(neta), h2(neta), rhs(neta,2), lhs(neta,2,neta,2), &
                  N1(neta), N2(neta), ipiv(2*neta), dg(neta), dh(neta), &
                  gbest(neta), hbest(neta) )

!.... make the grid and metrics

        if (stretch) then
          A = four + 0.4 * sqrt(R)
          ds = one / real(nxi-1)
          s(1)  = zero
          xi(1) = zero
          write(10,10) s(1), xi(1)
          do i = 2, nxi-1
            s(i) = real(i-1) * ds
            xi(i) = calc_xi( s(i), A, xi(i-1)+0.0001, 1.0e5)
            write(10,10) s(i), xi(i), s1(xi(i),A), s2(xi(i),A)
          end do
          s(nxi) = one
          xi(nxi) = infty
          write(10,10) s(nxi), xi(nxi)
        else
          xi_max = 1000.0
          ds = xi_max / real(nxi-1)
          do i = 1, nxi
            s(i) = real(i-1) * ds
            xi(i) = s(i)
            write(10,10) s(i), xi(i), one, zero
          end do
        end if
        
        dn = one / real(neta-1)
        do j = 1, neta
          N(j) = real(j-1) * dn
          if (j .ne. neta ) then
            eta(j) = ( five * N(j) / ( one - N(j) ) ) + Sqrt(R)
          else
            eta(j) = infty
          end if
          N1(j) = (one - N(j))**2 / 5.0
          N2(j) = two * (N(j) - one)**3 / 25.0
        end do

!.... compute the potential flow solution as the initial condition at the LE

        g = zero
        h = -Sqrt(R)

        g(1) = 0.860 ! good guess for R=10
        write(*,"('Enter guess for g(1) (For Re=10 g(1)=0.86 ==> ',$)")
        read(*,*) g(1)
        
        i = 1
        
100     continue
!=============================================================================!
!       M a i n   i t e r a t i o n
!=============================================================================!

        do iter = 1, niter

!.... satisfy the boundary conditions

        h(1) = -Sqrt(R)
        
        g(neta) = zero
        h(neta) = -( gc5 * h(neta-4)  + &
                     gc4 * h(neta-3)  + &
                     gc3 * h(neta-2)  + &
                     gc2 * h(neta-1)  ) / gc1

!.... compute first derivatives in computational space

        call grad1( 1, neta, g, g1, dN, -1, .false., .false. ) 
        call grad1( 1, neta, h, h1, dN, -1, .false., .false. ) 

!.... compute second derivatives in computational space

        call grad2( 1, neta, g, g2, dN, -1, .false., .false. ) 
        call grad2( 1, neta, h, h2, dN, -1, .false., .false. ) 

!.... transfer to physical coordinates

        g2(:) = g2 * N1**2 + g1 * N2
        g1(:) = g1 * N1
        
        h2(:) = h2 * N1**2 + h1 * N2
        h1(:) = h1 * N1

!=============================================================================!
!       F o r m   t h e   R H S  
!=============================================================================!

        rhs(:,1) = g2/eta**2 + (h/eta**2 + one/eta - four/eta**3) * g1 + &
                   ( -(h1+one)/eta**2 - two*(h/eta**3+one/eta**2) ) * g
        rhs(:,2) = h2 - g

!..... boundary conditions on the rhs

        rhs(1,1) = zero
        rhs(1,2) = zero

        rhs(neta,1) = zero
        rhs(neta,2) = -gc1 * h(neta  ) - gc2 * h(neta-1) &
                      -gc3 * h(neta-2) - gc4 * h(neta-3) &
                      -gc5 * h(neta-4)

!=============================================================================!
!       F o r m   t h e   L H S 
!=============================================================================!

        lhs = zero

!.... second derivatives
        
        b1 = dd1 / dN**2
        b2 = dd2 / dN**2
        b3 = dd3 / dN**2
        b4 = dd4 / dN**2
        b5 = dd5 / dN**2

        j = 1
        k = 1
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b1
        lhs(j,2,k,2) = N1(j)**2 * b1
        k = 2
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b2
        lhs(j,2,k,2) = N1(j)**2 * b2
        k = 3
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b3
        lhs(j,2,k,2) = N1(j)**2 * b3
        k = 4
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b4
        lhs(j,2,k,2) = N1(j)**2 * b4
        k = 5
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b5
        lhs(j,2,k,2) = N1(j)**2 * b5

        b1 = db1 / dN**2
        b2 = db2 / dN**2
        b3 = db3 / dN**2
        b4 = db4 / dN**2
        b5 = db5 / dN**2

        j = 2
        k = 1
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b1
        lhs(j,2,k,2) = N1(j)**2 * b1
        k = 2
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b2
        lhs(j,2,k,2) = N1(j)**2 * b2
        k = 3
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b3
        lhs(j,2,k,2) = N1(j)**2 * b3
        k = 4
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b4
        lhs(j,2,k,2) = N1(j)**2 * b4
        k = 5
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b5
        lhs(j,2,k,2) = N1(j)**2 * b5

        b1 =  da1 / dN**2
        b2 =  da2 / dN**2
        b3 =  da3 / dN**2
        b4 =  da4 / dN**2
        b5 =  da5 / dN**2

        do j = 3, neta-2
          k = j-2
          lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b1
          lhs(j,2,k,2) = N1(j)**2 * b1
          k = j-1
          lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b2
          lhs(j,2,k,2) = N1(j)**2 * b2
          k = j
          lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b3
          lhs(j,2,k,2) = N1(j)**2 * b3
          k = j+1
          lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b4
          lhs(j,2,k,2) = N1(j)**2 * b4
          k = j+2
          lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b5
          lhs(j,2,k,2) = N1(j)**2 * b5
        end do

        b1 = db5 / dN**2
        b2 = db4 / dN**2
        b3 = db3 / dN**2
        b4 = db2 / dN**2
        b5 = db1 / dN**2

        j = neta-1
        k = neta-4
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b1
        lhs(j,2,k,2) = N1(j)**2 * b1
        k = neta-3
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b2
        lhs(j,2,k,2) = N1(j)**2 * b2
        k = neta-2
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b3
        lhs(j,2,k,2) = N1(j)**2 * b3
        k = neta-1
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b4
        lhs(j,2,k,2) = N1(j)**2 * b4
        k = neta
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b5
        lhs(j,2,k,2) = N1(j)**2 * b5

        b1 = dd5 / dN**2
        b2 = dd4 / dN**2
        b3 = dd3 / dN**2
        b4 = dd2 / dN**2
        b5 = dd1 / dN**2

        j = neta
        k = neta-4
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b1
        lhs(j,2,k,2) = N1(j)**2 * b1
        k = neta-3
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b2
        lhs(j,2,k,2) = N1(j)**2 * b2
        k = neta-2
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b3
        lhs(j,2,k,2) = N1(j)**2 * b3
        k = neta-1
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b4
        lhs(j,2,k,2) = N1(j)**2 * b4
        k = neta
        lhs(j,1,k,1) = N1(j)**2/eta(j)**2 * b5
        lhs(j,2,k,2) = N1(j)**2 * b5

!.... first derivatives

        b1 = gc1 / dN
        b2 = gc2 / dN
        b3 = gc3 / dN
        b4 = gc4 / dN
        b5 = gc5 / dN

        j = 1
        k = 1
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b1
        k = 2
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b2
        k = 3
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b3
        k = 4
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b4
        k = 5
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b5
        
        b1 = gb1 / dN
        b2 = gb2 / dN
        b3 = gb3 / dN
        b4 = gb4 / dN
        b5 = gb5 / dN

        j = 2
        k = 1
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b1
        k = 2
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b2
        k = 3
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b3
        k = 4
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b4
        k = 5
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b5

        b1 = ga1  / dN
        b2 = ga2  / dN
        b3 = zero / dN
        b4 = ga3  / dN
        b5 = ga4  / dN

        do j = 3, neta-2
          k = j-2
          lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b1
          lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b1
          lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b1
          k = j-1
          lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b2
          lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b2
          lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b2
          k = j
          lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b3
          lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b3
          lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b3
          k = j+1
          lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b4
          lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b4
          lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b4
          k = j+2
          lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                        N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b5
          lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b5
          lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b5
        end do

        b1 = -gb5 / dN
        b2 = -gb4 / dN
        b3 = -gb3 / dN
        b4 = -gb2 / dN
        b5 = -gb1 / dN

        j = neta-1
        k = neta-4
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b1
        k = neta-3
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b2
        k = neta-2
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b3
        k = neta-1
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b4
        k = neta
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b5

        b1 = -gc5 / dN
        b2 = -gc4 / dN
        b3 = -gc3 / dN
        b4 = -gc2 / dN
        b5 = -gc1 / dN

        j = neta
        k = neta - 4
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b1
        k = neta-3
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b2
        k = neta-2
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b3
        k = neta-1
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b4
        k = neta
        lhs(j,1,k,1) = lhs(j,1,k,1) + ( N2(j)/eta(j)**2 + &
                       N1(j)*(h(j)/eta(j)**2+one/eta(j)-four/eta(j)**3) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) - g(j)/eta(j)**2 * N1(j) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + N2(j) * b5

!.... diagonal terms

        do j = 1, neta
          lhs(j,1,j,1) = lhs(j,1,j,1) - (h1(j)+one)/eta(j)**2 - &
                         two * (h(j)/eta(j)**3 + one/eta(j)**2)
          lhs(j,1,j,2) = lhs(j,1,j,2) + g1(j)/eta(j)**2 - two*g(j)/eta(j)**3
          lhs(j,2,j,1) = lhs(j,2,j,1) - one
        end do

!.... boundary conditions on LHS

        lhs(1,1,:,:) = zero
        lhs(1,1,1,1) = one

        lhs(1,2,:,:) = zero
        lhs(1,2,1,2) = one

        lhs(neta,1,:,:) = zero
        lhs(neta,1,neta,1) = one

        lhs(neta,2,:,:) = zero
        lhs(neta,2,neta  ,2) = -gc1
        lhs(neta,2,neta-1,2) = -gc2
        lhs(neta,2,neta-2,2) = -gc3
        lhs(neta,2,neta-3,2) = -gc4
        lhs(neta,2,neta-4,2) = -gc5

!=============================================================================!
!       S o l v e 
!=============================================================================!

        rhs = -rhs      ! move to rhs
        
#ifdef CRAY
        call SGESV( 2*neta, 1, lhs, 2*neta, ipiv, rhs, 2*neta, info )
#else
        call DGESV( 2*neta, 1, lhs, 2*neta, ipiv, rhs, 2*neta, info )
#endif

        if (info .lt. 0) then
          write(*,*) 'Illegal value at row ',abs(info)
        else if (info .gt. 0) then
          write(*,*) 'Matrix is singular at row ',abs(info)
        endif

!.... update the unknowns

        g = g + rhs(:,1)
        h = h + rhs(:,2)
        
        end do          ! iter
        
!.... compute the norm of the delta

        norm = zero
        do j = 1, neta
          norm = norm + rhs(j,1)**2 + rhs(j,2)**2
        end do
        norm = sqrt( norm / real(neta) )
        
!.... adjust wall vorticity and solve again if neccessary
        
        count = count + 1
        
        call grad1( 1, neta, h, h1, dN, -1, .false., .false. ) 
        h1(:) = h1 * N1

        write(*,"(i5,1x,i5,1x,4(1pe20.13,1x))") i, count, g(1), h1(1), &
                                                norm, xi(i)
        
        if ( abs(h1(1)+one) .le. eps2 ) goto 200
        
        gold2 = gold1
        gold1 = g(1)

        h1old2 = h1old1
        h1old1 = h1(1)
        
        if (gold2 .eq. zero) then
          g(1) = 1.01 * g(1)
          goto 100
        else
          if (h1old1-h1old2 .ne. zero) then
            g(1) = gold2 + (gold1-gold2)/(h1old1-h1old2)*(-one-h1old2)
            goto 100
          else
            write(*,*) 'Could not converge to requested tolerance...'
            goto 200
          end if
        end if
        
200     continue
        do j = 1, neta
          write(20+i,10) N(j), g(j), h(j)
        end do

!.... save the LE solution

        q(:,1,1) = g
        q(:,1,2) = h

        write(11,10) s(i), xi(i), q(1,i,1), q(1,i,2), h1(1)

!=============================================================================!
!       M a r c h   D o w n s t r e a m
!=============================================================================!

        do i = 2, nxi

5000    continue                ! redo with different scale

        count  = 0
        scale  = one
        time   = zero
        gold1  = zero
        h1old1 = zero

!.... predictor stage

        if (i.eq.2) then
          q(:,i,1) = q(:,i-1,1)
          q(:,i,2) = q(:,i-1,2)
        else if (i.eq.3) then
          q(:,i,1) = two * q(:,i-1,1) - q(:,i-2,1)
          q(:,i,2) = two * q(:,i-1,2) - q(:,i-2,2)
!         b1 = (xi(i)-xi(i-1))/(xi(i-2)-xi(i-1))
!         b2 = (xi(i)-xi(i-2))/(xi(i-1)-xi(i-2))
!         q(:,i,1) = b1 * q(:,i-2,1) + b2 * q(:,i-1,1)
!         q(:,i,2) = b1 * q(:,i-2,2) + b2 * q(:,i-1,2)
        else ! if (i.eq.4) then
          q(:,i,1) = three * q(:,i-1,1) - three * q(:,i-2,1) + q(:,i-3,1)
          q(:,i,2) = three * q(:,i-1,2) - three * q(:,i-2,2) + q(:,i-3,2)
!         b1 = (xi(i)-xi(i-1))*(xi(i)-xi(i-2))/(xi(i-3)-xi(i-1))/ &
!              (xi(i-3)-xi(i-2))
!         b2 = (xi(i)-xi(i-1))*(xi(i)-xi(i-3))/(xi(i-2)-xi(i-1))/ &
!              (xi(i-2)-xi(i-3))
!         b3 = (xi(i)-xi(i-2))*(xi(i)-xi(i-3))/(xi(i-1)-xi(i-2))/ &
!              (xi(i-1)-xi(i-3))
!         q(:,i,1) = b1 * q(:,i-3,1) + b2 * q(:,i-2,1) + b3 * q(:,i-1,1)
!         q(:,i,2) = b1 * q(:,i-3,2) + b2 * q(:,i-2,2) + b3 * q(:,i-1,2)
!       else
!         q(:,i,1) = 4.0 * q(:,i-1,1) - 6.0 * q(:,i-2,1) + &
!                    4.0 * q(:,i-3,1) - q(:,i-4,1)
!         q(:,i,2) = 4.0 * q(:,i-1,2) - 6.0 * q(:,i-2,2) + &
!                    4.0 * q(:,i-3,2) - q(:,i-4,2)
        end if

        if (gtmp.ne.zero) then
          q(1,i,1) = gtmp*q(1,i,1)
          write(*,*) 'Bad initial guess, trying ',q(1,i,1)
        end if
        
        gbest = q(:,i,1)
        hbest = q(:,i,2)
!=============================================================================!
!       M a i n   i t e r a t i o n
!=============================================================================!
1000    continue
!       if (count .gt. niter) goto 4000
        if (count .gt. niter) goto 2000
        
        norm = infty
        
        do iter = 1, niter      !.... Newton iterations

        if (midpnt) then        !.... midpoint rule
          alpha = pt5
          g = pt5 * ( q(:,i,1) + q(:,i-1,1) )
          h = pt5 * ( q(:,i,2) + q(:,i-1,2) )
          sm  = pt5 * ( s(i) + s(i-1) )
          if (stretch) then
            xim = calc_xi( sm, A, xi(i-1)+0.00001, xi(i) )
            s1m = s1( xim, A )
          else
            xim = sm
            s1m = one
          end if
        else                    !.... backward Euler
          alpha = one
          g = q(:,i,1)
          h = q(:,i,2)
          sm  = s(i)
          xim = xi(i)
          if (stretch) then
            s1m = s1( xim, A )
          else
            s1m = one
          end if
        end if
        
!.... compute the delta values

        dg = q(:,i,1) - q(:,i-1,1)
        dh = q(:,i,2) - q(:,i-1,2)
        
!.... satisfy the boundary conditions

        h(1) = -Sqrt(R)
        
        g(neta) = zero
        h(neta) = -( gc5 * h(neta-4)  + &
                     gc4 * h(neta-3)  + &
                     gc3 * h(neta-2)  + &
                     gc2 * h(neta-1)  ) / gc1

!.... compute first derivatives in computational space

        call grad1( 1, neta, g, g1, dN, -1, .false., .false. ) 
        call grad1( 1, neta, h, h1, dN, -1, .false., .false. ) 

!.... compute second derivatives in computational space

        call grad2( 1, neta, g, g2, dN, -1, .false., .false. ) 
        call grad2( 1, neta, h, h2, dN, -1, .false., .false. ) 

!.... transfer to physical coordinates

        g2(:) = g2 * N1**2 + g1 * N2
        g1(:) = g1 * N1
        
        h2(:) = h2 * N1**2 + h1 * N2
        h1(:) = h1 * N1

!=============================================================================!
!       F o r m   t h e   R H S  
!=============================================================================!

        jac = xim**2 + eta**2
        
        rhs(:,1) = -(xim*(h1+one) + four*xim/jac) * s1m/ds * dg + g2 + &
                    (h + eta + xim*s1m/ds*dh - four*eta/jac) * g1 + &
                    one/jac * ( (xim**2-eta**2)*(h1+one) - &
                    two*eta*(h + eta + xim*s1m/ds*dh) ) * g
        rhs(:,2) = h2 - g

        if (delt.ne.zero) then
          rhs(:,1) = -delt * rhs(:,1)
          rhs(:,2) = -delt * rhs(:,2)
        end if
        
!..... boundary conditions on the rhs

        rhs(1,1) = zero
        rhs(1,2) = zero

        rhs(neta,1) = zero
        rhs(neta,2) = -gc1 * h(neta  ) - gc2 * h(neta-1) &
                      -gc3 * h(neta-2) - gc4 * h(neta-3) &
                      -gc5 * h(neta-4)

!.... compute the norm of the residual

        normold = norm
        norm = zero
        do j = 1, neta
          norm = norm + rhs(j,1)**2 + rhs(j,2)**2
        end do
        norm = sqrt( norm / real(neta) )

!       write(*,"(i3,1x,1pe20.13)") iter, norm

        if (norm .lt. eps1) goto 3000
!       if (norm .gt. normold) then
!         write(*,*) 'Newton iteration failed to converge at i = ',i
!         call exit(1)
!       end if
        
!=============================================================================!
!       F o r m   t h e   L H S 
!=============================================================================!

        lhs = zero

!.... second derivatives
        
        b1 = dd1 / dN**2
        b2 = dd2 / dN**2
        b3 = dd3 / dN**2
        b4 = dd4 / dN**2
        b5 = dd5 / dN**2

        j = 1
        k = 1
        lhs(j,1,k,1) = alpha * N1(j)**2 * b1
        lhs(j,2,k,2) = alpha * N1(j)**2 * b1
        k = 2
        lhs(j,1,k,1) = alpha * N1(j)**2 * b2
        lhs(j,2,k,2) = alpha * N1(j)**2 * b2
        k = 3
        lhs(j,1,k,1) = alpha * N1(j)**2 * b3
        lhs(j,2,k,2) = alpha * N1(j)**2 * b3
        k = 4
        lhs(j,1,k,1) = alpha * N1(j)**2 * b4
        lhs(j,2,k,2) = alpha * N1(j)**2 * b4
        k = 5
        lhs(j,1,k,1) = alpha * N1(j)**2 * b5
        lhs(j,2,k,2) = alpha * N1(j)**2 * b5

        b1 = db1 / dN**2
        b2 = db2 / dN**2
        b3 = db3 / dN**2
        b4 = db4 / dN**2
        b5 = db5 / dN**2

        j = 2
        k = 1
        lhs(j,1,k,1) = alpha * N1(j)**2 * b1
        lhs(j,2,k,2) = alpha * N1(j)**2 * b1
        k = 2
        lhs(j,1,k,1) = alpha * N1(j)**2 * b2
        lhs(j,2,k,2) = alpha * N1(j)**2 * b2
        k = 3
        lhs(j,1,k,1) = alpha * N1(j)**2 * b3
        lhs(j,2,k,2) = alpha * N1(j)**2 * b3
        k = 4
        lhs(j,1,k,1) = alpha * N1(j)**2 * b4
        lhs(j,2,k,2) = alpha * N1(j)**2 * b4
        k = 5
        lhs(j,1,k,1) = alpha * N1(j)**2 * b5
        lhs(j,2,k,2) = alpha * N1(j)**2 * b5

        b1 =  da1 / dN**2
        b2 =  da2 / dN**2
        b3 =  da3 / dN**2
        b4 =  da4 / dN**2
        b5 =  da5 / dN**2

        do j = 3, neta-2
          k = j-2
          lhs(j,1,k,1) = alpha * N1(j)**2 * b1
          lhs(j,2,k,2) = alpha * N1(j)**2 * b1
          k = j-1
          lhs(j,1,k,1) = alpha * N1(j)**2 * b2
          lhs(j,2,k,2) = alpha * N1(j)**2 * b2
          k = j
          lhs(j,1,k,1) = alpha * N1(j)**2 * b3
          lhs(j,2,k,2) = alpha * N1(j)**2 * b3
          k = j+1
          lhs(j,1,k,1) = alpha * N1(j)**2 * b4
          lhs(j,2,k,2) = alpha * N1(j)**2 * b4
          k = j+2
          lhs(j,1,k,1) = alpha * N1(j)**2 * b5
          lhs(j,2,k,2) = alpha * N1(j)**2 * b5
        end do

        b1 = db5 / dN**2
        b2 = db4 / dN**2
        b3 = db3 / dN**2
        b4 = db2 / dN**2
        b5 = db1 / dN**2

        j = neta-1
        k = neta-4
        lhs(j,1,k,1) = alpha * N1(j)**2 * b1
        lhs(j,2,k,2) = alpha * N1(j)**2 * b1
        k = neta-3
        lhs(j,1,k,1) = alpha * N1(j)**2 * b2
        lhs(j,2,k,2) = alpha * N1(j)**2 * b2
        k = neta-2
        lhs(j,1,k,1) = alpha * N1(j)**2 * b3
        lhs(j,2,k,2) = alpha * N1(j)**2 * b3
        k = neta-1
        lhs(j,1,k,1) = alpha * N1(j)**2 * b4
        lhs(j,2,k,2) = alpha * N1(j)**2 * b4
        k = neta
        lhs(j,1,k,1) = alpha * N1(j)**2 * b5
        lhs(j,2,k,2) = alpha * N1(j)**2 * b5

        b1 = dd5 / dN**2
        b2 = dd4 / dN**2
        b3 = dd3 / dN**2
        b4 = dd2 / dN**2
        b5 = dd1 / dN**2

        j = neta
        k = neta-4
        lhs(j,1,k,1) = alpha * N1(j)**2 * b1
        lhs(j,2,k,2) = alpha * N1(j)**2 * b1
        k = neta-3
        lhs(j,1,k,1) = alpha * N1(j)**2 * b2
        lhs(j,2,k,2) = alpha * N1(j)**2 * b2
        k = neta-2
        lhs(j,1,k,1) = alpha * N1(j)**2 * b3
        lhs(j,2,k,2) = alpha * N1(j)**2 * b3
        k = neta-1
        lhs(j,1,k,1) = alpha * N1(j)**2 * b4
        lhs(j,2,k,2) = alpha * N1(j)**2 * b4
        k = neta
        lhs(j,1,k,1) = alpha * N1(j)**2 * b5
        lhs(j,2,k,2) = alpha * N1(j)**2 * b5

!.... first derivatives

        b1 = gc1 / dN
        b2 = gc2 / dN
        b3 = gc3 / dN
        b4 = gc4 / dN
        b5 = gc5 / dN

        j = 1
        k = 1
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b1
        k = 2
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b2
        k = 3
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b3
        k = 4
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b4
        k = 5
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b5
        
        b1 = gb1 / dN
        b2 = gb2 / dN
        b3 = gb3 / dN
        b4 = gb4 / dN
        b5 = gb5 / dN

        j = 2
        k = 1
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b1
        k = 2
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b2
        k = 3
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b3
        k = 4
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b4
        k = 5
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b5

        b1 = ga1  / dN
        b2 = ga2  / dN
        b3 = zero / dN
        b4 = ga3  / dN
        b5 = ga4  / dN

        do j = 3, neta-2
          k = j-2
          lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                          xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b1
          lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                          g(j) - xim*s1m/ds*dg(j) ) * b1
          lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b1
          k = j-1
          lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                          xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b2
          lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                          g(j) - xim*s1m/ds*dg(j) ) * b2
          lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b2
          k = j
          lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                          xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b3
          lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                          g(j) - xim*s1m/ds*dg(j) ) * b3
          lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b3
          k = j+1
          lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                          xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b4
          lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                          g(j) - xim*s1m/ds*dg(j) ) * b4
          lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b4
          k = j+2
          lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                          xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b5
          lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                          g(j) - xim*s1m/ds*dg(j) ) * b5
          lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b5
        end do

        b1 = -gb5 / dN
        b2 = -gb4 / dN
        b3 = -gb3 / dN
        b4 = -gb2 / dN
        b5 = -gb1 / dN

        j = neta-1
        k = neta-4
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b1
        k = neta-3
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b2
        k = neta-2
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b3
        k = neta-1
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b4
        k = neta
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b5

        b1 = -gc5 / dN
        b2 = -gc4 / dN
        b3 = -gc3 / dN
        b4 = -gc2 / dN
        b5 = -gc1 / dN

        j = neta
        k = neta-4
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b1
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b1
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b1
        k = neta-3
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b2
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b2
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b2
        k = neta-2
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b3
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b3
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b3
        k = neta-1
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b4
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b4
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b4
        k = neta
        lhs(j,1,k,1) = lhs(j,1,k,1) + alpha * ( N2(j) + N1(j)*( h(j) + eta(j) + &
                       xim*s1m/ds*dh(j) - four*eta(j)/jac(j) ) ) * b5
        lhs(j,1,k,2) = lhs(j,1,k,2) + alpha * N1(j) * ( (xim**2-eta(j)**2)/jac(j) * &
                       g(j) - xim*s1m/ds*dg(j) ) * b5
        lhs(j,2,k,2) = lhs(j,2,k,2) + alpha * N2(j) * b5

!.... diagonal terms

        do j = 1, neta
          lhs(j,1,j,1) = lhs(j,1,j,1) + ( alpha/jac(j)*( (xim**2-eta(j)**2)* &
                         (h1(j)+one) - two*eta(j)*(h(j)+eta(j)+xim*s1m/ds*dh(j)) ) &
                         -xim*s1m/ds*( (h1(j)+one) + four/jac(j) ) )
          lhs(j,1,j,2) = lhs(j,1,j,2) + ( (alpha + xim*s1m/ds)*g1(j) - &
                         two*eta(j)*g(j)/jac(j)*(alpha + xim*s1m/ds) )
          lhs(j,2,j,1) = lhs(j,2,j,1) - alpha
        end do

!.... divide the first equation by Jac

!       do k = 1, neta
!         lhs(:,1,k,1) = lhs(:,1,k,1) / Jac
!         lhs(:,1,k,2) = lhs(:,1,k,2) / Jac
!       end do
        
!.... add in the time term

        if (delt.ne.zero) then
          lhs = -delt * lhs
          do j = 1, neta
            lhs(j,1,j,1) = one + lhs(j,1,j,1)
            lhs(j,2,j,2) = one + lhs(j,2,j,2)
          end do
        end if
        
!.... boundary conditions on LHS

        lhs(1,1,:,:) = zero
        lhs(1,1,1,1) = one

        lhs(1,2,:,:) = zero
        lhs(1,2,1,2) = one

        lhs(neta,1,:,:) = zero
        lhs(neta,1,neta,1) = one

        lhs(neta,2,:,:) = zero
        lhs(neta,2,neta  ,2) = -gc1
        lhs(neta,2,neta-1,2) = -gc2
        lhs(neta,2,neta-2,2) = -gc3
        lhs(neta,2,neta-3,2) = -gc4
        lhs(neta,2,neta-4,2) = -gc5

!=============================================================================!
!       S o l v e 
!=============================================================================!

        rhs = -rhs      ! move to rhs
#ifdef CRAY        
        call SGESV( 2*neta, 1, lhs, 2*neta, ipiv, rhs, 2*neta, info )
#else
        call DGESV( 2*neta, 1, lhs, 2*neta, ipiv, rhs, 2*neta, info )
#endif
        if (info .lt. 0) then
          write(*,*) 'Illegal value at row ',abs(info)
          call exit(1)
        else if (info .gt. 0) then
          write(*,*) 'Matrix is singular at row ',abs(info)
          call exit(1)
        endif

!.... update the unknowns

        q(:,i,1) = q(:,i,1) + rhs(:,1)
        q(:,i,2) = q(:,i,2) + rhs(:,2)
        
        end do          ! iter
        
!.... adjust wall vorticity and solve again if neccessary
        
3000    continue

        count = count + 1
        
        if (norm.gt.eps1) then
          if (count.eq.1) then
            if (gtmp.eq.zero) then
              gtmp = 0.9991
            else
              gtmp = gtmp * 0.9991
            end if
            goto 5000
          else
            q(:,i,1) = gbest
            q(:,i,2) = hbest
            count = 0
            scale = scale * two
            if (scale.gt.2048) goto 4000        ! exit
            write(*,*) 'scale = ',scale
            goto 1000
          end if
        end if
        gtmp = zero

        g = q(:,i,1)
        h = q(:,i,2)

        gbest = g
        hbest = h
                
        call grad1( 1, neta, h, h1, dN, -1, .false., .false. ) 
        h1(:) = h1 * N1

        if ( count.eq.1 .or. (-one-h1(1))*(-one-h1old1).lt.zero ) then
          gold2 = gold1
          gold1 = g(1)
          h1old2 = h1old1
          h1old1 = h1(1)
        else if ( count .eq. 2 ) then
          gold2 = gold1
          gold1 = g(1)
          h1old2 = h1old1
          h1old1 = h1(1)
        else
          gold1 = g(1)
          h1old1 = h1(1)
        end if    

!       write(*,"(i5,1x,i5,1x,4(1pe20.13,1x))") i, count, g(1), h1(1), h1old1, h1old2
        write(*,"(i5,1x,i5,1x,4(1pe20.13,1x))") i, count, g(1), h1(1), norm, xi(i)
        if ( abs(h1(1)+one) .le. eps2 ) goto 2000
        
        if (gold2 .eq. zero .or. (-one-h1old1)*(-one-h1old2).gt.zero) then
          if (count.gt.1 .and. abs(-one-h1old1).gt.abs(-one-h1old2) ) then
            write(*,*) 'switching polarity'
            isign = -isign
          end if
          g(1) = (one + isign * (-one-h1(1))/scale ) * g(1)
          q(1,i,1) = g(1)
          goto 1000
        else
          if (h1old1-h1old2 .ne. zero) then
            g(1) = gold2 + (gold1-gold2)/(h1old1-h1old2)*(-one-h1old2)
            q(1,i,1) = g(1)
            goto 1000
          else
            write(*,*) 'Could not converge to requested tolerance...'
            goto 2000
          end if
        end if
        
2000    continue

        write(11,10) s(i), xi(i), q(1,i,1), q(1,i,2), h1(1)
        
!       do j = 1, neta
!         write(20+i,10) N(j), g(j), h(j)
!       end do

        end do  ! i

4000    continue

        if (i.eq.nxi+1) then
          nxi2 = nxi
        else
          nxi2 = i
        end if
        
!.... output the results in a plot3d file

        if (.false.) then
          open(90,file='pcyl.dat',form='unformatted')
          write(90) nxi, neta, 1
          write(90) (( s(i), i=1,nxi), j=1,neta), &
                    (( n(j), i=1,nxi), j=1,neta), &
                    (( zero, i=1,nxi), j=1,neta)
          close(90)
          
          open(90,file='flow.dat',form='unformatted')
          write(90) nxi, neta, 1
          write(90) zero, zero, zero, zero
          write(90) ((      one, i=1,nxi), j=1,neta), &
                    (( q(j,i,1), i=1,nxi), j=1,neta), &
                    (( q(j,i,2), i=1,nxi), j=1,neta), &
                    ((     zero, i=1,nxi), j=1,neta), &
                    ((      one, i=1,nxi), j=1,neta)
          close(90)
        end if
        
!.... allocate memory for the physical grid

        open (unit=90, file='grid.dat', form='unformatted', status='old')
        read(90) nx, ny, nz
        write(*,"('Nx = ',i4,', Ny = ',i4,', Nz = ',i4)") nx, ny, nz
        allocate ( xy(ny,nx,2), v(ny,nx,5) )
        read(90) (((xy(j,i,1), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((xy(j,i,2), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((      tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(90)

!.... Spline the results

        korder = 5      ! five is adequate
!       write(*,"('Enter korder ==> ',$)")
!       read(*,*) korder
        write(*,"('Enter Ma ==> ',$)")
        read(*,*) Ma

        allocate( xknot(nxi2+korder), yknot(neta+korder), bs(neta,nxi2,2) )
        
        call BSNAK( nxi2, s, korder, xknot)
        call BSNAK( neta, N, korder, yknot)
#ifdef USE_IMSL
        call BS2IN( neta, N, nxi2, s, q(:,1:nxi2,1), neta, &
                    korder, korder, yknot, xknot, bs(:,:,1) )
        call BS2IN( neta, N, nxi2, s, q(:,1:nxi2,2), neta, &
                    korder, korder, yknot, xknot, bs(:,:,2) )
#else
        write(*,*) "Need to implement B-spline routines"
        stop
#endif

!.... interpolate to the physical mesh

        do i = 1, nx
          write(*,*) 'i = ',i
          do j = 1, ny
            xl   = xy(j,i,1)
            yl   = xy(j,i,2)
            etal = Sqrt( R*(pt5-xl) + pt5*Sqrt( (two*R*xl-R)**2 + &
                   four * R**2 * yl**2 ) )
            xil  = R * yl / etal
            if (xil .eq. zero) then
              sl = zero
            else
              sl = one - A / xil * log( one + xil / A )
            end if
            Nl = (etal - Sqrt(R)) / (5.0 + etal - Sqrt(R))
            if (Nl .lt. zero) Nl = zero
            if (sl.le.s(nxi2) .and. Nl.le.one) then
#ifdef USE_IMSL
              gl  = BS2DR( 0, 0, Nl, sl, korder, korder, yknot, xknot, &
                           neta, nxi2, bs(:,:,1) )
              hl  = BS2DR( 0, 0, Nl, sl, korder, korder, yknot, xknot, &
                           neta, nxi2, bs(:,:,2) )
              g1l = BS2DR( 0, 1, Nl, sl, korder, korder, yknot, xknot, &
                           neta, nxi2, bs(:,:,1) )
              g2l = BS2DR( 1, 0, Nl, sl, korder, korder, yknot, xknot, &
                           neta, nxi2, bs(:,:,1) )
              h1l = BS2DR( 0, 1, Nl, sl, korder, korder, yknot, xknot, &
                           neta, nxi2, bs(:,:,2) )
              h2l = BS2DR( 1, 0, Nl, sl, korder, korder, yknot, xknot, &
                           neta, nxi2, bs(:,:,2) )
#else
              write(*,*) "Need to implement B-spline routines"
              stop
#endif
            end if

            fl  = hl + etal
            if (i.eq.1) then
              f1l = zero
            else
              f1l = s1(xil,A) * h1l
            end if
            f2l = (one-Nl)**2/5.0 * h2l + one
            
            jl = (xil**2 + etal**2)/R
            m1l =  xil/jl
            m2l =  etal/jl
            n1l = -etal/jl
            n2l =  xil/jl
            
            psil = xil * fl / R
            vorl = -xil * gl / (xil**2 + etal**2) / R
            ul   =  ( m2l * fl + xil * ( m2l * f1l + n2l * f2l ) ) / R
            vl   = -( m1l * fl + xil * ( m1l * f1l + n1l * f2l ) ) / R
            tl   = one - pt5 * gamma1 * Ma**2 * ( ul**2 + vl**2 - one )

!.... now switch to inviscid, potential flow to get pressure and density

            etal = sqrt( pt5 - xl + pt5 * sqrt( (two*xl-one)**2 + &
                    four * yl**2 ) )
            xil  = yl / etal
            pl   = one/(gamma*Ma**2) * ( one - pt5 * gamma1 * Ma**2 * ( &
                   ((etal**2 - etal + xil**2) / (xil**2 + etal**2))**2 + &
                   (xil / (xil**2 + etal**2))**2 - one ) ) ** (gamma/gamma1)
            rhol = gamma * Ma**2 * pl / tl
            
            v(j,i,1) = rhol
            v(j,i,2) = ul
            v(j,i,3) = vl
            v(j,i,4) = zero
            v(j,i,5) = tl
          end do
        end do

!.... write out a restart file
        
        open(90,file='par.dat',form='unformatted')
        write(90) 0, 0.0, nx, ny, nz, 5,  &
                  R, Ma, one, gamma, cv
        write(90) v
        close(90)

        stop
10      format(1p,8(e13.6,1x))
        end

!=============================================================================!
        function calc_xi( s, A, xi1, xi2 )
!=============================================================================!
        implicit none
        real :: calc_xi, s, A, xi1, xi2
        real, external :: func, rtflsp
        
        if (s .eq. 0.0) then
          calc_xi = 0.0
          return
        else if (s .eq. 1.0) then
          calc_xi = 1.0e99
          return
        end if
        
        calc_xi = rtflsp(func, xi1, xi2, 1.0d-8, s, A)

        return
        end

!=============================================================================!
        function func( xi, s, A )
!=============================================================================!
        implicit none
        real :: func, xi, s, A

        func = s - 1.0 + A / xi * log( 1.0 + xi / A)    
        return
        end

!=============================================================================!
        function s1( xi, A )
!=============================================================================!
        implicit none
        real :: s1, xi, A

        if (xi .eq. 0.0) then
          write(*,*) 'ERROR:  xi = 0 in s1'
          call exit(1)
        end if
        
        s1 = A/xi**2 * log( 1.0 + xi/A) - 1.0/(xi*(1.0 + xi/A))
                
        return
        end

!=============================================================================!
        function s2( xi, A )
!=============================================================================!
        implicit none
        real :: s2, xi, A

        if (xi .eq. 0.0) then
          write(*,*) 'ERROR:  xi = 0 in s2'
          call exit(1)
        end if

        s2 = 1.0/(A*xi*(1.0+xi/A)**2) + 2.0/(xi**2*(1.0+xi/A)) - &
             2.0*A*log(1.0+xi/A)/xi**3
                
        return
        end
