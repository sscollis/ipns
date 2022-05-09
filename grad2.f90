!=============================================================================!
	subroutine grad2( ndof, ny, v, g22v, dy, opty, yper, carp)
!
!  Take the second derivative of a 2-D field.
!  Updated to sixth-order accurate differencing on the interior with
!  the option for optimized fourth-order differencing.
!
!  Revised: 6-28-95
!
!=============================================================================!
	use diff
	implicit none
	
	integer :: ndof, ny, opty
	logical :: yper
	logical :: carp
	real    :: v(ny,ndof)
	real    :: g22v(ny,ndof)
	real    :: dy

	real, parameter :: zero = 0.0, one = 1.0, pt5 = 0.5, two = 2.0
	real :: dyinv, dysinv
	real :: a, b, c, w
	real :: gy1, gy2, gy3, gy4, gy5, gy6
	real :: dy1, dy2, dy3, dy4, dy5, dy6, dy7
	
	integer :: i, j, idof, isign
	real :: eps = 1.0e-12
!=============================================================================!

	dyinv  = one / dy
	dysinv = one / dy**2

!.... seven point stencil in y

	if (opty.eq.0) then
	  c = 1.0 / 60.0
	else if (opty.eq.-1) then
	  c = 0.0
	else
	  w = 2.0 * 3.1415926535897932385e+0 / 12.0
	  c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
	      (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
	end if
	
	a = (2.0 + 15.0 * c) / 3.0
	b = -(1.0 + 48.0 * c) / 12.0

	gy1 =  -c * dyinv
	gy2 =  -b * dyinv
	gy3 =  -a * dyinv
	gy4 =   a * dyinv
	gy5 =   b * dyinv
	gy6 =   c * dyinv

	if (opty.eq.0) then
	  c = 1.0 / 90.0
	else if (opty.eq.-1) then
	  c = 0.0
	else
	  w = 2.0 * 3.1415926535897932385e+0 / 12.0
	  c = -(3.0*w**2 - 16.0*Sin(w/2.0)**2 + Sin(w)**2) / &
               (12.0*(-15.0*Sin(w/2.0)**2 + 6.0*Sin(w)**2 -  &
               Sin(3.0*w/2.0)**2))
	end if
	
	a = (4.0 + 45.0 * c) / 3.0
	b = -(1.0 + 72.0 * c) / 12.0

	dy1 =  c * dysinv
	dy2 =  b * dysinv
	dy3 =  a * dysinv
	dy4 = -2.0 * ( a + b + c ) * dysinv
	dy5 =  a * dysinv
	dy6 =  b * dysinv
	dy7 =  c * dysinv

!=============================================================================!
!.... compute the second derivative in y
!=============================================================================!

	if (yper) then
	
	  g22v(1,:)       = ( dy1 * v(ny-3,:)   + &
				dy2 * v(ny-2,:)   + &
				dy3 * v(ny-1,:)   + &
				dy4 * v(1,:)      + &
				dy5 * v(2,:)	 + &
				dy6 * v(3,:)	 + &
				dy7 * v(4,:)      ) 
  
	  g22v(2,:)       = ( dy1 * v(ny-2,:)   + &
				dy2 * v(ny-1,:)   + &
				dy3 * v(1,:)      + &
				dy4 * v(2,:)      + &
				dy5 * v(3,:)      + &
				dy6 * v(4,:)      + &
				dy7 * v(5,:)      ) 
  
	  g22v(3,:)       = ( dy1 * v(ny-1,:)   + &
				dy2 * v(1,:)      + &
				dy3 * v(2,:)      + &
				dy4 * v(3,:)      + &
				dy5 * v(4,:)      + &
				dy6 * v(5,:)      + &
				dy7 * v(6,:)      ) 
  
	  g22v(ny-2,:)    = ( dy1 * v(ny-5,:)   + &
				dy2 * v(ny-4,:)   + &
				dy3 * v(ny-3,:)   + &
				dy4 * v(ny-2,:)   + &
				dy5 * v(ny-1,:)   + &
				dy6 * v(1,:)      + &
				dy7 * v(2,:)      ) 
  
	  g22v(ny-1,:)    = ( dy1 * v(ny-4,:)   + &
				dy2 * v(ny-3,:)   + &
				dy3 * v(ny-2,:)   + &
				dy4 * v(ny-1,:)   + &
				dy5 * v(1,:)      + &
				dy6 * v(2,:)      + &
				dy7 * v(3,:)   ) 
  
	  g22v(ny,:) = g22v(1,:)

	else
	
	  g22v(1,:)      =  ( dd1 * v(1,:) + &
				dd2 * v(2,:) + &
				dd3 * v(3,:) + &
				dd4 * v(4,:) + &
				dd5 * v(5,:) ) * dysinv

	  g22v(2,:)      =  ( db1 * v(1,:) + &
				db2 * v(2,:) + &
				db3 * v(3,:) + &
				db4 * v(4,:) + &
				db5 * v(5,:) ) * dysinv
  
	  g22v(3,:)       = ( da1 * v(1,:) + &
				da2 * v(2,:) + &
				da3 * v(3,:) + &
				da4 * v(4,:) + &
				da5 * v(5,:) ) * dysinv
  
	  g22v(ny-2,:)    = ( da1 * v(ny-4,:) + &
				da2 * v(ny-3,:) + &
				da3 * v(ny-2,:) + &
				da4 * v(ny-1,:) + &
				da5 * v(ny,:)   ) * dysinv
  
	  g22v(ny-1,:)   =  ( db1 * v(ny,:)   + &
				db2 * v(ny-1,:) + &
				db3 * v(ny-2,:) + &
				db4 * v(ny-3,:) + &
				db5 * v(ny-4,:) ) * dysinv
  
	  g22v(ny,:)     =  ( dd1 * v(ny,:)   + &
				dd2 * v(ny-1,:) + &
				dd3 * v(ny-2,:) + &
				dd4 * v(ny-3,:) + &
				dd5 * v(ny-4,:) ) * dysinv
			       
	endif

!.... interior

	g22v(4:ny-3,:) = ( dy1 * v(1:ny-6,:)   + &
	                     dy2 * v(2:ny-5,:)   + &
	                     dy3 * v(3:ny-4,:)   + &
	                     dy4 * v(4:ny-3,:)   + &
		             dy5 * v(5:ny-2,:)   + &
			     dy6 * v(6:ny-1,:)   + &
			     dy7 * v(7:ny  ,:)   ) 

	return
	end
