!=============================================================================!
	subroutine grad1( ndof, ny, v, g2v, dy, opty, yper, carp)
!
!  Take the gradient of a 2-D field in the y-direction.
!  Updated to sixth-order accurate differencing on the interior with
!  the option for optimized fourth-order differencing.
!
!=============================================================================!
	use diff
	implicit none
	
	integer :: ndof, ny, opty
	logical :: yper, carp
	real    :: v(ny,ndof), g2v(ny,ndof)
	real    :: dy
	
	real, parameter :: zero = 0.0, one = 1.0, pt5 = 0.5
	real dyinv
	real a, b, c, w
	real gy1, gy2, gy3, gy4, gy5, gy6
	
	integer :: i, j, idof, isign
	real :: eps = 1.0e-12
!=============================================================================!
	dyinv  = one / dy

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

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

	if (yper) then
	
	  g2v(1,:)       = ( gy1 * v(ny-3,:)	+ &
			     gy2 * v(ny-2,:)	+ &
			     gy3 * v(ny-1,:)	+ &
			     gy4 * v(2,:)	+ &
			     gy5 * v(3,:)	+ &
			     gy6 * v(4,:)  	) 
  
	  g2v(2,:)       = ( gy1 * v(ny-2,:)	+ &
			     gy2 * v(ny-1,:)	+ &
			     gy3 * v(1,:)	+ &
			     gy4 * v(3,:)	+ &
			     gy5 * v(4,:)	+ &
			     gy6 * v(5,:)	) 
  
	  g2v(3,:)       = ( gy1 * v(ny-1,:)	+ &
			     gy2 * v(1,:)	+ &
			     gy3 * v(2,:)	+ &
			     gy4 * v(4,:)	+ &
			     gy5 * v(5,:)	+ &
			     gy6 * v(6,:)	) 
  
	  g2v(ny-2,:)    = ( gy1 * v(ny-5,:)	+ &
			     gy2 * v(ny-4,:)	+ &
			     gy3 * v(ny-3,:)	+ &
			     gy4 * v(ny-1,:)	+ &
			     gy5 * v(1,:)	+ &
			     gy6 * v(2,:)	) 
  
	  g2v(ny-1,:)    = ( gy1 * v(ny-4,:)	+ &
			     gy2 * v(ny-3,:)	+ &
			     gy3 * v(ny-2,:)	+ &
			     gy4 * v(1,:)	+ &
			     gy5 * v(2,:)	+ &
			     gy6 * v(3,:)	) 
  
	  g2v(ny,:) = g2v(1,:)

	else

	  g2v(1,:)     = ( gc1 * v(1,:)	+ &
			   gc2 * v(2,:)  + &
			   gc3 * v(3,:)  + &
			   gc4 * v(4,:)  + &
			   gc5 * v(5,:)  ) * dyinv

	  g2v(2,:)     = ( gb1 * v(1,:)  + &
			   gb2 * v(2,:)  + &
			   gb3 * v(3,:)  + &
			   gb4 * v(4,:)  + &
			   gb5 * v(5,:)  ) * dyinv
  
	  g2v(3,:)     = ( ga1 * v(1,:)  + &
			   ga2 * v(2,:)  + &
			   ga3 * v(4,:)  + &
			   ga4 * v(5,:)  ) * dyinv
  
	  g2v(ny-2,:)  = ( ga1 * v(ny-4,:)  + &
			   ga2 * v(ny-3,:)  + &
			   ga3 * v(ny-1,:)  + &
			   ga4 * v(ny  ,:)  ) * dyinv
  
	  g2v(ny-1,:) = -( gb1 * v(ny  ,:)  + &
			   gb2 * v(ny-1,:)  + &
			   gb3 * v(ny-2,:)  + &
			   gb4 * v(ny-3,:)  + &
			   gb5 * v(ny-4,:)  ) * dyinv
  
	  g2v(ny,:)   = -( gc1 * v(ny  ,:)  + &
			   gc2 * v(ny-1,:)  + &
			   gc3 * v(ny-2,:)  + &
			   gc4 * v(ny-3,:)  + &
			   gc5 * v(ny-4,:)  ) * dyinv
			    
	end if
	
!.... interior

	g2v(4:ny-3,:) = ( gy1 * v(1:ny-6,:)	+ &
			  gy2 * v(2:ny-5,:)	+ &
			  gy3 * v(3:ny-4,:)	+ &
			  gy4 * v(5:ny-2,:)	+ &
			  gy5 * v(6:ny-1,:)	+ &
			  gy6 * v(7:ny  ,:)	) 

!.... Implement Carpenter's boundary stencil

	if (carp .and. (.not. yper) ) then
	  g2v(1,:)     = ( gg1 * v(1,:)	 + &
			   gg2 * v(2,:)  + &
			   gg3 * v(3,:)  + &
			   gg4 * v(4,:)  + &
			   gg5 * v(5,:)  + &
			   gg6 * v(6,:)  ) * dyinv
  
	  g2v(2,:)     = ( gh1 * v(1,:)  + &
			   gh2 * v(2,:)  + &
			   gh3 * v(3,:)  + &
			   gh4 * v(4,:)  + &
			   gh5 * v(5,:)  + &
			   gh6 * v(6,:)  ) * dyinv

	  g2v(3,:)     = ( gi1 * v(1,:)  + &
			   gi2 * v(2,:)  + &
			   gi3 * v(3,:)  + &
			   gi4 * v(4,:)  + &
			   gi5 * v(5,:)  + &
			   gi6 * v(6,:)  ) * dyinv

	  g2v(4,:)     = ( gj1 * v(1,:)  + &
			   gj2 * v(2,:)  + &
			   gj3 * v(3,:)  + &
			   gj4 * v(4,:)  + &
			   gj5 * v(5,:)  + &
			   gj6 * v(6,:)  ) * dyinv

	  g2v(ny-3,:) = -( gj1 * v(ny,:)    + &
			   gj2 * v(ny-1,:)  + &
			   gj3 * v(ny-2,:)  + &
			   gj4 * v(ny-3,:)  + &
			   gj5 * v(ny-4,:)  + &
			   gj6 * v(ny-5,:)  ) * dyinv

	  g2v(ny-2,:) = -( gi1 * v(ny,:)    + &
			   gi2 * v(ny-1,:)  + &
			   gi3 * v(ny-2,:)  + &
			   gi4 * v(ny-3,:)  + &
			   gi5 * v(ny-4,:)  + &
			   gi6 * v(ny-5,:)  ) * dyinv

	  g2v(ny-1,:) = -( gh1 * v(ny,:)    + &
			   gh2 * v(ny-1,:)  + &
			   gh3 * v(ny-2,:)  + &
			   gh4 * v(ny-3,:)  + &
			   gh5 * v(ny-4,:)  + &
			   gh6 * v(ny-5,:)  ) * dyinv

	  g2v(ny,:)   = -( gg1 * v(ny,:)    + &
			   gg2 * v(ny-1,:)  + &
			   gg3 * v(ny-2,:)  + &
			   gg4 * v(ny-3,:)  + &
			   gg5 * v(ny-4,:)  + &
			   gg6 * v(ny-5,:)  ) * dyinv
	end if

	return
	end
