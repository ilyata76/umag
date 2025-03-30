subroutine initial_magnetization(inmag_type,m0)

use magmod

implicit real*8 (a-h,o-z)

character*15 inmag_type

type (vector) m0 ! direction of initial magnetization 

type (vector) v,norm_vector

!---------------------------generate initial magnetization----------------
write(8,*)inmag_type

select case(trim(inmag_type))

case('random')


!-----------------------set random spin orientation------------------------------------

  do is=1,ns

    call random_number(w1)
	v%x=2.0*w1-1.0
    call random_number(w2)
	v%y=2.0*w2-1.0
    call random_number(w3)
    v%z=2.0*w3-1.0 

    sp(is)%m=norm_vector(v)

  end do

case('vector')

!-----------------------set spin orientation along Vector (vx,vy,vz)-----------------

  do is=1,ns

    sp(is)%m%x=m0%x
    sp(is)%m%y=m0%y   
    sp(is)%m%z=m0%z

  end do


case('head_to_head')

  direction=1.0

  do is=1,ns

    sp(is)%m%x=direction
    sp(is)%m%y=0.0d0   
    sp(is)%m%z=0.0d0     

    direction=-direction

  end do


case('vortex_cw')

  x0=sizex/2.0
  y0=sizey/2.0

  do is=1,ns 
      
    x1=sp(is)%r%x-x0
    y1=sp(is)%r%y-y0         
    r=dsqrt(x1*x1+y1*y1)
    sp(is)%m%x=y1/r
    sp(is)%m%y=-x1/r
    sp(is)%m%z=0.0
  end do

  print*,' initial magnetization - VORTEX_CW'
  write(8,*)' initial magnetization - VORTEX_CW' 

case('vortex_ccw')

  x0=sizex/2.0
  y0=sizey/2.0

  do is=1,ns 
      
    x1=sp(is)%r%x-x0
    y1=sp(is)%r%y-y0         
    r=dsqrt(x1*x1+y1*y1)
    sp(is)%m%x=-y1/r
    sp(is)%m%y=x1/r
    sp(is)%m%z=0.0
  end do

  print*,' initial magnetization - VORTEX_CCW'
  write(8,*)' initial magnetization - VORTEX_CCW' 
case default

  print*,' initial magnetization - nothing selected'
  write(8,*)' initial magnetization - nothing selected'


end select


return

end


