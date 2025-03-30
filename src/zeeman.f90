!-------------calculate ZEEMAN energy------------------------
!
!------------------------------------------------------------

subroutine zeeman(ez)

use magmod

implicit real*8 (a-h,o-z)


type (vector) vector_dot_scalar

type (vector) m , h_zeeman, vector_plus_vector, h_exch_bias

ez=0.0d0

do is=1,ns

  factor=sp(is)%cp%ms*sp(is)%cp%vol 

  m= vector_dot_scalar(sp(is)%m,factor)
  
!---------determine direction of external field

  if( extfield_geometry .eq. 'circular' )then
     rx=sp(is)%r%x 
     ry=sp(is)%r%y      
     r=dsqrt(rx**2+ry**2)
     if(r .ne. 0)then
       b%x=bmod*(ry/r)
       b%y=bmod*(-rx/r)
     else
       b%x=0.d0
       b%y=0.d0
     end if

      b%z=0.0 
!     write(8,*)'!!!!',extfield_type
!     write(8,*)' is=',is,' Mx,y,z=',sp(is)%m%x,sp(is)%m%y,sp(is)%m%z
!     write(8,*)'=======Bx,y,z=',b%x,b%y,b%z
  endif   

  ez=ez - vector_dot_vector(b,m)   !   -B*M

  h_zeeman=vector_dot_scalar(b,1/mu0)  !  H=B/mu0

  sp(is)%h = vector_plus_vector(sp(is)%h, h_zeeman)

!--------------add exchange bias----------------------
  h_exch_bias=vector_dot_scalar(sp(is)%cp%eb,1/mu0)

  sp(is)%h = vector_plus_vector(sp(is)%h, h_exch_bias)

end do

return
end