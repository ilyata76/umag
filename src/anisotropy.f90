!-----------compute Anisotropy contribution-----------

!
!------------------------------------------------------------

subroutine anisotropy(ea)
use magmod
implicit real*8 (a-h,o-z)

type (vector) ha

real*8 mi_ui

real*8 k1,ms

ea=0.0

do i=1,ns

  k1=sp(i)%cp%k1
  ms=sp(i)%cp%ms

  mi_ui=sp(i)%m%x*sp(i)%cp%u%x+sp(i)%m%y*sp(i)%cp%u%y+sp(i)%m%z*sp(i)%cp%u%z

  factor=2.0*k1*mi_ui/(mu0*ms)

  ha%x=factor*sp(i)%cp%u%x
  ha%y=factor*sp(i)%cp%u%y
  ha%z=factor*sp(i)%cp%u%z

  sp(i)%h%x=sp(i)%h%x+ha%x
  sp(i)%h%y=sp(i)%h%y+ha%y
  sp(i)%h%z=sp(i)%h%z+ha%z
  
  term=1.0d0-mi_ui**2

  ea=ea + k1*term 

end do

ea=ea*cvol

return
end