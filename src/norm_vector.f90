!-----------NORMALIZE VECTOR----------------------
!
!----------------------------------------------
function norm_vector(v)

use magmod

implicit real*8 (a-h,o-z)

type (vector) v,norm_vector

vv=dsqrt(v%x**2 + v%y**2 +v%z**2)

if(vv.ne.0)then
  norm_vector%x=v%x/vv
  norm_vector%y=v%y/vv
  norm_vector%z=v%z/vv
else
  norm_vector%x=0.0
  norm_vector%y=0.0
  norm_vector%z=0.0
end if

return 
end