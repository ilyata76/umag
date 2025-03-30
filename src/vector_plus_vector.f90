!-----------VECTOR SUM (V1+V2)----------------------
!
!----------------------------------------------
function vector_plus_vector(v1,v2)

use magmod

implicit real*8 (a-h,o-z)

type (vector) v1,v2,vector_plus_vector

vector_plus_vector%x=v1%x+v2%x
vector_plus_vector%y=v1%y+v2%y
vector_plus_vector%z=v1%z+v2%z

return
end