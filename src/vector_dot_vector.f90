!---------------VECTOR DOT VECTOR  (scalar product)

 function vector_dot_vector(v1,v2)

use magmod

implicit real*8 (a-h,o-z)

type (vector) v1,v2

vector_dot_vector=v1%x*v2%x+v1%y*v2%y+v1%z*v2%z

return
end 