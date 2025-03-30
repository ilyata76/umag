!---------VECTOR CROSS VECTOR (vector product)

function vector_cross_vector(v1,v2)

use magmod

implicit real*8 (a-h,o-z)

type (vector) v1,v2,vector_cross_vector

vector_cross_vector%x=v1%y*v2%z-v1%z*v2%y
vector_cross_vector%y=v1%z*v2%x-v1%x*v2%z
vector_cross_vector%z=v1%x*v2%y-v1%y*v2%x

return
end