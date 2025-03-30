!---------VECTOR DOT SCALAR---------------------

function vector_dot_scalar(v,scalar)

use magmod

implicit real*8 (a-h,o-z)


type (vector) v,vector_dot_scalar

vector_dot_scalar%x=v%x*scalar
vector_dot_scalar%y=v%y*scalar
vector_dot_scalar%z=v%z*scalar

return
end