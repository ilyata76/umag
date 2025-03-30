!---------MAKE VECTOR EQUAL TO ANOTHER--------------------

function vector_equal(v)

use magmod

type (vector) v,vector_equal 

vector_equal%x=v%x
vector_equal%y=v%y
vector_equal%z=v%z

return
end