!--------------Nxx(x,y,z)---for demag tensor calculation-----------------------
!
!----------------------------------------------------------------------------

real*8 function nxx(x,y,z,dx,dy,dz)

implicit real*8 (a-h,o-z)

term1=2.0d0*f(x,y,z,dx,dy,dz) 

term2=      f(x+dx,y,z,dx,dy,dz)   

term3=      f(x-dx,y,z,dx,dy,dz)

nxx= term1  - term2  - term3                   

return

end