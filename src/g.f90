!--------------G(x,y,z)---for demag tensor calculation-----------------------
!
!----------------------------------------------------------------------------

real*8 function g(x,y,z,dx,dy,dz)


implicit real*8 (a-h,o-z)


term1=g1(x,y,z,dx,dy,dz)  

term2=g1(x,y-dy,z,dx,dy,dz)

term3=g1(x,y,z-dz,dx,dy,dz) 

term4=g1(x,y-dy,z-dz,dx,dy,dz)

g=  term1  -  term2  -  term3   +  term4


return
end
