!--------------Nxy(x,y,z)---for demag tensor calculation-----------------------
!
!----------------------------------------------------------------------------

real*8 function nxy(x,y,z,dx,dy,dz)

implicit real*8 (a-h,o-z)

term1=g(x,y,z,dx,dy,dz)

term2=g(x-dx,y,z,dx,dy,dz) 

term3=g(x,y+dy,z,dx,dy,dz) 

term4=g(x-dx,y+dy,z,dx,dy,dz)

nxy=  term1   -   term2  -  term3  +  term4

return
end