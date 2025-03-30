!--------------F(x,y,z)---for demag tensor calculation-----------------------
!
!----------------------------------------------------------------------------

real*8 function f(x,y,z,dx,dy,dz)

implicit real*8 (a-h,o-z)

term1=f1(x  ,y+dy,z+dz,dx,dy,dz)

term2=f1(x  ,y   ,z+dz,dx,dy,dz) 

term3=f1(x  ,y+dy,z   ,dx,dy,dz) 

term4=f1(x  ,y   ,z   ,dx,dy,dz)

f=  term1 - term2 - term3 +  term4

return
end
