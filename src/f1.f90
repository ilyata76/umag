!--------------F1(X,Y,Z)----------for demag tensor--------------
!
!----------------------------------------------------------------

real*8 function f1(x,y,z,dx,dy,dz)

implicit real*8 (a-h,o-z)

f1=f2(x,y,z)-f2(x,y-dy,z)-f2(x,y,z-dz)+f2(x,y-dy,z-dz)

return
end