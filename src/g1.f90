!--------------G1(X,Y,Z)----------for demag tensor--------------
!
!----------------------------------------------------------------

real*8 function g1(x,y,z,dx,dy,dz)

implicit real*8 (a-h,o-z)

g1=g2(x+dx,y,z+dz)-g2(x+dx,y,z)-g2(x,y,z+dz)+g2(x,y,z)

return
end