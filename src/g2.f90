!--------------G2(X,Y,Z)----------for demag tensor--------------
!
!----------------------------------------------------------------

real*8 function g2(x,y,z)

implicit real*8 (a-h,o-z)


if(abs(x).le.1d-20)x=0.d0
if(abs(y).le.1d-20)y=0.d0
if(abs(z).le.1d-20)z=0.d0

g2=gun(x,y,z)-gun(x,y,0.d0)


return
end
