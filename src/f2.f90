!--------------F2(X,Y,Z)----------for demag tensor--------------
!
!----------------------------------------------------------------

real*8 function f2(x,y,z)

implicit real*8 (a-h,o-z)

if(abs(x).le.1d-20)x=0.d0
if(abs(y).le.1d-20)y=0.d0
if(abs(z).le.1d-20)z=0.d0

f2=fun(x,y,z)-fun(x,0.d0,z)-fun(x,y,0.d0)+fun(x,0.d0,0.d0)

return
end
