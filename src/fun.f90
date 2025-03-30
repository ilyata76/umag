!--------------f(x,y,z)----------for demag tensor--------------
!   FUN
!----------------------------------------------------------------

real*8 function fun(x,y,z)
implicit real*8 (a-h,o-z)


r=dsqrt(x**2+y**2+z**2)

znam1=dsqrt(x**2+z**2)
if(znam1.ne.0)then
  term1=(y/2.0)*(z**2-x**2)*ash(y/znam1) 
else
  term1=0.d0
end if

znam2=dsqrt(x**2+y**2)

if(znam2.ne.0.0)then
  term2=(z/2.0d0)*(y**2-x**2)*ash(z/znam2) 
else
  term2=0.0
end if

znam3=x*r

if(znam3.ne.0.0)then
  term3=x*y*z*atan(y*z/znam3)
else
  term3=0.0
end if

term4=(1.0d0/6.0d0)*(2.0d0*x**2-y**2-z**2)*r 

fun=   term1  +  term2  - term3  +  term4

return
end