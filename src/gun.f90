!--------------g(x,y,z)----------for demag tensor--------------
!   GUN
!----------------------------------------------------------------
real*8 function gun(x,y,z)

implicit real*8 (a-h,o-z)

r=dsqrt(x**2+y**2+z**2)

!--------------------------------TERM1----------------------------
znam1=dsqrt(x**2+y**2)

if(znam1.ne.0.d0)then
  term1=(x*y*z)*ash(z/znam1) 
else
  term1=0.d0
end if

!--------------------------------TERM2------------------------------------
znam2=dsqrt(y**2+z**2)

if(znam2.ne.0.d0)then
  term2=(y/6.0d0)*(3.0d0*z**2-y**2)*ash(x/znam2)
else
  term2=0.0d0
end if

!-------------------------------TERM3-------------------------------------
znam3=dsqrt(x**2+z**2)

if(znam3.ne.0.d0)then
  term3=(x/6.0d0)*(3.0d0*z**2-x**2)*ash(y/znam3)  
else
  term3=0.0d0
end if

!-------------------------------TERM4------------------------------------
znam4a=(z*r)

if(znam4a.ne.0.d0)then
  term4a=(z**3/6.0)*atan(x*y/znam4a)
else
  term4a=0.0d0
end if


znam4b=(y*r)

if(znam4b.ne.0.d0)then
  term4b=(z*y**2/2.0d0)*atan(x*z/znam4b)
else
  term4b=0.d0
end if

term4=term4a+term4b

!-------------------------------TERM5------------------------------------
znam5=(x*r)

if(znam5.ne.0.d0)then
  term5a=(z*x**2)*atan(y*z/znam5)
else
  term5a=0.d0
end if

term5b=x*y*r/3.0d0

term5=term5a+term5b

!---------------------------------GUN----------------------------------------
gun=    term1  + term2 + term3 - term4 - term5

return
end