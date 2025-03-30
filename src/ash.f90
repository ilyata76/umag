!--------------arcsh(x)----------for demag tensor--------------
!   ASH
!----------------------------------------------------------------

real*8 function ash(x)

implicit real*8 (a-h,o-z)

arg=x+dsqrt(1.d0+x**2)

ash=dlog(arg)


return
end