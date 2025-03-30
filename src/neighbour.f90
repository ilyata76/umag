!-------------Find 6 neighbours for each cell------------------------
!
!------------------------------------------------------------

subroutine neighbour()

use magmod
implicit real*8 (a-h,o-z)


do i=1,ns-1
  do j=i+1,ns

    xxij=sp(i)%r%x-sp(j)%r%x
    yyij=sp(i)%r%y-sp(j)%r%y
    zzij=sp(i)%r%z-sp(j)%r%z

    rrij=dsqrt(xxij**2+yyij**2+zzij**2)

!---------------determine 6 neighbors for each cell------------------------
    if(rrij.le.1.2*hx)then

      if(xxij.gt.0 .and. xxij .le. 1.1*hx)then
	    sp(i)%cp%xm=j
	    sp(j)%cp%xp=i	  
      end if

      if(xxij.lt.0 .and. xxij .ge. -1.1*hx)then
	    sp(i)%cp%xp=j
	    sp(j)%cp%xm=i	  
      end if

      if(yyij.gt.0 .and. yyij .le. 1.1*hy)then
	    sp(i)%cp%ym=j
	    sp(j)%cp%yp=i	  
      end if

      if(yyij.lt.0 .and. yyij .ge. -1.1*hy)then
	    sp(i)%cp%yp=j
	    sp(j)%cp%ym=i	  
      end if

      if(zzij.gt.0 .and. zzij .le. 1.1*hz)then
	    sp(i)%cp%zm=j
	    sp(j)%cp%zp=i	  
      end if

      if(zzij.lt.0 .and. zzij .ge. -1.1*hz)then
	    sp(i)%cp%zp=j
	    sp(j)%cp%zm=i	  
      end if	

    end if

  end do
end do


return
end
