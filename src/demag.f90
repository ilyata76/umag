!-------------calculate MAGNETOSTATIC energy-of DEMAGNETIZATION field-----------------------
!
!------------------------------------------------------------

subroutine demag(ed)

use magmod

implicit real*8 (a-h,o-z)

type (vector) mi,mj

real*8 ms_i,ms_j

ed=0.0


do i=1,ns
  do j=1,ns

    ms_i=sp(i)%cp%ms
    
    mi%x=sp(i)%m%x*ms_i
    mi%y=sp(i)%m%y*ms_i
    mi%z=sp(i)%m%z*ms_i

    ms_j=sp(j)%cp%ms

	mj%x=sp(j)%m%x*ms_j
	mj%y=sp(j)%m%y*ms_j
	mj%z=sp(j)%m%z*ms_j

!------------demag field acting on I dipole

	sp(i)%h%x=sp(i)%h%x-nn(i,j)%xx*mj%x-nn(i,j)%xy*mj%y-nn(i,j)%xz*mj%z
	sp(i)%h%y=sp(i)%h%y-nn(i,j)%yy*mj%y-nn(i,j)%yx*mj%x-nn(i,j)%yz*mj%z
	sp(i)%h%z=sp(i)%h%z-nn(i,j)%zz*mj%z-nn(i,j)%zx*mj%x-nn(i,j)%zy*mj%y

!---------contribution to demag energy -------------

    ed=ed + nn(i,j)%xx*mi%x*mj%x+nn(i,j)%xy*mi%x*mj%y+nn(i,j)%xz*mi%x*mj%z+  &
            nn(i,j)%yy*mi%y*mj%y+nn(i,j)%yx*mi%y*mj%x+nn(i,j)%yz*mi%y*mj%z+  &
            nn(i,j)%zz*mi%z*mj%z+nn(i,j)%zx*mi%z*mj%x+nn(i,j)%zy*mi%z*mj%y

  end do
end do


ed=mu0*ed*cvol/2.0d0

return
end