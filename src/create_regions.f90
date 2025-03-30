!---------creates magnetic regions ------------------------

subroutine create_regions()
use magmod
implicit real*8 (a-h,o-z)

character*60 title_regions
character*150 list_parameter

!----------------------------------------------------------------------

open(15,file='regions.dat')

read(15,'(a60)')title_regions
write(8,'(a60)')title_regions

read(15,'(a150)')list_parameter
write(8,'(a150)')list_parameter

xmax=-100000.0
xmin=100000.0

ymax=-100000.0
ymin=100000.0

zmax=-100000.0
zmin=100000.0

do ireg=1,nregion   ! maximal number of regions

  read(15,*,end=1)reg(ireg)%xmin,reg(ireg)%xmax,                               &
                  reg(ireg)%ymin,reg(ireg)%ymax,                               &
                  reg(ireg)%zmin,reg(ireg)%zmax,                               &
                  reg(ireg)%mt,reg(ireg)%ms,reg(ireg)%a,                       &
                  reg(ireg)%k1,reg(ireg)%u%x,reg(ireg)%u%y,reg(ireg)%u%z,      &
                  reg(ireg)%eb%x,reg(ireg)%eb%y,reg(ireg)%eb%z
!---------------find MAX and MIN-----------------------------
  if(reg(ireg)%xmin .lt. xmin) xmin=reg(ireg)%xmin
  if(reg(ireg)%xmax .gt. xmax) xmax=reg(ireg)%xmax

  if(reg(ireg)%ymin .lt. ymin) ymin=reg(ireg)%ymin
  if(reg(ireg)%ymax .gt. ymax) ymax=reg(ireg)%ymax

  if(reg(ireg)%zmin .lt. zmin) zmin=reg(ireg)%zmin
  if(reg(ireg)%zmax .gt. zmax) zmax=reg(ireg)%zmax

!--------------output borders---------------------------------

  write(8,'(2x,6(g10.3,1x),a10,1x,3(g10.3,1x),1x,3(f3.1,1x),3(f9.3,1x))')    &
            reg(ireg)%xmin,reg(ireg)%xmax,                                   &
            reg(ireg)%ymin,reg(ireg)%ymax,                                   &
            reg(ireg)%zmin,reg(ireg)%zmax,                                   &
            reg(ireg)%mt,reg(ireg)%ms,reg(ireg)%a,reg(ireg)%k1,              &
			reg(ireg)%u%x,reg(ireg)%u%y,reg(ireg)%u%z,                       &
            reg(ireg)%eb%x,reg(ireg)%eb%y,reg(ireg)%eb%z

end do

1 nreg=ireg-1

mesh_zone%xmin=xmin
mesh_zone%xmax=xmax

mesh_zone%ymin=ymin
mesh_zone%ymax=ymax

mesh_zone%zmin=zmin
mesh_zone%zmax=zmax

write(8,*)' MESH_ZONE ',xmin,xmax,ymin,ymax,zmin,zmax


write(8,*)'====================================='

close(15)

return
end