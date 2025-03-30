!---------creates mesh for calculations

subroutine create_mesh()

use magmod

implicit real*8 (a-h,o-z)

!----------------------------------------------------------------------

sizex=mesh_zone%xmax-mesh_zone%xmin
mesh_nx=sizex/hx+1

sizey=mesh_zone%ymax-mesh_zone%ymin
mesh_ny=sizey/hy+1

sizez=mesh_zone%zmax-mesh_zone%zmin
mesh_nz=sizez/hz+1

ns=0 ! initialization of spin counter

write(8,*)' sizeX,Y,Z=',sizex,sizey,sizez

write(8,*)' MESH nx,ny,nz=',mesh_nx,mesh_ny,mesh_nz

!-------------------generate rectangular mesh--------------------------------------

do iz=1,mesh_nz
   do iy=1,mesh_ny
      do ix=1,mesh_nx

        xsite=mesh_zone%xmin+hx/2.0+hx*(ix-1)
		ysite=mesh_zone%ymin+hy/2.0+hy*(iy-1)
		zsite=mesh_zone%zmin+hz/2.0+hz*(iz-1) 

        do ireg=1,nreg

          if(xsite .gt. reg(ireg)%xmin  .and.  xsite .lt.  reg(ireg)%xmax   .and. &
		     ysite .gt. reg(ireg)%ymin  .and.  ysite .lt.  reg(ireg)%ymax   .and. &
             zsite .gt. reg(ireg)%zmin  .and.  zsite .lt.  reg(ireg)%zmax   )then

             ns=ns+1

             sp(ns)%r%x=xsite
             sp(ns)%r%y=ysite
             sp(ns)%r%z=zsite  

             sp(ns)%cp%a=reg(ireg)%a 
             sp(ns)%cp%ms=reg(ireg)%ms
             sp(ns)%cp%vol=hx*hy*hz
		     sp(ns)%cp%k1=reg(ireg)%k1

     !--------------easy direction for anisotropy-----------   
		     sp(ns)%cp%u%x=reg(ireg)%u%x
		     sp(ns)%cp%u%y=reg(ireg)%u%y
		     sp(ns)%cp%u%z=reg(ireg)%u%z

     !--------------exchange bias field --------------------   
		     sp(ns)%cp%eb%x=reg(ireg)%eb%x
		     sp(ns)%cp%eb%y=reg(ireg)%eb%y
		     sp(ns)%cp%eb%z=reg(ireg)%eb%z

          end if   

        end do
    
      end do
   end do
end do


write(8,*)' ns=',ns

return
end