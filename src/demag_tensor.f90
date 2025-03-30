!-------------compute DEMAG tensor  (Nij)-----------------------
!
!------------------------------------------------------------

subroutine demag_tensor()

use magmod
implicit real*8 (a-h,o-z)

real*8 nxx,nxy

print*,'start demag tensor' 
write(8,*)'start demag tensor  !!!!' 

znam= 4.0d0*pi*cvol

!----------assign demag tensor to each I cell-----------------------

idone=0

print*,' total number of cells ns=',ns
ns_ns=ns*ns
print*,' total number of tensor elements  ns_ns=',ns_ns

iprint=ns_ns/100

!---------------------------------------------------------------------------

do i=1,ns

  xx=0.0d0
  yy=0.0d0
  zz=0.0d0

!------------------tensor Nij-------------------------------------------
!-----------------Nxx------------------------------------------------
    nn(i,i)%xx=nxx(xx,yy,zz,hx,hy,hz)/znam  
    
!-------------------Nxy-------------------------------------------------
    nn(i,i)%xy=nxy(xx,yy,zz,hx,hy,hz)/znam

!-------------------Nxz-------------------------------------------------
    nn(i,i)%xz=nxy(xx,zz,yy,hx,hz,hy)/znam

!-----------------Nyy------------------------------------------------
    nn(i,i)%yy=nxx(yy,xx,zz,hy,hx,hz)/znam  

!-------------------Nyx-------------------------------------------------
    nn(i,i)%yx=nn(i,i)%xy

!-------------------Nyz-------------------------------------------------
    nn(i,i)%yz=nxy(yy,zz,xx,hy,hz,hx)/znam

!-----------------Nzz------------------------------------------------
    nn(i,i)%zz=nxx(zz,yy,xx,hz,hy,hx)/znam  

!-------------------Nxz-------------------------------------------------
    nn(i,i)%zx=nn(i,i)%xz

!-------------------Nzy-------------------------------------------------
    nn(i,i)%zy=nn(i,i)%yz

    idone=idone+1 
     
    if(mod(idone,iprint).eq.0)print*,' calc diagonal idone=',idone,' tot num of tens elem ns_ns=',ns_ns

end do

!----------assign demag tensor to each I-J pair-----------------------

do i=1,ns-1
  do j=i+1,ns

    xx=sp(i)%r%x-sp(j)%r%x
    yy=sp(i)%r%y-sp(j)%r%y
    zz=sp(i)%r%z-sp(j)%r%z

!------------------tensor Nij-------------------------------------------
!-----------------Nxx------------------------------------------------
    nn(i,j)%xx=nxx(xx,yy,zz,hx,hy,hz)/znam  
    
!-------------------Nxy-------------------------------------------------
    nn(i,j)%xy=nxy(xx,yy,zz,hx,hy,hz)/znam

!-------------------Nxz-------------------------------------------------
    nn(i,j)%xz=nxy(xx,zz,yy,hx,hz,hy)/znam

!-----------------Nyy------------------------------------------------
    nn(i,j)%yy=nxx(yy,xx,zz,hy,hx,hz)/znam  

!-------------------Nyx-------------------------------------------------
    nn(i,j)%yx=nn(i,j)%xy

!-------------------Nyz-------------------------------------------------
    nn(i,j)%yz=nxy(yy,zz,xx,hy,hz,hx)/znam

!-----------------Nzz------------------------------------------------
    nn(i,j)%zz=nxx(zz,yy,xx,hz,hy,hx)/znam  

!-------------------Nxz-------------------------------------------------
    nn(i,j)%zx=nn(i,j)%xz

!-------------------Nzy-------------------------------------------------
    nn(i,j)%zy=nn(i,j)%yz

!------------tensor fot Nji------------------------------------------------------
    nn(j,i)%xx=nn(i,j)%xx
    nn(j,i)%xy=nn(i,j)%xy
    nn(j,i)%xz=nn(i,j)%xz
    nn(j,i)%yy=nn(i,j)%yy
    nn(j,i)%yx=nn(i,j)%yx
    nn(j,i)%yz=nn(i,j)%yz
    nn(j,i)%zz=nn(i,j)%zz
    nn(j,i)%zx=nn(i,j)%zx
    nn(j,i)%zy=nn(i,j)%zy

    idone=idone+1
    
    if(mod(idone,iprint).eq.0)print*,' calc others idone=',idone,' tot num of tens elem ns_ns=',ns_ns

  end do
end do

print*,' demag tensor done'
write(8,*)' demag tensor done'

return
end