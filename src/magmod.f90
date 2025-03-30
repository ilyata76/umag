!--------------MAGMOD---------------------------------
!
!   common module for MICROMAG project
!---------------------------------------------------------
module magmod
implicit none

!---------------basic physical constants--------------------
real*8, parameter :: pi=3.1415926d0
real*8  mu0  !  magnetic constant 
real*8 gamma ! gyromagnetic ratio
real*8 alfa !  damping parameter

!--------------------vector----------------------------------------------------------------------
type vector
  real*8 x,y,z
end type vector

!-------------------Regions--------------------------------------
type region
  real*8 xmin,xmax             ! region borders along X
  real*8 ymin,ymax             ! region borders along Y 
  real*8 zmin,zmax             ! region borders along Z
  character*10 mt              ! type of material
  real*8 ms                    ! magnetization saturation
  real*8 a                     ! exchange parameter
  real*8 k1                    ! anisotropy parameter K1
  type (vector) u              ! easy direction
  type (vector) eb             ! exchange bias field (T)
end type region

!------------------ COMPLETE SYSTEM ------------------------
type (region) mesh_zone

real*8 sizex,sizey,sizez

!----------------subdivision--------------------------------------
integer, parameter :: nregion=10

integer nreg
type (region) reg(nregion)

!----------------mesh--------------------------------------------
real*8 hx,hy,hz,cvol

!----------------PROPERTY OF MESH CELL----------------------------------------------------
type cell_property

  real*8 a  !  exchange constant per unite volume  
  real*8 ms ! saturation magnetization (a/m) per unite volume
  real*8 k1 !  uniaxial anisotropy constant K1 per unite volume
  real*8 vol ! cell volume
  type (vector) eb !  exchange bias field (B,T)

!--------6 neighbors around for exchange interaction-----------------

  integer xm,xp  !  neighbors along X axis (xm- x minus, xp - x plus)
  integer ym,yp  !  neighbors along Y axis (ym- y minus, yp - y plus)
  integer zm,zp  !  neighbors along Z axis (zm- z minus, zp - z plus)

!----------------easy direction-----------------------
  type (vector) u  

end type cell_property


!------------------magnitization vector at mesh centers-----------

type spin

  type (vector) r !  radius vector of mesh cell 

  type (vector) m !  orientation of the magnetization vector

  type (vector) h !  effective field on cell

  type (cell_property) cp

end type spin

!---------------------array of spins------------------------------

integer, parameter :: nsp=5100  ! size of array allocated for spins

type (spin) sp(nsp)

integer ns


!---------------Demag tensor-(Nij)--------------------------------

type demagtensor

  real*8 xx,xy,xz,yx,yy,yz,zx,zy,zz

end type demagtensor


type (demagtensor) nn(nsp,nsp)

!--------------external field----B-----------------------------------------------------

character(len=15) extfield_geometry

type (vector) b

real*8 bmod

!-------------IMPULSE B-------------------------------------------
integer nt_imp
real*8 ht_imp
real*8 bmt(300)


!-------------------parameters of numerical scheme------------------

integer nstop, iout

real*8 dt   ! (sec) time step for LLG evolve

real*8 max_torque

real*8 max_angle

end module magmod
