!-----------MICROMAG--7--------------------
!   modification of MICROMAG5 to work with multi region systems 
!   (each region in separate file)
!   project started on 5th october 2010
!----------------------------------------------
program micromag

use magmod

implicit real*8 (a-h,o-z)


type (vector) def ! direction of external field

type (vector) mag ! total magnetization of sample

type (vector) m0 ! direction of initial magnetization

type (vector) tangent ! tangent to the ring

type (spin) spw(nsp) ! magnetization of grid cell

!----------------initial parameters-----------------------

character(len=60) job_title, config_title, comments  

character(len=80),dimension(40)::text_row

character(len=60) line  

character(len=15),dimension(40) ::  wrd
real(8),dimension(40) :: xinp

real*8 m

real*8 impulse

character(len=15) initial_configuration
character(len=15) material_type,geometry_type,grid_type,extfield_type
character(len=15) inmag_type
character(len=15) numerics_type

character(len=15),dimension(10) :: file_layer

!=====    input of initial parameters==================================

!---------SAMPLE GEOMETRY-----------------------------------

 irow_start_geometry= 2
 ntext_geometry=2
 ndata_geometry=0
 nrows_geometry=ntext_geometry+ndata_geometry

 irow_finish_geometry=irow_start_geometry+nrows_geometry-1

!-----------GRID----------------------------------------
 irow_start_grid= irow_finish_geometry+2
  ntext_grid=1
  ndata_grid=3
  nrows_grid=ntext_grid+ndata_grid

 irow_finish_grid=irow_start_grid+nrows_grid-1

!-----------INMAG----------------------------------------
  irow_start_inmag= irow_finish_grid+2

  ntext_inmag=1
  ndata_inmag=3
  nrows_inmag=ntext_inmag+ndata_inmag

 irow_finish_inmag=irow_start_inmag+nrows_inmag-1


!-----------EXTFIELD----------------------------------------
 irow_start_extfield= irow_finish_inmag+2
  ntext_extfield=2
  ndata_extfield=6
  nrows_extfield=ntext_extfield+ndata_extfield

 irow_finish_extfield=irow_start_extfield+nrows_extfield-1

!-----------NUMERICS----------------------------------------
   irow_start_numerics= irow_finish_extfield+2
  ntext_numerics=1
  ndata_numerics=5
  nrows_numerics=ntext_numerics+ndata_numerics

 irow_finish_numerics=irow_start_numerics+nrows_numerics-1

 nrows_total=irow_finish_numerics

open(5,file='micromag.dat')

read(5,'(a60)')job_title

do irow=1,nrows_total

  read(5,'(a80)')text_row(irow)

!------------------GEOMETRY-----------------2-------------------
  if( irow .ge. irow_start_geometry .and. irow .le. irow_finish_geometry )then
     if(irow.le.irow_start_geometry+ntext_geometry-1)then
        read(text_row(irow),'(a60,a15)')line,wrd(irow)   
     else
        read(text_row(irow),'(a60,g15.5)')line,xinp(irow)   
     end if
  end if

!------------------GRID--------------------3----------------
  if( irow .ge. irow_start_grid .and. irow .le. irow_finish_grid )then
     if(irow .le. irow_start_grid+ntext_grid-1)then
        read(text_row(irow),'(a60,a15)')line,wrd(irow)   
     else
        read(text_row(irow),'(a60,g15.5)')line,xinp(irow)   
     end if
  end if

!------------------INMAG------------------4------------------
  if( irow .ge. irow_start_inmag .and. irow .le. irow_finish_inmag )then
     if(irow .le. irow_start_inmag+ntext_inmag-1)then
        read(text_row(irow),'(a60,a15)')line,wrd(irow)   
     else
        read(text_row(irow),'(a60,g15.5)')line,xinp(irow)   
     end if
  end if

!------------------EXTFIELD-----------------5-------------------
  if( irow .ge. irow_start_extfield  .and.  irow .le. irow_finish_extfield )then
     if(irow .le. irow_start_extfield+ntext_extfield-1)then
        read(text_row(irow),'(a60,a15)')line,wrd(irow)   
     else
        read(text_row(irow),'(a60,g15.5)')line,xinp(irow)   
     end if
  end if

!------------------NUMERICS--------------------6----------------
  if( irow .ge. irow_start_numerics .and. irow .le. irow_finish_numerics )then
     if(irow .le. irow_start_numerics+ntext_numerics-1)then
        read(text_row(irow),'(a60,a15)')line,wrd(irow)   
     else
        read(text_row(irow),'(a60,g15.5)')line,xinp(irow)   
     end if
  end if

end do

close(5)


!------------output parameters

open(8,file="micromag.rez")

write(8,'(a60)')job_title
do i=1,nrows_total
  write(8,'(a80)')text_row(i)
end do


!-------------assign initail parameters---------------------------------
!-----------------------GEOMETRY------------------------------------------
initial_configuration=wrd(irow_start_geometry)
geometry_type=wrd(irow_start_geometry+1)

write(8,'(a15)')geometry_type
!-------------------------GRID---------------------------------------------
grid_type=wrd(irow_start_grid)
hx=xinp(irow_start_grid+ntext_grid)
hy=xinp(irow_start_grid+ntext_grid+1)
hz=xinp(irow_start_grid+ntext_grid+2)

write(8,'(a15)')grid_type
write(8,*)hx,hy,hz
!------------------------INITIAL MAGNETIZATION------------------------------
inmag_type=wrd(irow_start_inmag)
m0%x=xinp(irow_start_inmag+ntext_inmag)
m0%y=xinp(irow_start_inmag+ntext_inmag+1)
m0%z=xinp(irow_start_inmag+ntext_inmag+2)

write(8,'(a15)')inmag_type
write(8,*)m0%x,m0%y,m0%z
!--------------------------EXTERNAL FIELD---------------------------------
extfield_type=wrd(irow_start_extfield)
extfield_geometry=wrd(irow_start_extfield+1)
def%x=xinp(irow_start_extfield+ntext_extfield)
def%y=xinp(irow_start_extfield+ntext_extfield+1)
def%z=xinp(irow_start_extfield+ntext_extfield+2)
bmin=xinp(irow_start_extfield+ntext_extfield+3)
hb=xinp(irow_start_extfield+ntext_extfield+4)
bmax=xinp(irow_start_extfield+ntext_extfield+5)

write(8,'(a15)')extfield_type
write(8,*)def%x,def%y,def%z
write(8,*)bmin,hb,bmax
!-----------------------PARAMETERS FOR NUMERICAL SCHEME----------------
numerics_type=wrd(irow_start_numerics)
nstop=xinp(irow_start_numerics+ntext_numerics)
iout=xinp(irow_start_numerics+ntext_numerics+1)
dt=xinp(irow_start_numerics+ntext_numerics+2)
alfa=xinp(irow_start_numerics+ntext_numerics+3)
torque_stop=xinp(irow_start_numerics+ntext_numerics+4)

write(8,'(a15)')numerics_type
write(8,*)nstop,iout,alfa,dt


!------------------------------------------------------------------------

mu0=4*pi*1d-7   !   magnetic constant

e=1.6e-19       !  charge of the electron
 
m=9e-31         !  mass of the electron 

gamma=mu0*(e/m)

cvol=hx*hy*hz

write(8,*)' gamma=',gamma

!gamma=2.21e5 !  gyromagnetic ratio from master thesis gamma=(e/m)*mu0 



!--------------------INITIAL CONFIGURATION---------------------------

select case(initial_configuration)

case('generator.3d')

!------------------create REGIONS -(subdivisions)------------------

  call create_regions()

!-----------------create Mesh----------------------------------------

  call create_mesh()

!---------------create initial magnetization-------------------------
  call initial_magnetization(inmag_type,m0)


!-----------------------------------------------MULTI LAYER-------------------
case('multi_layer')

  open(7,file='multi_layer.dat')
  read(7,'(a60)')line
  write(8,'(a60)')line
  do l=1,10
    read(7,'(a15)',end=2)file_layer(l)
  end do
2 number_of_layers=l-1
  write(8,*)' number of layers=',number_of_layers

  xmax=-100000
  xmin=1000000
  ymax=-100000
  ymin=1000000
  zmax=-100000
  zmin=1000000

  ns=0
  is=0

  do l=1,number_of_layers
    open(9,file=file_layer(l))
    !--------read MATERIAL PARAMETERS FROM file (REGION - l)
    write(8,*)'layer =',l
    read(9,'(a60)')config_title
    write(8,'(a60)')config_title
    read(9,'(a60)')comments 
    write(8,'(a60)')comments
    read(9,*)reg(l)%a ,reg(l)%ms,reg(l)%k1,reg(l)%u%x,reg(l)%u%y,reg(l)%u%z
    write(8,'(2x,3(g12.3,1x),3(f9.3,1x))')reg(l)%a ,reg(l)%ms,reg(l)%k1,reg(l)%u%x,reg(l)%u%y,reg(l)%u%z  

  
    do is_layer=1,10000
      
      read(9,*,end=3)sprx,spry,sprz,spmx,spmy,spmz

      is=is+1
      sp(is)%r%x=sprx
      sp(is)%r%y=spry
      sp(is)%r%z=sprz
      sp(is)%m%x=spmx
      sp(is)%m%y=spmy
      sp(is)%m%z=spmz
      
!-----------calculate size of box------------------------------

      if(sp(is)%r%x .gt. xmax)xmax=sp(is)%r%x
      if(sp(is)%r%x .lt. xmin)xmin=sp(is)%r%x   

      if(sp(is)%r%y .gt. ymax)ymax=sp(is)%r%y
      if(sp(is)%r%y .lt. ymin)ymin=sp(is)%r%y   
   
      if(sp(is)%r%z .gt. zmax)zmax=sp(is)%r%z
      if(sp(is)%r%z .lt. zmin)zmin=sp(is)%r%z   
   
   
!-------------assign material parameters------------------------
  
      sp(is)%cp%a=reg(l)%a 
      sp(is)%cp%ms=reg(l)%ms
      sp(is)%cp%vol=hx*hy*hz
      sp(is)%cp%k1=reg(l)%k1
  !--------------easy direction --------------------   
      sp(is)%cp%u%x=reg(l)%u%x
      sp(is)%cp%u%y=reg(l)%u%y
      sp(is)%cp%u%z=reg(l)%u%z

    end do

3 ns_layer=is_layer-1

    write(8,*)' input layer done ns_layer=',ns_layer 
    close(9)

    ns=ns+ns_layer

  end do
   
  sizex=xmax-xmin
  sizey=ymax-ymin
  sizez=zmax-zmin
  if(sizez.le.0)sizez=hz 
  
  write(8,*)' xmin,xmax=',xmin,xmax
  write(8,*)' ymin,ymax=',ymin,ymax  
  write(8,*)' zmin,zmax=',zmin,zmax  
  write(8,*)' sizex,sizey,sizez=',sizex,sizey,sizez 

  write(8,*)' initial config input done ns=',ns 

case default

  open(7,file=initial_configuration)

!--------read MATERIAL PARAMETERS FROM file (REGION 1)

  read(7,'(a60)')config_title
  write(8,'(a60)')config_title
  read(7,'(a60)')comments 
  write(8,'(a60)')comments
  read(7,*)reg(1)%a ,reg(1)%ms,reg(1)%k1,reg(1)%u%x,reg(1)%u%y,reg(1)%u%z
  write(8,'(2x,3(g12.3,1x),3(f9.3,1x))')reg(1)%a ,reg(1)%ms,reg(1)%k1,reg(1)%u%x,reg(1)%u%y,reg(1)%u%z  

  xmax=-100000
  xmin=1000000
  ymax=-100000
  ymin=1000000
  zmax=-100000
  zmin=1000000
  
  
  do is=1,10000

    read(7,*,end=1)sp(is)%r%x,sp(is)%r%y,sp(is)%r%z,sp(is)%m%x,sp(is)%m%y,sp(is)%m%z

!-----------calculate size of box------------------------------

   if(sp(is)%r%x .gt. xmax)xmax=sp(is)%r%x
   if(sp(is)%r%x .lt. xmin)xmin=sp(is)%r%x   

   if(sp(is)%r%y .gt. ymax)ymax=sp(is)%r%y
   if(sp(is)%r%y .lt. ymin)ymin=sp(is)%r%y   
   
   if(sp(is)%r%z .gt. zmax)zmax=sp(is)%r%z
   if(sp(is)%r%z .lt. zmin)zmin=sp(is)%r%z   
   
   
!-------------assign material parameters------------------------
    
    
    sp(is)%cp%a=reg(1)%a 
    sp(is)%cp%ms=reg(1)%ms
    sp(is)%cp%vol=hx*hy*hz
    sp(is)%cp%k1=reg(1)%k1
  !--------------easy direction --------------------   
    sp(is)%cp%u%x=reg(1)%u%x
	sp(is)%cp%u%y=reg(1)%u%y
	sp(is)%cp%u%z=reg(1)%u%z

  end do

1 ns=is-1

  write(8,*)' initial config input done ns=',ns 

  close(7)

  sizex=xmax-xmin
  sizey=ymax-ymin
  sizez=zmax-zmin
  if(sizez.le.0.0)sizez=hz 

  
  write(8,*)' xmin,xmax=',xmin,xmax
  write(8,*)' ymin,ymax=',ymin,ymax  
  write(8,*)' zmin,zmax=',zmin,zmax  
  write(8,*)' sizex,sizey,sizez=',sizex,sizey,sizez 

end select

!---------volume of the sample-----------------------------------------------

vol=sizex*sizey*sizez

!----------output initial magnetization vector field--------------------------

call output_magfield("magfield0.map  ")

!---------------------compute Demag tensor (Nij)---------------------

call demag_tensor()

!open(18,file='demag_tensor.dat')
!write(18,*)' -----------------DEMAG TENSOR----------------------------'
!do i=1,ns
!  do j=1,ns
!      write(18,*)' =============='
!      write(18,*)' i=',i,' j=',j
!      write(18,*)' NN xx,xy,xz=', nn(i,j)%xx,nn(i,j)%xy,nn(i,j)%xz
!      write(18,*)' NN yy,yx,yz=', nn(i,j)%yy,nn(i,j)%yx,nn(i,j)%yz
!      write(18,*)' NN zz,zx,zy=', nn(i,j)%zz,nn(i,j)%zx,nn(i,j)%zy
!      write(18,*)' j=',j,' i=',i
!      write(18,*)' NN xx,xy,xz=', nn(j,i)%xx,nn(j,i)%xy,nn(j,i)%xz
!      write(18,*)' NN yy,yx,yz=', nn(j,i)%yy,nn(j,i)%yx,nn(j,i)%yz
!      write(18,*)' NN zz,zx,zy=', nn(j,i)%zz,nn(j,i)%zx,nn(j,i)%zy
!  end do
!end do
!close(18)


call neighbour()

write(8,*)'--------------NEIGHBORS------------------'
!do i=1,ns
!  write(8,*)' i=',i,' jpx=',sp(i)%cp%xp, ' jmx=',sp(i)%cp%xm, &
!                    ' jpy=',sp(i)%cp%yp, ' jmy=',sp(i)%cp%ym, &
!		     ' jpz=',sp(i)%cp%zp, ' jmz=',sp(i)%cp%zm 
!end do

!--------------------------------open file for comparing with OOMMF---------------

if(extfield_type.eq.'fixed')open(16,file='compare.rez')

!--------------LOOP ON EXTERNAL FIELD-------------------------------

if(extfield_type.eq.'loop')open(12,file='loop.rez')

!------------open movie file-------------------------

open(11,file='magmovie.m')

write(11,'(a60)')'MAGMOVIE 3D test'
write(11,'(2x,3(g12.5,1x))')sizex,sizey,sizez

!-------output 3d movie-snapshot--istep=0---------------------------
istep=0
bmod=0.d0
if(extfield_type.eq.'fixed')then
  write(11,'(a1)')'c'
  write(11,'(2x,i6,2x,f10.3,2x,3(g12.5,1x))')istep,bmod,ez,ed,ee,ea,et
  do is=1,ns
    write(11,'(2x,3(g12.5,1x),3(f12.5,1x))')sp(is)%r%x,sp(is)%r%y,sp(is)%r%z,  &
                                            sp(is)%m%x,sp(is)%m%y,sp(is)%m%z
  end do
end if

!-------------------compute number of steps on external field------------------

select case(extfield_type)

case('loop')
  nb=2*(bmax-bmin)/hb
  dhb=hb  ! step for B
  bmod=bmin-dhb
case('impulse')
  open(19,file='impulse.dat')
  read(19,*)ht_imp
  do it=1,1000
    read(19,*,end=4)bmt(it)
  end do
4 nt_imp=it-1
  write(8,*)' IMPULSE data loaded nt_imp=',nt_imp
  nb=1
  open(20,file='impulse.rez')
case default !  fixed 
  nb=1
  dhb=hb  ! step for B
  bmod=bmin-dhb
end select


write(8,'(a15)')extfield_type
write(8,*)' nb=',nb

do ib=1,nb

  if(extfield_type.ne.'impulse')then
    bmod=bmod+dhb
    print*,' external field  B= ', bmod
    write(8,*)'!!!! external field  B= ', bmod
    if(bmod.ge.bmax)dhb=-hb
!----------SET LINEAR EXTERNAL FIELD-----------------------------
    if(extfield_type .ne. 'circular')then ! linear field 
      b%x=bmod*def%x
      b%y=bmod*def%y
      b%z=bmod*def%z
    end if
  end if

!----------------TIME LOOP--------------------------
  do istep=1,nstop
     
    t=dt*(istep-1)
  
    if(extfield_type.eq.'impulse')then
      bmod=impulse(t)
      if( extfield_geometry .eq. 'line')then
        b%x=bmod*def%x
        b%y=bmod*def%y
        b%z=bmod*def%z
      end if 
    end if 

!------------initially make the effective field equal to zero
    do is=1,ns
      sp(is)%h%x=0.0
      sp(is)%h%y=0.0
      sp(is)%h%z=0.0
    end do

!-----------------calculate zeeman contribution--------------
    call zeeman(ez)
    
!-----------------calculate demag contribution---------------
    call demag(ed)

!----------------calculate exchange contribution------------
    call exchange(ee)

!----------------calculate anisotropy contribution--------------
    call anisotropy(ea)

!---------------sum up all contributions--------------------------
    et=ez+ed+ee+ea

!--------------make LLG step---------------------------
    call llg_step()

!-----------control information output-------------------------
    if(mod(istep,iout).eq.0)then
      print*,' istep=',istep,' t=',t,' bmod=',bmod,' ez=',ez,' ed=',ed,' ee=',ee,' ea=',ea,' et=',et,      &
	  ' max_torque=',max_torque,' max_angle=',max_angle
      write(8,*)' istep=',istep,' t=',t,' bmod=',bmod,' ez=',ez,' ed=',ed,' ee=',ee,' ea=',ea,' et=',et,   &
	  ' max_torque=',max_torque,' max_angle=',max_angle

!------------output data per unit volume to compare with OOMMF------------------------
 
      if(vol.ne.0.0)then
        ed_1=ed/vol
        ee_1=ee/vol
        ea_1=ea/vol
        ez_1=ez/vol
        et_1=et/vol
!------------output for comparison with OOMMF---------------------------------------
        write(16,*)' istep=',istep,' ez_1=',ez_1,' ed_1=',ed_1,        &
		           ' ee_1=',ee_1,' ea_1=',ea_1,' et_1=',et_1 
      else
        write(16,*)' sizex=',sizex,' sizey=',sizey,' sizez=',sizez,' vol=',vol
      end if		    

!-------output 3d movie-snapshot--istep---------------------------
      if(extfield_type.ne.'loop')then
	write(11,'(a1)')'c'
        write(11,'(2x,i6,2x,f10.3,2x,3(g12.5,1x))')istep,bmod,ez,ed,ee,ea,et
        do is=1,ns
        write(11,'(2x,3(g12.5,1x),3(f12.5,1x))')sp(is)%r%x,sp(is)%r%y,sp(is)%r%z,  &
                                              sp(is)%m%x,sp(is)%m%y,sp(is)%m%z
        end do

        if( extfield_geometry .eq. 'circular')then
          tm_tangent=0.d0
          do is=1,ns
            radius= dsqrt( sp(is)%r%x**2 + sp(is)%r%y**2 )
            tangent%x= sp(is)%r%y/radius
            tangent%y=-sp(is)%r%x/radius   
            tm_tangent=tm_tangent+sp(is)%m%x*sp(is)%cp%ms*tangent%x + sp(is)%m%y*sp(is)%cp%ms*tangent%y    
          end do
          write(20,'(2x,i7,2x,g12.3,2x,f10.3,2x,g12.3)')istep,t,bmod,tm_tangent
        else
          tmx=0.0d0
	  tmy=0.0d0
	  tmz=0.0d0
          do is=1,ns
            tmx=tmx+sp(is)%m%x*sp(is)%cp%ms   !-----24/05/11 added multiply to Ms
            tmy=tmy+sp(is)%m%y*sp(is)%cp%ms
            tmz=tmz+sp(is)%m%z*sp(is)%cp%ms
          end do
          write(20,'( 2x,i7,2x,g12.3,2x,f10.3,2x,3(g12.3,1x) )')istep,t,bmod,tmx,tmy,tmz
        end if

      end if
    end if
!-----------------------------------------------------------------

    if(numerics_type.eq.'llg_min' .and. max_torque.le.torque_stop)then
      write(8,*)' minimization stoped  max_torque=',max_torque 
      exit 
    end if 

  end do  ! end of loop on time


  if(extfield_type.eq.'loop')then  !  output movie and loop

!-------------------compute total magnetization and output movie snapshot
    write(11,'(a1)')'c'
    write(11,'(2x,i6,2x,f10.3,2x,3(g12.5,1x))')istep,bmod,ez,ed,ee,ea,et

    if( extfield_geometry .eq. 'circular')then
       tm_tangent=0.d0
              
       do is=1,ns
        
         radius= dsqrt( sp(is)%r%x**2 + sp(is)%r%y**2 )

         tangent%x= sp(is)%r%y/radius

         tangent%y=-sp(is)%r%x/radius   
         
         tm_tangent=tm_tangent+sp(is)%m%x*sp(is)%cp%ms*tangent%x + sp(is)%m%y*sp(is)%cp%ms*tangent%y    

        write(11,'(2x,3(g12.5,1x),3(f12.5,1x))')sp(is)%r%x,sp(is)%r%y,sp(is)%r%z,  &
                                                sp(is)%m%x,sp(is)%m%y,sp(is)%m%z

      end do
   
      print*,'-----CIRCULAR---------- bmod=',bmod,' tm_tangent=',tm_tangent,' et=',et
      write(8,*)'--CIRCULAR---------- bmod=',bmod,' tm_tangent=',tm_tangent,' et=',et

      write(12,'(2x,f9.3,2x,3(g12.4,1x),1x,g12.5)')bmod,tm_tangent,et
       
    else  !  linear field
      tmx=0.0d0
	  tmy=0.0d0
	  tmz=0.0d0
      do is=1,ns
        tmx=tmx+sp(is)%m%x*sp(is)%cp%ms   !-----24/05/11 added multiply to Ms
        tmy=tmy+sp(is)%m%y*sp(is)%cp%ms
        tmz=tmz+sp(is)%m%z*sp(is)%cp%ms
        write(11,'(2x,3(g12.5,1x),3(f12.5,1x))')sp(is)%r%x,sp(is)%r%y,sp(is)%r%z,  &
                                                sp(is)%m%x,sp(is)%m%y,sp(is)%m%z
      end do
   
      print*,'----LINEAR----------- bmod=',bmod,' tmx=',tmx,' tmy=',tmy,' tmz=',tmz,' et=',et
      write(8,*)'--LINEAR------------- bmod=',bmod,' tmx=',tmx,' tmy=',tmy,' tmz=',tmz,' et=',et

      write(12,'(2x,f9.3,2x,3(g12.4,1x),1x,g12.5)')bmod,tmx,tmy,tmz,et
    end if

  end if
end do  !  end of loop on external field (B)

close(11)

close(12)


call output_magfield("magfield1.map  ")

!---------------------------------------------------------

close(16)
close(8)

stop

end
