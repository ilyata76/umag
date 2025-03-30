!-----------compute Exchange contribution-----------

!
!------------------------------------------------------------

subroutine exchange(ee)

use magmod
implicit real*8 (a-h,o-z)

type (vector) d2m_dx2,d2m_dy2,d2m_dz2,hei

type (vector) lmi

real*8 mi_lmi

real*8 ms

ee=0.0d0


max_angle=0.0d0

do i=1,ns

  a=sp(i)%cp%a
  ms=sp(i)%cp%ms

  factor=2.0*a/(mu0*ms)   !  factor 2 added 24/03/11

  sp_i_mx=sp(i)%m%x
  sp_i_my=sp(i)%m%y
  sp_i_mz=sp(i)%m%z

!-------------second derivative along X------D2M_DX2--------------------
!---------------neighbor from the right (x plus)---------------
  jp=sp(i)%cp%xp
!---------------neighbor from the left (x minus)---------------
  jm=sp(i)%cp%xm

!--------------FREE BORDER CONDITION-----------------------

  if(jp.eq.0)then
    sp_jp_mx=sp_i_mx  
    sp_jp_my=sp_i_my  
    sp_jp_mz=sp_i_mz  
  else
    sp_jp_mx=sp(jp)%m%x
    sp_jp_my=sp(jp)%m%y
    sp_jp_mz=sp(jp)%m%z
  end if

  if(jm.eq.0)then
    sp_jm_mx=sp_i_mx  
    sp_jm_my=sp_i_my  
    sp_jm_mz=sp_i_mz  
  else
    sp_jm_mx=sp(jm)%m%x
    sp_jm_my=sp(jm)%m%y  
    sp_jm_mz=sp(jm)%m%z
  end if

  d2m_dx2%x=(sp_jp_mx-2.0*sp_i_mx+sp_jm_mx)/(hx*hx) 
  d2m_dx2%y=(sp_jp_my-2.0*sp_i_my+sp_jm_my)/(hx*hx) 
  d2m_dx2%z=(sp_jp_mz-2.0*sp_i_mz+sp_jm_mz)/(hx*hx) 

!-calculate angle between Magnetization vectors I - JP and determine maximal angle

  spmi=dsqrt(sp_i_mx*sp_i_mx + sp_i_my*sp_i_my + sp_i_mz*sp_i_mz)
  
  spmjp=dsqrt(sp_jp_mx*sp_jp_mx + sp_jp_my*sp_jp_my + sp_jp_mz*sp_jp_mz)

  cos_i_jp=(sp_i_mx*sp_jp_mx+sp_i_my*sp_jp_my+sp_i_mz*sp_jp_mz)/(spmi*spmjp)

  if(cos_i_jp.gt.1.0d0)then
    angle=0.0d0
  else 
    if(cos_i_jp.lt.-1.0d0)then
	  angle=180.0d0
    else
      angle=acosd(cos_i_jp)
    end if 
  end if

  if(angle.gt.max_angle)max_angle=angle

!-calculate angle between Magnetization vectors I - JM and determine maximal angle

  spmjm=dsqrt(sp_jm_mx*sp_jm_mx + sp_jm_my*sp_jm_my + sp_jm_mz*sp_jm_mz)

  cos_i_jm=(sp_i_mx*sp_jm_mx+sp_i_my*sp_jm_my+sp_i_mz*sp_jm_mz)/(spmi*spmjm)

  if(cos_i_jm.gt.1.0d0)then
    angle=0.0d0
  else 
    if(cos_i_jm.lt.-1.0d0)then
	  angle=180.0d0
    else
      angle=acosd(cos_i_jm)
    end if 
  end if

  if(angle.gt.max_angle)max_angle=angle

!-------------second derivative along Y------D2M_DY2--------------------
!---------------neighbor from the right (y plus)---------------
  jp=sp(i)%cp%yp
!---------------neighbor from the left (y minus)---------------
  jm=sp(i)%cp%ym

!--------------FREE BORDER CONDITION-----------------------

  if(jp.eq.0)then
    sp_jp_mx=sp_i_mx  
    sp_jp_my=sp_i_my  
    sp_jp_mz=sp_i_mz  
  else
    sp_jp_mx=sp(jp)%m%x
    sp_jp_my=sp(jp)%m%y
    sp_jp_mz=sp(jp)%m%z
  end if

  if(jm.eq.0)then
    sp_jm_mx=sp_i_mx  
    sp_jm_my=sp_i_my  
    sp_jm_mz=sp_i_mz  
  else
    sp_jm_mx=sp(jm)%m%x
    sp_jm_my=sp(jm)%m%y  
    sp_jm_mz=sp(jm)%m%z
  end if

  d2m_dy2%x=(sp_jp_mx-2.0*sp_i_mx+sp_jm_mx)/(hy*hy) 
  d2m_dy2%y=(sp_jp_my-2.0*sp_i_my+sp_jm_my)/(hy*hy) 
  d2m_dy2%z=(sp_jp_mz-2.0*sp_i_mz+sp_jm_mz)/(hy*hy) 

!-calculate angle between Magnetization vectors I - JP and determine maximal angle

  spmjp=dsqrt(sp_jp_mx*sp_jp_mx + sp_jp_my*sp_jp_my + sp_jp_mz*sp_jp_mz)

  cos_i_jp=(sp_i_mx*sp_jp_mx+sp_i_my*sp_jp_my+sp_i_mz*sp_jp_mz)/(spmi*spmjp)

  if(cos_i_jp.gt.1.0d0)then
    angle=0.0d0
  else 
    if(cos_i_jp.lt.-1.0d0)then
	  angle=180.0d0
    else
      angle=acosd(cos_i_jp)
    end if 
  end if

  if(angle.gt.max_angle)max_angle=angle

!-calculate angle between Magnetization vectors I - JM and determine maximal angle

  spmjm=dsqrt(sp_jm_mx*sp_jm_mx + sp_jm_my*sp_jm_my + sp_jm_mz*sp_jm_mz)

  cos_i_jm=(sp_i_mx*sp_jm_mx+sp_i_my*sp_jm_my+sp_i_mz*sp_jm_mz)/(spmi*spmjm)

  if(cos_i_jm.gt.1.0d0)then
    angle=0.0d0
  else 
    if(cos_i_jm.lt.-1.0d0)then
	  angle=180.0d0
    else
      angle=acosd(cos_i_jm)
    end if 
  end if

  if(angle.gt.max_angle)max_angle=angle

!-------------second derivative along Z------D2M_DZ2--------------------
!---------------neighbor from the right (z plus)---------------
  jp=sp(i)%cp%zp
!---------------neighbor from the left (z minus)---------------
  jm=sp(i)%cp%zm

!--------------FREE BORDER CONDITION-----------------------

  if(jp.eq.0)then
    sp_jp_mx=sp_i_mx  
    sp_jp_my=sp_i_my  
    sp_jp_mz=sp_i_mz  
  else
    sp_jp_mx=sp(jp)%m%x
    sp_jp_my=sp(jp)%m%y
    sp_jp_mz=sp(jp)%m%z
  end if

  if(jm.eq.0)then
    sp_jm_mx=sp_i_mx  
    sp_jm_my=sp_i_my  
    sp_jm_mz=sp_i_mz  
  else
    sp_jm_mx=sp(jm)%m%x
    sp_jm_my=sp(jm)%m%y  
    sp_jm_mz=sp(jm)%m%z
  end if

  d2m_dz2%x=(sp_jp_mx-2.0*sp_i_mx+sp_jm_mx)/(hz*hz) 
  d2m_dz2%y=(sp_jp_my-2.0*sp_i_my+sp_jm_my)/(hz*hz) 
  d2m_dz2%z=(sp_jp_mz-2.0*sp_i_mz+sp_jm_mz)/(hz*hz) 


!-calculate angle between Magnetization vectors I - JP and determine maximal angle

  spmjp=dsqrt(sp_jp_mx*sp_jp_mx + sp_jp_my*sp_jp_my + sp_jp_mz*sp_jp_mz)

  cos_i_jp=(sp_i_mx*sp_jp_mx+sp_i_my*sp_jp_my+sp_i_mz*sp_jp_mz)/(spmi*spmjp)

  if(cos_i_jp.gt.1.0d0)then
    angle=0.0d0
  else 
    if(cos_i_jp.lt.-1.0d0)then
	  angle=180.0d0
    else
      angle=acosd(cos_i_jp)
    end if 
  end if

  if(angle.gt.max_angle)max_angle=angle

!-calculate angle between Magnetization vectors I - JM and determine maximal angle

  spmjm=dsqrt(sp_jm_mx*sp_jm_mx + sp_jm_my*sp_jm_my + sp_jm_mz*sp_jm_mz)

  cos_i_jm=(sp_i_mx*sp_jm_mx+sp_i_my*sp_jm_my+sp_i_mz*sp_jm_mz)/(spmi*spmjm)

  if(cos_i_jm.gt.1.0d0)then
    angle=0.0d0
  else 
    if(cos_i_jm.lt.-1.0d0)then
	  angle=180.0d0
    else
      angle=acosd(cos_i_jm)
    end if 
  end if

  if(angle.gt.max_angle)max_angle=angle

!-------------LAPLASIAN=sum of derivatives--------------------------------

  lmi%x=d2m_dx2%x+d2m_dy2%x+d2m_dz2%x
  lmi%y=d2m_dx2%y+d2m_dy2%y+d2m_dz2%y
  lmi%z=d2m_dx2%z+d2m_dy2%z+d2m_dz2%z

!-------------------------contribution to effective field---------------------

  hei%x=lmi%x*factor
  hei%y=lmi%y*factor
  hei%z=lmi%z*factor

  sp(i)%h%x=sp(i)%h%x + hei%x
  sp(i)%h%y=sp(i)%h%y + hei%y
  sp(i)%h%z=sp(i)%h%z + hei%z

  hei_mod=dsqrt(hei%x**2+hei%y**2+hei%z**2)

!-------------------------contribution to exchange energy------------
  mi_lmi=sp(i)%m%x*lmi%x+sp(i)%m%y*lmi%y+sp(i)%m%z*lmi%z 

  ee = ee - a*mi_lmi 

end do

ee=ee*cvol

return
end