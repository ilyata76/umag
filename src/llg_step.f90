!---------MAKE LLG STEP--------------------------------

subroutine llg_step()

use magmod

implicit real*8 (a-h,o-z)

type (vector) heff ! effective field 

type (vector) part1, part2, m, dm 
type (vector) m_x_heff, m_x_m_x_heff

max_torque=0.0d0

factor1=-dt*gamma

factor2=-dt*alfa*gamma

do is=1,ns

 !  effective field on given cell

  heff%x=sp(is)%h%x 
  heff%y=sp(is)%h%y 
  heff%z=sp(is)%h%z 
  
 !  magnetization on given cell

  m%x=sp(is)%m%x
  m%y=sp(is)%m%y
  m%z=sp(is)%m%z

  m_x_heff%x=m%y * heff%z - m%z * heff%y
  m_x_heff%y=m%z * heff%x - m%x * heff%z
  m_x_heff%z=m%x * heff%y - m%y * heff%x

  torque2=m_x_heff%x * m_x_heff%x + m_x_heff%y * m_x_heff%y + m_x_heff%z * m_x_heff%z

  abs_torque=dsqrt(torque2)

  if(abs_torque.gt.max_torque)max_torque=abs_torque  

  part1%x=m_x_heff%x*factor1
  part1%y=m_x_heff%y*factor1
  part1%z=m_x_heff%z*factor1

  m_x_m_x_heff%x=m%y * m_x_heff%z - m%z * m_x_heff%y
  m_x_m_x_heff%y=m%z * m_x_heff%x - m%x * m_x_heff%z
  m_x_m_x_heff%z=m%x * m_x_heff%y - m%y * m_x_heff%x

  part2%x=m_x_m_x_heff%x*factor2
  part2%y=m_x_m_x_heff%y*factor2
  part2%z=m_x_m_x_heff%z*factor2

  dm%x=part1%x+part2%x
  dm%y=part1%y+part2%y
  dm%z=part1%z+part2%z
 
  sp(is)%m%x=sp(is)%m%x+dm%x
  sp(is)%m%y=sp(is)%m%y+dm%y
  sp(is)%m%z=sp(is)%m%z+dm%z

!-----------normalize vector----------------------------------------
   
  ss=dsqrt(sp(is)%m%x**2 + sp(is)%m%y**2 +sp(is)%m%z**2)

  if(ss.ne.0.d0)then
    sp(is)%m%x=sp(is)%m%x/ss
    sp(is)%m%y=sp(is)%m%y/ss
    sp(is)%m%z=sp(is)%m%z/ss
  else
    sp(is)%m%x=0.0d0
    sp(is)%m%y=0.0d0
    sp(is)%m%z=0.0d0
  end if

end do

return
end
