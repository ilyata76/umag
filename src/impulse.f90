!--------------IMPULSE---(time dependence of external B field)-----------------------
!
!----------------------------------------------------------------------------

real*8 function impulse(t)

use magmod

implicit real*8 (a-h,o-z)

it=t/ht_imp+1

if(it.lt.nt_imp)then
  t_it=ht_imp*(it-1)
  bm=bmt(it)+(t-t_it)*(bmt(it+1)-bmt(it))/ht_imp
else
  bm=bmt(nt_imp)
end if

!write(8,*)' IMPULSE:  t=',t,' it=',it,' t_it=',t_it,' bmt(it)=',bmt(it),' bm=',bm

impulse=bm

return
end