!---------output Magnetization vector field------------


subroutine output_magfield(filename)

use magmod

implicit real*8 (a-h,o-z)

character*15 filename


open(10,file=filename)

do is=1,ns

  write(10,'(2x,3(g12.5,1x),1x,$)')sp(is)%r%x,sp(is)%r%y,sp(is)%r%z
  write(10,'(1x,3(g12.5,1x))')sp(is)%m%x,sp(is)%m%y,sp(is)%m%z

end do

close(10)

return 
end