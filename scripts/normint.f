      program normint
      implicit none
      real n,n0
c
      read(UNIT=5,FMT=*,END=10) n
      write(UNIT=6,FMT='(F12.7)') 1.0
      n0=n
c
      do 
      read(UNIT=5,FMT=*,END=10) n
      write(UNIT=6,FMT='(F12.7)') n/n0
      end do
c   
 10   close(5)
      end
