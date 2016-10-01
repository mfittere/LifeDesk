      program lossrate
      implicit none
      integer n,n0
c
      read(UNIT=5,FMT=*,END=10) n
      write(UNIT=6,FMT='(I6)') n
      n0=n
c
      do 
      read(UNIT=5,FMT=*,END=10) n
      write(UNIT=6,FMT='(I6)') n-n0
      n0=n
      end do
c   
 10   close(5)
      end
