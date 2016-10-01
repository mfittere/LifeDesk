      program speclumi
      implicit none
      real l,l0,i
c
      open(UNIT=3,FILE='intensity.txt')
      open(UNIT=4,FILE='lumi.txt')
c
      read(UNIT=4,FMT=*,END=10) l
      read(UNIT=3,FMT=*,END=10) i
      l0=l
      write(UNIT=6,FMT='(E14.6)') l/l0*i
      do 
      read(UNIT=4,FMT=*,END=10) l
      read(UNIT=3,FMT=*,END=10) i
      write(UNIT=6,FMT='(E14.6)') l/l0*i
      end do
c   
 10   close(3)
      close(4)
      end
