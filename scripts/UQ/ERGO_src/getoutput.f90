program main

implicit none

character(len=6) gst
character(len=200) rundir, mycase
double precision sse, mysse
integer i

do i=1,64
  write(gst, '(I6)'), 100000+i
  rundir = '/home/zdr/models/clm4_5_ornl/run/UQ/US-SPR_I1850CLM45/g' // gst(2:6)
  mycase = 'US-SPR_I1850CLM45'
  print*, rundir
  print*, sse(rundir, mycase)
end do


end program main
