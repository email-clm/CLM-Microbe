#!/bin/bash

source ~/.bashrc

# previous-month: 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
d=( 31 31 28 31 30 31 30 31 31 30 31 30 )

y0=1985
y1=2015
fileheader=../CLM1PT_data/

ymd=0
for y in {1985..2015}
do
  for m in {1..12}
  do

     if [ $m -lt 10 ]; then
       echo "File: $fileheader$y-0$m.nc"
       ncks -O -h --mk_rec_dim time $fileheader$y-0$m.nc -o tmp.nc
     else
       echo "File: $fileheader$y-$m.nc"
       ncks -O -h --mk_rec_dim time $fileheader$y-$m.nc -o tmp.nc
     fi

     if [ $y -eq $y0 ] && [ $m -eq 1 ]; then
       mv tmp.nc all_hourly0.nc
     else
       
       let ymd=ymd+${d[$m-1]}

       ncap2 -O -h -s "time=time+$ymd" tmp.nc -o tmp.nc
   
       ncrcat -O -h all_hourly0.nc tmp.nc -o all_hourly0.nc

     fi

     echo "--------------"


  done
done

rm -f tmp.nc
ncrename -O -d time,DTIME all_hourly0.nc -o all_hourly0.nc
ncrename -O -v time,DTIME all_hourly0.nc -o all_hourly0.nc

#scaling precipitation by 1.e7 and air pressure by 1.e-2, 
#so that if rounding to integer will not produce large error when scaling back later on. 
ncap2 -O -s "PRECTmms=PRECTmms*1.0e7;PSRF=PSRF*1.e-2" all_hourly0.nc all_hourly1.nc

# adding scaling factor
ncatted -O -a scale_factor,PRECTmms,a,d,"1.0e-7" all_hourly1.nc -o all_hourly2.nc
ncatted -O -a scale_factor,PSRF,a,d,"1.0e2" all_hourly2.nc -o all_hourly2.nc

ncks -O --fix_rec_dmn DTIME all_hourly2.nc -o all_hourly2.nc
ncpdq -a lat,lon,DTIME all_hourly2.nc all_hourly3.nc

# adding 2 constants
ncap2 -O -s "start_year=$y0" all_hourly3.nc -o all_hourly4.nc
ncap2 -O -s "end_year=$y1" all_hourly4.nc -o all_hourly4.nc

