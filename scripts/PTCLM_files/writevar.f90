program writevar

use netcdf
implicit none

integer i, j, k, ierr, npoints, ncid, nvars
integer year, month, dimid(4), varid(20)
integer numx, numy
character(len=200) filename_in, varnames_in(20)
character(len=4) yst, mst
double precision data_in(7, 1, 1, 1000), data_out(7,100,100,1000)
double precision lat_out(100,100), lon_out(100,100), zbot_out(100,100,1000)
double precision resx, resy, lat, lon, zbot, time_out(1000)

!Replace data in input meteorology netcdf file for CLM
!latitude, longitude, time, and met variables with a given datafile

!get filename
open(unit=8, file = 'data.in')
read(8,*) year, month
read(8,*) lat, lon, zbot
read(8,'(A200)') filename_in
read(8,*) nvars
read(8,*) varnames_in(1:nvars)
read(8,*) npoints
read(8,*) numx, numy, resx, resy

if (numx .gt. 100 .or. numy .gt. 100) then 
   print*, 'Error:  only up to 100x100 grid size is supported'
   return
end if

do i=1,npoints
   read(8,*) data_in(1:nvars,1,1,i)
end do
close(8)

write(yst,'(I4)') year
write(mst,'(I4)') 1000+month
ierr = nf90_create(trim(filename_in), nf90_clobber, ncid)
ierr = nf90_def_dim(ncid, 'time', npoints, dimid(3))
ierr = nf90_def_dim(ncid, 'lon', numx, dimid(1))
ierr = nf90_def_dim(ncid, 'lat', numy, dimid(2))
ierr = nf90_def_dim(ncid, 'scalar', 1, dimid(4))

ierr = nf90_def_var(ncid, 'time', nf90_double, dimid(3), varid(1))
ierr = nf90_put_att(ncid, varid(1), 'long_name', 'Time axis')
ierr = nf90_put_att(ncid, varid(1), 'units', 'days since ' // yst // &
     '-' // mst(3:4) // '-01 00:00:00')
ierr = nf90_put_att(ncid, varid(1), 'calendar', 'noleap')
ierr = nf90_def_var(ncid, 'LONGXY', nf90_double, dimid(1:2), varid(2))
ierr = nf90_put_att(ncid, varid(2), 'long_name', 'longitude')
ierr = nf90_put_att(ncid, varid(2), 'units', 'degrees E')
ierr = nf90_def_var(ncid, 'LATIXY', nf90_double, dimid(1:2), varid(3))
ierr = nf90_put_att(ncid, varid(3), 'long_name', 'latitude')
ierr = nf90_put_att(ncid, varid(3), 'units', 'degrees N')
ierr = nf90_def_var(ncid, 'ZBOT', nf90_double, dimid(1:3), varid(4))
ierr = nf90_put_att(ncid, varid(4), 'long_name', 'observational height')
ierr = nf90_put_att(ncid, varid(4), 'units', 'm')
ierr = nf90_def_var(ncid, 'EDGEW', nf90_double, dimid(4), varid(5))
ierr = nf90_put_att(ncid, varid(5), 'long_name', 'western edge in atmospheric data')
ierr = nf90_put_att(ncid, varid(5), 'units', 'degrees N')
ierr = nf90_def_var(ncid, 'EDGEE', nf90_double, dimid(4), varid(6))
ierr = nf90_put_att(ncid, varid(6), 'long_name', 'eastern edge in atmospheric data')
ierr = nf90_put_att(ncid, varid(6), 'units', 'degrees N')
ierr = nf90_def_var(ncid, 'EDGES', nf90_double, dimid(4), varid(7))
ierr = nf90_put_att(ncid, varid(7), 'long_name', 'southern edge in atmospheric data')
ierr = nf90_put_att(ncid, varid(7), 'units', 'degrees N')
ierr = nf90_def_var(ncid, 'EDGEN', nf90_double, dimid(4), varid(8))
ierr = nf90_put_att(ncid, varid(8), 'long_name', 'northern edge in atmospheric data')
ierr = nf90_put_att(ncid, varid(8), 'units', 'degrees N')
ierr = nf90_def_var(ncid, 'TBOT', nf90_double, dimid(1:3), varid(9))
ierr = nf90_put_att(ncid, varid(9), 'long_name', 'temperature at the lowest atm level (TBOT)')
ierr = nf90_put_att(ncid, varid(9), 'units', 'K')
ierr = nf90_def_var(ncid, 'RH', nf90_double, dimid(1:3), varid(14))
ierr = nf90_put_att(ncid, varid(14), 'long_name', 'relative humidity at the lowest atm level (RH)')
ierr = nf90_put_att(ncid, varid(14), 'units', '%')
ierr = nf90_def_var(ncid, 'WIND', nf90_double, dimid(1:3), varid(15))
ierr = nf90_put_att(ncid, varid(15), 'long_name', 'wind at the lowest atm level (WIND)')
ierr = nf90_put_att(ncid, varid(15), 'units', 'm/s')
ierr = nf90_def_var(ncid, 'FSDS', nf90_double, dimid(1:3), varid(12))
ierr = nf90_put_att(ncid, varid(12), 'long_name', 'incident shortwave (FSDS)')
ierr = nf90_put_att(ncid, varid(12), 'units', 'W/m2')
ierr = nf90_def_var(ncid, 'FLDS', nf90_double, dimid(1:3), varid(11))
ierr = nf90_put_att(ncid, varid(11), 'long_name', 'incident longwave (FLDS)')
ierr = nf90_put_att(ncid, varid(11), 'units', 'W/m2')
ierr = nf90_def_var(ncid, 'PSRF', nf90_double, dimid(1:3), varid(13))
ierr = nf90_put_att(ncid, varid(13), 'long_name', 'pressure at the lowest atm level (PSRF)')
ierr = nf90_put_att(ncid, varid(13), 'units', 'Pa')
ierr = nf90_def_var(ncid, 'PRECTmms', nf90_double, dimid(1:3), varid(10))
ierr = nf90_put_att(ncid, varid(10), 'long_name', 'precipitation (PRECTmms)')
ierr = nf90_put_att(ncid, varid(10), 'units', 'mm/s')
ierr = nf90_enddef(ncid)

do i=1, numx
   do j = 1, numy
      lon_out(i,j) = lon + resx*(i-1)
      lat_out(i,j) = lat + resy*(j-1)
      do k=1,npoints
         zbot_out(i,j,k) = zbot
         data_out(:,i,j,k) = data_in(:,1,1,k)
         time_out(k) = (k-1)/24.0  !+0.5/24.0
      end do
   end do
end do

ierr = nf90_put_var(ncid, varid(1), time_out(1:npoints))
ierr = nf90_put_var(ncid, varid(2), lon_out(1:numx, 1:numy))
ierr = nf90_put_var(ncid, varid(3), lat_out(1:numx, 1:numy))
ierr = nf90_put_var(ncid, varid(4), zbot_out(1:numx,1:numy,1:npoints))
ierr = nf90_put_var(ncid, varid(5), lon-resx/2)
ierr = nf90_put_var(ncid, varid(6), lon+(numx-1)*resx+resx/2)
ierr = nf90_put_var(ncid, varid(7), lat-resy/2)
ierr = nf90_put_var(ncid, varid(8), lat+(numy-1)*resy+resy/2)
do i = 1,nvars
   ierr = nf90_put_var(ncid, varid(8+i), data_out(i,1:numx,1:numy,1:npoints))
end do
ierr = nf90_close(ncid)

end program writevar
