module globals

integer NVAR_MAX, NT_MAX
parameter(NVAR_MAX = 10)
parameter(NT_MAX   = 2000)

end module globals


function sse(rundir, mycase) 

use globals
use netcdf
implicit none

character(len=200) filename, mycase, outfile, optdir, dummy, rundir
character(len=50) varnames(NVAR_MAX), thisvar
characteR(len=6) rst, est
character(len=5) yst, nst
character(len=3) mst
character(len=21) pname_temp
integer i, r, y ,m, n, v, n_years, n_vars, nav, varid
integer myid, ierr, numprocs, npy, n_parms, ens_size, start_year
integer npf, ind(NVAR_MAX), thisgroup, thisinst
integer :: RCODE, NCID, opt
integer mystart(3), mycount(3)
double precision val(NT_MAX)
double precision obs_in(5), sse_min
double precision vals_temp(NVAR_MAX, NT_MAX)
double precision vals_out(NVAR_MAX,NT_MAX)     !n_vars, size, n_tsteps
double precision obs(NVAR_MAX,NT_MAX), parms_data(3)
double precision sse_part(NVAR_MAX), sse
real tempvalt(2,365), tempvalp(17)
varnames(:)=''

!------------------- user input ------------------------------------

start_year = 2010    !model year to start analysis
n_years    = 2      !number of years to average over
npy        = 1      !number of output files per model year
npf        = 365      !number of model output timesteps per file 
nav        = 1     !number of output timesteps to average over

!variable names
varnames(1)='ZWT'
varnames(2)='H2OSFC'
varnames(3)='QDRAI'
varnames(4)='QH2OSFC'

!---------------end user input------------------------------------------

!figure out how many varnames and parameters
n_vars=0
do n=1,NVAR_MAX
   if (varnames(n) .ne. '') n_vars = n_vars+1
end do

vals_out(:,:)=0d0
ind(:)=1
do y=1,n_years         !loop through years
   do m=1,npy          !loop though months
      write(yst,'(I5)') 10000+start_year+y-1  !year 
      write(mst,'(I3)') 100+m                 !month
        
      !get varnames
      if (npy .eq. 12) then   !default monthly output 
         filename=trim(rundir) // '/' // trim(mycase) // '.clm2.h0.' // yst(2:5) // '-' // mst(2:3) // '.nc'      
      end if
      if (npy .eq. 1) then     !annual output 
         filename = trim(rundir) // '/' // trim(mycase) // '.clm2.h0.' // yst(2:5) // '-01-01-00000.nc' 
      end if
      
      RCODE = NF90_OPEN(filename, NF90_NOWRITE, NCID)
      !print*, filename, RCODE

      do v=1,n_vars
         RCODE = NF90_INQ_VARID(NCID, trim(varnames(v)), varid)  
         RCODE = NF90_GET_VAR(NCID, varid, tempvalt(1:2,1:npf))
         !print*, trim(filename)
         !print*, v, tempvalt(1:2,1)
         
         do i=1,npf
            !print*, i
            val(i)=tempvalt(1,i)*1d0
            if (v .eq. 2 .or. v .eq. 4) val(i)=tempvalt(2,i)*1d0
            vals_temp(v,ind(v))=val(i)
            ind(v)=ind(v)+1
         end do
      end do
      RCODE = NF90_CLOSE(NCID)
   end do
end do

!do time averaging
do i=1,n_years*(npy*npf)/nav
   do v=1,n_vars
      if (v .eq. 1) then !ZWT depth
         vals_out(v,i) = sum(vals_temp(v,(i-1)*nav+1:i*nav))/nav*-1+0.3
      else
         vals_out(v,i) = sum(vals_temp(v,(i-1)*nav+1:i*nav))/nav
      end if
   end do
end do

!load observations
obs(:,:)=-999. 
filename='/home/zdr/models/clm45microbe_spruce_2015Jan15/scripts/UQ/daily_water_table_new.csv'
open(unit=8, file=trim(filename))
read(8,*) dummy 
n=1
do i=1,10000
   read(8,*,end=20) obs_in(1:4)
   if (obs_in(3) .ne. -999 .and. obs_in(4) .ne. -999) then 
      !obs(1,132+i) = (obs_in(2))/1000.
      !obs(2,132+i) = (obs_in(3))/1000.
      obs(1,i) = obs_in(3)/1000.
      obs(2,i) = obs_in(4)/1000.
   end if
end do
20 continue
close(8)


!calculate sse
sse_part(:)=0d0
sse_min=99999
open(unit=8, file='obscomp.txt')
do i=1,n_years*npy*npf
   if (obs(1,i) .ne. -999 .and. obs(2,i) .ne. -999) then 
      !print*, i, vals_out(1,i), obs(1,i)
      sse_part(1)=sse_part(1)+(obs(1,i)-vals_out(1,i))**2
      sse_part(2)=sse_part(2)+(obs(2,i)-vals_out(1,i))**2
   end if
   write(8, '(I4,1x,3(f13.6,1x))'), i, vals_out(1,i), obs(1,i), obs(2,i)
end do
close(unit=8)

sse = sum(sse_part)

end function sse


