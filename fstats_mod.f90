!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  fstats_mod
!!  module for calculating various statistical quantities
!!
!!  20.12.2012
!!  Copyright 2012 Jakob Simeth
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fstats_mod
    implicit none
    private
    integer, parameter :: allocsize=1000 !allocate allocsize reals at once
    real, allocatable:: alldata(:), temp(:)
    integer :: npts = 0
    real :: mean_cached = 0.0, rho0_cached = 0.0
    logical :: m_is_cached = .false., rho0_is_cached = .false.

    !constants:
    integer, parameter,public :: ERROR_NAIVE=0, ERROR_BOOTSTRAP=1

    !interface:
    interface fstats_add
        module procedure add_real
        module procedure add_realarray
    end interface 
    public :: fstats_add, fstats_free, & 
        & fstats_mean, fstats_rho, fstats_acf, fstats_error, fstats_tau_int

    contains
        !! add values to container
        subroutine add_real(rnum)
            real, intent(in) :: rnum
            if(mod(npts,allocsize)==0) then
                if(.not.allocated(alldata)) then
                    !very first allocation here:
                    allocate(alldata(allocsize))
                else
                    !all others: copy to temp and then back.
                    allocate(temp(1:npts), source=alldata)
                    deallocate(alldata)
                    allocate(alldata(npts+allocsize))
                    alldata(1:npts) = temp
                    deallocate(temp)
                end if
            end if
            npts = npts + 1
            alldata(npts) = rnum
            m_is_cached = .false.
            rho0_is_cached = .false.
        end subroutine
        
        !! add vector of reals to container
        subroutine add_realarray(rnum)
            real, intent(in) :: rnum(:)
            integer :: newblocks, usedblocks

            !how many new blocks of data must be allocated?
            newblocks = size(rnum) / allocsize
            !if rest does not fit into already allocated space, allocate one
            !block more:
            if(.not. allocated(alldata) .or. size(alldata) - npts < mod(size(rnum),allocsize)) then
                newblocks = newblocks + 1
            end if

            !take new memory if needed:
            if(.not.allocated(alldata)) then
                !very first allocation here:
                allocate(alldata(newblocks*allocsize))
            else if (newblocks > 0) then
                !all others: copy to temp and then back.
                usedblocks=npts/allocsize+1
                allocate(temp(1:npts), source=alldata(1:npts))
                deallocate(alldata)
                allocate(alldata((usedblocks+newblocks)*allocsize))
                alldata(1:npts) = temp
                deallocate(temp)
            end if
            alldata(npts+1:npts+size(rnum)) = rnum
            npts = npts + size(rnum)
            m_is_cached = .false.
            rho0_is_cached = .false.
        end subroutine

        !!free memory (forget about data)
        subroutine fstats_free
            npts = 0
            m_is_cached = .false.
            rho0_is_cached = .false.
            deallocate(alldata)
        end subroutine

        !!simple arithmetic mean of saved data
        function fstats_mean()
            real :: sum, fstats_mean
            integer :: i
            sum = 0.0
            !calculate only if not cached
            if(.not. m_is_cached) then
                write(*,*) "npts=", npts, "sum=", sum
                do i=1,npts
                    sum = sum + alldata(i)
                end do
                mean_cached = sum / real(npts)
                m_is_cached = .true.
            end if
            fstats_mean = mean_cached
        end function


        !!all autocorrelation things here:
        !! rho0 = R(tau=0)
        function rho0()
            real :: rho0, m
            integer :: t
            if(.not. rho0_is_cached) then
                rho0_cached = 0.0
                m = fstats_mean()
                do t=1,npts
                    rho0_cached = rho0_cached + (alldata(t)-m)**2
                end do
                rho0_cached = rho0_cached / real(npts)
                rho0_is_cached = .true.
            end if
            rho0 = rho0_cached
        end function
        !! autocorrelation function rho(tau) := R(tau)/R(0)
        function fstats_rho(tau)
            integer, intent(in) :: tau
            integer :: nmax, i, t
            real :: fstats_rho, m
            fstats_rho = 0.0
            m = fstats_mean()
            do t=1,npts-tau
                fstats_rho = fstats_rho + (alldata(t)-m)*(alldata(t+tau)-m)
            end do
            fstats_rho = fstats_rho / (real(npts-tau)*rho0())
        end function
        !! complete autocorrelation function
        function fstats_acf(maxlength)
            real,allocatable :: fstats_acf(:)
            integer, optional :: maxlength
            integer :: maxlength_, i
            if(present(maxlength).and.maxlength <= npts) then
                maxlength_ = maxlength
            else
                maxlength_ = int(log(real(npts)))
            end if
            allocate(fstats_acf(maxlength_))
            do i=1,maxlength_
                fstats_acf(i) = fstats_rho(i)
            end do
        end function

        !! integrated auto correlation function, up to taumax
        function tau_int_max(taumax)
            integer, intent(in) :: taumax
            real :: tau_int_max
            integer :: i
            tau_int_max = 0.5
            do i=1,taumax
                tau_int_max = tau_int_max + fstats_rho(i)
            end do
        end function
        !! integrated auto corr fn, self-determined truncation of sum
        function fstats_tau_int()
            real :: fstats_tau_int, tau_int_old
            ! iteration stops, when tau_int
            ! doesn't change more than eps
            real, parameter :: eps = 0.01 
            !starting value for tau_int
            integer, parameter :: start_tau_max = 10 
            integer, parameter :: maxiter = 1000
            integer :: start_tau,i
            start_tau = min(start_tau_max, npts)
            tau_int_old = tau_int_max(start_tau)
            do i=1,maxiter
                !rule of thumb: truncate sum after 4*tau_int
                fstats_tau_int = tau_int_max(int(4*tau_int_old))
                if(abs((tau_int_old - fstats_tau_int)/(tau_int_old + fstats_tau_int)) < eps) then
                    exit
                end if
                tau_int_old = fstats_tau_int
            end do
        end function

        !!BINNING
        !
        function error_binning(binsize, verbose, error_method)
            !error_method may be one of the following:
            !ERROR_NAIVE
            !ERROR_BOOTSTRAP
            real :: error_binning, bin_mean, rn
            integer, intent(in) :: binsize
            integer, intent(in), optional :: error_method
            integer :: error_method_
            logical, intent(in),optional :: verbose
            logical :: verbose_t
            real, allocatable :: bins(:), bsdata(:)
            integer :: i, nbins, ind, j

            if(present(verbose)) then
                verbose_t = verbose
            else
                verbose_t = .false.
            end if
            if(present(error_method)) then
                error_method_ = error_method
            else
                error_method_ = ERROR_NAIVE
            end if

            nbins = npts/binsize

            allocate(bins(nbins))
            
            if(verbose_t) then
                write(6,*) "binning with binsize=", binsize, ", i.e. ", nbins, " bins" 
                write(6,*) "data covered:        ", real(nbins*binsize*100)/real(npts), "%"
            end if

            !bin, i.e. fill each bin with the mean of (binsize) data
            do i = 1, nbins
                bins(i) = sum(alldata((i-1)*binsize+1:i*binsize))/real(binsize)
            end do
            
            !as it is possible, that not all data is covered in bins, recompute
            !the bin-mean (and do not use the overall data-mean):
            bin_mean = sum(bins)/real(nbins)

            !binned error:
            select case(error_method_)
                case(ERROR_NAIVE)
                    if(verbose_t) then
                        write(6,*) "using naive error"
                    end if
                    error_binning=0.0
                    do i=1,nbins
                        error_binning=error_binning + (bins(i)-bin_mean)**2
                    end do
                    error_binning = sqrt(error_binning / real(nbins*(nbins-1)))
                case(ERROR_BOOTSTRAP)
                    if(verbose_t) then
                        write(6,*) "using bootstrap error"
                    end if
                    allocate(bsdata(10*nbins))
                    !sample data:
                    !draw 10*N samples, each of length N
                    do i=1,10*nbins
                        bsdata(i) = 0.0
                        do j=1,nbins
                            !randomly select a bin:
                            call random_number(rn)
                            ind=int(rn*nbins + 1) !round to next int
                            bsdata(i) = bsdata(i) + bins(ind)
                        end do
                        bsdata(i) = bsdata(i) / real(nbins)
                    end do
                    !calculate error:
                    error_binning=0.0
                    do i=1,10*nbins
                        error_binning=error_binning + (bsdata(i)-bin_mean)**2
                    end do
                    error_binning = sqrt(error_binning / real(10*nbins-1))
                    deallocate(bsdata)
            end select

            deallocate(bins)
        end function
        function error_binning_auto(error_method, iter)
            !find maximal error by binning
            !stops, when change is less than eps
            !iter is intent(out)! returns number of iteration steps
            real :: error_binning_auto
            integer, optional,intent(out) :: iter
            integer :: iter_
            integer, parameter :: startexp=1, maxsizeexp=10
            integer, intent(in), optional :: error_method
            real, parameter :: eps = 0.1
            real :: oldse, newse, maxse
            integer isize
            oldse = error_binning(2**startexp)
            maxse = oldse
            iter_ = 0
            do isize=startexp+1,maxsizeexp
                iter_ = iter_ + 1
                newse = error_binning(2**isize, .false., error_method)
                if(newse > maxse) maxse = newse
                if(abs(oldse-newse)/(oldse+newse) < eps) exit
            end do
            if(present(iter)) iter=iter_
            error_binning_auto = maxse
        end function
        !! sample variance := 1/(N-1) * sum(x - xbar)^2
        function var()
            real :: dat2, var, m
            integer :: i
            dat2 = 0.0
            m = fstats_mean()
            do i=1,npts
                dat2 = dat2 + (alldata(i)-m)**2
            end do
            var = dat2/real(npts-1)
        end function
        !! stderror := sqrt(var/N)
        function get_error_naive()
            real :: get_error_naive

            get_error_naive = sqrt(var()/real(npts))
        end function
        function fstats_error(error_method, binned)
            integer, intent(in), optional :: error_method
            integer :: err_method_
            logical, intent(in), optional :: binned
            logical :: bin_
            real :: fstats_error
            if(present(error_method)) then
                if(error_method /= ERROR_NAIVE .and. &
                    & error_method /= ERROR_BOOTSTRAP) then
                    write(0,*) "error: error_method must be " // &
                        & "either ERROR_NAIVE or ERROR_BOOTSTRAP"
                    fstats_error=0.0
                    return
                end if
                err_method_ = error_method
            else
                err_method_ = ERROR_NAIVE
            end if
            if(present(binned)) then
                bin_ = binned
            else
                bin_ = .false.
            end if

            if(bin_ .or. err_method_ == ERROR_BOOTSTRAP) then
                fstats_error = error_binning_auto(error_method)
            else
                fstats_error = get_error_naive()
            end if
        end function

end module
