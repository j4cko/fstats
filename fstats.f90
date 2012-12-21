!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  fstats
!!  cl-frontend for fstats_mod
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
program fstats
    use fstats_mod
    implicit none

    integer :: readunit = 56
    real :: summed=0.0, tmp,testarr(100)
    integer :: stat, i
    logical :: fexists
    character(3) :: fcanread
    character(240) :: fname

    
    if(command_argument_count() < 1) then
        !read from pipe / stdin
        readunit = 5
    else
        !first argument is filename:
        call get_command_argument(1, fname)
        !check if file exists, readable, ...
        inquire(file=fname, exist=fexists)
        if(.not.fexists) then
            write(0, *) "File ", trim(adjustl(fname)), " does not exist."
            stop
        end if
        !open file:
        open(file=fname, unit=readunit, status='old', iostat=stat)
        inquire(file=fname, read=fcanread)
        if(fcanread /= 'YES') then
            write(0, *) "File ", trim(adjustl(fname)), " cannot be read."
            stop
        end if
        if(stat /= 0) then
            write(0, *) "File ", trim(adjustl(fname)), " cannot be opened for reading"
            stop
        end if
    end if

    !read in:
    do
        read(readunit,*,iostat=stat) tmp
        call fstats_add(tmp)
        if(stat /= 0) then
            exit
        end if
    end do
    if(readunit /= 5) then
        close(readunit)
    end if
    write(*,*) "mean:        ", fstats_mean()
    write(*,*) "naive error: ", fstats_error()
    write(*,*) "tau_int:     ", fstats_tau_int()
    write(*,*) "se(binning): ", fstats_error(ERROR_NAIVE, .true.)

    write(*,*) "binned bootstrap"
    do i=1,20
        write(*,*) i, fstats_error(ERROR_BOOTSTRAP, .true.)
    end do
    call fstats_free
end program
