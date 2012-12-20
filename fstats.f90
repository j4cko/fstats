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

    integer :: readunit = 5
    real :: summed=0.0, tmp,testarr(100)
    integer :: stat, i ,j

   ! do
   !     read(readunit,*,iostat=stat) tmp
   !     if(stat /= 0) then
   !         exit
   !     end if
   !     call add_real(tmp)
   ! end do
    do i=0,9
        do j=1,100
            testarr(j)=real(i*100+j)
        end do
        call fstats_add(testarr)
    end do
    call fstats_add(19203.234)
    
    write(*,*) "mean:        ", mean()
    write(*,*) "naive error: ", se_naive()
    !write(*,*) "tau_int:     ", tau_int()
    !write(*,*) "se(binning): ", se_binning(i)
    !write(*,*) "        (after ", i, " steps, i.e. binsize=", 2**i, ")"

    write(*,*) "binned bootstrap"
    do i=1,20
        write(*,*) i, error_binning(i,.false.,ERROR_BOOTSTRAP)
    end do
    call fstats_free
end program
