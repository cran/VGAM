! Test F90 subroutines; cannot be translated by ratfor
! 20080225
! Author: T. W. Yee 

!module VGAMf90  ! =============================================

subroutine VGAM_F90_fill9(vec, veclen, ansvec)
implicit none
! 20080225

integer :: veclen
double precision :: vec(veclen), ansvec(veclen)
double precision, allocatable :: workspace1(:)

! Local variables
integer :: iii

allocate(workspace1(veclen))
do iii = 1, veclen
    workspace1(iii) = iii
    ansvec(iii) = vec(iii) + workspace1(iii)
end do
deallocate(workspace1)

end subroutine VGAM_F90_fill9


!end module VGAMf90  ! =========================================


