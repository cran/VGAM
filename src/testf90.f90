


subroutine vgamf90fill9(vec, veclen, ansvec)
implicit none

integer          :: veclen
double precision :: vec(veclen), ansvec(veclen)
double precision, allocatable :: workspace1(:)

integer :: iii

allocate(workspace1(veclen))
do iii = 1, veclen
    workspace1(iii) = iii
    ansvec(iii) = vec(iii) + workspace1(iii)
end do
deallocate(workspace1)

end subroutine vgamf90fill9








subroutine vgamf90mux34(xmat, Dmat, nrowx, ncolx, symmetric, ansvec)
implicit none


integer          :: nrowx, ncolx, symmetric
double precision :: xmat(nrowx,ncolx), Dmat(ncolx,ncolx), ansvec(nrowx)

integer :: iii, jay, kay

if(ncolx .eq. 1) then
    do iii = 1, nrowx
        ansvec(iii) = Dmat(1,1) * xmat(iii, 1)**2
    end do
    return
end if

if(symmetric .eq. 1) then
    do iii = 1, nrowx
        ansvec(iii) = 0.0d0
        do jay = 1, ncolx
            ansvec(iii) = ansvec(iii) + Dmat(jay,jay) * xmat(iii, jay)**2
        end do
        if(ncolx .gt. 1) then
            do jay = 1, ncolx
                do kay = jay+1, ncolx
                    ansvec(iii) = ansvec(iii) + 2.0 * Dmat(jay,kay) * &
                                  xmat(iii, jay) * xmat(iii, kay)
                end do
            end do
        end if
    end do
else
    do iii = 1, nrowx
        ansvec(iii) = 0.0d0
        do jay = 1, ncolx
            do kay = 1, ncolx
                ansvec(iii) = ansvec(iii) + &
                              Dmat(jay,kay) * xmat(iii, jay) * xmat(iii, kay)
            end do
        end do
    end do
end if

return
end subroutine vgamf90mux34



