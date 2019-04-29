module global
    implicit none
    integer :: ntraj, nstep, ist, ied, nb, natom, ind(3), tmax, ng = 1
    real*8, parameter :: au2a = 1 / 1.889725989 ! from au to angstrom
    real*8, parameter :: au2fs = 2.418884254d-2 ! from au to fs
    real*8 :: cv0, dt, invnb
    real*8, allocatable :: t(:)
    character(len=200) :: rootdir, date*30, inpnm*10
    logical :: tdiff, tproj, tdist 

    contains
    subroutine stopgm(msgp1, msgp2)
        implicit none
    
        character(len=*), intent(in) :: msgp1
        character(len=*), optional, intent(in) :: msgp2
        character(len=50) :: msg
    
        if (present(msgp2)) then
            write(msg, '(3a)') trim(msgp1), trim(msgp2), '.'
        else
            write(msg, '(2a)') trim(msgp1), '.'
        end if
        ! for some unknown reasons, ifort does not stop + arg
#if defined(__INTEL_COMPILER)
        write(*, '(a)') msg
        stop
#else
        stop msg
#endif
    end subroutine

    subroutine shuffle(a, n)
        implicit none

        integer, intent(in) :: n
        integer, intent(out) :: a(n)
        integer :: i, randpos, tmp
        real(kind=8) :: r

        do i = 1, n
            a(i) = i
        end do
        do i = n, 2, -1
            call random_number(r)
            randpos = int(r * i) + 1
            tmp = a(randpos)
            a(randpos) = a(i)
            a(i) = tmp
        end do

        end subroutine shuffle

end module global

module cv
    implicit none

    contains
    ! TODO: there must be some redundant codes.
    subroutine get_xy_xz(r, xy, xz)
        implicit none

        real*8, intent(in) :: r(3, 3)
        real*8, intent(out) :: xy(3), xz(3)

        xy = r(2, :) - r(1, :)
        xz = r(3, :) - r(1, :)
        
        return
    end subroutine 

    subroutine get_xy_yz(r, xy, yz)
        implicit none

        real*8, intent(in) :: r(3, 3)
        real*8, intent(out) :: xy(3), yz(3)

        xy = r(2, :) - r(1, :)
        yz = r(3, :) - r(2, :)
        
        return
    end subroutine 

    subroutine proj_dist(xy, xz, s)
        implicit none
    
        real*8, intent(in) :: xy(3), xz(3)
        real*8, intent(out) :: s
        real*8 :: lxz, dotxyxz
    
        lxz = norm2(xz)
        dotxyxz = dot_product(xy, xz)
        s = dotxyxz / lxz
    
    end subroutine
    
    subroutine proj_diff(xy, xz, s)
        implicit none
    
        real*8, intent(in) :: xy(3), xz(3)
        real*8, intent(out) :: s
        real*8 :: lxz, dotxyxz
    
        lxz = norm2(xz)
        dotxyxz = dot_product(xy, xz)
        s = dotxyxz / lxz * 2 - lxz
    
    end subroutine
    
    subroutine diff(xy, yz, s)
        implicit none
    
        real*8, intent(in) :: xy(3), yz(3)
        real*8, intent(out) :: s
        real*8 :: lxy, lyz
    
        lxy = norm2(xy)
        lyz = norm2(yz)
        s = lxy - lyz
    
    end subroutine
end module cv

module deriv
    implicit none

    contains
    ! here presented are all analytical derivatives
    ! TODO: consider to replace this with efficient numerical way
    subroutine d_diff(xy, yz, dr)
        implicit none
        
        real*8, intent(in) :: xy(3), yz(3)
        real*8, intent(out) :: dr(3, 3)
        real*8 :: lxy, lyz

        lxy = norm2(xy)
        lyz = norm2(yz)
        dr(1, 1) = -xy(1) / lxy
        dr(1, 2) = -xy(2) / lxy
        dr(1, 3) = -xy(3) / lxy
        dr(3, 1) = -yz(1) / lyz
        dr(3, 2) = -yz(2) / lyz
        dr(3, 3) = -yz(3) / lyz
        dr(2, 1) = -dr(1, 1) - dr(3, 1)
        dr(2, 2) = -dr(1, 2) - dr(3, 2)
        dr(2, 3) = -dr(1, 3) - dr(3, 3)

        return
    end subroutine
      
    subroutine d_proj_diff(xy, xz, dr)
        implicit none

        real*8, intent(in) :: xy(3), xz(3)
        real*8, intent(out) :: dr(3, 3)
        real*8 :: dxyxz, lxz, invr, invr2, t1, t2
    
        
        lxz = norm2(xz)
        invr = 1 / lxz
        dxyxz = dot_product(xy, xz)
        invr2 = invr * invr
        invr = invr * 2
        t1 = dxyxz * invr2 - 0.5
        t2 = t1 + 1
        dr(1, 1) = invr * (t1 * xz(1) - xy(1))
        dr(1, 2) = invr * (t1 * xz(2) - xy(2))
        dr(1, 3) = invr * (t1 * xz(3) - xy(3))
        dr(2, 1) = xz(1) * invr
        dr(2, 2) = xz(2) * invr
        dr(2, 3) = xz(3) * invr
        dr(3, 1) = invr * (xy(1) - t2 * xz(1))
        dr(3, 2) = invr * (xy(2) - t2 * xz(2))
        dr(3, 3) = invr * (xy(3) - t2 * xz(3))
    
        return
    end subroutine

    subroutine d_proj_dist(xy, xz, dr)
        implicit none

        real*8, intent(in) :: xy(3), xz(3)
        real*8, intent(out) :: dr(3, 3)
        real*8 :: dxyxz, lxz, invr, invr3
    
        lxz = norm2(xz)
        invr = 1 / lxz
        invr3 = invr ** 3
        dxyxz = dot_product(xy, xz)
        dr(1, 1) = -(xy(1) + xz(1)) * invr + dxyxz * xz(1) * invr3
        dr(1, 2) = -(xy(2) + xz(2)) * invr + dxyxz * xz(2) * invr3
        dr(1, 3) = -(xy(3) + xz(3)) * invr + dxyxz * xz(3) * invr3
        dr(2, 1) = xz(1) * invr
        dr(2, 2) = xz(2) * invr
        dr(2, 3) = xz(3) * invr
        dr(3, 1) = xy(1) * invr - dxyxz * xz(1) * invr3
        dr(3, 2) = xy(2) * invr - dxyxz * xz(2) * invr3
        dr(3, 3) = xy(3) * invr - dxyxz * xz(3) * invr3
        return
    end subroutine
end module deriv
