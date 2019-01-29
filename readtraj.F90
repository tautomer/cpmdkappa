subroutine proc_traj(kappa)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use cv 
    use deriv
    implicit none
    
    abstract interface
        subroutine sub_vec(r, r1, r2)
            implicit none
            real*8, intent(in) :: r(3, 3)
            real*8, intent(out) :: r1(3), r2(3)
        end subroutine
        subroutine sub_deriv(r1, r2, dr)
            implicit none
            real*8, intent(in) :: r1(3), r2(3)
            real*8, intent(out) :: dr(3, 3)
        end subroutine
        subroutine sub_cv(r1, r2, s)
            implicit none
            real*8, intent(in) :: r1(3), r2(3)
            real*8, intent(out) :: s
        end subroutine
    end interface

    real*8, intent(out) :: kappa(nstep)
    integer :: i, j, tid, flag, summ(2, 2)
    real*8 :: v0, vsum, vsum_old, h(nstep)
    procedure(sub_vec), pointer :: getr1r2 => null()
    procedure(sub_deriv), pointer :: derivative => null()
    procedure(sub_cv), pointer :: colvar => null()

    ! assign pointers to correct subroutines
    if (tproj) then
        getr1r2 => get_xy_xz
        if (tdist) then
            derivative => d_proj_dist
            colvar => proj_dist
        else
            derivative => d_proj_diff
            colvar => proj_diff
        end if
    else
        getr1r2 => get_xy_yz
        derivative => d_diff
        colvar => diff
    end if

    vsum = 0
    vsum_old = 0
    kappa = 0
    j = 0
    summ = 0
#if defined(_OPENMP)
    !$omp parallel do &
    !$omp private(i, tid, flag, v0, h, vsum_old) &
    !$omp reduction(+: kappa, vsum, summ)
#endif
    do i = ist, ied
#if defined(_OPENMP)
        tid = omp_get_thread_num()
#else
        tid = 1
#endif
        vsum_old = vsum
        call read_init(i, v0, flag, tid, getr1r2, derivative)
        if (flag /= 0) then
            cycle
        end if
        if (v0 <= 0) then
            h(1) = 1d0
            vsum = vsum + v0
        else
            h(1) = 0d0
        end if
        call read_traj(i, v0, h, flag, tid, getr1r2, colvar)
        if (flag > 0) then
            vsum = vsum_old
            cycle
        end if
        if (v0 <= 0) then
            if (flag == 0) then
                summ(1, 1) = summ(1, 1) + 1
            else
                summ(1, 2) = summ(1, 2) + 1
            end if
        else 
            if (flag == 0) then
                summ(2, 1) = summ(2, 1) + 1
            else
                summ(2, 2) = summ(2, 2) + 1
            end if
        end if
        kappa = kappa + v0 * h      
        j = j + 1
        write(*, "(a,i0,a)") "Data from trajectory ", i, " collected"
    end do
#if defined(_OPENMP)
    !$omp end parallel do
#endif
    kappa = kappa / vsum
    write(*, "(i0,a,i0,a)") j, " out of ", ntraj, " trajectores processed"
    write(*, "(a,i0)") "number of forward trajectores without recrossing:  ", \
                       summ(1, 1)
    write(*, "(a,i0)") "number of forward trajectores with recrossing:     ", \
                       summ(1, 2)
    write(*, "(a,i0)") "number of backward trajectores without recrossing: ", \
                       summ(2, 1)
    write(*, "(a,i0)") "number of backward trajectores with recrossing:    ", \
                       summ(2, 2)
    ntraj = j

    return
end subroutine

subroutine read_init(idx, v0, flag, tid, getr1r2, derivative)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use cv 
    use deriv
    implicit none

    integer, intent(in) :: idx, tid
    integer, intent(out) :: flag
    real*8, intent(out) :: v0
    external :: getr1r2, derivative, colvar
    integer :: i, j, pu, vu, ioerr, ioerr2
    real*8 :: pr(3, 3, nb), pv(3, 3, nb), r(3, 3), v(3, 3)
    real*8 :: r1(3), r2(3), dr(3, 3), f
    character(len=*), parameter :: cvp = "/cv_pos", cvv = "/cv_vel" 
    character(len=25) :: cvpos, cvvel

    flag = 0
    write(cvpos, "(a,i0,a)") "./", idx, cvp
    write(cvvel, "(a,i0,a)") "./", idx, cvv
#if defined(_OPENMP)
    pu = 10 + tid
    vu = 10 + tmax + tid
#else
    pu = 21
    vu = 22
#endif
    ! TODO: put reading part into a universial subroutine
    open(unit=pu, file=cvpos, action='read', iostat=ioerr)
    open(unit=vu, file=cvvel, action='read', iostat=ioerr2)
    if (ioerr /= 0 .or. ioerr2 /= 0) then
        write(*, "(a,i0)") "Cannot find initial velocities or poitions file &
        &from trajectory ", idx
        flag = 1
        return
    end if

    do j = 1, nb
        do i = 1, 3
            read(pu, *, iostat=ioerr) pr(i, :, j)
            read(vu, *, iostat=ioerr2) pv(i, :, j)
            if (ioerr /= 0 .or. ioerr2 /= 0) then
                flag = 1
                write(*, "(a,i0)") "Error in initial velocities or poitions &
                &file from trajectory ", i
                close(pu)
                close(vu)
                return
            end if
        end do
    end do
    close(pu)
    close(vu)
    
    r = sum(pr, 3) * invnb
    v = sum(pv, 3) * invnb
    call getr1r2(r, r1, r2)
    call derivative(r1, r2, dr)

    v0 = 0
    do i = 1, 3
        do j = 1, 3
            v0 = v0 + dr(i, j) * v(i, j)
        end do
    end do
    
    return
end subroutine

subroutine read_traj(idx, v0, h, flag, tid, getr1r2, colvar)
    use global
    use cv 
    use deriv
    implicit none 

    integer, intent(in) :: idx, tid
    integer, intent(out) :: flag
    real*8, intent(in) :: v0
    real*8, intent(inout) :: h(nstep)
    external :: getr1r2, colvar
    integer :: i, j, k, tmp, tu, ou, ioerr
    real*8 :: pr(natom, 3, nb), r(3, 3), r1(3), r2(3), f
    character(len=*), parameter :: traj = "/TRAJECTORY"
    character(len=25) :: ftraj, fcv

    flag = 0
    write(ftraj, "(a,i0,a)") "./", idx, traj
    write(fcv, "(a,i0,a)") "./cvs/traj_", idx, '.dat'
#if defined(_OPENMP)
    tu = 10 + tid
    ou = 10 + tmax + tid
#else
    tu = 23
    ou = 24
#endif
    open(unit=tu, file=ftraj, action='read', iostat=ioerr)
    if (ioerr /= 0) then
        flag = 1
        write(*, "(a,i0)") "Cannot find TRAJECTORY file from trajectory ", i
        return
    end if
    open(unit=ou, file=fcv)
    if (h(1) == 1d0) then
        write(ou, '(a)') '#forward'
    else
        write(ou, '(a)') '#backward'
    end if
    write(ou, '(f7.3,2f11.7)') t(1), cv0 * au2a, v0 * h(1)

    do i = 2, nstep
        do k = 1, nb
            do j = 1, natom
                read(tu, *, iostat=ioerr) tmp, pr(j, :, k)
                if (ioerr /= 0) then 
                    flag = 1
                    write(*, "(a,i0)") "Error in reading TRAJECTORY file from&
                    & trajectory ", idx
                    close(tu)
                    close(ou, status='delete')
                    return
                end if
            end do
        end do
        r = sum(pr(ind, :, :), 3) * invnb
        call getr1r2(r, r1, r2)
        call colvar(r1, r2, f)
        if (f <= cv0) then
            h(i) = 1d0
        else
            h(i) = 0d0
        end if
        if (h(i-1) - h(i) /= 0) then
            write(ou, '(a)') '#recrossing over here'
            flag = -1
        end if
        write(ou,'(f7.3,2f11.7)') t(i), f * au2a, v0 * h(i)
    end do
    close(tu)
    close(ou)

    return
end subroutine
