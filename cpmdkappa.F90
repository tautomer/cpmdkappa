Program CalcKappa
    use global
    implicit none

    real*8, allocatable :: kappa(:)

    call fdate(date)
    call read_conf()
    allocate(kappa(nstep))
    call proc_traj(kappa)
    call write_output(kappa)

end program

subroutine write_output(kappa)
    use global
    implicit none

    real*8, intent(in) :: kappa(ng,nstep)
    integer :: i, j
    real*8 :: avg(nstep), stddev(nstep), stepsum, stepmean
    character(len=15) :: outfmt

    open(unit=25, file="kappa.dat")
    write(25, "(a,i0)") "# number of beads: ", nb
    write(25, "(a,i0)") "# number of trajectories: ", ntraj
    write(25, "(2a)") "# program started on ", date
#if defined(_OPENMP)
    write(25, "(a,i8)") &
        "# the number of threads available: ", tmax
#else
    write(25, "(a)" ) "# OMP disabled"
#endif

    avg = kappa(1, :)
    stddev = 0
    do i = 2, ng
        do j = 1, nstep
            stepsum = kappa(i, j) - avg(j)
            stepmean = stepsum / ng
            avg(j) = avg(j) + stepmean
            stddev(j) = stddev(j) + stepmean * stepsum * (i - 1)
        end do
    end do
    stddev = sqrt(stddev) / ng 
    write(outfmt, '(a,i0,a)') '(f7.3,', ng+2, 'f11.7)'
    do i = 1, nstep
        write(25, outfmt) dt*(i-1), avg(i), stddev(i), kappa(:,i)
    end do
    call fdate(date)
    write(25, "(2a)") "# program ended on ", date
    close(25)

end subroutine
