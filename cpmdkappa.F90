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

    real*8, intent(in) :: kappa(nstep)
    integer :: i

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
    do i = 1, nstep
        write(25, '(f7.3,f11.7)') t(i), kappa(i)
    end do
    call fdate(date)
    write(25, "(2a)") "# program ended on ", date
    close(25)

end subroutine
