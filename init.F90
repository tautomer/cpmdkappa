subroutine read_conf()
    implicit none

    logical :: ex
    character(len=30) :: config

    call get_command_argument(1, config)
    open(unit=20, file=config, action="read")
    call parse_conf()

#if defined(__INTEL_COMPILER)
    inquire(directory="./cvs", exist=ex)
#else
    inquire(file="./cvs", exist=ex)
#endif
    if (ex) then
        write(*, '(a)') "Dir cv already exists. Cleaning it up."
        call system("rm -f cvs/*")
    else
        call system("mkdir -p cvs")
    end if
    call system("rm -f exit")
    call read_cpmd_inp()
end subroutine

subroutine parse_conf()
    implicit none

    integer :: line = 0, ioerr = 0, pos
    character(len=200) :: buffer, label*20

    call def_val()
    do while(ioerr == 0)
        read(20, '(a)', iostat=ioerr) buffer
        if(ioerr == 0) then
            line = line + 1
            pos = scan(buffer, ' ')
            label = buffer(1: pos)
            buffer = buffer(pos+1: )
            call read_buffer(label, buffer, line)
        end if
    end do
    close(20)
    call fix_unspecd()

end subroutine

subroutine def_val()
    use global
    implicit none

    tdiff = .false.
    tdist = .false.
    tproj = .false.
    ist = 0
    ied = 0
    grp = 1
end subroutine

subroutine read_buffer(label, buffer, line)
    use global
    implicit none

    integer :: ioerr
    integer, intent(in) :: line
    character(len=*), intent(in) :: buffer, label
    select case (label)
    ! general settings
    ! need fix twham & tbl & tdiff & tproj
    case ("diff")
        tdiff = .true.
    case ("proj")
        tproj = .true.
        select case(buffer)
        case ("diff")
            tdiff = .true.
        case ("dist")
            tdist = .true.
        case default
            call stopgm("Wrong projected collective variable type.")
        end select
    case ("cv")
        read(buffer, *, iostat=ioerr) cv0
        if (ioerr > 0) then
            call stopgm("Intial value of cv should be a real number.")
        else if (ioerr < 0) then
            call stopgm("Intial value of cv not specified.")
        end if
        cv0 = cv0 / au2a
    case ("directory")
        read(buffer, *, iostat=ioerr) ist, ied
        if (ioerr > 0) then
            call stopgm("Number of windows should be an integer.")
        else if (ioerr < 0) then
            call stopgm("Number of windows not sepcified.")
        end if
        ntraj = ied - ist + 1
        !allocate()
    case ("input")
        read(buffer, *, iostat=ioerr) inpnm
        if (ioerr > 0) then
            call stopgm("Number of windows should be an string.")
        else if (ioerr < 0) then
            call stopgm("CPMD input filename to be read not provided.")
        end if
    case ("group")
        read(buffer, *, iostat=ioerr) grp
        if (ioerr > 0) call stopgm("Number of groups should be an integer.")
    case default
        write(*, "(a, i0)") "Skipping invalid label at line ", line
    end select
    return

end subroutine

subroutine fix_unspecd()
    use global
    implicit none

    if (.not.tproj .and. .not.tdiff) call stopgm("No collective variable &
    &specified")
    if (ied < ist) call stopgm("Wrong Start and end directory indices")
    ! TODO: how to handle undefined cv0?

end subroutine

subroutine read_cpmd_inp()
#if defined(__INTEL_COMPILER)
    use ifport, only: system
#endif
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    implicit none

    integer :: i, msg, ioerr
    character(len=*), parameter :: getmol2 = "; grep -A 1 'FLUX SIDE' $inp |&
    & tail -n 1 > tmpin; m=$(grep -c INTEG $inp); if [[ $m -eq 0 ]]; then&
    & nb=1; else nb=$(sed -n '/TROT/{n;p;}' $inp); fi; geo=$(ls GEOMETRY*.xyz&
    &| head -n 1); [[ -z $geo ]] && exit 1; nat=$(wc -l $geo|cut -d' ' -f1);&
    & echo $nb $(($nat-2)) >> tmpin; dt=$(sed -n '/TIMES/{n;p;}' $inp); nst=&
    &$(sed -n '/MAXS/{n;p;}' $inp); echo $dt $nst >> tmpin"

    character(len=500) :: getmol, path*20

    call getcwd(rootdir)
    write(path, '(a,i0)') './', ist
    call chdir(path)
    write(getmol, '(3a)') 'inp=', trim(inpnm), getmol2
    msg = system(getmol)
    if (msg /= 0) call stopgm('CPMD input file not found.')
    open(unit=11, file='tmpin')
    read(11, *, iostat=ioerr) ind
    if (ioerr < 0) call stopgm('Incomplete input file')
    read(11, *, iostat=ioerr) nb, natom
    if (ioerr < 0) call stopgm('Incomplete input file')
    read(11, *, iostat=ioerr) dt, nstep
    if (ioerr < 0) call stopgm('Incomplete input file')
    close(11, status='delete')
    ! need t = 0, so one more step
    nstep = nstep + 1
    allocate(t(nstep))
    dt = dt * au2fs
    invnb = 1d0 / nb
#if defined(_OPENMP)
    tmax = omp_get_max_threads()
    !$omp parallel do &
    !$omp private(i) &
    !$omp shared(t, dt)
#else
    tmax = 1
#endif
    do i = 1, nstep
        t(i) = dt * (i - 1)
    end do
#if defined(_OPENMP)
    !$omp end parallel do
#endif
    ! TODO: shuffling all folders and split into groups
    call chdir(rootdir)
end subroutine
