program cz
    use base
    use Math, only : inter2d
    use atmo
    use io_mod
    Implicit None
    Real, dimension(:, :), allocatable :: ua, va, p
    Real, dimension(:, :), allocatable :: sst, divm, uam, vam, uom, vom, wem, TO
    
    Integer :: i, j
    call param_init
    
    call ocean_init(sst, divm, uam, vam, uom, vom, wem, TO)

    call atmo_init(ua, va, p)
    
    call read_bd(sst, divm, uam, vam, uom, vom, wem)

    call read_vars(TO)

    !print*, shape(sst)
    !print*, lon_a, lat_a
    call atmo_solve(ua, va, p, sst, divm, TO)


    
end
