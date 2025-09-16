program main
    use math
    use FFT_routine
    use matrix
    use netcdf
    Implicit none
    Integer, parameter :: nx=128, ny=64
    Real :: L=2, epsilon=0.1, lam=Pi/180, theta=PI/180
    Complex :: omega_k(nx)
    Integer :: ncid, uid
    Real :: lon(nx), lat(ny)
    Real, dimension(nx, ny) :: uwnd_1
    Real, dimension(ny, nx) :: uwnd, vwnd, P, Q
    Complex, dimension(ny, nx) :: Qm, Um, Vm, Pm
    Complex :: W(3*ny, 3*ny), b(3*ny, 1), x(3*ny, 1), prime, tes(nx)
    integer :: i, j, k
 
    lon = [(0.165*i-6, i=0, nx-1)]
    lat = [(0.16*i-5, i=0, ny-1)]
    print*, lon
    do i = 1, ny
        do j = 1, nx
            if(abs(lon(j)) < L) then
                Q(i, j) = cos(lon(j)*PI/2/L)*exp(-1.0/4*lat(i)**2)*16
            endif
        enddo
    enddo
    omega_k = [(cmplx(0, 1.0)*2.0*PI*(i - nx/2)/nx/lam, i=0, nx-1)]
    !call check(nf90_open('uwnd.mon.mean.nc', nf90_nowrite, ncid))
    !call check(nf90_inq_varid(ncid, 'uwnd', uid))
    !call check(nf90_get_var(ncid, uid, uwnd_1, start=(/1, 1, 3, 1/), count=(/nx, ny, 1, 1/)))
    !uwnd = transpose(uwnd_1)
    !dudy = matmul(deri(ny, 2.5, 1), uwind)
    open(10, file='qq.grd', form='binary')
    write(10) ((Q(i, j), j=1, nx), i=1, ny)
    Qm = fftshift(fft(Q, axio=2), axio=2)
    Do i = 1, 15
        j = nx/2 + i + 1 
        prime = omega_k(j)
        W(1:ny, 1:ny) = get_e(ny)*epsilon
        W(1:ny, ny+1:2*ny) = -1.0/2*diag(lat)
        W(1:ny, 2*ny+1:3*ny) = get_e(ny)*prime
        b(1:ny, 1) = 0

        W(ny+1:2*ny, 1:ny) = 1.0/2*diag(lat)
        W(ny+1:2*ny, ny+1:2*ny) = get_e(ny)*epsilon
        W(ny+1:2*ny, 2*ny+1:3*ny) = deri(ny, theta, 1)
        b(ny+1:2*ny, 1) = 0

        W(2*ny+1:3*ny, 1:ny) = get_e(ny)*prime
        W(2*ny+1:3*ny, ny+1:2*ny) = deri(ny, theta, 1)
        W(2*ny+1:3*ny, 2*ny+1:3*ny) = get_e(ny)*epsilon
        b(2*ny+1:3*ny, 1) = -Qm(:, j)

        !print*, b(2*ny+1:3*ny, 1)        
        x = solve_c(W, b)
        Um(:, j) = x(1:ny, 1)
        Vm(:, j) = x(ny+1:2*ny, 1)
        Pm(:, j) = x(2*ny+1:3*ny, 1)
    Enddo
    if(mod(nx, 2) .eq. 0) then
        Um(:, 2:nx/2) = conjg(Um(:, nx:nx/2+2:-1))
        Vm(:, 2:nx/2) = conjg(Vm(:, nx:nx/2+2:-1))
        Pm(:, 2:nx/2) = conjg(Pm(:, nx:nx/2+2:-1))
    endif
    uwnd = real(ifft(ifftshift(Um, axio=2), axio=2))
    vwnd = real(ifft(ifftshift(Vm, axio=2), axio=2))
    P    = real(ifft(ifftshift(Pm, axio=2), axio=2))
    open(20, file='pp.grd', form='binary')
    write(20) ((P(i, j), j=1, nx), i=1, ny)
    open(21, file='uwnd.grd', form='binary')
    write(21) ((uwnd(i, j), j=1, nx), i=1, ny)
    open(22, file='vwnd.grd', form='binary')
    write(22) ((vwnd(i, j), j=1, nx), i=1, ny)

    contains
        subroutine check(status)
            Integer :: status
            if(status /= 0)then
                print*, nf90_strerror(status)
                stop
            endif
        end subroutine check
end
