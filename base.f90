Module base
    use math
    Implicit None
    Public
    Real, parameter :: omega=7.292e-5, a=6.4e+6, g=9.8

    Real :: theta_a, lambda_a, cen_lon_a, cen_lat_a, dy_a, &
            dx_a, beta0, epsilon, alpha, beta, Ca 
    Real, dimension(:), allocatable :: lon_a, lat_a, x_a, y_a
    Complex, dimension(:), allocatable :: omega_k
    Integer :: nx_a, ny_a, file_type, wave_number
    
    Real :: theta_o, lambda_o, cen_lon_o, cen_lat_o, dy_o, &
            dx_o , dt, c    
    Real, dimension(:), allocatable :: lon_o, lat_o, x_o, y_o
    Integer :: nx_o, ny_o

    Integer :: wb, eb, nb, sb


    namelist /param_atmo/ nx_a, ny_a, theta_a, lambda_a, cen_lon_a,&
                    cen_lat_a, beta0, epsilon, alpha, beta, Ca, file_type,&
                    wave_number
    namelist /param_ocean/ nx_o, ny_o, theta_o, lambda_o, cen_lon_o,&
                    cen_lat_o, dt, c

    contains
        subroutine param_init()
            Implicit None
            Integer :: i
            Real :: scale_da, scale_eps, scale_ta
            Real :: scale_do, scale_to
            
            open(20, file='namelist', form='formatted', status='old', access='sequential')
            read(20, nml=param_atmo)
            rewind(20)
            read(20, nml=param_ocean)
            close(20)
            dx_a = a*pi/180*theta_a
            dy_a = a*pi/180*lambda_a
            allocate(lon_a(nx_a), lat_a(ny_a), x_a(nx_a), y_a(ny_a))
            if(mod(nx_a, 2) .ne. 0) then
                print*, 'FFT nx should be power of 2'
                stop
            endif
            x_a = [(dx_a*i, i=-nx_a/2, nx_a/2-1)]
            lon_a = [(cen_lon_a + theta_a*i, i=-nx_a/2, nx_a/2-1)]
            if(mod(ny_a, 2) .eq. 0) then
                !y_a = [(dy_a*i, i=-ny_a/2, ny_a/2-1)]
                lat_a = [(cen_lat_a+lambda_a*i, i=-ny_a/2, ny_a/2-1)]
            else
                !y_a = [(dy_a*i, i=-ny_a/2, ny_a/2)]
                lat_a = [(cen_lat_a+lambda_a*i, i=-ny_a/2, ny_a/2)]
            endif
            y_a = lat_a*PI/180*a

            dx_o = a*pi/180*theta_o
            dy_o = a*pi/180*lambda_o
            allocate(lon_o(nx_o), lat_o(ny_o), x_o(nx_o), y_o(ny_o))
            if(mod(nx_o, 2) .eq. 0) then
                x_o = [(dx_o*i, i=-nx_o/2, nx_o/2-1)]
                lon_o = [(cen_lon_o + theta_o*i, i=-nx_o/2, nx_o/2-1)]
            else
                x_o = [(dx_o*i, i=-nx_o/2, nx_o/2)]
                lon_o = [(cen_lon_o + theta_o*i, i=-nx_o/2, nx_o/2)]
            endif
            if(mod(ny_o, 2) .eq. 0) then
                y_o = [(dy_o*i, i=-ny_o/2, ny_o/2-1)]
                lat_o = [(cen_lat_o+lambda_o*i, i=-ny_o/2, ny_o/2-1)]
            else
                y_o = [(dy_o*i, i=-ny_o/2, ny_o/2)]
                lat_o = [(cen_lat_o+lambda_o*i, i=-ny_o/2, ny_o/2)]
            endif

            do i = 1, nx_a
                if(lon_a(i) > lon_o(1)) then
                    wb = i
                    exit
                endif
            enddo
            do i = 1, nx_a
                if(lon_a(i) > lon_o(nx_o)) then
                    eb = i-1
                    exit
                endif
            enddo
            do i = 1, ny_a
                if(lat_a(i) > lat_o(1)) then
                    sb = i
                    exit
                endif
            enddo
            do i = 1, ny_a
                if(lat_a(i) > lat_o(ny_o))then 
                    nb = i-1
                    exit
                endif
            enddo
            allocate(omega_k(nx_a))
            !omega_k = [(cmplx(0, 1.0)*2.0*PI*(i - nx_a/2)/nx_a/dx_a, i=0, nx_a-1)]
            !omega_k = [(cmplx(0, 1.0)*2.0*PI*(i - nx_a/2)/nx_a/(theta_a*PI/180.0), i=0, nx_a-1)]
            omega_k = [(cmplx(0, 1.0)*2.0*PI*(i - nx_a/2)/nx_a, i=0, nx_a-1)]

            ! atmophereic None dimensional variables
            scale_da = sqrt(Ca/(2.0*beta0))
            scale_eps = sqrt(2.0*beta0*Ca)
            scale_ta = 1.0/scale_eps

            dx_a = dx_a/scale_da
            dy_a = dy_a/scale_da

            x_a = x_a/scale_da
            y_a = y_a/scale_da

            epsilon = epsilon/scale_eps
            omega_k = omega_k/dx_a


            ! oceanic None dimensional variables
            scale_do = sqrt(c/beta0)
            scale_to = 1.0/sqrt(c*beta0)

            dx_o = dx_o/scale_do
            dy_o = dy_o/scale_do

            x_o = x_o/scale_do
            y_o = y_o/scale_do

            dt = dt*365.25/12*24*3600/scale_to
            print*, dt
            

        end subroutine param_init

        subroutine atmo_init(ua, va, p)
            Implicit None
            Real, dimension(:, :), allocatable :: ua, va, p, Qs, Q1

            allocate(ua(ny_a, nx_a))
            allocate(va, p, Qs, Q1, mold=ua)

            ua = 0
            va = 0
            p = 0
        end subroutine atmo_init

        subroutine ocean_init(sst, divm, uam, vam, uom, vom, wem, TO)
            Implicit None
            Real, dimension(:, :), allocatable :: sst, divm, uam, vam, uom, vom, wem, TO

            allocate(sst(nx_o, ny_o))
            allocate(divm, uam, vam, uom, vom, wem, TO, mold=sst)

        end subroutine ocean_init
End Module base




