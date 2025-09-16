Module atmo
    use base
    use FFT_routine
    use matrix
    contains
        subroutine gill_solve(ua, va, p, Qsum)
            Implicit None
            Real, dimension(ny_a, nx_a) :: ua, va, p, Qsum
            
            Complex, dimension(ny_a, nx_a) :: Um, Vm, Pm, Qm
            Complex, dimension(3*ny_a, 3*ny_a) :: W
            Complex, dimension(3*ny_a, 1) :: x, b
            Complex :: prime
            Integer :: i, j, nx, ny

            nx = nx_a
            ny = ny_a
            Qm = fftshift(fft(Qsum, axio=2), axio=2)
            !Do i = 2, 2
            Do i = 1, wave_number
                j = nx/2 + i + 1
                prime = omega_k(j)/(theta_a*PI/180.0)
                W(1:ny, 1:ny) = get_e(ny)*epsilon
                W(1:ny, ny+1:2*ny) = -beta0*diag(a*lat_a*PI/180)
                W(1:ny, 2*ny+1:3*ny) = cmplx(0, 1)*diag(aimag(prime)/(a*cos(lat_a*PI/180)))
                b(1:ny, 1) = 0

                W(ny+1:2*ny, 1:ny) =  beta0*diag(a*lat_a*PI/180)
                W(ny+1:2*ny, ny+1:2*ny) = get_e(ny)*epsilon
                W(ny+1:2*ny, 2*ny+1:3*ny) = deri(ny, lambda_a*PI/180, 1)/a
                b(ny+1:2*ny, 1) = 0

                W(2*ny+1:3*ny, 1:ny) = cmplx(0, 1)*diag(Ca**2*aimag(prime)/a*cos(lat_a*PI/180))
                W(2*ny+1:3*ny, ny+1:2*ny) = Ca**2*deri(ny, lambda_a*PI/180, 1)/a
                W(2*ny+1:3*ny, 2*ny+1:3*ny) = get_e(ny)*epsilon
                b(2*ny+1:3*ny, 1) = Qm(:, j)

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
            ua = real(ifft(ifftshift(Um, axio=2), axio=2))
            va = real(ifft(ifftshift(Vm, axio=2), axio=2))
            p  = real(ifft(ifftshift(Pm, axio=2), axio=2))
        end subroutine gill_solve

        subroutine gill_solve_no_dim(ua, va, p, Qsum)
            Implicit None
            Real, dimension(ny_a, nx_a) :: ua, va, p, Qsum

            Complex, dimension(ny_a, nx_a) :: Um, Vm, Pm, Qm
            Complex, dimension(3*ny_a, 3*ny_a) :: W
            Complex, dimension(3*ny_a, 1) :: x, b
            Complex :: prime

            Integer :: i, j, nx, ny

            nx = nx_a
            ny = ny_a
            Qm = fftshift(fft(Qsum, axio=2), axio=2)
            !print*, lat_a
            !Do i = 2, 2
            epsilon = 0.3
            Do i = 1, wave_number
                j = nx/2 + i + 1
                prime = omega_k(j)

                W(1:ny, 1:ny) = get_e(ny)*epsilon
                W(1:ny, ny+1:2*ny) = -1.0/2*diag(y_a)
                W(1:ny, 2*ny+1:3*ny) = get_e(ny)*prime
                b(1:ny, 1) = 0

                W(ny+1:2*ny, 1:ny) = 1.0/2*diag(y_a)
                W(ny+1:2*ny, ny+1:2*ny) = get_e(ny)*epsilon
                W(ny+1:2*ny, 2*ny+1:3*ny) = deri(ny, dy_a, 1)
                b(ny+1:2*ny, 1) = 0

                W(2*ny+1:3*ny, 1:ny) = get_e(ny)*prime
                W(2*ny+1:3*ny, ny+1:2*ny) = deri(ny, dy_a, 1)
                W(2*ny+1:3*ny, 2*ny+1:3*ny) = get_e(ny)*epsilon
                b(2*ny+1:3*ny, 1) = Qm(:, j)

                x = solve_c(W, b)
                Um(:, j) = x(1:ny, 1)
                Vm(:, j) = x(ny+1:2*ny, 1)
                Pm(:, j) = x(2*ny+1:3*ny, 1)
            Enddo
            !print*, 'dx', dx_a
            !print*, 'y', lat, 'eps', epsilon/scale_eps
            if(mod(nx, 2) .eq. 0) then
                Um(:, 2:nx/2) = conjg(Um(:, nx:nx/2+2:-1))
                Vm(:, 2:nx/2) = conjg(Vm(:, nx:nx/2+2:-1))
                Pm(:, 2:nx/2) = conjg(Pm(:, nx:nx/2+2:-1))
            endif
            ua = real(ifft(ifftshift(Um, axio=2), axio=2))
            va = real(ifft(ifftshift(Vm, axio=2), axio=2))
            p  = real(ifft(ifftshift(Pm, axio=2), axio=2))
        end subroutine gill_solve_no_dim


        subroutine atmo_solve(ua, va, p, sst, divm, TO)
            Implicit None
            Real, dimension(ny_a, nx_a) :: ua, va, p, sst_a, Qs, Q1, Qsum, divm_a, div_a, div_aa, TO_a
            Real, dimension(nx_a, ny_a) :: sst_w, div_w, TO_w
            Real, dimension(nx_o, ny_o) :: sst, divm, TO

            Integer :: i, j, nx, ny
            Real :: sst_ave
            Real :: lon(nx_a), lat(ny_a)

            Real :: ctlo=0.1, max_div                 ! threshold value for maximum differences in wind divergence 
            Integer :: IMAX=6, it=0

            nx = nx_a
            ny = ny_a
            sst_w = 0
            div_w = 0
            TO_w = 0

            sst_w(wb:eb, sb:nb) = inter2d(sst, lon_o, lat_o, lon_a(wb:eb), lat_a(sb:nb))
            div_w(wb:eb, sb:nb) = inter2d(divm, lon_o, lat_o, lon_a(wb:eb), lat_a(sb:nb))
            TO_w(wb:eb, sb:nb) = inter2d(TO, lon_o, lat_o, lon_a(wb:eb), lat_a(sb:nb))

            sst_a = transpose(sst_w)
            divm_a = transpose(div_w)
            TO_a = transpose(TO_w)

            lon = (lon_a-140.0)/20
            lat = lat_a/12
            do i = 1, ny
            do j = 1, nx
            if(lon(j) < 3 .and. lon(j) > -3 .and. lat(i) < 3 .and. lat(i) > -3) then
                Qs(i, j) = cos(PI*(lon(j)+0)/2/3)*cos(PI*(lat(i)-0)/2/3)
                
                !Qs(i, j) = 0.5
            endif
            if(lat(i) > 0.5 .and. lat(j) < 1.5) then
                divm_a(i, j) = -6.0
            else
                divm_a(i, j) = 2.0
            endif
            enddo
            enddo

            !Qs = 1.6*TO_a*exp(0.06*(sst_a-29.8))

            open(10, file='qq.grd', form='binary')
            open(20, file='pp.grd', form='binary')
            open(21, file='uwnd.grd', form='binary')
            open(22, file='vwnd.grd', form='binary')
            
            !divm_a = 0
            div_aa = 0
            Do 
            div_a = div_aa
            do i = 1, nx_a
                do j = 1, ny_a
                    if(divm_a(j, i) .gt. 0) then
                        if(divm_a(j, i) + div_a(j, i) .gt. 0) then
                            Qsum(j, i) = Qs(j, i)
                        else
                            Qsum(j, i) = Qs(j, i) - beta*(divm_a(j, i) + div_a(j, i))
                        endif
                    else
                        if(divm_a(j, i) + div_a(j, i) .gt. 0) then
                            Qsum(j, i) = Qs(j, i) + beta*divm_a(j, i)
                        else
                            Qsum(j, i) = Qs(j, i) - beta*div_a(j, i)
                        endif
                    endif
                enddo
            enddo

            call gill_solve_no_dim(ua, va, p, -Qsum)
            div_aa = gradient(ua, d=dx_a, axio=2) + gradient(va, d=dy_a, axio=1)
            max_div = maxval(abs(div_a(sb:nb, wb:eb) - div_aa(sb:nb, wb:eb)))
            if(max_div .lt. ctlo) then
                print*, 'Num of iterations:', it
                exit
            endif
            it = it+1
            print*, 'div:', max_div

            if(it > IMAX) then
                print*, 'div:', max_div
                print*, 'Solution cannot converge.....'
                stop
            endif
            Enddo
            write(10) ((Qsum(i, j), j=1, nx), i=1, ny)
            write(20) ((p(i, j), j=1, nx), i=1, ny)
            write(21) ((ua(i, j), j=1, nx), i=1, ny)
            write(22) ((va(i, j), j=1, nx), i=1, ny)


        end subroutine atmo_solve
End Module atmo
