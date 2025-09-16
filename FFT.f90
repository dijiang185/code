Module FFT_routine
    use math
    Interface fft
        module procedure fft_1
        module procedure fft_2
    End Interface fft
    
    Interface ifft
        module procedure ifft_1
        module procedure ifft_2
    End Interface ifft

    Interface fftshift
        module procedure fftshift_1
        module procedure fftshift_2
    End Interface fftshift

    Interface ifftshift
        module procedure ifftshift_1
        module procedure ifftshift_2
    End Interface ifftshift
    !Real, parameter :: PI=3.1415926
    contains
        Function swap(n)
            Implicit None
            Integer :: n
            Integer, dimension(0:n-1) :: swap

            Integer :: i

            swap(0) = 0
            Do i = 1, n-1
                swap(i) = swap(i/2)/2 + mod(i, 2)*(n)/2
            Enddo
        End Function swap

        Function fft_1(x)
            Implicit None
            Real, dimension(:) :: x
            Complex, allocatable, dimension(:) :: fft_1, w
            Integer, allocatable, dimension(:) :: ind

            Integer :: m, n
            Integer :: i, j

            n = ubound(x, 1)
            allocate(ind(n), fft_1(0:n-1), w(n))
            ind = swap(n) + 1

            Do i = 0, n/2-1
                m = 1
                do j = 1, n
                    w(j) = x(ind(j))
                enddo
                do while(2**m .lt. n)
                    do j = 1, n/2**m
                        w(j) = w(2*j-1) + exp(cmplx(0, -2*pi/2**m*i))*w(2*j)
                    enddo
                    m = m + 1
                enddo
                fft_1(i) = w(1) + exp(cmplx(0, -2*pi/n*i))*w(2)
                fft_1(i+n/2) = w(1) - exp(cmplx(0, -2*pi/n*i))*w(2)
            Enddo
        End Function fft_1

        Function fft_2(x, axio)
            Implicit None
            Real, dimension(:, :) :: x
            Integer, optional :: axio
            Complex, dimension(:, :), allocatable :: fft_2
            
            Integer :: axios
            Integer :: m, n, i, j
            m = ubound(x, 1)
            n = ubound(x, 2)
            allocate(fft_2(m, n))
            if(present(axio)) then
                axios = axio
            else
                axios = 1
            endif

            If(axios .eq. 1)then
                do i = 1, n
                    fft_2(:, i) = fft_1(x(:, i))
                enddo
            Else
                do i = 1, m
                    fft_2(i, :) = fft_1(x(i, :))
                enddo
            Endif
        End Function fft_2

        Function ifft_1(x)
            Implicit None
            Complex, dimension(:) :: x
            Integer :: m, n
            Complex, allocatable, dimension(:) :: ifft_1, w
            Integer, allocatable, dimension(:) ::  ind

            Integer :: i, j
            n = ubound(x, 1)
            allocate(ind(n), ifft_1(0:n-1), w(n))
            ind = swap(n) + 1

            Do i = 0, n/2-1
                m = 1
                do j = 1, n
                    w(j) = x(ind(j))
                enddo
                do while(2**m .lt. n)
                    do j = 1, n/2**m
                        w(j) = w(2*j-1) + exp(cmplx(0, 2*pi/2**m*i))*w(2*j)
                    enddo
                    m = m + 1
                enddo
                ifft_1(i) = w(1) + exp(cmplx(0, 2*pi/n*i))*w(2)
                ifft_1(i+n/2) = w(1) - exp(cmplx(0, 2*pi/n*i))*w(2)
            Enddo
            ifft_1 = ifft_1/n
        End Function ifft_1

        Function ifft_2(x, axio)
            Implicit None
            Complex, dimension(:, :) :: x
            Integer, optional :: axio
            Complex, dimension(:, :), allocatable :: ifft_2
            
            Integer :: axios
            Integer :: m, n, i, j
            m = ubound(x, 1)
            n = ubound(x, 2)
            allocate(ifft_2(m, n))
            if(present(axio)) then
                axios = axio
            else
                axios = 1
            endif

            If(axios .eq. 1)then
                do i = 1, n
                    ifft_2(:, i) = ifft_1(x(:, i))
                enddo
            Else
                do i = 1, m
                    ifft_2(i, :) = ifft_1(x(i, :))
                enddo
            Endif
        End Function ifft_2

        Function fftshift_1(x) result(ans)
            Implicit None
            Complex, dimension(:) :: x
            Complex, dimension(:), allocatable :: ans, w

            Integer :: n, mid

            n = ubound(x, 1)
            mid = n/2
            allocate(ans(n))
            if(mod(n, 2) .eq. 0) then
                ans(mid+1:n) = x(1:mid)
                ans(1:mid) = x(mid+1:n)
            else
                ans(1:mid) = x(mid+2:n)
                ans(mid+1:n) = x(1:mid+1)
            endif
        End Function fftshift_1

        Function fftshift_2(x, axio) result(ans)
            Implicit None
            Complex, dimension(:, :) :: x
            Integer, optional :: axio
            Complex, dimension(:, :), allocatable :: ans

            Integer :: axios, m, n, i
            
            m = ubound(x, 1)
            n = ubound(x, 2)
            allocate(ans(m, n))
            if(present(axio)) then
                axios = axio
            else
                axios = 2
            endif

            if(axios .eq. 1) then
                do i = 1, n
                    ans(:, i) = fftshift_1(x(:, i))
                enddo
            else
                do i = 1, m
                    ans(i, :) = fftshift_1(x(i, :))
                enddo
            endif
        End Function fftshift_2


        Function ifftshift_1(x) result(ans)
            Implicit None
            Complex, dimension(:) :: x
            Complex, dimension(:), allocatable :: ans, w

            Integer :: n, mid

            n = ubound(x, 1)
            mid = n/2
            allocate(ans(n))
            if(mod(n, 2) .eq. 0) then
                ans(mid+1:n) = x(1:mid)
                ans(1:mid) = x(mid+1:n)
            else
                ans(1:mid+1) = x(mid+1:n)
                ans(mid+2:n) = x(1:mid)
            endif
        End Function ifftshift_1

        Function ifftshift_2(x, axio) result(ans)
            Implicit None
            Complex, dimension(:, :) :: x
            Integer, optional :: axio
            Complex, dimension(:, :), allocatable :: ans

            Integer :: axios, m, n, i
            
            m = ubound(x, 1)
            n = ubound(x, 2)
            allocate(ans(m, n))
            if(present(axio)) then
                axios = axio
            else
                axios = 2
            endif

            if(axios .eq. 1) then
                do i = 1, n
                    ans(:, i) = ifftshift_1(x(:, i))
                enddo
            else
                do i = 1, m
                    ans(i, :) = ifftshift_1(x(i, :))
                enddo
            endif
        End Function ifftshift_2


End Module FFT_routine
