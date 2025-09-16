Module Math
        Implicit None
        Interface gradient
                module procedure gradient_1d
                module procedure gradient_2d
                module procedure gradient_3d
                module procedure gradient_4d
        End interface gradient
        Real, parameter :: PI=3.1415926
        contains
                Function gradient_1d(x, d) result(ans)
                        Implicit None
                        Real, optional :: d
                        Real :: dd
                        Real, dimension(:) :: x
                        Real, allocatable, dimension(:) :: ans
                        Integer :: n
                        Integer :: i

                        n = ubound(x, 1)
                        allocate(ans(n))
                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1
                        Endif
                        If(dd == 0. .or. dd < 1e-6) Then
                               ans = -9999
                        Else 
                                Do i = 1, n
                                        If (i == 1) Then
                                                ans(i) = (x(i+1) - x(i))/dd
                                        Elseif(i == n) Then
                                                ans(i) = (x(i) - x(i-1))/dd
                                        Else
                                                ans(i) = (x(i+1) - x(i-1))/(2*dd)
                                        Endif
                                Enddo
                        Endif
                End Function gradient_1d
                
                Function gradient_2d(x, axio, d) result(ans)
                        Implicit None
                        Integer  :: d1, d2
                        Integer, optional :: axio
                        Integer :: axios
                        Real, optional :: d
                        Real :: dd
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: ans
                        Integer :: i

                        d1 = ubound(x, 1)
                        d2 = ubound(x, 2)
                        allocate(ans(d1, d2))
                        If(present(axio)) Then
                                axios = axio
                        Else
                                axios = 1
                        Endif
                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1.0
                        Endif
            
                        If(axios /= 1 .and. axios /= 2) Then
                                Print*, 'axios must bs 1 or 2'
                                stop
                        Endif
                        If(axios == 1) Then
                                Do i = 1, d2
                                        ans(:, i) = gradient_1d(x(:, i), d)
                                Enddo
                        Else
                                Do i = 1, d1
                                        ans(i, :) = gradient_1d(x(i, :), d)
                                Enddo
                        Endif


                End Function gradient_2d

                Function gradient_3d(x, axio, d) result(ans)
                        Implicit None
                        Real, dimension(:, :, :) :: x
                        Real, allocatable, dimension(:, :, :) :: ans
                        Integer, optional :: axio
                        Real, optional :: d
                        Integer :: d1, d2, d3, axios
                        Real :: dd
                        Integer :: i, j

                        If(present(axio)) Then
                                axios = axio
                        Else
                                axios = 1
                        Endif

                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1.0
                        Endif

                        d1 = ubound(x, 1)
                        d2 = ubound(x, 2)
                        d3 = ubound(x, 3)
                        allocate(ans(d1, d2, d3))

                        If(axios == 1) Then
                                Do i = 1, d2
                                        do j = 1, d3
                                                ans(:, i, j) = gradient_1d(x(:, i, j), dd)
                                        Enddo
                                Enddo
                        Elseif(axios == 2) Then
                                Do i = 1, d1
                                        do j = 1, d3
                                                ans(i, :, j) = gradient_1d(x(i, :, j), dd)
                                        Enddo
                                Enddo
                        Else
                                Do i = 1, d1
                                        do j = 1, d2
                                                ans(i, j, :) = gradient_1d(x(i, j, :), dd)
                                        Enddo
                                Enddo
                        Endif
                End Function gradient_3d
                
                Function gradient_4d(x, axio, d) result(ans)
                        Implicit None
                        Real, dimension(:, :, :, :) :: x
                        Real, allocatable, dimension(:, :, :, :) :: ans
                        Integer, optional :: axio
                        Real, optional :: d
                        Integer :: i, j, k, axios
                        Integer :: d1, d2, d3, d4
                        Real :: dd

                        If(present(axio)) Then
                                axios = axio
                        Else
                                axios = 1
                        Endif

                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1.0
                        Endif

                        d1 = ubound(x, 1)
                        d2 = ubound(x, 2)
                        d3 = ubound(x, 3)
                        d4 = ubound(x, 4)
                        allocate(ans(d1, d2, d3, d4))
                        If(axios == 1) Then
                                Do i = 1, d2
                                        do j = 1, d3
                                                Do k = 1, d4
                                                        ans(:, i, j, k) = gradient_1d(x(:, i, j, k), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Elseif(axio == 2) Then
                                Do i = 1, d1
                                        do j = 1, d3
                                                Do k  =1, d4
                                                        ans(i, :, j, k) = gradient_1d(x(i, :, j, k), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Elseif(axio == 3) Then
                                Do i = 1, d1
                                        do j = 1, d2
                                                Do k = 1, d4
                                                        ans(i, j, :, k) = gradient_1d(x(i, j, :, k), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Else
                                Do i = 1, d1
                                        do j = 1, d2
                                                Do k = 1, d3
                                                        ans(i, j, k, :) = gradient_1d(x(i, j, k, :), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Endif
                End Function gradient_4d 
                
                Function cov(x, y) result(ans)
                        Implicit None
                        Real, dimension(:) :: x, y
                        Real :: ans, Ex, Ey
                        Integer :: n

                        n = ubound(x, 1)
                        If(n /= ubound(y, 1)) Then
                                Print*, 'The length of x and y must be same.'
                                stop
                        Endif
                        Ex = sum(x)/n
                        Ey = sum(y)/n
                        ans = sum((x - Ex)*(x - Ey))/(n - 1)
                End Function cov

                Function corr(x, y) result(ans)
                        Implicit None
                        Real, dimension(:) :: x, y
                        Real :: ans
                        
                        ans = cov(x, y)/sqrt(cov(x, x)*cov(y, y))
                End Function corr

                Function regression(x, y) result(ans)
                        Implicit None
                        Real, dimension(:) :: y, x
                        Real, dimension(2) :: ans
                        Real :: Ex, Ey
                        Integer :: n

                        n = ubound(x, 1)
                        Ex = sum(x)/n
                        Ey = sum(y)/n

                        ans(1) = (sum(x*y)/n - Ex*Ey)/(sum(x**2)/n - Ex**2)
                        ans(2) = Ey - ans(1)*Ex
                End Function regression

                Function legendre_(n, u) result(P)
                        Implicit None
                        Integer :: n
                        Real :: u
                        Real, dimension(0:n, 0:1) :: P

                        Integer :: i, j

                        P(0, 0) = 1
                        P(0, 1) = 0
                        P(1, 0) = sqrt(3.0)*u
                        P(1, 1) = sqrt(3.0)
                        Do i = 1, n-1
                                P(i+1, 0) = sqrt((2.0*i + 1)*(2.0*i + 3))/(i + 1)*u*P(i, 0) - &
                                          1.0*i/(i + 1)*sqrt((2.0*i + 3)/(2.0*i - 1))*P(i-1, 0)
                                P(i+1, 1) = 1.0*(i + 1)/(1 - u**2) * (sqrt((2.0*i + 3)/(2*i + 1))*P(i, 0) - u*P(i+1, 0))
                        Enddo
                        !print*, p(:, 0)
                End Function legendre_

                Function legendre(n, u) result(P)
                        Implicit None
                        Integer :: n
                        Real :: u
                        Real, dimension(0:n+1, 0:n+1) :: P

                        Integer :: i, j
                        Real :: f, g, c, d, e

                        P = 0_8
                        P(0, 0) = 1_8
                        P(1, 0) = sqrt(3.0_8)*u
                        Do i = 1, n
                                P(i+1, 0) = sqrt((2.0*i + 1)*(2.0*i + 3))/(i + 1)*u*P(i, 0) - &
                                          1.0*i/(i+1)*sqrt((2.0*i + 3)/(2.0*i - 1))*P(i-1, 0)
                        Enddo

                        
                        P(1, 1) = sqrt(3.0_8*(1 - u**2)/2)
                        Do i = 1, n
                                f = sqrt((2.0*i + 1)*(2*i + 3)/i/(i + 2))
                                g = sqrt((i - 1)*(i + 1)*(2.0*i + 3)/i/(i + 2)/(2*i - 1))
                                P(i+1, 1) = f*u*P(i, 1) - g*P(i-1, 1)
                        Enddo
                        
                        Do j = 2, n+1
                                do i = j, n+1
                                        c = sqrt((2.0*i + 1)*(i + j - 1)*(i + j - 3)/(2.0*i - 3)/(i + j)/(i + j - 2))
                                        d = sqrt((2.0*i + 1)*(i + j - 1)*(i - j + 1)/(2.0*i - 1)/(i + j)/(i + j - 2))
                                        e = sqrt((2.0*i + 1)*(i - j)/(2.0*i - 1)/(i + j))
                                        P(i, j) = c*P(i-2, j-2) - d*u*P(i-1, j-2) + e*u*P(i-1, j)
                                enddo
                        Enddo
                End Function legendre

                Function legendre_deri(P, mu, M, Jg) result(P_)
                        Implicit None
                        Integer :: M, Jg
                        Real, dimension(0:M+1, 0:M+1, Jg) :: P
                        Real, dimension(Jg) :: mu
                        Real, allocatable, dimension(:, :, :) :: P_

                        Integer :: i, j, k
                        Real :: d1, d2

                        allocate(P_(0:M+1, 0:M, Jg))
                        P_ = 0
                        Do k = 1, Jg
                                do i = 0, M
                                        Do j = i, M
                                                if(j /= 0) then
                                                        !P_(j, i, k) = (j - i + 1)/2.0*P(j-1, i+1, k) - &
                                                        !              (j - i + 1)*P(j+1, i, k)
                                                        d1 = sqrt((1.0*(j+1)**2 - i**2)/(4.0*(j+1)**2 - 1))
                                                        d2 = sqrt((1.0*j**2 - i**2)/(4.0*j**2 - 1))
                                                        P_(j, i, k) = (j*d1*P(j+1, i, k) - (j + 1)*d2*P(j-1, i, k))/(mu(k)**2-1)
                                                        
                                                 endif
                                        Enddo
                                enddo
                        Enddo
                        
                End Function legendre_deri

                Function gauss_integ(n) result(gauss)
                        Implicit None
                        Integer :: n
                        Real, dimension(2, n) :: gauss

                        Integer :: i, j
                        Real, dimension(0:n, 2) :: P
                        Real :: x, epsilon=1e-4
                        
                        Do i = 1, n
                                gauss(1, i) = sin((i - 0.5)*PI/n - PI/2)
                                do
                                        P = legendre_(n, gauss(1, i))
                                        gauss(1, i) = gauss(1, i) - P(n, 1)/P(n, 2)
                                        If(abs((P(n, 1)/P(n, 2))/(gauss(1, i))) <= epsilon .or. abs(gauss(1, i)) == 0) then
                                                exit
                                        Endif
                                enddo
                                gauss(2, i) = 2*(2.0*n - 1)*(1 - gauss(1, i)**2)/(n*P(n-1, 1))**2
                        Enddo

                End Function gauss_integ

                Function inter2d(f, x, y, x1, y1) result(f1)
                        Implicit None
                        Real, dimension(:, :) :: f 
                        Real, dimension(:) :: x, y, x1, y1
                        Real, allocatable, dimension(:, :) :: f1
                        Integer, allocatable, dimension(:) :: indm, indn
                        Integer :: m, n, m1, n1
                        Integer :: i, j, ind
                        Real :: w1, w2
                        
                        n = ubound(x, 1)
                        m = ubound(y, 1)
                        n1 = ubound(x1, 1)
                        m1 = ubound(y1, 1)
                        allocate(f1(n, m), indm(m1), indn(n1))
                        
                        ind = 1
                        Do i = 1, m1
                                do
                                        If(y1(i) > y(ind)) Then
                                                if(ind == m) then
                                                        indm(i) = ind - 1
                                                        exit
                                                elseif(y1(i) < y(ind)) then
                                                        indm(i) = ind
                                                else
                                                        ind = ind + 1
                                                        
                                                endif
                                        Else
                                                indm(i) = ind
                                                exit
                                        Endif
                                enddo
                        Enddo
                        ind = 1
                        Do i = 1, n1
                                do
                                        If(x1(i) > x(ind)) Then
                                                if(ind == n) then
                                                        indn(i) = ind - 1
                                                        exit
                                                elseif(x1(i) < x(ind)) then
                                                        indn(i) = ind
                                                else
                                                        ind = ind + 1

                                                endif
                                        Else
                                                indn(i) = ind
                                                exit
                                        Endif
                                enddo
                        Enddo
                         
                        Do i = 1, n1
                                do j = 1, m1
                                        w1 = linear_inter(x(indn(i)), x(indn(i)+1), x1(i), f(indn(i), indm(j)), f(indn(i)+1,indm(j)))
                                        w2 = linear_inter(x(indn(i)), x(indn(i)+1), x1(i), f(indn(i), indm(j)+1), f(indn(i)+1,indm(j)+1))
                                        f1(i, j) = linear_inter(y(indm(j)), y(indm(j)+1), y1(j), w1, w2)
                                enddo
                        Enddo


                        contains
                              Function linear_inter(x1, x2, x, y1, y2) result(ret)
                                        Implicit None
                                        Real :: x1, x2, y1, y2, x, ret

                                        ret = (x2 - x)/(x2 - x1)*y1 + (x - x1)/(x2 - x1)*y2
                              End Function linear_inter
                        End Function inter2d 


End Module Math


