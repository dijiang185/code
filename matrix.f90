module matrix

        contains
                Subroutine write_matrix(x)
                        Implicit None
                        Complex, dimension(:, :) :: x
                        Integer :: i, n

                        n = ubound(x, 1)
                        print*, 'print matrix'
                        do i = 1, n
                                print*, x(i, :)
                        enddo
                End Subroutine write_matrix
                Function sign(x)
                        Implicit None
                        Real :: x
                        Real :: sign
                        if(x .lt. 0) then
                                sign = -1.0
                        elseif(x .eq. 0) then 
                                sign = 0
                        else
                                sign = 1.0
                        endif
                End Function sign
                Function get_e(n) result(E)
                        Implicit None
                        Integer :: n
                        Real :: E(n, n)

                        Integer :: i
                        E = 0
                        Do i = 1, n
                                E(i, i) = 1
                        Enddo
                End Function get_e
                
                Function diag(x)
                    Implicit None
                    Real, dimension(:) :: x
                    Real, allocatable, dimension(:, :) ::  diag

                    Integer :: n, i

                    n = ubound(x, 1)
                    allocate(diag(n, n))
                    diag = 0
                    do i = 1, n
                        diag(i, i) = x(i)
                    enddo
                End Function diag

                Function deri(n, d, axio)
                    Implicit None
                    Integer :: n, axio
                    Real :: d
                    Real, dimension(n, n) :: deri
                    
                    Integer :: i
                    deri = 0
                    deri(1, 1:2) = [-1.0, 1.0]/d
                    deri(n, n-1:n) = [-1.0, 1.0]/d
                    do i = 2, n-1
                        deri(i, i-1:i+1) = [-1.0, 0.0, 1.0]/2/d
                    enddo
                    if(axio .eq. 2) deri = transpose(deri)
                End Function deri
                
                Function eigen_main(x) result(ans)
                        Implicit None
                        Real, dimension(:, :) :: x

                        Real, allocatable, dimension(:) :: ans
                        Integer :: n, i, j
                        Real, allocatable, dimension(:) :: u, v

                        n = ubound(x, 1)
                        allocate(ans(n+1), u(n), v(n))
                        u = 1
                        v = 1
                        Do i = 1, 40
                                v = matmul(x, u)
                                u = v/maxval(v)
                                print*, i, u, maxval(v)
                        Enddo
                        ans(1) = maxval(v)
                        ans(2:n+1) = u
                End Function eigen_main

                Subroutine eigen_sym(x, eigenvalue, eigenvector)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: eigenvector
                        Real, allocatable, dimension(:) :: eigenvalue
                        
                        Real, allocatable, dimension(:, :, :) :: w
                        Real, allocatable, dimension(:, :) :: U
                        Integer :: n, i, j, p, q, k = 0
                        Real :: maxvalue, o, t, c, s
                        Real :: erro=1e-8

                        n = ubound(x, 1)
                        allocate(eigenvalue(n), eigenvector(n, n), w(n, n, 2))
                        eigenvector = get_e(n)
                        w(:, :, 1) = x
                        Do
                        k = k+1
                        maxvalue = abs(w(1, 2, 1))
                        p = 1
                        q = 2
                        do i = 1, n
                                do j = i+1, n
                                if(abs(w(i, j, 1)) .gt. maxvalue) then
                                        p = i
                                        q = j
                                        maxvalue = abs(w(i, j, 1))
                                endif
                                enddo
                        enddo
                        if(maxvalue .le. erro) exit
                        o = (w(p, p, 1) - w(q, q, 1))/(2*w(p, q, 1))
                        t = sign(o)*(sqrt(o**2 + 1) - abs(o))
                        c = 1.0/sqrt(1 + t**2)
                        s = t*c
                        U = get_e(n)
                        U(p, p) = c
                        U(q, q) = c
                        U(p, q) = -s
                        U(q, p) = s

                        w(p, q, 2) = 0
                        w(q, p, 2) = 0
                        w(p, p, 2) = c**2*w(p, p, 1) + s**2*w(q, q, 1) + 2*c*s*w(p, q, 1)
                        w(q, q, 2) = s**2*w(p, p, 1) + c**2*w(q, q, 1) - 2*c*s*w(p, q, 1)


                        do i = 1, n
                                do j = 1, n
                                if(i .ne. p .and. i .ne. q .and. j .ne. p .and. j .ne. q) then
                                        w(i, j, 2) = w(i, j, 1)
                                endif
                                if(i .eq. p .and. j .ne. p .and. j .ne. q) then
                                        w(i, j, 2) = c*w(p, j, 1) + s*w(q, j, 1)
                                endif
                                if(i .eq. q .and. j .ne. p .and. j .ne. q) then
                                        w(i, j, 2) = -s*w(p, j, 1) + c*w(q, j, 1)
                                endif
                                if(j .eq. p .and. i .ne. p .and. i .ne. q) then
                                        w(i, j, 2) = c*w(p, i, 1) + s*w(q, i, 1)
                                endif
                                if(j .eq. q .and. i .ne. p .and. i .ne. q) then
                                        w(i, j, 2) = -s*w(p, i, 1) + c*w(q, i, 1)
                                endif
                                enddo
                        enddo 
                        w(:, :, 1) = w(:, :, 2)
                        eigenvector = matmul(eigenvector, U)
                        Enddo
                        do i = 1, n
                                eigenvalue(i) = w(i, i, 1)
                        enddo
                End Subroutine eigen_sym

                Function Householder(x) result(H)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: u
                        Real, allocatable, dimension(:, :) :: H, x_1
                        Real :: beta, k
                        Integer :: n

                        n = ubound(x, 1)
                        allocate(H(n, n), u(n, 1), x_1(n, 1))
                        x_1 = x/maxval(abs(x))
                        k = sign(x_1(1, 1))*sqrt(sum(x_1**2))
                        u = x_1
                        u(1, 1) = u(1, 1) + k
                        beta = k*(k + x_1(1, 1))
                        H = get_e(n) - matmul(u, transpose(u))/beta
                End Function Householder

                Subroutine QR_resolve(x, Q, R)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: Q, R, H, D

                        Integer :: n, i, j

                        n = ubound(x, 1)
                        allocate(Q(n, n), R(n, n), H(n, n), D(n, n))
                        Q = get_e(n)
                        R = x
                        do i = 1, n-1
                                H = get_e(n)
                                H(i:n, i:n) = Householder(reshape(R(i:n, i), (/n-i+1, 1/)))
                                R = matmul(H, R)
                                Q = matmul(H, Q)
                        enddo
                        D = get_e(n)
                        do i = 1, n
                                D(i, i) = sign(R(i, i))
                        enddo
                        Q = matmul(transpose(Q), D)
                        R = matmul(D, R)
                End Subroutine QR_resolve

                Subroutine Heisenberg(x, H, U)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: H, U, u_1

                        Integer :: n, i, j

                        n = ubound(x, 1)
                        allocate(H(n, n), U(n, n), u_1(n, n))
                        H = x
                        U = get_e(n)
                        do i = 1, n-2
                                u_1 = get_e(n)
                                u_1(i+1:n, i+1:n) = Householder(reshape(H(i+1:n, i), (/n-i, 1/)))
                                H = matmul(matmul(u_1, H), u_1)
                                U = matmul(u_1, U)
                        enddo
                End Subroutine Heisenberg

                Subroutine QR(x, eig_val, eig_vec)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: eig_vec, Q, R, A
                        Real, allocatable, dimension(:) :: eig_val

                        Integer :: n, i, j
                        Real :: s, s1, erro=1e-6

                        n = ubound(x, 1)
                        allocate(eig_val(n), eig_vec(n, n), A(n, n))
                        A = x
                        eig_vec = get_e(n)
                        s1 = 9999
                        Do 
                        call QR_resolve(A, Q, R)
                        A = matmul(R, Q)
                        eig_vec = matmul(eig_vec, Q)
                        s = sum(abs(A))
                        do i = 1, n
                                s = s - abs(A(i, i))
                        enddo
                        if(s .lt. erro .or. abs(s1-s) .lt. erro) exit
                        s1 = s
                        deallocate(Q, R)
                        Enddo
                        !eig_vec = transpose(eig_vec)
                        do i = 1, n
                        eig_val(i) = A(i, i)
                        enddo
                End Subroutine QR
                Function GE(x)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: GE, xx
                        Real, allocatable, dimension(:) :: w
                        Integer :: n, i, j, ind
                        Real :: max_val

                        n = ubound(x, 1)
                        allocate(GE(n, n), xx(n, 2*n), w(2*n))
                        xx(:, 1:n) = x
                        xx(:, n+1:2*n) = get_e(n)
                        Do i = 1, n-1
                                max_val = abs(xx(i, i))
                                ind = i
                                do j = i+1, n
                                if(abs(xx(j, i)) .gt. max_val) then
                                        max_val = abs(xx(j, i))
                                        ind = j
                                endif
                                enddo
                                if(max_val .lt. 1e-6) stop
                                if(ind .ne. i) then
                                        w = xx(i, :)
                                        xx(i, :) = xx(ind, :)
                                        xx(ind, :) = w
                                endif
                                xx(i, :) = xx(i, :)/xx(i, i)
                                do j = i+1, n
                                        xx(j, :) = xx(j, :) - xx(i, :)*xx(j, i)
                                enddo
                        Enddo
                        Do i = n, 2, -1
                                max_val = abs(xx(i, i))
                                ind = i
                                do j = i-1, 1, -1
                                if(abs(xx(j, i)) .gt. max_val) then
                                        max_val = abs(xx(j, i))
                                        ind = j
                                endif
                                enddo
                                if(ind .ne. i) then
                                        w = xx(i, :)
                                        xx(i, :) = xx(ind, :)
                                        xx(ind, :) = w
                                endif
                                xx(i, :) = xx(i, :)/xx(i, i)
                                do j = i-1, 1, -1
                                        xx(j, :) = xx(j, :) - xx(i, :)*xx(j, i)
                                enddo
                        Enddo
                        GE = xx(:, n+1:2*n)
                End Function GE
                
                Function solve_c(A, b) result(x)
                    Implicit None
                    Complex, dimension(:, :) :: A, b
                    Complex, dimension(:, :), allocatable :: x, xx

                    Integer :: n, i, j
                    Integer, dimension(1) :: ind
                    Complex, allocatable, dimension(:) :: w
                    Real :: max_val
                    
                    n = ubound(A, 1)
                    if(n .ne. ubound(b, 1)) then
                        print*, 'shape of A', n, ubound(A, 2)
                        print*, 'shape of b', ubound(b, 1), ubound(A, 2)
                        stop
                    endif
                    allocate(x(n, 1), xx(n, n+1), w(n+1))
                    xx(:, 1:n) = A
                    xx(:, n+1) = b(:, 1)
                    do i = 1, n-1
                        ind = maxloc(abs(xx(i:n, i))) + i - 1
                        w = xx(i, :)
                        xx(i, :) = xx(ind(1), :)
                        xx(ind(1), :) = w

                        xx(i, i:n+1) = xx(i, i:n+1)*conjg(xx(i, i))/abs(xx(i, i))**2
                        do j = i+1, n
                            xx(j, i:n+1) = xx(j, i:n+1) - xx(i, i:n+1)*xx(j, i)
                        enddo
                    !call write_matrix(xx(:, 1:n))
                    enddo
                    do i = n, 2, -1
                        xx(i, i:n+1) = xx(i, i:n+1)*conjg(xx(i, i))/abs(xx(i, i))**2
                        do j = i-1, 1, -1
                            xx(j, i:n+1) = xx(j, i:n+1) - xx(i, i:n+1)*xx(j, i)
                        enddo
                    enddo
                    !call write_matrix(xx(:, 1:n))
                    x(:, 1) = xx(:, n+1)

                End Function solve_c
                    

End module matrix
