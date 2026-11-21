!****************************************************************************
!  x5pxq_solver.f90
!
!  四精度（REAL(128)）五次方程复根求解核心。
!  从命令行版 x5pxq.f90 抽取，供对话框界面复用。
!****************************************************************************

module x5pxq_solver
    use, intrinsic :: iso_fortran_env, only: wp => real128
    implicit none

    private
    public :: wp, zero_thresh, solve_quintic, root_errors, fmt_wp, root_text
    public :: poly_text, fmt_coef, eval_expr, expr_report

    real(wp), parameter :: zero_thresh = 1.0e-20_wp
    real(wp), parameter :: coef_eps    = 1.0e-25_wp

contains

    ! 求 5 个复根，并返回全部根中的最大相对残差
    subroutine solve_quintic(a, rhs, roots, maxrel)
        real(wp), intent(in)     :: a(6), rhs
        complex(wp), intent(out) :: roots(5)
        real(wp), intent(out)    :: maxrel
        complex(wp) :: b(6), b0(6), temp(5)
        real(wp) :: bound
        integer  :: i, deg

        ! 归一化为首一多项式，避免 a5 过大/过小导致的上溢、下溢
        b(1) = cmplx(1.0_wp, 0.0_wp, wp)
        b(2) = cmplx(a(2) / a(1), 0.0_wp, wp)
        b(3) = cmplx(a(3) / a(1), 0.0_wp, wp)
        b(4) = cmplx(a(4) / a(1), 0.0_wp, wp)
        b(5) = cmplx(a(5) / a(1), 0.0_wp, wp)
        b(6) = cmplx((a(6) - rhs) / a(1), 0.0_wp, wp)

        b0    = b
        bound = root_bound(b, 5)
        deg   = 5

        do i=1,5
            roots(i) = newton_root(b(1:deg+1), deg, bound)
            call deflate(b(1:deg+1), deg, roots(i), temp(1:deg))
            b(1:deg) = temp(1:deg)
            deg = deg - 1
        end do

        ! 对原始五次多项式精修，消除逐次降次的精度损失
        do i=1,5
            roots(i) = polish_root(roots(i), b0, 5)
        end do

        ! 实系数多项式的复根必须共轭成对
        call symmetrize_conjugates(roots)

        maxrel = 0.0_wp
        do i=1,5
            maxrel = max(maxrel, rel_residual(a, rhs, roots(i)))
        end do
    end subroutine solve_quintic

    function poly(x, c, deg) result(res)
        complex(wp), intent(in) :: x
        complex(wp), intent(in) :: c(:)
        integer, intent(in)     :: deg
        complex(wp) :: res
        integer :: j

        res = c(1)
        do j=2, deg+1
            res = res * x + c(j)
        end do
    end function poly

    function poly_deriv(x, c, deg) result(res)
        complex(wp), intent(in) :: x
        complex(wp), intent(in) :: c(:)
        integer, intent(in)     :: deg
        complex(wp) :: res
        integer :: j

        res = c(1) * deg
        do j=2, deg
            res = res * x + c(j) * (deg - j + 1)
        end do
    end function poly_deriv

    ! 根的模长上界（Fujiwara 型），用于给出与系数尺度匹配的牛顿初值
    function root_bound(c, deg) result(bound)
        complex(wp), intent(in) :: c(:)
        integer, intent(in)     :: deg
        real(wp) :: bound, t
        integer  :: i

        bound = 0.0_wp
        do i=2, deg+1
            t = abs(real(c(i), wp))
            t = t ** (1.0_wp / real(i-1, wp))
            bound = max(bound, t)
        end do
        bound = bound + 1.0_wp
        if (bound < 1.0_wp) bound = 1.0_wp
    end function root_bound

    ! 牛顿迭代：多初值 + 收敛判定 + 残差择优
    function newton_root(c, deg, bound) result(root)
        complex(wp), intent(in) :: c(:)
        integer, intent(in)     :: deg
        real(wp), intent(in)    :: bound
        complex(wp) :: root, x, best, dx, f, df
        real(wp) :: best_res, res, r0, theta
        real(wp), parameter :: tol = 1.0e-32_wp
        real(wp), parameter :: golden_angle = 2.399963229728653_wp
        integer, parameter  :: maxit = 600, nseeds = 8
        integer :: j, k
        logical :: ok

        best     = cmplx(0.0_wp, 0.0_wp, wp)
        best_res = huge(1.0_wp)

        do k=1, nseeds
            theta = golden_angle * real(k, wp)
            r0    = bound * (0.15_wp + 0.12_wp * real(k-1, wp))
            x     = cmplx(r0 * cos(theta), r0 * sin(theta), wp)
            ok    = .false.

            do j=1, maxit
                f  = poly(x, c, deg)
                df = poly_deriv(x, c, deg)
                if (abs(df) == 0.0_wp) exit
                dx = f / df
                x  = x - dx
                if (.not. (abs(x) < huge(1.0_wp))) exit
                if (abs(x) > 4.0_wp * bound + 10.0_wp) exit
                if (abs(dx) <= tol * (1.0_wp + abs(x))) then
                    ok = .true.
                    exit
                end if
            end do

            res = abs(poly(x, c, deg))
            if (res < huge(1.0_wp)) then
                if (res < best_res) then
                    best_res = res
                    best     = x
                end if
            end if
            if (ok) exit
        end do

        root = best
    end function newton_root

    subroutine deflate(c, deg, r, c_new)
        complex(wp), intent(in)  :: c(:)
        integer, intent(in)      :: deg
        complex(wp), intent(in)  :: r
        complex(wp), intent(out) :: c_new(:)
        integer :: i

        c_new(1) = c(1)
        do i=2, deg
            c_new(i) = c(i) + r * c_new(i-1)
        end do
    end subroutine deflate

    ! 用原始多项式做牛顿精修
    function polish_root(x0, c, deg) result(x)
        complex(wp), intent(in) :: x0
        complex(wp), intent(in) :: c(:)
        integer, intent(in)     :: deg
        complex(wp) :: x, xn, dx, f, df
        real(wp), parameter :: tol = 1.0e-32_wp
        integer, parameter  :: maxit = 80
        integer :: j

        x = x0
        do j=1, maxit
            f  = poly(x, c, deg)
            df = poly_deriv(x, c, deg)
            if (abs(df) == 0.0_wp) exit
            dx = f / df
            xn = x - dx
            if (.not. (abs(xn) < huge(1.0_wp))) exit
            if (abs(poly(xn, c, deg)) > abs(f)) exit
            x = xn
            if (abs(dx) <= tol * (1.0_wp + abs(x))) exit
        end do
    end function polish_root

    subroutine symmetrize_conjugates(roots)
        complex(wp), intent(inout) :: roots(5)
        complex(wp) :: z
        real(wp) :: d, bestd
        integer  :: i, j, bestj
        logical  :: used(5)

        used = .false.
        do i=1,5
            if (used(i)) cycle
            if (abs(aimag(roots(i))) <= zero_thresh) then
                roots(i) = cmplx(real(roots(i), wp), 0.0_wp, wp)
                used(i)  = .true.
                cycle
            end if

            bestj = 0
            bestd = huge(1.0_wp)
            do j=i+1,5
                if (used(j)) cycle
                if (abs(aimag(roots(j))) <= zero_thresh) cycle
                d = abs(roots(i) - conjg(roots(j)))
                if (d < bestd) then
                    bestd = d
                    bestj = j
                end if
            end do

            if (bestj > 0 .and. bestd <= 1.0e-10_wp * (1.0_wp + abs(roots(i)))) then
                z            = (roots(i) + conjg(roots(bestj))) * 0.5_wp
                roots(i)     = z
                roots(bestj) = conjg(z)
                used(i)      = .true.
                used(bestj)  = .true.
            else
                used(i) = .true.
            end if
        end do
    end subroutine symmetrize_conjugates

    ! 根代入原方程后的实部/虚部误差
    subroutine root_errors(a, rhs, x, err_re, err_im)
        real(wp), intent(in)    :: a(6), rhs
        complex(wp), intent(in) :: x
        real(wp), intent(out)   :: err_re, err_im
        complex(wp) :: c(6), f_val

        c(1) = cmplx(a(1), 0.0_wp, wp)
        c(2) = cmplx(a(2), 0.0_wp, wp)
        c(3) = cmplx(a(3), 0.0_wp, wp)
        c(4) = cmplx(a(4), 0.0_wp, wp)
        c(5) = cmplx(a(5), 0.0_wp, wp)
        c(6) = cmplx(a(6) - rhs, 0.0_wp, wp)

        f_val  = poly(x, c, 5)
        err_re = abs(real(f_val, wp))
        err_im = abs(aimag(f_val))
    end subroutine root_errors

    ! 相对残差 |p(x)| / sum(|a_i|*|x|^i)
    function rel_residual(a, rhs, x) result(rel)
        real(wp), intent(in)    :: a(6), rhs
        complex(wp), intent(in) :: x
        real(wp) :: rel, denom, av
        complex(wp) :: c(6)

        c(1) = cmplx(a(1), 0.0_wp, wp)
        c(2) = cmplx(a(2), 0.0_wp, wp)
        c(3) = cmplx(a(3), 0.0_wp, wp)
        c(4) = cmplx(a(4), 0.0_wp, wp)
        c(5) = cmplx(a(5), 0.0_wp, wp)
        c(6) = cmplx(a(6) - rhs, 0.0_wp, wp)

        av    = abs(x)
        denom = abs(a(1))
        denom = denom * av + abs(a(2))
        denom = denom * av + abs(a(3))
        denom = denom * av + abs(a(4))
        denom = denom * av + abs(a(5))
        denom = denom * av + abs(a(6) - rhs)

        rel = abs(poly(x, c, 5)) / max(denom, tiny(1.0_wp))
    end function rel_residual

    ! 复数根的显示文本：a、a+bi、bi，按零阈值处理掉数值噪声
    function root_text(z) result(s)
        complex(wp), intent(in) :: z
        character(len=100) :: s
        real(wp) :: re, im

        re = real(z, wp)
        im = aimag(z)
        if (abs(re) < zero_thresh) re = 0.0_wp
        if (abs(im) < zero_thresh) im = 0.0_wp

        if (abs(im) < zero_thresh) then
            s = fmt_wp(re)
        else if (abs(re) < zero_thresh) then
            s = trim(fmt_wp(im)) // "i"
        else if (im > 0.0_wp) then
            s = trim(fmt_wp(re)) // " + " // trim(fmt_wp(im)) // "i"
        else
            s = trim(fmt_wp(re)) // " - " // trim(fmt_wp(abs(im))) // "i"
        end if
    end function root_text

    ! 实数格式化：常规量程用定点 20 位小数，其余用科学计数（24 位有效数字）
    ! 注意：必须去掉格式串产生的多余空格，否则拼成复数时会插入大段空白
    function fmt_wp(v) result(s)
        real(wp), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=48) :: t

        if (abs(v) >= 1.0e-4_wp .and. abs(v) < 1.0e3_wp) then
            write(t,'(f26.20)') v
        else
            write(t,'(es36.24e3)') v
        end if
        s = trim(adjustl(t))
    end function fmt_wp

    ! ------------------------------------------------------------------
    !  系数格式：整数样值直接写整数，其余保留四精度数字并去掉尾随 0
    ! ------------------------------------------------------------------
    function fmt_coef(v) result(s)
        real(wp), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=48) :: t
        integer :: k

        if (abs(v) < coef_eps) then
            s = "0"
            return
        end if
        if (abs(v) < 1.0e15_wp) then
            if (abs(v - real(nint(v), wp)) <= 1.0e-24_wp * max(1.0_wp, abs(v))) then
                write(t,'(i0)') nint(v)
                s = trim(t)
                return
            end if
        end if

        s = fmt_wp(v)
        if (index(s,'E') == 0 .and. index(s,'e') == 0) then
            if (index(s,'.') > 0) then
                k = len_trim(s)
                do while (k > 1)
                    if (s(k:k) /= '0') exit
                    k = k - 1
                end do
                if (s(k:k) == '.') k = k - 1
                s = s(1:k)
            end if
        end if
    end function fmt_coef

    ! ------------------------------------------------------------------
    !  把实际填入的系数写成方程文本（常量项写成“(A0 - RHS)”）
    !  例：10000*x^5 - 15*x^4 + 85*x^3 - 225*x^2 + 274*x + (-120 - 0) = 0
    ! ------------------------------------------------------------------
    function poly_text(a, rhs) result(s)
        real(wp), intent(in) :: a(6), rhs
        character(len=:), allocatable :: s, term
        integer :: i, n

        s = monomial_term(a(1), 5)
        do i = 2, 5
            n = 6 - i                       ! a(2)->x^4 … a(5)->x
            if (abs(a(i)) < coef_eps) cycle
            if (a(i) > 0.0_wp) then
                s = s // " + " // monomial_term(a(i), n)
            else
                s = s // " - " // monomial_term(abs(a(i)), n)
            end if
        end do
        if (abs(a(6)) >= coef_eps .or. abs(rhs) >= coef_eps) then
            s = s // " + (" // fmt_coef(a(6)) // " - " // fmt_coef(rhs) // ")"
        end if
        s = s // " = 0"
    end function poly_text

    ! 单项式文本：系数 1 / -1 时省略系数，次数 1 时写作 x
    function monomial_term(c, n) result(t)
        real(wp), intent(in) :: c
        integer,  intent(in) :: n
        character(len=:), allocatable :: t, pname

        if (n <= 1) then
            pname = "x"
        else
            pname = "x^" // char(48+n)
        end if

        if (abs(abs(c) - 1.0_wp) < coef_eps) then
            t = pname
        else
            t = fmt_coef(c) // "*" // pname
        end if
    end function monomial_term

    ! ------------------------------------------------------------------
    !  精确代数式求值（四精度）：
    !    + - * / ** 、括号、一元正负号
    !    常数 pi、e；函数 sqrt abs exp log log10 sin cos tan asin acos atan sinh cosh tanh
    ! ------------------------------------------------------------------
    subroutine eval_expr(str, v, ok)
        character(*), intent(in)  :: str
        real(wp),     intent(out) :: v
        logical,      intent(out) :: ok
        integer :: pos
        logical :: good

        pos  = 1
        good = .true.
        v    = parse_sum()
        call skip_blanks()
        if (pos <= len_trim(str)) good = .false.
        if (.not. good) v = 0.0_wp
        ok = good

    contains

        subroutine skip_blanks()
            do while (pos <= len(str))
                if (str(pos:pos) /= ' ' .and. str(pos:pos) /= char(9)) exit
                pos = pos + 1
            end do
        end subroutine skip_blanks

        logical function is_digit(c)
            character(1), intent(in) :: c
            is_digit = (c >= '0' .and. c <= '9')
        end function is_digit

        recursive function parse_sum() result(r)
            real(wp) :: r, t
            character(1) :: c
            r = parse_prod()
            do
                call skip_blanks()
                if (pos > len(str)) exit
                c = str(pos:pos)
                if (c /= '+' .and. c /= '-') exit
                pos = pos + 1
                t = parse_prod()
                if (c == '+') then
                    r = r + t
                else
                    r = r - t
                end if
            end do
        end function parse_sum

        recursive function parse_prod() result(r)
            real(wp) :: r, t
            character(1) :: c
            r = parse_power()
            do
                call skip_blanks()
                if (pos > len(str)) exit
                c = str(pos:pos)
                if (c == '*') then
                    if (pos + 1 <= len(str)) then
                        if (str(pos+1:pos+1) == '*') exit      ! 交给幂运算处理
                    end if
                    pos = pos + 1
                    r = r * parse_power()
                else if (c == '/') then
                    pos = pos + 1
                    t = parse_power()
                    if (t == 0.0_wp) then
                        good = .false.
                        r = 0.0_wp
                    else
                        r = r / t
                    end if
                else
                    exit
                end if
            end do
        end function parse_prod

        recursive function parse_power() result(r)
            real(wp) :: r, e
            character(2) :: cc
            r = parse_unary()
            call skip_blanks()
            if (pos + 1 <= len(str)) then
                cc = str(pos:pos+1)
                if (cc == '**') then
                    pos = pos + 2
                    e = parse_power()
                    if (r < 0.0_wp .and. abs(e - real(nint(e), wp)) > 1.0e-15_wp) then
                        good = .false.
                        r = 0.0_wp
                    else
                        r = r**e
                    end if
                end if
            end if
        end function parse_power

        recursive function parse_unary() result(r)
            character(1) :: c
            real(wp) :: r
            call skip_blanks()
            r = 0.0_wp
            if (pos > len(str)) then
                good = .false.
                return
            end if
            c = str(pos:pos)
            if (c == '+') then
                pos = pos + 1
                r = parse_unary()
            else if (c == '-') then
                pos = pos + 1
                r = -parse_unary()
            else
                r = parse_atom()
            end if
        end function parse_unary

        recursive function parse_atom() result(r)
            real(wp) :: r, arg
            character(24) :: name
            character(1) :: c
            integer :: p1, p2, ios

            call skip_blanks()
            r = 0.0_wp
            if (pos > len(str)) then
                good = .false.
                return
            end if

            c = str(pos:pos)

            if (c == '(') then
                pos = pos + 1
                r = parse_sum()
                call skip_blanks()
                if (pos <= len(str)) then
                    if (str(pos:pos) == ')') then
                        pos = pos + 1
                    else
                        good = .false.
                    end if
                else
                    good = .false.
                end if
                return
            end if

            ! 数值：数字串 + 可选小数点 + 可选指数
            if (is_digit(c) .or. c == '.') then
                p1 = pos
                do while (pos <= len(str))
                    if (is_digit(str(pos:pos)) .or. str(pos:pos) == '.') then
                        pos = pos + 1
                    else
                        exit
                    end if
                end do
                if (pos <= len(str)) then
                    if (str(pos:pos) == 'e' .or. str(pos:pos) == 'E') then
                        p2 = pos
                        pos = pos + 1
                        if (pos <= len(str)) then
                            if (str(pos:pos) == '+' .or. str(pos:pos) == '-') pos = pos + 1
                        end if
                        do while (pos <= len(str))
                            if (is_digit(str(pos:pos))) then
                                pos = pos + 1
                            else
                                exit
                            end if
                        end do
                        if (pos == p2 + 1) good = .false.   ! “1e” 这种残缺指数
                    end if
                end if
                read(str(p1:pos-1), *, iostat=ios) r
                if (ios /= 0) good = .false.
                return
            end if

            ! 名称：常数或函数
            if ((c >= 'a' .and. c <= 'z') .or. (c >= 'A' .and. c <= 'Z')) then
                p1 = pos
                do while (pos <= len(str))
                    c = str(pos:pos)
                    if ((c >= 'a' .and. c <= 'z') .or. (c >= 'A' .and. c <= 'Z') .or. is_digit(c)) then
                        pos = pos + 1
                    else
                        exit
                    end if
                end do
                name = str(p1:pos-1)
                call to_lower(name)

                call skip_blanks()
                if (pos <= len(str)) then
                    if (str(pos:pos) == '(') then
                        pos = pos + 1
                        arg = parse_sum()
                        call skip_blanks()
                        if (pos <= len(str)) then
                            if (str(pos:pos) == ')') then
                                pos = pos + 1
                            else
                                good = .false.
                            end if
                        else
                            good = .false.
                        end if
                        r = apply_func(name, arg, good)
                        return
                    end if
                end if

                select case (trim(name))
                case ("pi")
                    r = acos(-1.0_wp)
                case ("e")
                    r = exp(1.0_wp)
                case default
                    good = .false.
                    r = 0.0_wp
                end select
                return
            end if

            good = .false.
        end function parse_atom

        subroutine to_lower(s)
            character(*), intent(inout) :: s
            integer :: k, ic
            do k = 1, len(s)
                ic = ichar(s(k:k))
                if (ic >= 65 .and. ic <= 90) s(k:k) = char(ic + 32)
            end do
        end subroutine to_lower

        recursive function apply_func(fname, x, okf) result(r)
            character(*), intent(in) :: fname
            real(wp), intent(in)     :: x
            logical, intent(inout)   :: okf
            real(wp) :: r

            select case (trim(fname))
            case ("sqrt");  r = sqrt(x)
            case ("abs");   r = abs(x)
            case ("exp");   r = exp(x)
            case ("log");   r = log(x)
            case ("log10"); r = log10(x)
            case ("sin");   r = sin(x)
            case ("cos");   r = cos(x)
            case ("tan");   r = tan(x)
            case ("asin");  r = asin(x)
            case ("acos");  r = acos(x)
            case ("atan");  r = atan(x)
            case ("sinh");  r = sinh(x)
            case ("cosh");  r = cosh(x)
            case ("tanh");  r = tanh(x)
            case default
                okf = .false.
                r = 0.0_wp
            end select
        end function apply_func

    end subroutine eval_expr

    ! ------------------------------------------------------------------
    !  把表达式求值结果整理成对话框里显示的文本
    ! ------------------------------------------------------------------
    function expr_report(str, v, ok) result(s)
        character(*), intent(in) :: str
        real(wp),     intent(in) :: v
        logical,      intent(in) :: ok
        character(len=900) :: s
        character(len=80)  :: line
        character(2)       :: nl

        nl = achar(13) // achar(10)

        if (.not. ok) then
            s = "错误：无法解析表达式 """ // trim(str) // """。" // nl // &
                "支持：+ - * / **、括号、( )；常数 pi、e；" // nl // &
                "函数 sqrt abs exp log log10 sin cos tan asin acos atan sinh cosh tanh。" // nl // &
                "示例：1/3、sqrt(2)、pi/4、(1+sqrt(5))/2、2**0.5"
            return
        end if

        s = "表达式：" // trim(str) // nl // nl

        if (abs(v) >= 1.0e-4_wp .and. abs(v) < 1.0e10_wp) then
            write(line,'(f42.34)') v
            s = trim(s) // "定点表示：" // nl // "  " // trim(adjustl(line)) // nl
        end if

        write(line,'(es44.34e4)') v
        s = trim(s) // "科学计数法：" // nl // "  " // trim(adjustl(line)) // nl // nl
        s = trim(s) // "说明：REAL(128) 四精度的机器精度约 1.93E-34，约 33–34 位十进制有效数字。" // nl // &
                      "上面的数字就是该代数式在四精度下的最佳二进制近似。"
    end function expr_report

end module x5pxq_solver
