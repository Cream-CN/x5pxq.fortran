program quintic_complex_solver
    use, intrinsic :: iso_fortran_env, only: wp => real128
    implicit none
    real(wp)    :: a(6), rhs
    complex(wp) :: roots(5)
    integer     :: i, ierr
    real(wp)    :: maxrel, err_re, err_im
    real(wp), parameter :: zero_thresh = 1.0e-20_wp  ! 零根判断阈值（与目标精度一致）
    real(wp), parameter :: err_limit = 1.0e-20_wp    ! 误差阈值（真正参与判定）

    write(*,*) "=== High Precision Quintic Equation Solver (10^-20) ==="
    write(*,*) "Equation: A5*x^5 + A4*x^4 + A3*x^3 + A2*x^2 + A1*x + A0 = RHS"
    write(*,*) "Input 7 parameters: A5 A4 A3 A2 A1 A0 RHS"
    read(*,*,iostat=ierr) a(1),a(2),a(3),a(4),a(5),a(6),rhs
    if (ierr /= 0) error stop "Invalid input! Please enter 7 numbers."
    if (a(1) == 0.0_wp) error stop "Invalid input! A5 must be non-zero (a quintic is required)."

    call find_all_roots(a, rhs, roots)

    write(*,*) "=========================================================="
    write(*,*) "Results (Format: a+bi | Real/Imaginary Errors)"
    write(*,*) "=========================================================="

    maxrel = 0.0_wp
    do i=1,5
        call print_root_format(roots(i), i, zero_thresh)
        call calc_separate_error(a, rhs, roots(i), err_re, err_im)
        call print_separate_error(roots(i), err_re, err_im, zero_thresh)
        maxrel = max(maxrel, rel_residual(a, rhs, roots(i)))
        write(*,*)
    end do

    write(*,*) "=========================================================="
    write(*,'(a,es12.3,a,es10.2,a)') " Max relative residual = ", maxrel, "  (limit ", err_limit, ")"
    if (maxrel <= err_limit) then
        write(*,*) " Precision check: PASS"
    else
        write(*,*) " Precision check: FAIL (relative residual exceeds 1e-20)"
    end if

contains

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

! 根的模长上界：Fujiwara 型界，用于给出与系数尺度匹配的牛顿初值
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

! 单个根的牛顿迭代：多初值 + 收敛判定 + 残差择优，避免发散产生垃圾根
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

! 用原始多项式对根做牛顿精修，消除逐次降次带来的精度损失
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
        if (abs(poly(xn, c, deg)) > abs(f)) exit   ! 残差不再下降就停止，避免震荡
        x = xn
        if (abs(dx) <= tol * (1.0_wp + abs(x))) exit
    end do
end function polish_root

subroutine find_all_roots(a, rhs, roots)
    real(wp), intent(in)     :: a(6), rhs
    complex(wp), intent(out) :: roots(5)
    complex(wp) :: b(6), b0(6), temp(5)
    real(wp) :: bound
    integer  :: i, deg

    ! 归一化为首一多项式：避免 a5 很小或很大时的上溢/下溢
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

    ! 精修：对原始五次多项式再收敛一次
    do i=1,5
        roots(i) = polish_root(roots(i), b0, 5)
    end do

    ! 实系数多项式的复根必须共轭成对，抑制数值噪声造成的不对称
    call symmetrize_conjugates(roots)
end subroutine find_all_roots

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
            z        = (roots(i) + conjg(roots(bestj))) * 0.5_wp
            roots(i)      = z
            roots(bestj)  = conjg(z)
            used(i)  = .true.
            used(bestj) = .true.
        else
            used(i) = .true.
        end if
    end do
end subroutine symmetrize_conjugates

! 实数输出：常规量程用定点 20 位小数，超量程自动切科学计数，杜绝 *****
function fmt_real(v) result(s)
    real(wp), intent(in) :: v
    character(len=40) :: s
    character(len=40) :: t

    if (abs(v) >= 1.0e-4_wp .and. abs(v) < 1.0e3_wp) then
        write(t,'(f26.20)') v
    else
        write(t,'(es36.24e3)') v
    end if
    s = adjustl(t)
end function fmt_real

subroutine print_root_format(x, idx, zero_tol)
    complex(wp), intent(in) :: x
    integer, intent(in)     :: idx
    real(wp), intent(in)    :: zero_tol
    real(wp) :: re, im

    re = real(x, wp)
    im = aimag(x)
    if (abs(re) < zero_tol) re = 0.0_wp
    if (abs(im) < zero_tol) im = 0.0_wp

    write(*,'(a,i2,a)',advance='no') "Root ",idx,": "

    if (abs(im) < zero_tol) then
        write(*,'(a)') trim(fmt_real(re))
    else if (abs(re) < zero_tol) then
        write(*,'(a,a)') trim(fmt_real(im)), "i"
    else if (im > 0.0_wp) then
        write(*,'(a,a,a,a)') trim(fmt_real(re)), " + ", trim(fmt_real(im)), "i"
    else
        write(*,'(a,a,a,a)') trim(fmt_real(re)), " - ", trim(fmt_real(abs(im))), "i"
    end if
end subroutine print_root_format

subroutine calc_separate_error(a, rhs, x, err_re, err_im)
    real(wp), intent(in)    :: a(6), rhs
    complex(wp), intent(in) :: x
    real(wp), intent(out)   :: err_re, err_im
    complex(wp)             :: c(6), f_val

    c(1) = cmplx(a(1), 0.0_wp, wp)
    c(2) = cmplx(a(2), 0.0_wp, wp)
    c(3) = cmplx(a(3), 0.0_wp, wp)
    c(4) = cmplx(a(4), 0.0_wp, wp)
    c(5) = cmplx(a(5), 0.0_wp, wp)
    c(6) = cmplx(a(6) - rhs, 0.0_wp, wp)

    f_val  = poly(x, c, 5)
    err_re = abs(real(f_val, wp))
    err_im = abs(aimag(f_val))
end subroutine calc_separate_error

! 相对残差 |p(x)| / sum(|a_i|*|x|^i)，避免大系数下绝对残差天然偏大导致误判
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

subroutine print_separate_error(x, err_re, err_im, zero_tol)
    complex(wp), intent(in) :: x
    real(wp), intent(in)    :: err_re, err_im, zero_tol
    real(wp) :: re, im

    re = real(x, wp)
    im = aimag(x)

    if (abs(im) < zero_tol) then
        write(*,'(a)',advance='no') "  Real error: "
        write(*,'(es32.20e3)') err_re
    else if (abs(re) < zero_tol) then
        write(*,'(a)',advance='no') "  Imag error: "
        write(*,'(es32.20e3)') err_im
    else
        write(*,'(a)',advance='no') "  Real error: "
        write(*,'(es32.20e3)',advance='no') err_re
        write(*,'(a)',advance='no') "  Imag error: "
        write(*,'(es32.20e3)') err_im
    end if
end subroutine print_separate_error

end program quintic_complex_solver
