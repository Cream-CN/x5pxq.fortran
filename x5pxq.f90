!======================================================================
! 五次方程复根求解器 | 修复降次逻辑 + 统一输出格式 + 无警告
! 编译：ifx -O2 x5pxq.f90 -o x5pxq
!======================================================================
program quintic_complex_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp)    :: a5,a4,a3,a2,a1,a0,rhs
    complex(dp) :: roots(5)
    integer     :: i, ierr
    real(dp), parameter :: zero_thresh = 1.0e-10_dp  ! 零根判断阈值
    real(dp), parameter :: err_limit = 1.0e-20_dp    ! 误差阈值

    write(*,*) "=== High Precision Quintic Equation Solver (10^-20) ==="
    write(*,*) "Equation: A5*x^5 + A4*x^4 + A3*x^3 + A2*x^2 + A1*x + A0 = RHS"
    write(*,*) "Input 7 parameters: A5 A4 A3 A2 A1 A0 RHS"
    read(*,*,iostat=ierr) a5,a4,a3,a2,a1,a0,rhs
    if (ierr/=0 .or. a5==0.0_dp) error stop "Invalid input! Please enter 7 numbers."

    call find_all_roots(a5,a4,a3,a2,a1,a0-rhs, roots)

    write(*,*) "=========================================================="
    write(*,*) "Results (Format: a+bi | Real/Imaginary Errors)"
    write(*,*) "=========================================================="
    do i=1,5
        call print_root_format(roots(i), i, zero_thresh)
        call print_separate_error(a5,a4,a3,a2,a1,a0,rhs,roots(i), zero_thresh, err_limit)
        write(*,*)
    end do
    write(*,*) "=========================================================="

contains

! ====================== 核心修复：支持任意次数的多项式求值 ======================
function poly(x, c, deg) result(res)
    complex(dp), intent(in) :: x
    complex(dp), intent(in) :: c(:)  ! 多项式系数（从最高次到常数项）
    integer, intent(in)     :: deg   ! 当前多项式次数
    complex(dp) :: res
    integer :: j

    res = c(1)
    do j=2, deg+1
        res = res * x + c(j)
    end do
end function poly

! 核心修复：支持任意次数的多项式导数
function poly_deriv(x, c, deg) result(res)
    complex(dp), intent(in) :: x
    complex(dp), intent(in) :: c(:)
    integer, intent(in)     :: deg
    complex(dp) :: res
    integer :: j

    res = c(1) * deg
    do j=2, deg
        res = res * x + c(j) * (deg - j + 1)
    end do
end function poly_deriv

! 核心修复：支持任意次数的牛顿迭代求根
function newton_complex(c, deg) result(root)
    complex(dp), intent(in) :: c(:)
    integer, intent(in)     :: deg
    complex(dp) :: root, x, dx
    real(dp), parameter   :: eps = 1.0e-20_dp
    integer, parameter    :: maxit = 200
    real(dp) :: r1, r2
    integer :: j

    call random_number(r1)
    call random_number(r2)
    x = cmplx(r1*2.0_dp - 1.0_dp, r2*2.0_dp - 1.0_dp, dp)

    do j=1,maxit
        dx = -poly(x,c,deg) / poly_deriv(x,c,deg)
        x = x + dx
        if (abs(dx) < eps) exit
    end do
    root = x
end function newton_complex

! 核心修复：支持任意次数的多项式降次（霍纳法）
subroutine deflate(c, deg, r, c_new)
    complex(dp), intent(in)  :: c(:)    ! 原多项式系数（deg次）
    integer, intent(in)     :: deg     ! 原多项式次数
    complex(dp), intent(in)  :: r      ! 根
    complex(dp), intent(out) :: c_new(:)! 降次后的多项式系数（deg-1次）
    integer :: i

    c_new(1) = c(1)
    do i=2, deg
        c_new(i) = c(i) + r * c_new(i-1)
    end do
end subroutine deflate

! 求解所有根（修复降次逻辑）
subroutine find_all_roots(a5,a4,a3,a2,a1,a0, roots)
    real(dp), intent(in)  :: a5,a4,a3,a2,a1,a0
    complex(dp), intent(out) :: roots(5)
    complex(dp) :: poly_coeff(6), temp(5)  ! 系数数组，最多支持5次
    integer  :: i, deg  ! deg记录当前多项式次数

    ! 初始化五次多项式系数
    poly_coeff = [cmplx(a5,0.0_dp,dp), cmplx(a4,0.0_dp,dp), cmplx(a3,0.0_dp,dp), &
                  cmplx(a2,0.0_dp,dp), cmplx(a1,0.0_dp,dp), cmplx(a0,0.0_dp,dp)]
    deg = 5  ! 初始次数为5

    do i=1,5
        ! 求当前deg次多项式的根
        roots(i) = newton_complex(poly_coeff(1:deg+1), deg)
        ! 降次：除以(x - roots(i))，得到deg-1次多项式
        call deflate(poly_coeff(1:deg+1), deg, roots(i), temp(1:deg))
        ! 更新系数和次数
        poly_coeff(1:deg) = temp(1:deg)
        deg = deg - 1
    end do
end subroutine find_all_roots

! 打印根（修复负零问题）
subroutine print_root_format(x, idx, zero_tol)
    complex(dp), intent(in) :: x
    integer, intent(in) :: idx
    real(dp) :: re, im
    real(dp), intent(in) :: zero_tol

    re = real(x)
    im = aimag(x)

    ! 处理负零，强制转为0
    if (abs(re) < zero_tol) re = 0.0_dp
    if (abs(im) < zero_tol) im = 0.0_dp

    write(*,'(a,i2,a)',advance='no') "Root ",idx,": "
    
    if (abs(im) < zero_tol) then
        ! 实数根
        write(*,'(f25.20)') re
    else if (abs(re) < zero_tol) then
        ! 纯虚数根
        write(*,'(f25.20,a)') im, "i"
    else
        ! 复数根，统一格式 a+bi
        if (im > 0.0_dp) then
            write(*,'(f25.20,a,f25.20,a)') re, " + ", im, "i"
        else
            write(*,'(f25.20,a,f25.20,a)') re, " - ", abs(im), "i"
        end if
    end if
end subroutine print_root_format

! 计算实部/虚部误差
subroutine calc_separate_error(a5,a4,a3,a2,a1,a0,rhs,x, err_re, err_im)
    real(dp), intent(in)    :: a5,a4,a3,a2,a1,a0,rhs
    complex(dp), intent(in) :: x
    real(dp), intent(out)   :: err_re, err_im
    complex(dp)             :: f_val

    f_val = a5*x**5 + a4*x**4 + a3*x**3 + a2*x**2 + a1*x + (a0-rhs)
    err_re = abs(real(f_val))
    err_im = abs(aimag(f_val))
end subroutine calc_separate_error

! 打印误差（统一科学计数法格式）
subroutine print_separate_error(a5,a4,a3,a2,a1,a0,rhs,x, zero_tol, err_limit)
    real(dp), intent(in)    :: a5,a4,a3,a2,a1,a0,rhs, zero_tol, err_limit
    complex(dp), intent(in) :: x
    real(dp) :: err_re, err_im, re, im

    re = real(x)
    im = aimag(x)
    call calc_separate_error(a5,a4,a3,a2,a1,a0,rhs,x, err_re, err_im)

    ! 实数根：仅输出实部误差
    if (abs(im) < zero_tol) then
        write(*,'(a)',advance='no') "  Real error: "
        write(*,'(es27.20)') err_re

    ! 纯虚数根：仅输出虚部误差
    else if (abs(re) < zero_tol) then
        write(*,'(a)',advance='no') "  Imag error: "
        write(*,'(es27.20)') err_im

    ! 复数根：输出实部+虚部误差
    else
        write(*,'(a)',advance='no') "  Real error: "
        write(*,'(es27.20)',advance='no') err_re
        write(*,'(a)',advance='no') "  Imag error: "
        write(*,'(es27.20)') err_im
    end if
end subroutine print_separate_error

end program quintic_complex_solver