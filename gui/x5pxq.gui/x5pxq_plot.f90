!****************************************************************************
!  x5pxq_plot.f90
!
!  五次函数图像渲染模块。
!
!  在主对话框右侧的子窗口（静态控件 IDC_PLOT）内自绘：
!    1) 首一化五次多项式 y = (p(x)-RHS)/A5 的实函数曲线；
!    2) 5 个根中实根在 x 轴上的位置（编号与结果区 Root n 一致）；
!    3) 复根在右上角复平面小图中的位置——与曲线在同一个子窗口内。
!
!  绘制全部使用 GDI：Polyline 画曲线，FillRect 画根标记，TextOut 写刻度。
!****************************************************************************

module x5pxq_plot
    use ifwinty
    use user32
    use gdi32
    use x5pxq_solver, only: wp, zero_thresh

    implicit none
    private
    public :: plot_attach, plot_set_data, plot_refresh

    integer(INT_PTR), parameter :: GWLP_LOCAL = -4
    integer, parameter :: ncurve = 1201          ! 曲线采样点数

    ! 当前绘图数据（由“求解”回调写入）
    real(wp), save :: plot_a(6) = 0.0_wp
    real(wp), save :: plot_rhs = 0.0_wp
    real(wp), save :: plot_maxrel = 0.0_wp
    complex(wp), save :: plot_roots(5) = (0.0_wp, 0.0_wp)
    logical, save :: plot_valid = .false.

    integer(HANDLE), save :: plot_hwnd = 0
    integer(INT_PTR), save :: plot_oldproc = 0

contains

!****************************************************************************
!  把绘图子窗口挂到自己的窗口过程上（子类化静态控件）
!****************************************************************************
    subroutine plot_attach(hdlg)
        integer(HANDLE), intent(in) :: hdlg

        include 'resource.fd'

        plot_hwnd = GetDlgItem(hdlg, int(IDC_PLOT, SINT))
        if (plot_hwnd == 0) return
        plot_oldproc = SetWindowLongPtr(plot_hwnd, GWLP_LOCAL, loc(PlotWndProc))
    end subroutine plot_attach

!****************************************************************************
!  保存一次求解结果：系数、右端项、5 个根与最大相对残差
!****************************************************************************
    subroutine plot_set_data(a, rhs, roots, maxrel)
        real(wp), intent(in)    :: a(6), rhs, maxrel
        complex(wp), intent(in) :: roots(5)

        plot_a      = a
        plot_rhs    = rhs
        plot_roots  = roots
        plot_maxrel = maxrel
        plot_valid  = .true.
    end subroutine plot_set_data

!****************************************************************************
!  请求重绘（求解后调用）
!****************************************************************************
    subroutine plot_refresh()
        type(T_RECT) :: rc
        integer(BOOL) :: ok

        if (plot_hwnd == 0) return
        ok = GetClientRect(plot_hwnd, rc)
        ok = InvalidateRect(plot_hwnd, rc, .TRUE.)
        ok = UpdateWindow(plot_hwnd)
    end subroutine plot_refresh

!****************************************************************************
!  子窗口过程：只处理 WM_PAINT / WM_ERASEBKGND，其余交给原静态控件过程
!****************************************************************************
    integer(INT_PTR) function PlotWndProc(hwnd, msg, wparam, lparam)
!DEC$ ATTRIBUTES STDCALL, ALIAS : 'PLOTWNDPROC' :: PlotWndProc

        integer(HANDLE)  :: hwnd
        integer(UINT)    :: msg
        integer(UINT_PTR):: wparam
        integer(LONG_PTR):: lparam
        integer(INT_PTR) :: res

        if (msg == WM_PAINT) then
            call paint_plot(hwnd)
            PlotWndProc = 0
        else if (msg == WM_ERASEBKGND) then
            PlotWndProc = 1
        else
            res = CallWindowProc(plot_oldproc, hwnd, msg, wparam, lparam)
            PlotWndProc = res
        endif
    end function PlotWndProc

    subroutine paint_plot(hwnd)
        integer(HANDLE) :: hwnd, hdc
        type(T_PAINTSTRUCT) :: ps
        integer(BOOL) :: r

        hdc = BeginPaint(hwnd, ps)
        call draw_scene(hdc, hwnd)
        r = EndPaint(hwnd, ps)
    end subroutine paint_plot

!****************************************************************************
!  全部绘制工作
!****************************************************************************
    subroutine draw_scene(hdc, hwnd)
        integer(HANDLE), intent(in) :: hdc, hwnd

        type(T_RECT) :: rc
        integer(SINT) :: w, h, left, right, top, bot, i, xpix, ypix
        integer(SINT) :: nreal, ncomp
        real(wp) :: b(6), xlo, xhi, ylo, yhi, span, x, y, xt, yt
        real(wp) :: stepx, stepy, yv, rmax, sc, re, im
        type(T_POINT) :: pts(ncurve)
        integer(HANDLE) :: hpen, hold
        integer(BOOL) :: r
        character(len=32)  :: tick
        character(len=120) :: status
        character(len=6)   :: idx

        r = GetClientRect(hwnd, rc)
        w = int(rc%right, SINT)
        h = int(rc%bottom, SINT)

        call solid_rect(hdc, 0, 0, w, h, rgb(255, 255, 255))

        left  = 62
        right = w - 16
        top   = 34
        bot   = h - 44

        if (right - left < 90 .or. bot - top < 70) then
            call put_text(hdc, 8, 8, "绘图区域太小", rgb(120, 120, 120), &
                          .false., TA_LEFT)
            return
        endif

        if (.not. plot_valid) then
            call put_text(hdc, left, (top + bot) / 2, &
                          "点击“求解”后，这里显示五次函数曲线与 5 个根的位置", &
                          rgb(120, 120, 120), .false., TA_LEFT)
            return
        endif

        ! ---- 首一化系数：b1=1，b6=(A0-RHS)/A5 ----
        b(1) = 1.0_wp
        do i = 2, 5
            b(i) = plot_a(i) / plot_a(1)
        end do
        b(6) = (plot_a(6) - plot_rhs) / plot_a(1)

        ! ---- x 视窗：包含全部根的实部，并留 8% 余量 ----
        nreal = 0
        ncomp = 0
        xlo = huge(1.0_wp)
        xhi = -huge(1.0_wp)
        do i = 1, 5
            if (abs(aimag(plot_roots(i))) <= zero_thresh) then
                nreal = nreal + 1
            else
                ncomp = ncomp + 1
            endif
            re = real(plot_roots(i), wp)
            if (abs(re) < 1.0e30_wp) then
                xlo = min(xlo, re)
                xhi = max(xhi, re)
            endif
        end do
        if (xlo >= xhi) then
            xlo = -1.0_wp
            xhi =  1.0_wp
        endif
        span = xhi - xlo
        if (span <= 0.0_wp) span = 1.0_wp
        xlo = xlo - 0.08_wp * span
        xhi = xhi + 0.08_wp * span
        span = xhi - xlo
        if (xlo > 0.0_wp .and. xlo <= 0.30_wp * span) xlo = 0.0_wp
        if (xhi < 0.0_wp .and. -xhi <= 0.30_wp * span) xhi = 0.0_wp

        ! ---- y 视窗：由导数的实零点（曲线拐点）决定，避免两端陡壁压扁图形 ----
        call calc_yrange(b, xlo, xhi, ylo, yhi)

        ! ---- 边框与网格 ----
        call frame_rect(hdc, left, top, right, bot, rgb(150, 150, 150), 1)

        stepx = nice_step(xhi - xlo, 5)
        xt = real(ceiling(xlo / stepx), wp) * stepx
        do while (xt <= xhi + 1.0e-10_wp * max(1.0_wp, abs(xhi)))
            xpix = mapx(xt)
            call stroke_line(hdc, xpix, top, xpix, bot, rgb(228, 228, 228), 1)
            tick = fmt_tick(xt)
            call put_text(hdc, xpix, bot + 6, tick, rgb(90, 90, 90), &
                          .false., TA_CENTER)
            xt = xt + stepx
        end do

        stepy = nice_step(yhi - ylo, 4)
        yt = real(ceiling(ylo / stepy), wp) * stepy
        do while (yt <= yhi + 1.0e-10_wp * max(1.0_wp, abs(yhi)))
            ypix = mapy(yt)
            call stroke_line(hdc, left, ypix, right, ypix, rgb(228, 228, 228), 1)
            tick = fmt_tick(yt)
            call put_text(hdc, left - 6, ypix - 7, tick, rgb(90, 90, 90), &
                          .false., TA_RIGHT)
            yt = yt + stepy
        end do

        ! ---- 坐标轴 ----
        if (ylo <= 0.0_wp .and. yhi >= 0.0_wp) then
            ypix = mapy(0.0_wp)
            call stroke_line(hdc, left, ypix, right, ypix, rgb(70, 70, 70), 1)
        endif
        if (xlo <= 0.0_wp .and. xhi >= 0.0_wp) then
            xpix = mapx(0.0_wp)
            call stroke_line(hdc, xpix, top, xpix, bot, rgb(70, 70, 70), 1)
        endif

        ! ---- 曲线 y = p(x)/A5 ----
        do i = 1, ncurve
            x = xlo + (xhi - xlo) * real(i - 1, wp) / real(ncurve - 1, wp)
            y = pval(x, b)
            if (y > 1.0e30_wp)  y = 1.0e30_wp
            if (y < -1.0e30_wp) y = -1.0e30_wp
            ypix = mapy(y)
            if (ypix < top) ypix = top
            if (ypix > bot) ypix = bot
            pts(i)%x = int(mapx(x), LONG)
            pts(i)%y = int(ypix, LONG)
        end do

        hpen = CreatePen(PS_SOLID, 2, rgb(31, 119, 180))
        hold = SelectObject(hdc, hpen)
        r = Polyline(hdc, pts(1), int(ncurve, SINT))
        hold = SelectObject(hdc, hold)
        r = DeleteObject(hpen)

        ! ---- 实根标记：红方块 + 编号 ----
        do i = 1, 5
            if (abs(aimag(plot_roots(i))) > zero_thresh) cycle
            re = real(plot_roots(i), wp)
            xpix = mapx(re)
            ypix = mapy(0.0_wp)
            call solid_rect(hdc, xpix - 4, ypix - 4, xpix + 4, ypix + 4, &
                            rgb(214, 39, 40))
            write(idx, '(i1)') i
            if (mod(i, 2) == 1) then
                call put_text(hdc, xpix, ypix - 24, trim(idx), &
                              rgb(214, 39, 40), .true., TA_CENTER)
            else
                call put_text(hdc, xpix, ypix + 8, trim(idx), &
                              rgb(214, 39, 40), .true., TA_CENTER)
            endif
        end do

        ! ---- 复根：右上角复平面小图（仍在同一子窗口内） ----
        if (ncomp > 0) then
            call draw_complex_map(hdc, left, right, top, bot)
        endif

        ! ---- 标题、图例与状态行 ----
        call put_text(hdc, left, 6, "曲线 y = (p(x) - RHS)/A5（零点即根）", &
                      rgb(40, 40, 40), .true., TA_LEFT)
        call solid_rect(hdc, right - 176, 10, right - 168, 18, rgb(214, 39, 40))
        call put_text(hdc, right - 162, 7, "实根（编号 = Root n）", &
                      rgb(90, 90, 90), .false., TA_LEFT)

        write(status, '(a,i0,a,i0,a,es9.2)') &
            "实根 ", nreal, " 个，复根 ", ncomp, " 个，最大相对残差 ", &
            real(plot_maxrel, 8)
        call put_text(hdc, left, h - 20, trim(status), &
                      rgb(100, 100, 100), .false., TA_LEFT)

    contains

        integer(SINT) function mapx(v)
            real(wp), intent(in) :: v
            mapx = nint(real(left, wp) + (v - xlo) / (xhi - xlo) &
                        * real(right - left, wp), SINT)
        end function mapx

        integer(SINT) function mapy(v)
            real(wp), intent(in) :: v
            mapy = nint(real(bot, wp) - (v - ylo) / (yhi - ylo) &
                        * real(bot - top, wp), SINT)
        end function mapy

    end subroutine draw_scene

!****************************************************************************
!  复平面小图：复根（实根已在曲线 x 轴上标出）
!****************************************************************************
    subroutine draw_complex_map(hdc, left, right, top, bot)
        integer(HANDLE), intent(in) :: hdc
        integer(SINT), intent(in) :: left, right, top, bot
        integer(SINT) :: bx1, bx2, by1, by2, ccx, ccy, i, xpix, ypix
        integer(SINT) :: ptop, pbot
        real(wp) :: rmax, sc, re, im
        character(len=6) :: idx

        bx2 = right - 4
        bx1 = bx2 - 150
        by1 = top + 4
        by2 = by1 + 116
        if (bx1 <= left + 40 .or. by2 >= bot - 8) return

        call solid_rect(hdc, bx1, by1, bx2, by2, rgb(255, 255, 255))

        ptop = by1 + 20
        pbot = by2 - 8
        ccx  = (bx1 + bx2) / 2
        ccy  = (ptop + pbot) / 2

        rmax = 0.0_wp
        do i = 1, 5
            if (abs(aimag(plot_roots(i))) <= zero_thresh) cycle
            re = abs(plot_roots(i))
            if (re < 1.0e30_wp) rmax = max(rmax, re)
        end do
        if (rmax <= 0.0_wp) rmax = 1.0_wp
        sc = 0.40_wp * real(min(bx2 - bx1 - 16, pbot - ptop), wp) / rmax

        call stroke_line(hdc, bx1 + 8, ccy, bx2 - 8, ccy, rgb(190, 190, 190), 1)
        call stroke_line(hdc, ccx, ptop, ccx, pbot, rgb(190, 190, 190), 1)

        do i = 1, 5
            if (abs(aimag(plot_roots(i))) <= zero_thresh) cycle
            re = real(plot_roots(i), wp)
            im = aimag(plot_roots(i))
            xpix = ccx + int(nint(re * sc), SINT)
            ypix = ccy - int(nint(im * sc), SINT)
            call solid_rect(hdc, xpix - 4, ypix - 4, xpix + 4, ypix + 4, &
                            rgb(123, 63, 160))
            write(idx, '(i1)') i
            call put_text(hdc, xpix + 6, ypix - 16, trim(idx), &
                          rgb(80, 80, 80), .false., TA_LEFT)
        end do

        call frame_rect(hdc, bx1, by1, bx2, by2, rgb(150, 150, 150), 1)
        call put_text(hdc, bx1 + 6, by1 + 3, "复根（复平面）", &
                      rgb(90, 90, 90), .false., TA_LEFT)
    end subroutine draw_complex_map

!****************************************************************************
!  由导数实零点估计曲线的纵向范围
!****************************************************************************
    subroutine calc_yrange(b, xlo, xhi, ylo, yhi)
        real(wp), intent(in)  :: b(6), xlo, xhi
        real(wp), intent(out) :: ylo, yhi

        integer, parameter :: nd = 401
        real(wp) :: xa, xb, da, db, xm, dm, ym, yc(8)
        integer :: k, nc, it

        nc = 0
        xa = xlo
        da = dval(xa, b)
        do k = 2, nd
            xb = xlo + (xhi - xlo) * real(k - 1, wp) / real(nd - 1, wp)
            db = dval(xb, b)
            if (da == 0.0_wp) then
                nc = nc + 1
                yc(min(nc, 8)) = pval(xa, b)
            else if (db /= 0.0_wp .and. da * db < 0.0_wp) then
                do it = 1, 60
                    xm = 0.5_wp * (xa + xb)
                    dm = dval(xm, b)
                    if (da * dm <= 0.0_wp) then
                        xb = xm
                        db = dm
                    else
                        xa = xm
                        da = dm
                    endif
                end do
                nc = nc + 1
                yc(min(nc, 8)) = pval(0.5_wp * (xa + xb), b)
            endif
            xa = xb
            da = db
        end do

        ylo = 0.0_wp
        yhi = 0.0_wp
        do k = 1, min(nc, 8)
            ylo = min(ylo, yc(k))
            yhi = max(yhi, yc(k))
        end do
        if (nc == 0) then
            ylo = min(0.0_wp, pval(xlo, b), pval(xhi, b))
            yhi = max(0.0_wp, pval(xlo, b), pval(xhi, b))
        endif
        if (yhi <= ylo) then
            yhi = ylo + 1.0_wp
            ylo = ylo - 1.0_wp
        endif
        ym  = yhi - ylo
        ylo = ylo - 0.10_wp * ym
        yhi = yhi + 0.10_wp * ym
    end subroutine calc_yrange

    pure real(wp) function pval(x, b)
        real(wp), intent(in) :: x, b(6)
        integer :: j

        pval = b(1)
        do j = 2, 6
            pval = pval * x + b(j)
        end do
    end function pval

    pure real(wp) function dval(x, b)
        real(wp), intent(in) :: x, b(6)
        integer :: j

        dval = 5.0_wp * b(1)
        do j = 2, 5
            dval = dval * x + real(6 - j, wp) * b(j)
        end do
    end function dval

    pure real(wp) function nice_step(span, target)
        real(wp), intent(in) :: span
        integer, intent(in)  :: target
        real(wp) :: raw, p10, f
        integer  :: e

        raw = span / real(max(target, 1), wp)
        if (raw <= 0.0_wp) then
            nice_step = 1.0_wp
            return
        endif
        e   = int(floor(log10(raw)))
        p10 = 10.0_wp ** real(e, wp)
        f   = raw / p10
        if (f < 1.5_wp) then
            nice_step = p10
        else if (f < 3.0_wp) then
            nice_step = 2.0_wp * p10
        else if (f < 7.0_wp) then
            nice_step = 5.0_wp * p10
        else
            nice_step = 10.0_wp * p10
        endif
    end function nice_step

    function fmt_tick(v) result(s)
        real(wp), intent(in) :: v
        character(len=32) :: s, t
        real(wp) :: a
        integer :: k

        a = abs(v)
        if (a >= 1.0e4_wp .or. (a > 0.0_wp .and. a < 1.0e-3_wp)) then
            write(t, '(es9.2)') real(v, 8)
            s = trim(adjustl(t))
            return
        endif

        write(t, '(f12.4)') real(v, 8)
        s = adjustl(t)
        k = len_trim(s)
        do while (k > 1)
            if (s(k:k) /= '0') exit
            k = k - 1
        end do
        if (k > 1) then
            if (s(k:k) == '.') k = k - 1
        endif
        s = s(1:k)
        if (len_trim(s) == 0) s = "0"
    end function fmt_tick

!****************************************************************************
!  GDI 小工具
!****************************************************************************
    pure integer(ULONG) function rgb(r, g, b)
        integer, intent(in) :: r, g, b
        rgb = int(r, ULONG) + 256_ULONG * int(g, ULONG) &
            + 65536_ULONG * int(b, ULONG)
    end function rgb

    subroutine solid_rect(hdc, x1, y1, x2, y2, color)
        integer(HANDLE), intent(in) :: hdc
        integer(SINT), intent(in)   :: x1, y1, x2, y2
        integer(ULONG), intent(in)  :: color
        type(T_RECT) :: rc
        integer(HANDLE) :: hbr
        integer(SINT) :: r

        rc%left   = x1
        rc%top    = y1
        rc%right  = x2
        rc%bottom = y2
        hbr = CreateSolidBrush(color)
        r = FillRect(hdc, rc, hbr)
        r = DeleteObject(hbr)
    end subroutine solid_rect

    subroutine frame_rect(hdc, x1, y1, x2, y2, color, width)
        integer(HANDLE), intent(in) :: hdc
        integer(SINT), intent(in)   :: x1, y1, x2, y2, width
        integer(ULONG), intent(in)  :: color
        integer(HANDLE) :: hpen, hbr, holdpen, holdbr
        integer(BOOL) :: r

        hpen = CreatePen(PS_SOLID, width, color)
        hbr  = GetStockObject(NULL_BRUSH)
        holdpen = SelectObject(hdc, hpen)
        holdbr  = SelectObject(hdc, hbr)
        r = MSFWIN$Rectangle(hdc, x1, y1, x2, y2)
        holdpen = SelectObject(hdc, holdpen)
        holdbr  = SelectObject(hdc, holdbr)
        r = DeleteObject(hpen)
    end subroutine frame_rect

    subroutine stroke_line(hdc, x1, y1, x2, y2, color, width)
        integer(HANDLE), intent(in) :: hdc
        integer(SINT), intent(in)   :: x1, y1, x2, y2, width
        integer(ULONG), intent(in)  :: color
        type(T_POINT) :: pt
        integer(HANDLE) :: hpen, hold
        integer(BOOL) :: r

        hpen = CreatePen(PS_SOLID, width, color)
        hold = SelectObject(hdc, hpen)
        r = MoveToEx(hdc, x1, y1, pt)
        r = MSFWIN$LineTo(hdc, x2, y2)
        hold = SelectObject(hdc, hold)
        r = DeleteObject(hpen)
    end subroutine stroke_line

    subroutine put_text(hdc, x, y, s, color, bold, align)
        integer(HANDLE), intent(in) :: hdc
        integer(SINT), intent(in)   :: x, y, align
        character(*), intent(in)    :: s
        integer(ULONG), intent(in)  :: color
        logical, intent(in)         :: bold
        integer(HANDLE) :: hfont, hold
        integer(SINT) :: oldmode
        integer(BOOL) :: r
        character(len=40) :: face

        if (len_trim(s) == 0) return
        face = "Microsoft YaHei UI" // char(0)
        hfont = CreateFont(-14, 0, 0, 0, merge(FW_BOLD, FW_NORMAL, bold), &
                           0, 0, 0, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, &
                           CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH, &
                           face)
        hold = SelectObject(hdc, hfont)
        oldmode = SetBkMode(hdc, TRANSPARENT)
        r = MSFWIN$SetTextColor(hdc, color)
        r = SetTextAlign(hdc, int(align, UINT))
        r = TextOut(hdc, x, y, s, len_trim(s))
        r = SetTextAlign(hdc, TA_LEFT)
        hold = SelectObject(hdc, hold)
        r = DeleteObject(hfont)
    end subroutine put_text

end module x5pxq_plot
