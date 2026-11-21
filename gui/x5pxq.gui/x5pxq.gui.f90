!  x5pxq.gui.f90
!
!  五次方程求根器 —— Intel Fortran 对话框（IFLOGM）界面
!
!  FUNCTIONS:
!   WinMain()       - 入口；创建对话框、处理消息循环
!   x5pxqguiSub()   - 对话框回调（初始化 / 销毁）
!   x5pxqguiSolve() - “求解”按钮回调：读入系数、调用四精度求解器、显示结果
!   x5pxqguiExit()  - “退出”按钮回调

!****************************************************************************
!
!  FUNCTION: WinMain( hInstance, hPrevInstance, lpszCmdLine, nCmdShow )
!
!  PURPOSE:  Entry point for the application
!
!  COMMENTS: Displays the main window and processes the message loop
!
!****************************************************************************

function WinMain( hInstance, hPrevInstance, lpszCmdLine, nCmdShow )
!DEC$ IF DEFINED(_X86_)
!DEC$ ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WinMain
!DEC$ ELSE
!DEC$ ATTRIBUTES STDCALL, ALIAS : 'WinMain' :: WinMain
!DEC$ ENDIF

    use user32
    use kernel32
    use iflogm
    use x5pxqguiGlobals
    use x5pxq_solver, only: wp, solve_quintic
    use x5pxq_plot, only: plot_attach, plot_set_data, plot_refresh

    implicit none

    integer(SINT) :: WinMain
    integer(HANDLE) hInstance
    integer(HANDLE) hPrevInstance
    integer(LPWSTR) lpszCmdLine
    integer(SINT)   nCmdShow

    include 'resource.fd'

    external x5pxqguiSub
    external x5pxqguiSolve
    external x5pxqguiExit
    external x5pxqguiExpr
    external x5pxqguiMethod
    external x5pxqguiAbout
    external DlgSubProc

    ! Variables
    type (T_MSG)            mesg
    integer*4               ret
    integer(LRESULT)        lret
    integer(BOOL)           bret
    integer(SINT)           iret
    logical                 lg
    integer(HANDLE)         hMenu
    integer, parameter :: GWLP_WNDPROC_LOCAL = -4

    ! 菜单在运行时挂到对话框上：把 MENU 写进对话框模板时 IFLOGM 建出的窗口
    ! 会保持隐藏，改为先显示对话框再 SetMenu 更稳妥。
    interface
      function LoadMenuOrdinal( hInstance, menuID ) result( hm )
        !DEC$ ATTRIBUTES DEFAULT, STDCALL, DECORATE, ALIAS : 'LoadMenuA' :: LoadMenuOrdinal
        use ifwinty
        integer(HANDLE) :: hInstance
        integer(HANDLE) :: menuID
        integer(HANDLE) :: hm
      end function LoadMenuOrdinal
    end interface

    ghInstance = hInstance
    ghModule = GetModuleHandle(NULL)
    ghwndMain = NULL

    lret = DlgInit(IDD_x5pxqgui_DIALOG, gdlg)
    if (lret == FALSE) goto 99999
    lret = DlgSetSub(gdlg, IDD_x5pxqgui_DIALOG, x5pxqguiSub)
    lret = DlgSetSub(gdlg, IDM_SOLVE, x5pxqguiSolve)
    lret = DlgSetSub(gdlg, IDM_EXIT, x5pxqguiExit)
    if (lret == FALSE) goto 99999

    ! 默认示例：x^5 - 15x^4 + 85x^3 - 225x^2 + 274x - 120 = 0（根 1,2,3,4,5）
    lg = DlgSetChar(gdlg, IDC_EDIT_A5, "1")
    lg = DlgSetChar(gdlg, IDC_EDIT_A4, "-15")
    lg = DlgSetChar(gdlg, IDC_EDIT_A3, "85")
    lg = DlgSetChar(gdlg, IDC_EDIT_A2, "-225")
    lg = DlgSetChar(gdlg, IDC_EDIT_A1, "274")
    lg = DlgSetChar(gdlg, IDC_EDIT_A0, "-120")
    lg = DlgSetChar(gdlg, IDC_EDIT_RHS, "0")
    lg = DlgSetChar(gdlg, IDC_RESULT, "点击“求解”按钮计算结果。")

    lret = DlgModeless(gdlg, nCmdShow)
    if (lret == FALSE) goto 99999

    ! ------------------------------------------------------------------
    !  菜单“工具 / 帮助”的菜单项不是对话框控件，IFLOGM 的 DlgSetSub 只认
    !  控件 ID（菜单命令的 lParam 为 0，IFLOGM 会直接丢弃），所以这里接管
    !  对话框窗口过程：只截获菜单 WM_COMMAND，其余消息原样转发给 IFLOGM。
    ! ------------------------------------------------------------------
    gOldDlgProc = SetWindowLongPtr(gdlg%hwnd, GWLP_WNDPROC_LOCAL, loc(DlgSubProc))

    ! 给主对话框挂上“工具 / 帮助”菜单
    hMenu = LoadMenuOrdinal(ghInstance, int(IDR_MAINMENU, HANDLE))
    if (hMenu /= NULL .and. hMenu /= 0_HANDLE) then
      lret = SetMenu(gdlg%hwnd, hMenu)
      lret = DrawMenuBar(gdlg%hwnd)
    end if

    ! ------------------------------------------------------------------
    !  把右侧静态控件子类化为绘图子窗口，并先用默认示例画出曲线和根，
    !  让窗口一打开就能看到函数图像与 5 个根的位置。
    ! ------------------------------------------------------------------
    call plot_attach(gdlg%hwnd)
    block
      real(wp)    :: a0(6), maxrel0
      complex(wp) :: roots0(5)
      a0 = [1.0_wp, -15.0_wp, 85.0_wp, -225.0_wp, 274.0_wp, -120.0_wp]
      call solve_quintic(a0, 0.0_wp, roots0, maxrel0)
      call plot_set_data(a0, 0.0_wp, roots0, maxrel0)
      call plot_refresh()
    end block

    ! Read and process messsages
    do while( GetMessage (mesg, NULL, 0, 0) ) 
       if ( DlgIsDlgMessage(mesg) .EQV. FALSE ) then
           bret = TranslateMessage( mesg )
           lret  = DispatchMessage( mesg )
       end if
    end do
    call DlgUninit(gdlg)

    WinMain = mesg.wParam
    return

99999 &

    iret = MessageBox(ghwndMain, "Error initializing application x5pxqgui"C, &
                     "Error"C, MB_OK)
    WinMain = 0

end 

!****************************************************************************
!
!  FUNCTION: DlgSubProc ( hwnd, msg, wparam, lparam )
!
!  PURPOSE:  主对话框的窗口过程。只截获菜单命令（工具/帮助），
!            其余消息全部转交 IFLOGM 原来的窗口过程处理。
!
!****************************************************************************

integer(INT_PTR) function DlgSubProc( hwnd, msg, wparam, lparam )
!DEC$ ATTRIBUTES STDCALL, ALIAS : 'DLGSUBPROC' :: DlgSubProc

  use user32
  use iflogm
  use x5pxqguiGlobals

  implicit none

  include 'resource.fd'

  integer(HANDLE)  hwnd
  integer(UINT)    msg
  integer(fWPARAM) wparam
  integer(fLPARAM) lparam
  integer(INT_PTR) cid

  external x5pxqguiExpr
  external x5pxqguiMethod
  external x5pxqguiAbout

  if (msg == WM_COMMAND) then
    cid = iand(int(wparam, INT_PTR), int(z'FFFF', INT_PTR))
    select case (int(cid))
    case (int(IDM_EXPR_TO_R128))
      call x5pxqguiExpr(gdlg, IDM_EXPR_TO_R128, dlg_clicked)
      DlgSubProc = 0
      return
    case (int(IDM_HELP_METHOD))
      call x5pxqguiMethod(gdlg, IDM_HELP_METHOD, dlg_clicked)
      DlgSubProc = 0
      return
    case (int(IDM_HELP_ABOUT))
      call x5pxqguiAbout(gdlg, IDM_HELP_ABOUT, dlg_clicked)
      DlgSubProc = 0
      return
    end select
  end if

  DlgSubProc = CallWindowProc(gOldDlgProc, hwnd, msg, wparam, lparam)
end function DlgSubProc

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiSub ( dlg, id, callbacktype )
!
!  PURPOSE:  Dialog box callback for destroy
!
!****************************************************************************

SUBROUTINE x5pxqguiSub( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiSub

  use user32
  use iflogm

  implicit none

  type (dialog) dlg
  integer id, callbacktype

  if (callbacktype == dlg_destroy) then
    call PostQuitMessage(0)
  endif

  END SUBROUTINE x5pxqguiSub

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiSolve ( dlg, id, callbacktype )
!
!  PURPOSE:  “求解”按钮：读取 7 个输入，调用四精度求解器，输出 5 个根
!
!****************************************************************************

SUBROUTINE x5pxqguiSolve( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiSolve

  use iflogm
  use x5pxq_solver
  use x5pxq_plot, only: plot_set_data, plot_refresh

  implicit none

  include 'resource.fd'

  type (dialog) dlg
  integer id, callbacktype

  ! 局部变量
  character(len=6) :: labels(7) = [character(len=6) :: "A5","A4","A3","A2","A1","A0","RHS"]
  integer, parameter :: edit_ids(7) = [IDC_EDIT_A5, IDC_EDIT_A4, IDC_EDIT_A3, &
                                       IDC_EDIT_A2, IDC_EDIT_A1, IDC_EDIT_A0, IDC_EDIT_RHS]
  character(len=64) :: text
  real(wp)          :: vals(7)
  real(wp)          :: a(6), rhs, maxrel, err_re, err_im
  complex(wp)       :: roots(5)
  character(len=4096) :: report
  character(len=200)  :: line
  integer :: i, j, ios
  logical :: lg
  real(wp), parameter :: err_limit = 1.0e-20_wp

  if (callbacktype /= dlg_clicked) return

  ! ---- 读取 7 个输入框 ----
  do j=1,7
    text = ""
    lg = DlgGetChar(dlg, edit_ids(j), text)
    if (len_trim(text) == 0) then
      report = "输入不完整：" // trim(labels(j)) // " 为空，请填写数字。"
      call set_result(dlg, report)
      return
    end if
    read(text, *, iostat=ios) vals(j)
    if (ios /= 0) then
      report = "输入无效：" // trim(labels(j)) // " = """ // trim(text) // """ 不是有效数字。"
      call set_result(dlg, report)
      return
    end if
  end do

  a(1:6) = vals(1:6)
  rhs    = vals(7)

  if (a(1) == 0.0_wp) then
    report = "输入无效：A5 不能为 0，否则方程不是五次。"
    call set_result(dlg, report)
    return
  end if

  ! ---- 求解 ----
  call solve_quintic(a, rhs, roots, maxrel)

  ! ---- 组织输出文本 ----
  report = ""
  report = trim(report) // "方程：" // poly_text(a, rhs) // char(13) // char(10)
  report = trim(report) // "精度：REAL(128) 四精度，残差阈值 1e-20" // char(13) // char(10)
  report = trim(report) // char(13) // char(10)

  do i=1,5
    call root_errors(a, rhs, roots(i), err_re, err_im)
    write(line,'(a,i1,a,a)') "Root ", i, ":  ", trim(root_text(roots(i)))
    report = trim(report) // trim(line) // char(13) // char(10)

    write(line,'(a,es11.3,a,es11.3)') "        实部误差 ", err_re, "    虚部误差 ", err_im
    report = trim(report) // trim(line) // char(13) // char(10)
  end do

  report = trim(report) // char(13) // char(10)
  write(line,'(a,es11.3,a,es9.2,a)') "最大相对残差 = ", maxrel, "   (阈值 ", err_limit, ")"
  report = trim(report) // trim(line) // char(13) // char(10)

  if (maxrel <= err_limit) then
    write(line,'(a)') "精度判定：PASS —— 全部残差均小于 1e-20"
  else
    write(line,'(a)') "精度判定：FAIL —— 存在残差超过 1e-20 的根"
  end if
  report = trim(report) // trim(line)

  call set_result(dlg, report)

  ! ---- 把本次系数与根送到右侧子窗口，重绘曲线与根的位置 ----
  call plot_set_data(a, rhs, roots, maxrel)
  call plot_refresh()

contains

  !--------------------------------------------------------------------
  !  把文本写进结果框，并把插入点复位到开头（否则编辑框会自动
  !  滚动到底部，最上面的“方程/精度/Root 1”几行会被挡住）
  !--------------------------------------------------------------------
  subroutine set_result(d, msgtext)
    use iflogm

    type (dialog), intent(inout) :: d
    character(*), intent(in)     :: msgtext
    logical :: l
    integer :: ir
    integer, parameter :: EM_SETSEL_MSG = 177   ! EM_SETSEL = 0x00B1

    l  = DlgSetChar(d, IDC_RESULT, msgtext)
    call DlgFlush(d, .true.)
    ir = DlgSendCtrlMessage(d, IDC_RESULT, EM_SETSEL_MSG, 0, 0)
  end subroutine set_result

  END SUBROUTINE x5pxqguiSolve

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiExit ( dlg, id, callbacktype )
!
!  PURPOSE:  “退出”按钮：结束消息循环
!
!****************************************************************************

SUBROUTINE x5pxqguiExit( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiExit

  use user32
  use iflogm

  implicit none

  type (dialog) dlg
  integer id, callbacktype

  if (callbacktype == dlg_clicked) then
    call PostQuitMessage(0)
  endif

  END SUBROUTINE x5pxqguiExit

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiExpr ( dlg, id, callbacktype )
!
!  PURPOSE:  菜单“工具 → 精确代数式转 REAL(128)”: 打开表达式转换对话框
!
!****************************************************************************

SUBROUTINE x5pxqguiExpr( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiExpr

  use user32
  use iflogm

  implicit none

  include 'resource.fd'

  type (dialog) dlg
  integer id, callbacktype
  type (dialog) ed
  logical lg

  external x5pxqguiExprConv
  external x5pxqguiClose

  if (callbacktype /= dlg_clicked) return

  lg = DlgInit(IDD_EXPR, ed)
  if (.not. lg) return

  lg = DlgSetSub(ed, IDM_EXPR_CONVERT, x5pxqguiExprConv)
  lg = DlgSetSub(ed, IDM_EXPR_CLOSE,   x5pxqguiClose)
  lg = DlgSetChar(ed, IDC_EXPR_IN,  "1/3")
  lg = DlgSetChar(ed, IDC_EXPR_OUT, "在上面的输入框里写下代数式，然后点击“转换”。")
  ! 立即落盘，否则 DlgGetChar 会返回缓存里的初始值而不是用户后来改的内容
  call DlgFlush(ed, .true.)
  lg = DlgModal(ed)
  call DlgUninit(ed)

  END SUBROUTINE x5pxqguiExpr

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiExprConv ( dlg, id, callbacktype )
!
!  PURPOSE:  表达式转换对话框里的“转换”按钮：四精度求值并显示
!
!****************************************************************************

SUBROUTINE x5pxqguiExprConv( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiExprConv

  use iflogm
  use x5pxq_solver

  implicit none

  include 'resource.fd'

  type (dialog) dlg
  integer id, callbacktype
  character(len=600) :: expr
  character(len=900) :: outtext
  real(wp) :: v
  logical  :: ok, lg

  if (callbacktype /= dlg_clicked) return

  expr = ""
  lg = DlgGetChar(dlg, IDC_EXPR_IN, expr)
  if (len_trim(expr) == 0) then
    lg = DlgSetChar(dlg, IDC_EXPR_OUT, "请输入代数式，例如 1/3 或 sqrt(2)。")
    return
  end if

  call eval_expr(trim(expr), v, ok)
  outtext = expr_report(trim(expr), v, ok)

  lg = DlgSetChar(dlg, IDC_EXPR_OUT, outtext)
  call DlgFlush(dlg, .true.)

  END SUBROUTINE x5pxqguiExprConv

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiMethod ( dlg, id, callbacktype )
!
!  PURPOSE:  菜单“帮助 → 五次方程近似法”: 显示算法说明
!
!****************************************************************************

SUBROUTINE x5pxqguiMethod( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiMethod

  use user32
  use iflogm

  implicit none

  include 'resource.fd'

  type (dialog) dlg
  integer id, callbacktype
  type (dialog) md
  logical lg
  character(len=2000) :: txt
  character(2) :: nl

  external x5pxqguiClose

  if (callbacktype /= dlg_clicked) return

  nl = achar(13) // achar(10)

  txt = "算法说明：五次方程的近似求解" // nl // nl
  txt = trim(txt) // "1) 尺度归一化" // nl // &
    "   先把方程除以 A5，化成首一多项式 x^5 + b4 x^4 + … + b0 = 0，" // nl // &
    "   这样求根过程与 A5 的量级无关，避免大系数造成的上溢/下溢。" // nl // nl // &
    "2) 根界与初值" // nl // &
    "   用 Fujiwara 型根界 R = 2*max(|b_i|^(1/(5-i))) 估计所有根的模上界，" // nl // &
    "   牛顿迭代的初值随机取在半径 R 的圆内（8 个黄金角种子，确定性可复现）。" // nl // nl // &
    "3) Newton 迭代 + 收敛判定" // nl // &
    "   对每个种子做复数 Newton 迭代 x ← x - p(x)/p'(x)，" // nl // &
    "   直到 |Δx| 相对步长小于 1e-33 或迭代次数用尽，再按残差择优。" // nl // nl // &
    "4) 对原始方程精修（polish）" // nl // &
    "   逐次降次会把误差累积到后面的根，所以五个根都求出之后，" // nl // &
    "   再对**原始**五次多项式逐个精修若干步，把残差压回四精度极限。" // nl // nl
  txt = trim(txt) // "5) 共轭配对" // nl // &
    "   实系数多项式的复根必然成对出现；对数值上接近共轭的两个根强制取平均，" // nl // &
    "   让虚部严格反号，避免出现 1e-30 量级的假虚部。" // nl // nl // &
    "6) 精度判定" // nl // &
    "   残差用相对残差 |p(x)| / Σ|a_i||x|^i，与阈值 1e-20 比较；" // nl // &
    "   REAL(128) 的机器精度约为 1.93E-34，所以正常情况下残差能到 1e-31～1e-35。" // nl // nl // &
    "注意：重根（例如 (x-1)^5）受噪声放大影响，只能精确到 eps^(1/5) ≈ 1e-7 量级，" // nl // &
    "      这是重根问题的固有限制，不是算法缺陷。"

  lg = DlgInit(IDD_METHOD, md)
  if (.not. lg) return
  lg = DlgSetSub(md, 1, x5pxqguiClose)      ! IDOK = 1
  lg = DlgSetChar(md, IDC_METHOD_TEXT, txt)
  call DlgFlush(md, .true.)
  lg = DlgModal(md)
  call DlgUninit(md)

  END SUBROUTINE x5pxqguiMethod

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiAbout ( dlg, id, callbacktype )
!
!  PURPOSE:  菜单“帮助 → 关于”: 显示版本信息
!
!****************************************************************************

SUBROUTINE x5pxqguiAbout( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiAbout

  use user32
  use iflogm

  implicit none

  include 'resource.fd'

  type (dialog) dlg
  integer id, callbacktype
  type (dialog) ad
  logical lg

  external x5pxqguiClose

  if (callbacktype /= dlg_clicked) return

  lg = DlgInit(IDD_ABOUT, ad)
  if (.not. lg) return
  lg = DlgSetSub(ad, 1, x5pxqguiClose)      ! IDOK = 1
  lg = DlgModal(ad)
  call DlgUninit(ad)

  END SUBROUTINE x5pxqguiAbout

!****************************************************************************
!
!  SUBROUTINE: x5pxqguiClose ( dlg, id, callbacktype )
!
!  PURPOSE:  模态小对话框的“确定/关闭”按钮：结束模态循环
!
!****************************************************************************

SUBROUTINE x5pxqguiClose( dlg, id, callbacktype )
!DEC$ ATTRIBUTES DEFAULT :: x5pxqguiClose

  use iflogm

  implicit none

  type (dialog) dlg
  integer id, callbacktype

  if (callbacktype == dlg_clicked) then
    call DlgExit(dlg)
  end if

  END SUBROUTINE x5pxqguiClose
