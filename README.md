# x5pxq.fortran

五次方程复根求解器（Intel ifx，四精度 REAL(128)，残差可达 1e-30 量级）。

## 编译

需先加载 Visual Studio 与 oneAPI 环境（ifx 位于 `D:\Program Files (x86)\Intel\oneAPI`）：

    call "D:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
    call "D:\Program Files (x86)\Intel\oneAPI\setvars.bat"
    ifx /O2 x5pxq.f90 /exe:build\x5pxq.exe

## 运行

输入顺序为 A5 A4 A3 A2 A1 A0 RHS，即 A5*x^5 + A4*x^4 + A3*x^3 + A2*x^2 + A1*x + A0 = RHS：

    echo 1 -15 85 -225 274 -120 0 | build\x5pxq.exe

每个根按 a+bi 输出并给出实部/虚部误差，末尾给出最大相对残差与 PASS/FAIL 判定（阈值 1e-20）。

## GUI 函数图像

Visual Fortran 工程 `gui\x5pxq.gui\x5pxq.gui.vfproj` 的左侧输入系数，右侧新增了一个绘图子窗口（`IDC_PLOT`），在同一窗口内显示：

- 首一化曲线 `y = (p(x) - RHS)/A5`（零点与原方程相同，纵向刻度不被首项系数压扁）；
- 实根在 x 轴上的位置（红色方块，编号对应结果区 `Root n`）；
- 复根在右上角复平面小图中的位置（紫色方块）。

窗口打开时即画出默认示例 `(x-1)(x-2)(x-3)(x-4)(x-5)`；每次“求解 (Solve)”后即时重绘。绘图实现在 `gui\x5pxq.gui\x5pxq_plot.f90`，全部使用 GDI（`Polyline` / `FillRect` / `TextOut`），无需第三方图形库。

## 回归测试

    python build\verify_test.py
