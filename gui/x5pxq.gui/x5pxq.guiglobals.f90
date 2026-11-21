!****************************************************************************
!  Global data, parameters, and structures 
!****************************************************************************

module x5pxqguiGlobals
use iflogm
use ifwinty

implicit none

!  Parameters

integer*4, parameter, public :: SIZEOFAPPNAME = 100

!  Global data

integer(HANDLE)		ghInstance
integer(HANDLE)		ghModule
integer(HANDLE)		ghwndMain
type (dialog) gdlg
integer(INT_PTR)	gOldDlgProc


end module
