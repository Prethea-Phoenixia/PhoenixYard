from tkinter import Tk
from IB import IB
from misc import center, loadfont, resolvepath
from ctypes import windll
import platform
import multiprocessing

from matplotlib import font_manager

# import matplotlib.font_manager as font_manager

if __name__ == "__main__":
    multiprocessing.freeze_support()
    # this tells windows that our program will handle scaling ourselves
    winRelease = platform.release()
    if winRelease in ("8", "10"):
        windll.shcore.SetProcessDpiAwareness(1)
    elif winRelease in ("7", "Vista"):
        windll.user32.SetProcessDPIAware()
    else:
        print("Unknown release: ", winRelease, ", skipping DPI handling")

    # this allows us to set our own taskbar icon
    # "mycompany.myproduct.subproduct.version"
    myappid = "Phoenix.Internal Ballistics.Solver.043"  # arbitrary string
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    root = Tk()
    root.iconbitmap(resolvepath("ui/logo.ico"))

    loadfont(resolvepath("ui/sarasa-mono-sc-regular.ttf"), False, True)

    font_manager.fontManager.addfont(
        resolvepath("ui/sarasa-mono-sc-regular.ttf")
    )

    dpi = root.winfo_fpixels("1i")

    # Tk was originally developed for a dpi of 72
    # root.tk.call("tk", "scaling", "-displayof", ".", dpi / 72.0)
    scale = 1.0 * dpi / 72.0
    root.tk.call("tk", "scaling", scale)

    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.12"))

    root.option_add("*tearOff", False)

    root.title("PIBS v0.4.3")

    ibPanel = IB(root, dpi, scale)

    center(root)

    # print(font.families())

    # root.minsize(root.winfo_width(), root.winfo_height())
    root.mainloop()
