from tkinter import Tk, Frame, Menu
from tkinter import ttk
from misc import resolvepath
from misc import loadfont
from multiprocessing import freeze_support

from ctypes import windll
import platform


class IB(Frame):
    def __init__(self, parent, menubar):
        ttk.Frame.__init__(self, parent)
        pass


if __name__ == "__main__":
    freeze_support()

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
    myappid = "Phoenix.External.Ballistics.Solver.0.1"  # arbitrary string
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    root = Tk()
    root.iconbitmap(resolvepath("ui/logo.ico"))

    loadfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"), False, True)

    dpi = root.winfo_fpixels("1i")
    root.tk.call("tk", "scaling", dpi / 72.0)

    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.12"))

    root.option_add("*tearOff", False)
    root.title("PEBS v0.1")
    menubar = Menu(root)
    root.config(menu=menubar)

    ibFrame = IB(root, menubar)
    ibFrame.pack(expand=1, fill="both", side="left")

    root.minsize(root.winfo_width(), root.winfo_height())  # set minimum size
    root.state("zoomed")  # maximize window

    root.mainloop()
