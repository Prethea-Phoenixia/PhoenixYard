from tkinter import *
from IB import IB
from misc import center, loadfont, resolvepath
from ctypes import windll
import platform


if __name__ == "__main__":
    # this tells windows that our program will handle scaling ourselves
    winRelease = platform.release()
    if winRelease in ("8", "10"):
        windll.shcore.SetProcessDpiAwareness(1)
    elif winRelease in ("7", "Vista"):
        windll.user32.SetProcessDPIAware()
    else:
        print("Unknown release: ", release, ", skipping DPI handling")

    # this allows us to set our own taskbar icon
    myappid = "mycompany.myproduct.subproduct.version"  # arbitrary string
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    root = Tk()
    root.iconbitmap(resolvepath("ui/logo.ico"))
    # one must supply the entire path
    loadfont(resolvepath("ui/Hack-Regular.ttf"))

    dpi = root.winfo_fpixels("1i")

    # Tk was originally developed for a dpi of 72
    # root.tk.call("tk", "scaling", "-displayof", ".", dpi / 72.0)
    root.tk.call("tk", "scaling", dpi / 72.0)

    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.12"))

    root.option_add("*Font", "Hack 8")
    root.option_add("*tearOff", FALSE)

    root.title("Phoenix's Internal Ballistics Solver v0.3")

    ibPanel = IB(root, dpi)

    center(root)

    root.minsize(root.winfo_width(), root.winfo_height())
    root.mainloop()
