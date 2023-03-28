from tkinter import *
from tkinter import ttk
from gun import *
import os
import sys

_prefix = {
    "y": 1e-24,  # yocto
    "z": 1e-21,  # zepto
    "a": 1e-18,  # atto
    "f": 1e-15,  # femto
    "p": 1e-12,  # pico
    "n": 1e-9,  # nano
    "u": 1e-6,  # micro
    "m": 1e-3,  # mili
    "": 1,  # unit
    "k": 1e3,  # kilo
    "M": 1e6,  # mega
    "G": 1e9,  # giga
    "T": 1e12,  # tera
    "P": 1e15,  # peta
    "E": 1e18,  # exa
    "Z": 1e21,  # zetta
    "Y": 1e24,  # yotta
}

from ctypes import windll, byref, create_unicode_buffer, create_string_buffer


FR_PRIVATE = 0x10
FR_NOT_ENUM = 0x20


def resolvepath(path):
    if getattr(sys, "frozen", False):
        # If the 'frozen' flag is set, we are in bundled-app mode!
        resolved_path = os.path.abspath(os.path.join(sys._MEIPASS, path))
    else:
        # Normal development mode. Use os.getcwd() or __file__ as appropriate in your case...
        resolved_path = os.path.abspath(os.path.join(os.getcwd(), path))

    return resolved_path


def loadfont(fontpath, private=True, enumerable=False):
    # Makes fonts located in file `fontpath` available to the font system.
    # `private`     if True, other processes cannot see this font, and this
    #               font will be unloaded when the process dies
    # `enumerable`  if True, this font will appear when enumerating fonts
    #
    # See https://msdn.microsoft.com/en-us/library/dd183327(VS.85).aspx

    # This function was taken from
    # https://github.com/ifwe/digsby/blob/f5fe00244744aa131e07f09348d10563f3d8fa99/digsby/src/gui/native/win/winfonts.py#L15
    # This function is written for Python 2.x. For 3.x, you
    # have to convert the isinstance checks to bytes and str
    # For 3.x, you have to convert the isinstance checks to bytes and str

    if isinstance(fontpath, bytes):
        pathbuf = create_string_buffer(fontpath)
        AddFontResourceEx = windll.gdi32.AddFontResourceExA

    elif isinstance(fontpath, str):
        pathbuf = create_unicode_buffer(fontpath)
        AddFontResourceEx = windll.gdi32.AddFontResourceExW
    else:
        raise TypeError("fontpath must be of type str or unicode")

    flags = (FR_PRIVATE if private else 0) | (
        FR_NOT_ENUM if not enumerable else 0
    )
    numFontsAdded = AddFontResourceEx(byref(pathbuf), flags, 0)
    return bool(numFontsAdded)


def toSI(v, dec=4):
    for prefix, magnitude in zip(_prefix.keys(), _prefix.values()):
        if 1 <= (v / magnitude) < 1e3:
            vstr = "{:#.{:}g}".format(v / magnitude, dec)
            return vstr + " " * (dec + 1 - len(vstr) + vstr.find(".")) + prefix
    if v == 0:
        return "{:#.{:}g}".format(v, dec)
    else:
        raise ValueError(v + " not possible to assign a SI prefix")


def validateNN(inp):
    """
    validate an input if it results in:
    - result >=0
    - result is empty
    in the latter case, the empty field will be handled by tracing
    change in variable.
    """
    if inp == "":
        return True
    try:
        if float(inp) >= 0:
            return True
        else:
            return False
    except ValueError:
        return False


def validatePI(inp):
    """
    validate an input such that the result is a positive integer"""

    if inp == "":
        return True  # we will catch this by filling the default value
    try:
        return float(inp).is_integer() and float(inp) > 0 and "." not in inp
    except ValueError:
        return False


def formatFloatInput(event):
    v = event.widget.get()
    if v == "":
        event.widget.insert(0, event.widget.default)
    else:
        event.widget.delete(0, END)
        event.widget.insert(0, float(v))


def formatIntInput(event):
    v = event.widget.get()
    if v == "":
        event.widget.insert(0, event.widget.default)
    else:
        event.widget.delete(0, END)
        event.widget.insert(0, int(v))


def dot_aligned(matrix):
    transposed = []

    for seq in zip(*matrix):
        snums = []
        for n in seq:
            try:
                snums.append(toSI(float(n)))
            except ValueError:
                snums.append(n)
        dots = [s.find(".") for s in snums]
        m = max(dots)
        transposed.append(tuple(" " * (m - d) + s for s, d in zip(snums, dots)))

    return tuple(zip(*transposed))


class IB(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        parent.title("Phoenix's Internal Ballistics Solver v0.1.1")

        self.compositions = GrainComp.readFile(
            resolvepath("data/propellants.csv")
        )
        self.geometries = {i.desc: i for i in Geometry}

        self.prop = None
        self.gun = None
        self.tableData = []
        self.errorLst = []
        self.geomError = False

        self.propOptions = tuple(self.compositions.keys())
        self.geoOptions = tuple(self.geometries.keys())
        self.domainOptions = ("time", "length")

        self.geoLocked = IntVar()

        parent.columnconfigure(0, weight=0)
        parent.columnconfigure(1, weight=1)
        parent.columnconfigure(2, weight=0)

        parent.rowconfigure(0, weight=1)
        parent.rowconfigure(1, weight=0)

        self.addSpecFrm(parent)
        self.addParFrm(parent)
        self.addErrFrm(parent)
        self.addTblFrm(parent)
        self.addOpsFrm(parent)
        parent.bind("<space>", self.calculate)

    def calculate(self, event=None):
        # force an immediate redraw after calculation
        compo = self.compositions[self.dropProp.get()]
        geom = self.geometries[self.dropGeom.get()]

        self.tableData = []

        # self.updateSpec(compo)

        self.va.set(self.va.default)
        self.be.set(self.be.default)
        self.te.set(self.te.default)

        try:
            self.prop = Propellant(
                compo,
                geom,
                float(self.webmm.get()) / 1000,
                float(self.permm.get()) / 1000,
                float(self.grlR.get()),
            )

        except Exception as e:
            self.prop = None
            self.errorLst.append("Exception when defining propellant:")
            self.errorLst.append(str(e))

        try:
            chamberVolume = (
                float(self.chgkg.get())
                / self.prop.rho_p
                / self.prop.maxLF
                / float(self.ldf.get())
                * 100
            )
            self.gun = Gun(
                float(self.calmm.get()) / 1000,
                float(self.shtkg.get()),
                self.prop,
                float(self.chgkg.get()),
                chamberVolume,
                float(self.stpMPa.get()) * 1e6,
                float(self.tblmm.get()) / 1000,
                float(self.clr.get()),
            )

            self.va.set(round(self.gun.v_j, 1))

        except Exception as e:
            self.gun = None
            self.errorLst.append("Exception when defining guns:")
            self.errorLst.append(str(e))

        if self.gun is not None:
            try:
                self.tableData = self.gun.integrate(
                    steps=int(self.steps.get()),
                    dom=self.dropOptn.get(),
                    tol=10 ** -(int(self.accExp.get())),
                )
                i = tuple(i[0] for i in self.tableData).index("SHOT EXIT")
                te = (self.tableData[i][4] / self.gun.v_j) ** 2
                self.te.set(round(te * 100, 1))
                self.be.set(round(te / self.gun.phi * 100, 1))
            except Exception as e:
                self.errorLst.append("Exception while solving:")
                self.errorLst.append(str(e))

        self.tv.delete(*self.tv.get_children())
        self.tableData = dot_aligned(self.tableData)

        for row in self.tableData:
            self.tv.insert("", "end", values=row, tags=(row[0], "monospace"))

        self.updateError()

    def updateError(self):
        self.errorText.delete("1.0", "end")
        if self.geomError:
            self.errorLst.append("Invalid geometry")

        for line in self.errorLst:
            self.errorText.insert("end", line + "\n")
        self.errorLst = []

    def addSpecFrm(self, parent):
        specFrm = ttk.LabelFrame(parent, text="Design Summary")
        specFrm.grid(row=0, column=0, rowspan=2, sticky="nsew")
        i = 0
        self.va, _, i = self.add12Disp(specFrm, i, "Asymptotic Velocity", "m/s")
        self.te, _, i = self.add12Disp(specFrm, i, "Thermal Efficiency", "%")
        self.be, _, i = self.add12Disp(specFrm, i, "Ballistic Efficiency", "%")

    def updateSpec(self, *inp):
        self.specs.config(state="normal")
        compo = self.compositions[self.dropProp.get()]
        self.specs.delete("1.0", "end")
        for line in compo.desc.split(","):
            self.specs.insert("end", line + "\n")
        self.specs.config(state="disabled")

        return True

    def addParFrm(self, parent):
        parFrm = ttk.LabelFrame(parent, text="Parameters")
        parFrm.grid(row=0, column=2, sticky="nsew")
        parFrm.columnconfigure(0, weight=3)
        parFrm.columnconfigure(1, weight=3)
        parFrm.columnconfigure(2, weight=1)

        # validation
        validationNN = parent.register(validateNN)

        i = 0
        self.calmm, _, i = self.add3Input(
            parFrm, i, "Caliber", "mm", "0.0", validationNN
        )
        self.tblmm, _, i = self.add3Input(
            parFrm, i, "Tube Length", "mm", "0.0", validationNN
        )
        self.shtkg, _, i = self.add3Input(
            parFrm, i, "Shot Mass", "kg", "0.0", validationNN
        )
        self.chgkg, _, i = self.add3Input(
            parFrm, i, "Charge Mass", "kg", "0.0", validationNN
        )

        ttk.Label(parFrm, text="Propellant").grid(
            row=4, column=0, padx=5, pady=5, sticky="nsew"
        )

        specVal = parent.register(self.updateSpec)
        # Create Dropdown menu
        # dp = StringVar()
        self.dropProp = ttk.Combobox(
            parFrm,
            values=self.propOptions,
            # textvariable=dp,
            state="readonly",
            validate="focusin",
            validatecommand=(specVal, "%P"),
            justify="center",
        )
        self.dropProp.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropProp.current(0)
        self.dropProp.grid(
            row=4, column=1, columnspan=2, sticky="nsew", padx=5, pady=5
        )

        self.dropProp.configure(width=20)

        specScroll = ttk.Scrollbar(parFrm, orient="vertical")
        specScroll.grid(row=5, column=2, sticky="nsew", padx=5, pady=5)
        self.specs = Text(
            parFrm,
            wrap=WORD,
            height=4,
            width=20,  # yscrollcommand=specScroll
        )
        self.specs.grid(row=5, column=1, sticky="nsew", padx=5, pady=5)
        self.updateSpec()

        ttk.Label(parFrm, text="Grain Shape").grid(
            row=6,
            column=0,
            sticky="nsew",
            padx=5,
            pady=5,
        )
        # Create Dropdown menu
        self.dropGeom = ttk.Combobox(
            parFrm, values=self.geoOptions, state="readonly", justify="center"
        )
        self.dropGeom.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropGeom.current(0)
        self.dropGeom.grid(
            row=6, column=1, columnspan=2, sticky="nsew", padx=5, pady=5
        )
        self.dropGeom.configure(width=20)

        parFrm.rowconfigure(4, weight=0)
        parFrm.rowconfigure(5, weight=0)
        parFrm.rowconfigure(6, weight=0)
        i += 3

        self.webmm, _, i = self.add3Input(
            parFrm, i, "Web Thickness", "mm", "0.0", validationNN
        )
        self.permm, self.perw, i = self.add3Input(
            parFrm, i, "Perf Diameter", "mm", "0.0", validationNN
        )
        self.grlR, _, i = self.add3Input(
            parFrm, i, "Grain L/D", "", "2.25", validationNN
        )

        self.ldf, _, i = self.add3Input(
            parFrm, i, "Load Factor", "%", "50.0", validationNN
        )
        self.clr, _, i = self.add3Input(
            parFrm, i, "Chamber L. Ratio", "", "1.1", validationNN
        )
        self.stpMPa, _, i = self.add3Input(
            parFrm, i, "Start Pressure", "MPa", "30", validationNN
        )

    def addOpsFrm(self, parent):
        opFrm = ttk.LabelFrame(parent, text="Options")
        opFrm.grid(row=1, column=2, sticky="nsew")

        opFrm.columnconfigure(0, weight=1)
        opFrm.columnconfigure(1, weight=1)
        opFrm.columnconfigure(1, weight=1)

        validationNN = parent.register(validateNN)
        i = 0
        self.webR, webRw, i = self.add2Input(
            opFrm, i, 0, "W.Th. ", "2.0", validation=validationNN
        )
        self.perR, perRw, i = self.add2Input(
            opFrm,
            i,
            0,
            "P.Dia. ",
            "1.0",
            validation=validationNN,
        )

        ratioEntrys = (webRw, perRw)
        directEntrys = (self.perw,)

        def configEntry():
            if self.geoLocked.get() == 0:
                for en in ratioEntrys:
                    en.config(state="disabled")
                for en in directEntrys:
                    en.config(state="normal")
            else:
                for en in ratioEntrys:
                    en.config(state="normal")
                for en in directEntrys:
                    en.config(state="disabled")

            self.webcallback(None, None, None)

        configEntry()

        self.webmm.trace_add("write", self.webcallback)
        self.webR.trace_add("write", self.webcallback)
        self.perR.trace_add("write", self.webcallback)

        ttk.Checkbutton(
            opFrm,
            text="Lock Geometry",
            variable=self.geoLocked,
            onvalue=1,
            offvalue=0,
            command=configEntry,
        ).grid(row=0, column=2, rowspan=2, sticky="nsew", padx=5, pady=5)

        ttk.Label(opFrm, text="∫ in").grid(
            row=2, column=0, sticky="nsew", padx=5, pady=5
        )

        self.dropOptn = ttk.Combobox(
            opFrm, values=self.domainOptions, state="readonly", justify="center"
        )
        self.dropOptn.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropOptn.current(0)
        self.dropOptn.grid(
            row=2, column=1, columnspan=1, sticky="nsew", padx=5, pady=5
        )

        self.dropOptn.configure(width=5)

        ttk.Label(opFrm, text="domain").grid(
            row=2, column=2, sticky="nsew", padx=5, pady=5
        )

        validationPI = parent.register(validatePI)

        self.steps, _, _ = self.add3Input(
            opFrm,
            3,
            "Show",
            "steps",
            "50",
            validation=validationPI,
            formatter=formatIntInput,
        )

        self.accExp, _, _ = self.add2Input(
            opFrm,
            4,
            0,
            "-log10(ε) =",
            default="5",
            validation=validationPI,
            formatter=formatIntInput,
            color="red",
        )

        ttk.Button(opFrm, text="Calculate", command=self.calculate).grid(
            row=5, column=0, columnspan=3, sticky="nsew", padx=5, pady=5
        )

        opFrm.rowconfigure(5, weight=0)

    def addErrFrm(self, parent):
        errorFrm = ttk.LabelFrame(parent, text="Errors")
        errorFrm.grid(row=1, column=1, sticky="nsew")
        errorFrm.columnconfigure(0, weight=1)
        errorFrm.rowconfigure(0, weight=1)

        errScroll = ttk.Scrollbar(errorFrm, orient="vertical")
        errScroll.grid(row=0, column=1, sticky="nsew")
        self.errorText = Text(
            errorFrm,
            yscrollcommand=errScroll.set,
            wrap=WORD,
            height=0,
            width=75,
        )
        self.errorText.grid(row=0, column=0, sticky="nsew")

    def addTblFrm(self, parent):
        columnList = [
            "Event",
            "Time/s",
            "Travel/m",
            "Burnup/1",
            "Velocity/ms^-1",
            "Pressure/Pa",
        ]
        tblFrm = ttk.LabelFrame(parent, text="Result Table")
        tblFrm.grid(row=0, column=1, sticky="nsew")
        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(0, weight=1)

        self.tv = ttk.Treeview(tblFrm, selectmode="browse")
        self.tv.grid(row=0, column=0, sticky="nsew")

        self.tv["columns"] = columnList
        self.tv["show"] = "headings"
        self.tv.tag_configure("PEAK PRESSURE", foreground="orange")
        self.tv.tag_configure("BURNOUT", foreground="red")
        self.tv.tag_configure("FRACTURE", foreground="brown")
        # self.tv.tag_configure("monospace", font=("TkFixedFont", 9))
        self.tv.tag_configure("monospace", font=("Hack", 8))

        for column in columnList:  # foreach column
            self.tv.heading(
                column, text=column
            )  # let the column heading = column name
            self.tv.column(column, stretch=1, width=0, anchor="w")

        vertscroll = ttk.Scrollbar(
            tblFrm, orient="vertical"
        )  # create a scrollbar
        vertscroll.configure(command=self.tv.yview)  # make it vertical
        vertscroll.grid(row=0, column=1, sticky="nsew")

        horzscroll = ttk.Scrollbar(tblFrm, orient="horizontal")
        horzscroll.configure(command=self.tv.xview)
        horzscroll.grid(row=1, column=0, sticky="nsew")
        self.tv.configure(
            yscrollcommand=vertscroll.set, xscrollcommand=horzscroll.set
        )  # assign the scrollbar to the Treeview Widget

    def add2Input(
        self,
        parent,
        rowIndex,
        colIndex,
        labelText,
        default="1.0",
        validation=None,
        entryWidth=5,
        formatter=formatFloatInput,
        color=None,
    ):
        ttk.Label(parent, text=labelText).grid(
            row=rowIndex, column=colIndex, sticky="nsew", padx=5, pady=5
        )
        parent.rowconfigure(rowIndex, weight=0)
        e = StringVar(parent)
        e.set(default)
        en = ttk.Entry(
            parent,
            textvariable=e,
            validate="key",
            validatecommand=(validation, "%P"),
            width=entryWidth,
            foreground=color,
            justify="center",
        )
        en.default = default
        en.grid(
            row=rowIndex, column=colIndex + 1, sticky="nsew", padx=5, pady=5
        )
        en.bind("<FocusOut>", formatter)
        return e, en, rowIndex + 1

    def add3Input(
        self,
        parent,
        rowIndex,
        labelText,
        unitText,
        default="0.0",
        validation=None,
        entryWidth=5,
        formatter=formatFloatInput,
    ):
        ttk.Label(parent, text=labelText).grid(
            row=rowIndex, column=0, sticky="nsew", padx=5, pady=5
        )
        parent.rowconfigure(rowIndex, weight=0)
        e = StringVar(parent)
        e.set(default)
        en = ttk.Entry(
            parent,
            textvariable=e,
            validate="key",
            validatecommand=(validation, "%P"),
            width=entryWidth,
            justify="center",
        )
        en.default = default
        en.grid(row=rowIndex, column=1, sticky="nsew", padx=5, pady=5)
        en.bind("<FocusOut>", formatter)
        ttk.Label(parent, text=unitText).grid(
            row=rowIndex, column=2, sticky="nsew", padx=5, pady=5
        )
        return e, en, rowIndex + 1

    def add12Disp(
        self,
        parent,
        rowIndex,
        labelText,
        unitText,
        default="0.0",
        entryWidth=10,
    ):
        ttk.Label(parent, text=labelText).grid(
            row=rowIndex, column=0, columnspan=2, sticky="nsew", padx=5, pady=5
        )
        e = StringVar(parent)
        e.default = default
        e.set(default)
        parent.rowconfigure(rowIndex, weight=0)
        en = ttk.Entry(
            parent, textvariable=e, width=entryWidth, state="disabled"
        )
        en.grid(row=rowIndex + 1, column=0, sticky="nsew", padx=5, pady=5)
        ttk.Label(parent, text=unitText).grid(
            row=rowIndex + 1, column=1, sticky="nsew", padx=5, pady=5
        )
        return e, en, rowIndex + 2

    def webcallback(self, var, index, mode):
        if any(
            (
                self.webmm.get() == "",
                self.webR.get() == "",
                self.perR.get() == "",
            )
        ):
            pass
        else:
            if self.geoLocked.get() == 1:
                try:
                    self.geomError = False
                    self.permm.set(
                        float(self.webmm.get())
                        / float(self.webR.get())
                        * float(self.perR.get())
                    )
                except ZeroDivisionError:
                    self.geomError = True


if __name__ == "__main__":
    # high dpi scaling
    windll.shcore.SetProcessDpiAwareness(1)
    root = Tk()
    # one must supply the entire path
    loadfont(resolvepath("ui/Hack-Regular.ttf"))
    # tksvg.load(root)
    dpi = root.winfo_fpixels("1i")

    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.12"))
    root.tk.call("package", "require", "awdark")
    style = ttk.Style(root)
    style.theme_use("awdark")
    # ensure that the treeview rows are roughly the same height
    # regardless of dpi. on Windows, default is Segoe UI at 9 points
    style.configure("Treeview", rowheight=round(12 * dpi / 72))
    style.configure("Treeview.Heading", font=("Hack", 9))
    style.configure("TButton", font=("Hack", 9))
    style.configure("TLabelframe.Label", font=("Hack", 10))
    root.option_add("*Font", "Hack 9")
    ibPanel = IB(root)
    root.mainloop()
