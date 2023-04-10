from tkinter import *
from tkinter import ttk
import traceback
import tksvg
from gun import *
import os
import sys

DEBUG = True

_prefix = {
    "y": 1e-24,  # yocto
    "z": 1e-21,  # zepto
    "a": 1e-18,  # atto
    "f": 1e-15,  # femto
    "p": 1e-12,  # pico
    "n": 1e-9,  # nano
    "μ": 1e-6,  # micro
    "m": 1e-3,  # mili
    " ": 1,  # unit
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


def toSI(v, dec=4, unit=None):
    for prefix, magnitude in zip(_prefix.keys(), _prefix.values()):
        if 1 <= (v / magnitude) < 1e3:
            vstr = "{:#.{:}g}".format(v / magnitude, dec)
            return (
                vstr
                + " " * (dec + 1 - len(vstr) + vstr.find("."))
                + prefix
                + (unit if unit is not None else "")
            )
    if v == 0:
        return (
            "{:#.{:}g}".format(v, dec)
            + "  "
            + (unit if unit is not None else "")
        )
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


def dot_aligned(matrix, units):
    transposed = []

    for seq, unit in zip(zip(*matrix), units):
        snums = []
        for n in seq:
            try:
                snums.append(toSI(float(n), unit=unit))
            except ValueError:
                snums.append(n)
        dots = [s.find(".") for s in snums]
        m = max(dots)
        transposed.append(tuple(" " * (m - d) + s for s, d in zip(snums, dots)))

    return tuple(zip(*transposed))


class ToolTip(object):
    """solution found in
    https://stackoverflow.com/questions/20399243/display-message-when-hovering-over-something-with-mouse-cursor-in-python
    """

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        """Display text in tooltip window"""
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx()  # + self.widget.winfo_width()
        y = y + cy + self.widget.winfo_rooty() + self.widget.winfo_height()
        """ initalize the tooltip window to the lower right corner of the widget"""
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(
            tw,
            text=self.text,
            justify=LEFT,
            background="#ffffe0",
            relief=SOLID,
            borderwidth=1,
            font=("tahoma", "8", "normal"),
        )
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()


def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)

    def enter(event):
        toolTip.showtip(text)

    def leave(event):
        toolTip.hidetip()

    widget.bind("<Enter>", enter)
    widget.bind("<Leave>", leave)


class IB(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.compositions = GrainComp.readFile(
            resolvepath("data/propellants.csv")
        )
        self.geometries = {i.desc: i for i in Geometry}

        self.prop = None
        self.gun = None
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
        self.addExcFrm(parent)
        self.addTblFrm(parent)
        self.addOpsFrm(parent)

        parent.bind("<Return>", self.calculate)

    def calculate(self, event=None):
        # force an immediate redraw after calculation
        compo = self.compositions[self.dropProp.get()]
        geom = self.geometries[self.dropGeom.get()]

        self.tableData = []

        self.resetSpec()

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
            if DEBUG:
                self.errorLst.append("".join(traceback.format_exception(e)))
            else:
                self.errorLst.append(str(e))

        try:
            chamberVolume = (
                float(self.chgkg.get())
                / self.prop.rho_p
                / self.prop.maxLF
                / float(self.ldf.get())
                * 100
            )
            self.cv.set(toSI(chamberVolume))
            if DEBUG:
                print(
                    *(
                        float(self.calmm.get()) / 1000,
                        float(self.shtkg.get()),
                        self.prop,
                        float(self.chgkg.get()),
                        chamberVolume,
                        float(self.stpMPa.get()) * 1e6,
                        float(self.tblmm.get()) / 1000,
                        float(self.clr.get()),
                    ),
                    sep=","
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

            self.va.set(toSI(self.gun.v_j))

            t_err, l_err, psi_e, v_err, p_err = self.gun.getErr(
                10 ** -(int(self.accExp.get()))
            )

            self.terr.set(toSI(t_err))
            self.lerr.set(toSI(l_err))
            self.psierr.set(toSI(psi_e))
            self.verr.set(toSI(v_err))
            self.perr.set(toSI(p_err))

        except Exception as e:
            self.gun = None
            self.errorLst.append("Exception when defining guns:")
            if DEBUG:
                self.errorLst.append("".join(traceback.format_exception(e)))
            else:
                self.errorLst.append(str(e))

        if self.gun is not None:
            try:
                self.tableData = self.gun.integrate(
                    steps=int(self.steps.get()),
                    dom=self.dropOptn.get(),
                    tol=10 ** -(int(self.accExp.get())),
                )
                i = tuple(i[0] for i in self.tableData).index("SHOT EXIT")
                vg = self.tableData[i][4]
                te, be = self.gun.getEff(vg)
                self.te.set(round(te * 100, 1))
                self.be.set(round(te / self.gun.phi * 100, 1))
            except Exception as e:
                self.errorLst.append("Exception while solving numerically:")
                if DEBUG:
                    self.errorLst.append("".join(traceback.format_exception(e)))
                else:
                    self.errorLst.append(str(e))

        self.tv.delete(*self.tv.get_children())
        self.tableData = dot_aligned(
            self.tableData, units=(None, "s", "m", None, "m/s", "Pa")
        )

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
        specFrm.columnconfigure(0, weight=1)
        i = 0
        self.va, _, i = self.add12Disp(
            specFrm, i, "Asymptotic Vel.", "m/s", justify="right"
        )
        self.te, _, i = self.add12Disp(specFrm, i, "Thermal Eff.", "%")
        self.be, _, i = self.add12Disp(specFrm, i, "Ballistic Eff.", "%")
        self.cv, _, i = self.add12Disp(
            specFrm, i, "Chamber Volume", "m^3", justify="right"
        )

        specFrm.rowconfigure(i, weight=1)
        errFrm = ttk.LabelFrame(specFrm, text="Maximum Error")
        errFrm.grid(row=i, column=0, columnspan=2, sticky="sew")
        errFrm.columnconfigure(0, weight=1)

        j = 0
        self.terr, _, j = self.add12Disp(
            errFrm, j, "Time", "s", justify="right"
        )
        self.lerr, _, j = self.add12Disp(
            errFrm, j, "Travel", "m", justify="right"
        )
        self.psierr, _, j = self.add12Disp(
            errFrm, j, "Burnup", "", justify="right"
        )
        self.verr, _, j = self.add12Disp(
            errFrm, j, "Velocity", "m/s", justify="right"
        )
        self.perr, _, j = self.add12Disp(
            errFrm, j, "Pressure", "Pa", justify="right"
        )

        self.specFrmDisps = (
            self.va,
            self.te,
            self.be,
            self.cv,
            self.terr,
            self.lerr,
            self.psierr,
            self.verr,
            self.perr,
        )

    def resetSpec(self):
        for specDisp in self.specFrmDisps:
            specDisp.set(specDisp.default)

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

        specVal = parent.register(self.updateSpec)

        parFrm.rowconfigure(i, weight=1)

        propFrm = ttk.LabelFrame(parFrm, text="Propellant")
        propFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        self.dropProp = ttk.Combobox(
            propFrm,
            values=self.propOptions,
            # textvariable=dp,
            state="readonly",
            validate="focusin",
            validatecommand=(specVal, "%P"),
            justify="center",
        )
        self.dropProp.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropProp.current(0)
        self.dropProp.grid(row=0, column=0, columnspan=2, sticky="nsew", pady=2)
        self.dropProp.configure(width=10)

        specScroll = ttk.Scrollbar(propFrm, orient="vertical")
        specScroll.grid(
            row=1,
            column=1,
            sticky="nsew",
            pady=2,
        )
        self.specs = Text(
            propFrm,
            wrap=WORD,
            height=4,
            width=10,
            yscrollcommand=specScroll.set,
        )
        self.specs.grid(row=1, column=0, sticky="nsew", pady=2)
        specScroll.config(command=self.specs.yview)

        propFrm.rowconfigure(1, weight=1)
        propFrm.columnconfigure(0, weight=1)

        self.updateSpec()

        i += 1

        ttk.Label(parFrm, text="Grain Shape").grid(
            row=i,
            column=0,
            sticky="nsew",
            padx=2,
            pady=2,
        )

        i += 1

        # Create Dropdown menu
        self.dropGeom = ttk.Combobox(
            parFrm, values=self.geoOptions, state="readonly", justify="center"
        )
        self.dropGeom.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropGeom.current(0)
        self.dropGeom.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        self.dropGeom.configure(width=10)

        i += 1

        self.webmm, _, i = self.add3Input(
            parFrm, i, "Web Thickness", "mm", "0.0", validationNN
        )
        self.permm, self.perw, i = self.add3Input(
            parFrm, i, "Perf Diameter", "mm", "0.0", validationNN
        )

        grlRtext = "\n".join(
            ("Length to diameter ratio of the", "propellant grain.")
        )

        self.grlR, _, i = self.add3Input(
            parFrm, i, "Grain L/D", "", "3.0", validationNN, infotext=grlRtext
        )

        ldftext = "\n".join(
            (
                "Percentage of chamber volume filled by",
                "the outlines of the grain. Value of",
                "0-100 % (not inclusive) are supported,",
                "although it is highly unlikely that",
                "values greater than 60% could be",
                "achieved in practice, due to packing",
                "behaviour when loading the cartridge,",
                "with realistic grain sizes and geometry.",
                "A high value is also undesirable for",
                "causing excessive peak pressure as well",
                "as moving the pressure spike closer to",
                "the breech.",
            )
        )

        self.ldf, _, i = self.add3Input(
            parFrm,
            i,
            "Load Factor",
            "%",
            "50.0",
            validationNN,
            infotext=ldftext,
        )

        clrtext = "\n".join(
            (
                "Chamber length ratio is the ratio between",
                "the length of reduced chamber (dividing",
                "the chamber volume with barrel cross",
                "section) and the actual chamber.",
            )
        )

        self.clr, _, i = self.add3Input(
            parFrm, i, "Chamber L.R.", "", "1.1", validationNN, infotext=clrtext
        )

        stpText = "\n".join(
            (
                "Peak pressure that initially resists",
                "the movement of shot. For rifled weapons",
                "this is predominantely caused by the",
                "shot or band being squeezed in to contact",
                "with the rifling. For large caliber guns",
                "25-30MPa is common. For rifles, 35-45MPa",
                "is common.",
                "For smoothbore weapons, this is instead",
                "caused by cartridge holding on to the",
                "projectile, greatly varying depending",
                "on desired rate of fire and shot weight.",
            )
        )

        self.stpMPa, _, i = self.add3Input(
            parFrm,
            i,
            "Start Pressure",
            "MPa",
            "10",
            validationNN,
            infotext=stpText,
        )

    def addOpsFrm(self, parent):
        opFrm = ttk.LabelFrame(parent, text="Options")
        opFrm.grid(row=1, column=2, sticky="nsew")

        opFrm.columnconfigure(1, weight=1)

        validationNN = parent.register(validateNN)
        i = 0
        self.webR, webRw, i = self.add2Input(
            opFrm, i, 0, "W.Th.", "2.0", validation=validationNN
        )
        self.perR, perRw, i = self.add2Input(
            opFrm,
            i,
            0,
            "P.Dia.",
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
            text="Lock\nGeom",
            variable=self.geoLocked,
            onvalue=1,
            offvalue=0,
            command=configEntry,
        ).grid(row=0, column=2, rowspan=2, sticky="nsew", padx=2, pady=2)

        ttk.Label(opFrm, text="∫ in").grid(
            row=2, column=0, sticky="nsew", padx=2, pady=2
        )

        self.dropOptn = ttk.Combobox(
            opFrm, values=self.domainOptions, state="readonly", justify="center"
        )
        self.dropOptn.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropOptn.current(0)
        self.dropOptn.grid(
            row=2, column=1, columnspan=1, sticky="nsew", padx=2, pady=2
        )
        self.dropOptn.configure(width=5)

        ttk.Label(opFrm, text="domain").grid(
            row=2, column=2, sticky="nsew", padx=2, pady=2
        )

        validationPI = parent.register(validatePI)

        self.steps, _, _ = self.add3Input(
            opFrm,
            3,
            "Show",
            "steps",
            "10",
            validation=validationPI,
            formatter=formatIntInput,
        )

        self.accExp, _, _ = self.add2Input(
            opFrm,
            4,
            0,
            "-log10(ε) ",
            default="5",
            validation=validationPI,
            formatter=formatIntInput,
            color="red",
        )

        ttk.Button(
            opFrm, text="Calculate", command=self.calculate, underline=0
        ).grid(row=5, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        opFrm.rowconfigure(5, weight=0)

    def addExcFrm(self, parent):
        errorFrm = ttk.LabelFrame(parent, text="Exceptions")
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
            "Time",
            "Travel",
            "Burnup",
            "Velocity",
            "Pressure",
        ]
        tblFrm = ttk.LabelFrame(parent, text="Result Table")
        tblFrm.grid(row=0, column=1, sticky="nsew")

        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(1, weight=1)

        # configure the numerical

        self.tv = ttk.Treeview(tblFrm, selectmode="browse")
        self.tv.grid(row=1, column=0, sticky="nsew")

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
        vertscroll.grid(row=1, column=1, sticky="nsew")

        horzscroll = ttk.Scrollbar(tblFrm, orient="horizontal")
        horzscroll.configure(command=self.tv.xview)
        horzscroll.grid(row=2, column=0, sticky="nsew")
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
            row=rowIndex, column=colIndex, sticky="nsew", padx=2, pady=2
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
            row=rowIndex, column=colIndex + 1, sticky="nsew", padx=2, pady=2
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
        entryWidth=10,
        formatter=formatFloatInput,
        infotext=None,
    ):
        lb = ttk.Label(parent, text=labelText)
        lb.grid(row=rowIndex, column=0, sticky="nsew", padx=2, pady=2)
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
        en.grid(row=rowIndex, column=1, sticky="nsew", padx=2, pady=2)
        en.bind("<FocusOut>", formatter)
        ttk.Label(parent, text=unitText).grid(
            row=rowIndex, column=2, sticky="nsew", padx=2, pady=2
        )
        if infotext is not None:
            CreateToolTip(lb, infotext)
        return e, en, rowIndex + 1

    def add12Disp(
        self,
        parent,
        rowIndex,
        labelText,
        unitText,
        default="0.0",
        entryWidth=5,
        justify="center",
    ):
        ttk.Label(parent, text=labelText).grid(
            row=rowIndex, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )
        e = StringVar(parent)
        e.default = default
        e.set(default)
        parent.rowconfigure(rowIndex, weight=0)
        en = ttk.Entry(
            parent,
            textvariable=e,
            width=entryWidth,
            state="disabled",
            justify=justify,
        )
        en.grid(row=rowIndex + 1, column=0, sticky="nsew", padx=2, pady=2)
        ttk.Label(parent, text=unitText).grid(
            row=rowIndex + 1, column=1, sticky="nsew", padx=2, pady=2
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
    tksvg.load(root)
    dpi = root.winfo_fpixels("1i")

    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    # root.tk.call("package", "require", "awdark")

    style = ttk.Style(root)
    style.theme_use("awdark")

    # ensure that the treeview rows are roughly the same height
    # regardless of dpi. on Windows, default is Segoe UI at 9 points
    style.configure("Treeview", rowheight=round(12 * dpi / 72))
    style.configure("Treeview.Heading", font=("Hack", 8))
    style.configure("TButton", font=("Hack", 9, "bold"))
    style.configure("TLabelframe.Label", font=("Hack", 9, "bold"))
    style.configure("TNotebook.Tab", font=("Hack", 8))
    root.option_add("*Font", "Hack 8")
    # root.option_add("*TCombobox*Listbox*Font", "Hack 8")

    root.title("Phoenix's Internal Ballistics Solver v0.2")

    ibPanel = IB(root)
    root.mainloop()
