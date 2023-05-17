from tkinter import *
from tkinter import ttk
import tkinter.font as tkFont
import traceback
from gun import *
import os
import sys
import ctypes
from math import ceil
from labellines import labelLine, labelLines

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk,
)

# Implement the default Matplotlib key bindings.
import matplotlib.pyplot as mpl
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

GEOM_CONTEXT = {
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.titlesize": 10,
    "axes.edgecolor": "white",
    "axes.facecolor": "#33393b",
    "figure.facecolor": "#33393b",
    "text.color": "white",
    "xtick.color": "white",
    "ytick.color": "white",
    "lines.markersize": 2,
}

FIG_CONTEXT = {
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 12,
    "lines.linewidth": 1,
    "axes.edgecolor": "white",
    "axes.facecolor": "#33393b",
    "axes.labelcolor": "white",
    "figure.facecolor": "#33393b",
    "text.color": "white",
    "xtick.color": "white",
    "font.weight": "bold",
    "ytick.color": "white",
    "lines.markersize": 4,
}
chgText = " ".join(
    (
        "Mass of propellant charge to be used. Chamber",
        "volume is determined from this and load fraction.",
    )
)
vinfText = " ".join(
    (
        "Velocity the shot would achieve if",
        "the barrel is extended to infinite length.",
    )
)
teffText = " ".join(
    (
        "Thermal efficiency of the gun system, i.e.",
        "the amount of work done to both gas and",
        "projectile over potential from propellant",
    )
)
beffText = " ".join(
    (
        "Ballistic efficiency of the gun system, i.e.",
        "the amount of work done the projectile over",
        "chemical potential energy of propellant.",
    )
)

specsText = " ".join(
    (
        "Specify the propellant composition used.\n",
        "Pure Nitrocellulose is plasticized with",
        "Dinitrotoluene to reduce burn rate and lower",
        "flame temperature, forming single based",
        "propellant.\n",
        "Double based propellant is formed when",
        "Nitroglycerin is used as gelatinizer instead.",
        "While more energetic, it also burns hotter",
        "and erodes barrel more.\n",
        "Triple base propellants contains Nitroguanidine,",
        "with higher energy content",
        "while keeping flame temperature low.",
        "However, it is mechanically the weakest.",
    )
)
geomPlotTxt = " ".join(
    (
        "Plot of σ, or burn surface area (unitless) with",
        "respect to Z, or the linear burnt ratio. A upward",
        "slope on this graph indicate progressive burning",
        "while a downward slope indicate regressive burning.",
        "Discontinuity indicate fracture of mutli-perforated",
        "propellant.\n",
        "Regressive burning puts the peak pressure point closer",
        "to the breech, leading to higher peak pressure, and",
        "vice versa for progressive burning.\n",
        "Note, the treatment for post fracture",
        "burn behaviour is not entirely rigorous.",
    )
)
ldftext = " ".join(
    (
        "Percentage of chamber volume filled by",
        "the outlines of the grain, also known",
        "as 'packing density'. Value of 0-100 % ",
        "(not inclusive) are supported, although",
        "it is highly unlikely that values greater",
        "than 60% could be achieved in practice,",
        "due to packing behaviour when loading the",
        "cartridg.\n",
        "A high value is also undesirable for",
        "causing excessive peak pressure as well",
        "as moving the pressure spike closer to",
        "the breech.\n",
        "Only internal voids are accounted for,",
        "thus this is not simply the usually quoted",
        "'load fraction' which is by weight.",
    )
)
arcText = " ".join(
    (
        "Specify the geometry of propellant grain",
        "using dimensions. Also known as the web.",
        "All other geometries are scaled by this.\n",
        "Burning radially outward along the arc is",
        "the shortest path, and therefore is the",
        "primary geometrical determinant of burn rapidity.\n",
        "In theory micrometer level precision",
        "is possible. In practice, tolerance for industrial",
        "bulk production is in the range of 1μm - 0.15mm",
        " - 1mm depending on sources.\n",
        "Arc thickness is generally found close to 1mm",
        "for small to intermediate calibers.",
    )
)

pDiaRText = " ".join(
    (
        "Specify diameter of perforation over arc width.\n",
        "Perforations are formed by protrusions in the",
        "copper casting die. Standard multi-perf grain tends",
        "to come in the range of 0.5-1, whereas longer, single",
        "perf, or tubular grains come to around 1.33 for the same.",
    )
)

diaText = " ".join(("Specify the diameter of the propellant grain.",))

perfLRtext = " ".join(
    (
        "Specify length to diameter ratio of the grain, for",
        "mutli perforated grains this is usually in the range of"
        " 1.82 to 3.57. A higher value facilitate progressive",
        "burning. If this is too low its possible for",
        "propellant to exhibit regressive burning behaviour,",
        "as the primary mode will be axial burning instead.",
        "Longer grains tends to suffer from premature fracture",
        "due to mechanical stress.",
    )
)

cylLRtext = " ".join(
    (
        "Specify length to diameter ratio of the grain.",
        "Cylindrical or tubular propellant can be made rather long",
    )
)

rodRtext = " ".join(
    (
        "Specify the length to width ratio of propellant rod or flake.",
        "Can be quite long for tape like propellant.",
    )
)

widthText = " ".join(
    (
        "Specify the width of propellant rod or flake.\n",
        "This value is used to scale all other dimension",
        "for rod or flake like propellants. Does not have",
        "to be the smallest dimension.",
    )
)
heightRtext = " ".join(
    ("Specify the height to width ratio of propellant rod or flake.",)
)


tolText = " ".join(
    (
        "The maximum relative error, ε, that is allowed in the integrands",
        "for each component. Some components may have significantly less",
        "error than specified here, shown under each entry in the table.",
    )
)
stpText = " ".join(
    (
        "Peak pressure that initially resists",
        "the movement of shot. This is made up of",
        "the rifling resisting the drive band or",
        "shell body, the static friction of the",
        "barrel, and the cartridge case holding onto",
        "the shot (greatly varying between 0.25 - 15MPa",
        "depending on caliber and desired RoF)\n",
        "For large caliber guns 25 - 30MPa is common.",
        "For rifles, reports range from 35 to 45MPa.",
    )
)
clrtext = " ".join(
    (
        "Chamber length ratio is the ratio between",
        "the length of reduced chamber (dividing",
        "the chamber volume with barrel cross",
        "section) to the actual chamber. Accounts",
        "for the necking of the cartridge.",
    )
)
dgctext = " ".join(
    (
        "Drag coefficient, or the pressure induced by barrel",
        "friction divided by the shot base pressure. Currently",
        "2%-7% is reported for rifled weapons, with smaller",
        "calibers on the higher end and larger caliber shots",
        "with driving band on the lower end.",
    )
)
sampTxt = " ".join(
    (
        "Samples are taken equidistantly along specified domain.",
        "Sampling is done after the system has been solved and thus",
        "does not influence the accuracy of calculation in anyway.",
    )
)

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
    # "c": 1e-2,  # centi
    # "d": 1e-1,  # deci
    " ": 1,  # unit
    # "da": 1e1, # deca
    # "h": 1e2,  # hecto
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


def toSI(v, dec=4, unit=None, useSN=False):
    if v >= 0:
        positive = True
    else:
        positive = False
    v = abs(v)
    for prefix, magnitude, nextMagnitude in zip(
        _prefix.keys(),
        tuple(_prefix.values())[:-1],
        tuple(_prefix.values())[1:],
    ):
        if 1 <= (v / magnitude) < (nextMagnitude / magnitude):
            # allows handling of non-uniformly log10 spaced prefixes
            vstr = "{:#.{:}g}".format(v / magnitude, dec)
            return (
                (" " if positive else "-")
                + vstr
                + " " * (dec + 1 - len(vstr) + vstr.find("."))
                + (
                    (
                        "E{:<3}".format(round(log(magnitude, 10)))
                        if magnitude != 1
                        else "    "  # 4 SPACES!
                    )
                    if useSN
                    else prefix
                )
                + (unit if unit is not None else "")
            )
    if v == 0:
        return (
            (" " if positive else "-")
            + "{:#.{:}g}".format(v, dec)
            + ("     " if useSN else "  ")
            + (unit if unit is not None else "")
        )
    else:  # return a result in SI as a last resort
        closest = log(v, 10) // 3
        magnitude = 10 ** (closest * 3)
        vstr = "{:#.{:}g}".format(v / magnitude, dec)
        return (
            (" " if positive else "-")
            + vstr
            + " " * (dec + 1 - len(vstr) + vstr.find("."))
            + ("E{:<3}".format(round(log(magnitude, 10))))
            + (unit if unit is not None else "")
        )


def validateNN(inp):
    """
    validate an input if it results in:
    - result >=0
    - result is empty
    in the latter case, the empty field will be handled by tracing
    change in variable.
    """
    if inp == "" or inp == ".":
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
    if v == "" or v == ".":
        event.widget.delete(0, END)
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


def dot_aligned(matrix, units, useSN):
    transposed = []

    for seq, unit, isSN in zip(zip(*matrix), units, useSN):
        snums = []
        for n in seq:
            try:
                snums.append(toSI(float(n), unit=unit, useSN=isSN))
            except ValueError:
                snums.append(n)
        dots = [s.find(".") for s in snums]
        m = max(dots)
        transposed.append(tuple(" " * (m - d) + s for s, d in zip(snums, dots)))

    return tuple(zip(*transposed))


class ToolTip(object):
    """solution found in
    https://stackoverflow.com/questions/20399243/display-message
    -when-hovering-over-something-with-mouse-cursor-in-python
    """

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        """Display text in tooltip window"""
        if isinstance(text, StringVar):
            self.text = text.get()
        else:
            self.text = text
        if self.tipwindow or not self.text:
            return

        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        root = self.widget.winfo_toplevel()

        t_Font = tkFont.Font(family="hack", size=8)
        # we use a fixed width font so any char will do
        columnWidth = 40
        width, height = t_Font.measure("m"), t_Font.metrics("linespace")
        x, y, cx, cy = self.widget.bbox("insert")
        rx, ry, crx, cry = root.bbox()
        # bouding box coordinate is in regard to origin of widget/window

        """ initalize the tooltip window to the lower right corner of the widget"""

        if (
            x + self.widget.winfo_rootx()
            > root.winfo_rootx() + 0.5 * root.winfo_width()
        ):
            x = x + self.widget.winfo_rootx() - width * (columnWidth + 1)
            y = y + self.widget.winfo_rooty()
        else:
            x = x + self.widget.winfo_rootx() + self.widget.winfo_width()
            y = y + self.widget.winfo_rooty()

        margin = (
            y
            + ceil(t_Font.measure(self.text) / (width * columnWidth) + 1)
            * height
            - ry
            - root.winfo_rooty()
            - root.winfo_height()
        )  # ensure that tooltip does not overrun the main window

        if margin > 0:
            y -= margin

        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(
            tw,
            text=self.text,
            justify=LEFT,
            background="#ffffe0",
            wraplength=width * columnWidth,
            relief=SOLID,
            borderwidth=1,
            font=t_Font,
        )
        label.config(width=columnWidth)  # characters
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.after(25, lambda: tw.destroy())


def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)

    def enter(event):
        toolTip.showtip(text)

    def leave(event):
        toolTip.hidetip()

    widget.bind("<Enter>", enter)
    widget.bind("<Leave>", leave)


class IB(Frame):
    def __init__(self, parent, dpi):
        Frame.__init__(self, parent)

        self.dpi = dpi

        self.compositions = GrainComp.readFile(
            resolvepath("data/propellants.csv")
        )
        self.geometries = GEOMETRIES
        self.prop = None
        self.gun = None
        self.errorLst = []

        self.propOptions = tuple(self.compositions.keys())
        self.geoOptions = tuple(self.geometries.keys())
        self.domainOptions = (DOMAIN_TIME, DOMAIN_LENG)

        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)

        self.addRightFrm(parent)
        self.addErrFrm(parent)
        self.addPlotFrm(parent)
        self.addParFrm(parent)
        self.addTblFrm(parent)

        parent.update_idletasks()
        self.addGeomPlot()

        parent.update_idletasks()
        self.addFigPlot()

        self.updateSpec(None, None, None)
        self.updateGeom(None, None, None)

        # parent.bind("<Return>", self.calculate)

        parent.bind("<Configure>", self.resizeFigPlot)

    def calculate(self, event=None):
        self.intgRecord = []

        # force an immediate redraw after calculation
        compo = self.compositions[self.dropProp.get()]
        # lookup dictionary using the string key to get
        # the requisite object
        geom = self.geometries[self.dropGeom.get()]

        if self.prop is None:
            return

        try:
            chamberVolume = (
                float(self.chgkg.get())
                / self.prop.rho_p
                / self.prop.maxLF
                / float(self.ldf.get())
                * 100
            )
            self.cv.set(toSI(chamberVolume, useSN=True).strip())
            self.ld.set(round(self.prop.maxLF * float(self.ldf.get()), 1))

            if DEBUG:
                print(
                    *(
                        float(self.calmm.get()) * 1e-3,
                        float(self.shtkg.get()),
                        self.prop,
                        float(self.chgkg.get()),
                        float(self.arcmm.get()) * 1e-3,
                        chamberVolume,
                        float(self.stpMPa.get()) * 1e6,
                        float(self.tblmm.get()) * 1e-3,
                        float(self.clr.get()),
                    ),
                    sep=","
                )
            self.gun = Gun(
                caliber=float(self.calmm.get()) * 1e-3,
                shotMass=float(self.shtkg.get()),
                propellant=self.prop,
                grainSize=float(self.arcmm.get()) * 1e-3,
                chargeMass=float(self.chgkg.get()),
                chamberVolume=chamberVolume,
                startPressure=float(self.stpMPa.get()) * 1e6,
                lengthGun=float(self.tblmm.get()) * 1e-3,
                chamberExpansion=float(self.clr.get()),
                dragCoe=float(self.dgc.get()) * 1e-2,
            )

            self.va.set(toSI(self.gun.v_j))

        except Exception as e:
            self.gun = None
            self.errorLst.append("Exception when defining guns:")
            if DEBUG:
                self.errorLst.append("".join(traceback.format_exception(e)))
            else:
                self.errorLst.append(str(e))

        if self.gun is not None:
            try:
                self.tableData, self.errorData = self.gun.integrate(
                    steps=int(self.steps.get()),
                    dom=self.dropOptn.get(),
                    tol=10 ** -(int(self.accExp.get())),
                    record=self.intgRecord,
                )

                i = tuple(i[0] for i in self.tableData).index("SHOT EXIT")
                vg = self.tableData[i][4]
                te, be = self.gun.getEff(vg)
                self.te.set(round(te * 100, 1))
                self.be.set(round(te / self.gun.phi * 100, 1))
            except Exception as e:
                self.gun = None
                self.tableData = []
                self.errorData = []
                self.errorLst.append("Exception while solving gun system:")
                if DEBUG:
                    self.errorLst.append("".join(traceback.format_exception(e)))
                else:
                    self.errorLst.append(str(e))

        self.tv.delete(*self.tv.get_children())
        useSN = (False, False, False, True, False, False, True)
        units = (None, "s", "m", None, "m/s", "Pa", "K")
        tableData = dot_aligned(
            self.tableData,
            units=units,
            useSN=useSN,
        )
        errorData = dot_aligned(self.errorData, units=units, useSN=useSN)
        # negErr, posErr = arrErr(self.errorData, units=units, useSN=useSN)
        i = 0
        for row, erow in zip(tableData, errorData):
            self.tv.insert(
                "", "end", str(i), values=row, tags=(row[0], "monospace")
            )
            self.tv.insert(
                str(i),
                "end",
                str(i + 1),
                values=tuple("±" + e if e != erow[0] else e for e in erow),
                tags="error",
            )

            self.tv.move(str(i + 1), str(i), "end")

            i += 2

        self.updateError()
        self.updateFigPlot()

    def addRightFrm(self, parent):
        rightFrm = ttk.Frame(parent)
        rightFrm.grid(row=0, column=2, rowspan=3, sticky="nsew")
        rightFrm.columnconfigure(0, weight=1)
        rightFrm.rowconfigure(0, weight=1)
        specFrm = ttk.LabelFrame(
            rightFrm, text="Design Summary", style="SubLabelFrame.TLabelframe"
        )
        specFrm.grid(row=0, column=0, sticky="nsew")
        specFrm.columnconfigure(0, weight=1)
        i = 0

        self.va, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText="Asymptotic Vel.",
            unitText="m/s",
            justify="right",
            infotext=vinfText,
        )

        self.te, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText="Thermal Eff.",
            unitText="%",
            infotext=teffText,
        )

        self.be, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText="Ballistic Eff.",
            unitText="%",
            infotext=beffText,
        )
        self.cv, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText="Chamber Volume",
            unitText="m³",
            justify="right",
        )
        self.ld, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText="Loading Density",
            unitText="%",
        )

        opFrm = ttk.LabelFrame(
            rightFrm, text="Forward Calc.", style="SubLabelFrame.TLabelframe"
        )
        opFrm.grid(row=1, column=0, sticky="nsew")

        # opFrm.columnconfigure(0, weight=1)
        opFrm.columnconfigure(1, weight=1)

        validationNN = parent.register(validateNN)
        validationPI = parent.register(validatePI)

        i = 0

        sampleFrm = ttk.LabelFrame(
            opFrm, text="Sampling", style="SubLabelFrame.TLabelframe"
        )
        sampleFrm.grid(
            row=i, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )

        j = 0

        sampleFrm.columnconfigure(0, weight=1)
        sampleFrm.columnconfigure(1, weight=1)
        self.dropOptn = ttk.Combobox(
            sampleFrm,
            values=self.domainOptions,
            state="readonly",
            justify="center",
        )

        self.dropOptn.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropOptn.current(0)
        self.dropOptn.grid(
            row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )
        # self.dropOptn.configure(width=0)

        j += 1

        self.steps, _, j = self.add2Input(
            parent=sampleFrm,
            rowIndex=j,
            colIndex=0,
            labelText="Steps",
            default="25",
            validation=validationNN,
            formatter=formatIntInput,
            reverse=True,
            anchor="center",
        )

        CreateToolTip(sampleFrm, sampTxt)

        i += 1

        self.accExp, _, i = self.add2Input(
            parent=opFrm,
            rowIndex=i,
            colIndex=0,
            labelText="-log10(ε)",
            default="3",
            validation=validationPI,
            formatter=formatIntInput,
            color="red",
            infotext=tolText,
        )

        calButton = ttk.Button(
            opFrm,
            text="Calculate",
            command=self.calculate,  # underline=0
        )
        calButton.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        opFrm.rowconfigure(i, weight=1)
        CreateToolTip(calButton, "Integrate system using RKF7(8) integrator")

    def addErrFrm(self, parent):
        errorFrm = ttk.LabelFrame(parent, text="Exceptions")
        errorFrm.grid(row=2, column=0, sticky="nsew")
        errorFrm.columnconfigure(0, weight=1)
        errorFrm.rowconfigure(0, weight=1)

        errScroll = ttk.Scrollbar(errorFrm, orient="vertical")
        errScroll.grid(row=0, column=1, sticky="nsew")
        self.errorText = Text(
            errorFrm,
            yscrollcommand=errScroll.set,
            wrap=WORD,
            height=6,
            width=0,
        )
        self.errorText.grid(row=0, column=0, sticky="nsew")

    def addParFrm(self, parent):
        parFrm = ttk.LabelFrame(parent, text="Parameters")
        parFrm.grid(row=0, column=1, rowspan=3, sticky="nsew")
        parFrm.columnconfigure(0, weight=1)
        # parFrm.columnconfigure(1, weight=1)
        # parFrm.columnconfigure(2, weight=1)

        # validation
        validationNN = parent.register(validateNN)

        i = 0
        self.calmm, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Caliber",
            unitText="mm",
            default="50.0",
            validation=validationNN,
        )
        self.tblmm, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Tube Length",
            unitText="mm",
            default="3500.0",
            validation=validationNN,
        )
        self.shtkg, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Shot Mass",
            unitText="kg",
            default="1.0",
            validation=validationNN,
        )
        self.chgkg, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Charge Mass",
            unitText="kg",
            default="1.0",
            validation=validationNN,
            infotext=chgText,
        )

        parFrm.rowconfigure(i, weight=1)

        propFrm = ttk.LabelFrame(
            parFrm, text="Propellant", style="SubLabelFrame.TLabelframe"
        )
        propFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        propFrm.rowconfigure(1, weight=1)
        propFrm.columnconfigure(0, weight=1)
        # propFrm.columnconfigure(1, weight=1)

        self.currProp = StringVar()
        self.dropProp = ttk.Combobox(
            propFrm,
            textvariable=self.currProp,
            values=self.propOptions,
            state="readonly",
            justify="center",
        )
        self.dropProp.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropProp.current(0)
        # self.dropProp.configure(width=10)
        self.dropProp.grid(row=0, column=0, columnspan=2, sticky="nsew", pady=2)

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
            height=5,
            width=0,
            yscrollcommand=specScroll.set,
        )
        self.specs.grid(row=1, column=0, sticky="nsew", pady=2)
        specScroll.config(command=self.specs.yview)
        CreateToolTip(propFrm, specsText)

        i += 1

        grainFrm = ttk.LabelFrame(
            parFrm, text="Grain Geometry", style="SubLabelFrame.TLabelframe"
        )
        grainFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        # grainFrm.rowconfigure(1, weight=0)
        grainFrm.columnconfigure(0, weight=1)
        # grainFrm.columnconfigure(1, weight=1)

        j = 0

        # geomVal = parent.register(self.updateGeom)

        # Create Dropdown menu
        self.currGeom = StringVar()

        self.dropGeom = ttk.Combobox(
            grainFrm,
            textvariable=self.currGeom,
            values=self.geoOptions,
            state="readonly",
            justify="center",
        )

        self.dropGeom.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropGeom.current(0)
        self.dropGeom.configure(width=30)

        self.dropGeom.grid(
            row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        self.lengthPrimaryAs = StringVar()
        self.lengthPrimaryTip = StringVar()

        j += 1
        self.arcmm, _, j = self.add3Input(
            parent=grainFrm,
            rowIndex=j,
            labelText=self.lengthPrimaryAs,
            # "Arc Thickness",
            unitText="mm",
            default="1.0",
            validation=validationNN,
            infotext=self.lengthPrimaryTip,
        )

        self.ratioAs = StringVar()
        self.lengthSecondaryTip = StringVar()
        self.grdR, self.grdRw, j = self.add3Input(
            parent=grainFrm,
            rowIndex=j,
            labelText=self.ratioAs,
            unitText="",
            default="1.0",
            validation=validationNN,
            infotext=self.lengthSecondaryTip,
        )

        self.lengthRatioAs = StringVar()
        self.lengthRatioTip = StringVar()

        self.grlR, self.grlRw, j = self.add3Input(
            parent=grainFrm,
            rowIndex=j,
            labelText=self.lengthRatioAs,
            unitText="",
            default="2.5",
            validation=validationNN,
            infotext=self.lengthRatioTip,
        )

        geomPlotFrm = ttk.LabelFrame(
            grainFrm, text="σ(Z)", style="SubLabelFrame.TLabelframe"
        )
        geomPlotFrm.grid(
            row=j,
            column=0,
            columnspan=3,
            sticky="nsew",
            padx=2,
            pady=2,
        )
        geomPlotFrm.rowconfigure(0, weight=1)
        geomPlotFrm.columnconfigure(0, weight=1)
        CreateToolTip(geomPlotFrm, geomPlotTxt)

        i += 1

        self.ldf, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Load Factor",
            unitText="%",
            default="50.0",
            validation=validationNN,
            infotext=ldftext,
        )

        self.clr, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Chamber L.R.",
            unitText="",
            default="1.1",
            validation=validationNN,
            infotext=clrtext,
        )

        self.dgc, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Drag coefficient",
            unitText="%",
            default="5.0",
            validation=validationNN,
            infotext=dgctext,
        )

        self.stpMPa, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText="Start Pressure",
            unitText="MPa",
            default="10",
            validation=validationNN,
            infotext=stpText,
        )

        self.currProp.trace_add("write", self.updateSpec)
        self.currGeom.trace_add("write", self.updateGeom)

        self.grdR.trace_add("write", self.callback)
        self.grlR.trace_add("write", self.callback)
        self.arcmm.trace_add("write", self.callback)

        self.parFrm = parFrm
        self.geomPlotFrm = geomPlotFrm

    def addGeomPlot(self):
        parFrm = self.parFrm
        geomPlotFrm = self.geomPlotFrm
        _, _, width, _ = parFrm.bbox("insert")
        width -= 12  # paddings

        """
        geomPlotFrm.config(width=width, height=width)
        # we lock the frame the plot is put in
        geomPlotFrm.grid_propagate(False)
        # and lock it there
        
        If we set the global dpi awareness for this window
        to false (0), this will result in a correctly sized
        window, and everything would be fine, ironically.
        (in this case, dpi = dpi)

        width = (
            geomPlotFrm.winfo_width() - 2
        )  # in pixels, -2 to account for label frame border thickness
        """
        dpi = self.dpi
        with mpl.rc_context(GEOM_CONTEXT):
            fig = Figure(
                figsize=(width / dpi, width / dpi), dpi=96, layout="constrained"
            )
            self.geomAx = fig.add_subplot()

            self.geomCanvas = FigureCanvasTkAgg(fig, master=geomPlotFrm)
            self.geomCanvas.get_tk_widget().grid(
                row=0, column=0, padx=0, pady=0, sticky="ne"
            )

    def addPlotFrm(self, parent):
        plotFrm = ttk.LabelFrame(parent, text="Plot")
        plotFrm.grid(row=0, column=0, sticky="nsew")
        plotFrm.columnconfigure(0, weight=1)
        plotFrm.rowconfigure(0, weight=1)

        self.plotFrm = plotFrm

    def addFigPlot(self):
        plotFrm = self.plotFrm
        width = (
            plotFrm.winfo_width() - 2
        )  # in pixels, -2 to account for label frame border thickness
        height = plotFrm.winfo_height() - 2

        dpi = self.dpi

        with mpl.rc_context(FIG_CONTEXT):
            fig = Figure(
                figsize=(width / dpi, height / dpi),
                dpi=96,
                layout="constrained",
            )
            # fig.subplots_adjust(right=0.75)

            axes = fig.add_subplot()

            ax = axes
            axP = ax.twinx()
            axv = ax.twinx()

            ax.yaxis.tick_right()
            ax.set(xlabel="Domain")
            axv.spines.right.set_position(("axes", 1.0 + 40 * dpi / 96 / width))

            # fig.tight_layout(pad=1)

            self.ax = ax
            self.axP = axP
            self.axv = axv
            self.fig = fig

            self.pltCanvas = FigureCanvasTkAgg(fig, master=plotFrm)
            self.pltCanvas.get_tk_widget().grid(
                row=0, column=0, padx=0, pady=0, sticky="nsew"
            )

    def resizeFigPlot(self, event):
        """
        width = (
            plotFrm.winfo_width() - 2
        )  # in pixels, -2 to account for label frame border thickness
        height = plotFrm.winfo_height() - 2
        """
        plotFrm = self.plotFrm
        _, _, width, height = plotFrm.bbox("insert")
        dpi = self.dpi

        with mpl.rc_context(FIG_CONTEXT):
            # print(width, height)
            # print(width / dpi, height / dpi)
            self.fig.set_size_inches(width / dpi, height / dpi)
            self.axv.spines.right.set_position(
                ("axes", 1 + 40 * dpi / 96 / width)
            )

    def updateFigPlot(self):
        with mpl.rc_context(FIG_CONTEXT):
            gun = self.gun

            self.ax.cla()
            self.axP.cla()
            self.axv.cla()

            dpi = self.dpi
            size = self.fig.get_size_inches() * self.fig.dpi
            self.axv.spines.right.set_position(
                ("axes", 1 + 40 * dpi / 96 / size[0])
            )
            self.ax.set(xlabel="Domain")

            if gun is not None:
                xs = []
                vs = []
                Ps = []
                psis = []
                dom = self.dropOptn.get()

                for i, (t, (l, psi, v, p)) in enumerate(self.intgRecord):
                    if dom == DOMAIN_TIME:
                        xs.append(t * 1000)
                    elif dom == DOMAIN_LENG:
                        xs.append(l)
                    vs.append(v)
                    Ps.append(p / 1e6)
                    psis.append(psi)

                self.axv.scatter(xs, vs, color="tab:blue", marker="s", s=8)

                self.axP.scatter(xs, Ps, color="tab:green", marker="s", s=8)

                self.ax.scatter(xs, psis, color="tab:red", marker="s", s=8)

                xPeak = 0
                for i, (tag, t, l, psi, v, p, T) in enumerate(self.tableData):
                    if dom == DOMAIN_TIME:
                        x = t * 1000
                    elif dom == DOMAIN_LENG:
                        x = l
                    xs.append(x)
                    if tag == POINT_PEAK:
                        xPeak = x
                    vs.append(v)
                    Ps.append(p / 1e6)
                    psis.append(psi)

                self.axP.spines.right.set_position(("data", xPeak))

                self.ax.set_xlim(left=0, right=xs[-1])
                self.ax.set_ylim(bottom=0, top=1.025)

                self.axP.set(ylim=(0, max(Ps) * 1.05))
                self.axv.set(ylim=(0, max(vs) * 1.05))

                (xs, vs, Ps, psis) = zip(
                    *sorted(zip(xs, vs, Ps, psis), key=lambda line: line[0])
                )

                (pv,) = self.axv.plot(
                    xs,
                    vs,
                    "tab:blue",
                    label="Shot Velocity\nm/s",
                    marker=".",
                    alpha=0.75,
                )
                (pP,) = self.axP.plot(
                    xs,
                    Ps,
                    "tab:green",
                    label="Avg. Pressure\nMPa",
                    marker=".",
                    alpha=0.75,
                )
                (ref,) = self.ax.plot(
                    (0, xs[-1]),
                    (1, 1),
                    "tab:red",
                    alpha=0.5,
                    linestyle="dashed",
                    label="Burnout",
                )
                (ppsi,) = self.ax.plot(
                    xs,
                    psis,
                    "tab:red",
                    label="Volume Burnup",
                    marker=".",
                    alpha=0.75,
                )

                tkw = dict(size=4, width=1.5)

                self.ax.yaxis.tick_right()
                self.ax.tick_params(axis="y", colors=ppsi.get_color(), **tkw)
                self.axv.tick_params(axis="y", colors=pv.get_color(), **tkw)
                self.axP.tick_params(axis="y", colors=pP.get_color(), **tkw)
                self.ax.tick_params(axis="x", **tkw)

                self.axP.yaxis.set_ticks(self.axP.get_yticks()[1:-1:])

                if dom == DOMAIN_TIME:
                    self.ax.set(xlabel="Time - ms")
                elif dom == DOMAIN_LENG:
                    self.ax.set(xlabel="Length - m")

                for ax, xvals in zip(
                    (
                        self.axP,
                        self.ax,
                        self.axv,
                    ),
                    (
                        (0.5 * xs[-1] + 0.5 * xPeak, xs[-1]),
                        (0, xPeak),
                        (xPeak, 0.5 * xs[-1] + 0.5 * xPeak),
                    ),
                ):
                    labelLines(
                        ax.get_lines(),
                        align=True,
                        zorder=3.5,
                        color="white",
                        xvals=xvals,
                    )

                # labelLines(, zorder=2.5, color="white")

            self.pltCanvas.draw()

    def addTblFrm(self, parent):
        columnList = [
            "Event",
            "Time",
            "Travel",
            "Burnup",
            "Velocity",
            "Avg. Pressure",
            "Avg. Temperature",
        ]
        tblFrm = ttk.LabelFrame(parent, text="Result Table")
        tblFrm.grid(row=1, column=0, sticky="nsew")

        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(1, weight=1)

        # configure the numerical

        self.tv = ttk.Treeview(
            tblFrm, selectmode="browse", height=10
        )  # this set the nbr. of values
        self.tv.grid(row=1, column=0, sticky="nsew")

        self.tv["columns"] = columnList
        self.tv["show"] = "headings"
        self.tv.tag_configure(POINT_PEAK, foreground="orange")
        self.tv.tag_configure(POINT_BURNOUT, foreground="red")
        self.tv.tag_configure(POINT_FRACTURE, foreground="brown")
        # self.tv.tag_configure("monospace", font=("TkFixedFont", 9))

        t_Font = tkFont.Font(family="hack", size=8)

        self.tv.tag_configure("monospace", font=t_Font)
        self.tv.tag_configure("error", font=("hack", 7), foreground="grey")

        # we use a fixed width font so any char will do
        width, height = t_Font.measure("m"), t_Font.metrics("linespace")

        for column in columnList:  # foreach column
            self.tv.heading(
                column, text=column
            )  # let the column heading = column name
            self.tv.column(
                column,
                stretch=True,  # will adjust to window resizing
                width=width * 16,
                minwidth=width * 16,
                anchor="e",
            )

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

    def updateSpec(self, var, index, mode):
        self.specs.config(state="normal")
        compo = self.compositions[self.dropProp.get()]
        self.specs.delete("1.0", "end")
        for line in compo.desc.split(","):
            self.specs.insert("end", line + "\n")
        self.specs.config(state="disabled")
        # this updates the specification description

        self.callback(None, None, None)

        return True

    def updateGeom(self, var, index, mode):
        geom = self.geometries[self.dropGeom.get()]
        if geom == SimpleGeometry.SPHERE:
            self.grdRw.config(state="disabled")
            self.grlRw.config(state="disabled")

        elif geom == SimpleGeometry.CYLINDER:
            self.grdRw.config(state="disabled")
            self.grlRw.config(state="normal")

        else:
            self.grdRw.config(state="normal")
            self.grlRw.config(state="normal")

        if geom == SimpleGeometry.SPHERE:
            self.lengthPrimaryAs.set("Diameter")

            self.lengthRatioAs.set("")
            self.ratioAs.set("")

            self.lengthPrimaryTip.set(diaText)
            self.lengthRatioTip.set("")
            self.lengthSecondaryTip.set("")

        elif geom == SimpleGeometry.ROD:
            self.lengthPrimaryAs.set("Width")
            self.lengthRatioAs.set("Length / Width")
            self.ratioAs.set("Height / Width")

            self.lengthPrimaryTip.set(widthText)
            self.lengthRatioTip.set(rodRtext)
            self.lengthSecondaryTip.set(heightRtext)

        elif geom == SimpleGeometry.CYLINDER:
            self.lengthPrimaryAs.set("Diameter")

            self.lengthRatioAs.set("Length / Diameter")
            self.ratioAs.set("")

            self.lengthPrimaryTip.set(diaText)
            self.lengthRatioTip.set(cylLRtext)
            self.lengthSecondaryTip.set("")

        else:
            self.lengthPrimaryAs.set("Arc Thickness")
            self.lengthRatioAs.set("Length / Diameter")
            self.ratioAs.set("Perf.Dia. / A.Th.")

            self.lengthPrimaryTip.set(arcText)
            self.lengthRatioTip.set(perfLRtext)
            self.lengthSecondaryTip.set(pDiaRText)

        self.callback(None, None, None)

        return True

    def updateGeomPlot(self):
        with mpl.rc_context(GEOM_CONTEXT):
            N = 50
            prop = self.prop
            self.geomAx.cla()
            if prop is not None:
                xs = [i / (N - 1) * prop.Z_b for i in range(N)]
                ys = [prop.f_sigma_Z(x) for x in xs]
                xs.append(xs[-1])
                ys.append(0)
                self.geomAx.plot(
                    xs,
                    ys,
                    color="#215d9c",
                )
                self.geomAx.grid(
                    which="major", color="grey", linestyle="dotted"
                )
                # self.geomAx.minorticks_on()
                self.geomAx.set_xlim(left=0, right=prop.Z_b)
                self.geomAx.xaxis.set_ticks(
                    [i * 0.2 for i in range(ceil(prop.Z_b / 0.2) + 1)]
                )
                self.geomAx.set_ylim(bottom=0, top=max(ys))
                self.geomAx.yaxis.set_ticks(
                    [i * 0.2 for i in range(ceil(max(ys) / 0.2) + 1)]
                )

            self.geomCanvas.draw()

    def updateError(self):
        self.errorText.delete("1.0", "end")
        for line in self.errorLst:
            self.errorText.insert("end", line + "\n")
        self.errorLst = []

    def callback(self, var, index, mode):
        """
        updates the propellant object on write to the ratio entry fields
        and, on changing the propellant or geometrical specification.

        Double calling is due value validation, no workaround has been found
        at this time!
        """
        # print("called", var, index, mode)

        geom = self.geometries[self.dropGeom.get()]
        compo = self.compositions[self.dropProp.get()]

        try:
            self.prop = Propellant(
                compo,
                geom,
                float(self.grdR.get()),
                float(self.grlR.get()),
            )
            self.updateGeomPlot()

        except Exception as e:
            self.prop = None
            if DEBUG:
                self.errorLst.append("".join(traceback.format_exception(e)))
            else:
                self.errorLst.append(str(e))

        self.updateError()

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
        infotext=None,
        anchor="w",
        reverse=False,
    ):
        if isinstance(labelText, StringVar):
            lb = ttk.Label(parent, textvariable=labelText, anchor=anchor)
        else:
            lb = ttk.Label(parent, text=labelText, anchor=anchor)

        lb.grid(
            row=rowIndex,
            column=colIndex + (1 if reverse else 0),
            sticky="nsew",
            padx=2,
            pady=2,
        )
        if infotext is not None:
            CreateToolTip(lb, infotext)
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
            row=rowIndex,
            column=colIndex + (0 if reverse else 1),
            sticky="nsew",
            padx=2,
            pady=2,
        )
        en.bind("<FocusOut>", formatter)
        return e, en, rowIndex + 1

    def add3Input(
        self,
        parent,
        rowIndex,
        labelText="",
        unitText="",
        default="0.0",
        validation=None,
        entryWidth=10,
        formatter=formatFloatInput,
        color=None,
        infotext=None,
        colIndex=0,
    ):
        e, en, _ = self.add2Input(
            parent,
            rowIndex,
            colIndex,
            labelText,
            default,
            validation,
            entryWidth,
            formatter,
            color,
            infotext,
        )
        ttk.Label(parent, text=unitText).grid(
            row=rowIndex, column=colIndex + 2, sticky="nsew", padx=2, pady=2
        )
        return e, en, rowIndex + 1

    def add12Disp(
        self,
        parent,
        rowIndex,
        colIndex=0,
        labelText="",
        unitText="",
        default="0.0",
        entryWidth=5,
        justify="center",
        infotext=None,
    ):
        lb = ttk.Label(parent, text=labelText)
        lb.grid(
            row=rowIndex,
            column=colIndex,
            columnspan=2,
            sticky="nsew",
            padx=2,
            pady=2,
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
        en.grid(
            row=rowIndex + 1, column=colIndex, sticky="nsew", padx=2, pady=2
        )
        ttk.Label(parent, text=unitText).grid(
            row=rowIndex + 1, column=colIndex + 1, sticky="nsew", padx=2, pady=2
        )
        if infotext is not None:
            CreateToolTip(lb, infotext)

        return e, en, rowIndex + 2


def center(win):
    """
    centers a tkinter window
    :param win: the main window or Toplevel window to center
    """
    win.update_idletasks()
    width = win.winfo_width()
    frm_width = win.winfo_rootx() - win.winfo_x()
    win_width = width + 2 * frm_width
    height = win.winfo_height()
    titlebar_height = win.winfo_rooty() - win.winfo_y()
    win_height = height + titlebar_height + frm_width
    x = win.winfo_screenwidth() // 2 - win_width // 2
    y = win.winfo_screenheight() // 2 - win_height // 2
    win.geometry("{}x{}+{}+{}".format(width, height, x, y))
    win.deiconify()


if __name__ == "__main__":
    # this tells windows that our program will handle scaling ourselves
    windll.shcore.SetProcessDpiAwareness(1)
    root = Tk()
    # one must supply the entire path
    loadfont(resolvepath("ui/Hack-Regular.ttf"))
    # dpi = round(root.winfo_fpixels("1i") / 96) * 96
    dpi = root.winfo_fpixels("1i")
    """
    # this will return the "windows scaling factor"
    # wsf = ctypes.windll.user32.GetDpiForWindow(root.winfo_id()) / 96
    MM_TO_IN = 0.0393700787
    dc = ctypes.windll.user32.GetDC(root.winfo_id())
    # the monitor's physical width in inches
    mw = ctypes.windll.gdi32.GetDeviceCaps(dc, 4) * MM_TO_IN
    # the monitor's physical height in inches
    mh = ctypes.windll.gdi32.GetDeviceCaps(dc, 6) * MM_TO_IN
    # horizontal resolution
    dw = ctypes.windll.gdi32.GetDeviceCaps(dc, 8)
    # vertical resolution
    dh = ctypes.windll.gdi32.GetDeviceCaps(dc, 10)

    # get the "true, physical" DPI
    hdpi, vdpi = dw / mw, dh / mh
    # get the "diagonal" physical dpi
    ddpi = (hdpi**2 + vdpi**2) ** 0.5
    """

    # Tk was originally developed for a dpi of 72
    root.tk.call("tk", "scaling", "-displayof", ".", dpi / 72.0)

    mpl.rc("figure", dpi=dpi)

    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.12"))
    # root.tk.call("package", "require", "awdark")

    style = ttk.Style(root)
    style.theme_use("awdark")

    # ensure that the treeview rows are roughly the same height
    # regardless of dpi. on Windows, default is Segoe UI at 9 points
    # so the default row height should be around 12

    style.configure("Treeview", rowheight=round(12 * dpi / 72.0))

    style.configure("Treeview.Heading", font=("Hack", 8))
    style.configure("TButton", font=("Hack", 10, "bold"))
    style.configure("TLabelframe.Label", font=("Hack", 10, "bold"))
    style.configure("SubLabelFrame.TLabelframe.Label", font=("Hack", 9))
    style.configure("TNotebook.Tab", font=("Hack", 10))
    root.option_add("*Font", "Hack 8")

    # root.option_add("*TCombobox*Listbox*Font", "Hack 8")

    root.title("Phoenix's Internal Ballistics Solver v0.3")

    ibPanel = IB(root, dpi)

    center(root)
    root.minsize(root.winfo_width(), root.winfo_height())
    root.mainloop()
