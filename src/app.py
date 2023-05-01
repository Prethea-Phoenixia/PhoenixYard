from tkinter import *
from tkinter import ttk
import tkinter.font as tkFont
import traceback
from gun import *
import os
import sys
from math import ceil

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
        magnitude = 10**closest
        vstr = "{:#.{:}g}".format(v / magnitude, dec)
        return (
            (" " if positive else "-")
            + "{:#.{:}g}".format(v / magnitude, dec)
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


def validateGZ(inp):
    """
    validate an input such that the result is greater than zero."""

    if inp == "" or inp == ".":
        return True  # we will catch this by filling the default value
    try:
        return float(inp) > 0
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

        if (
            y + t_Font.measure(self.text) / (width * columnWidth) * height
            > ry + root.winfo_rooty() + root.winfo_height()
        ):
            print("this")

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
    def __init__(self, parent):
        Frame.__init__(self, parent)

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

        self.wrapup()

        # parent.bind("<Return>", self.calculate)

    def calculate(self, event=None):
        # force an immediate redraw after calculation
        compo = self.compositions[self.dropProp.get()]
        # lookup dictionary using the string key to get
        # the requisite object
        geom = self.geometries[self.dropGeom.get()]

        self.tableData = []
        self.errorData = []

        # self.resetSpec()

        try:
            self.prop = Propellant(
                compo,
                geom,
                float(self.arcmm.get()) / 1000,
                float(self.permm.get()) / float(self.arcmm.get()),
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
            self.cv.set(toSI(chamberVolume, useSN=True).strip())
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
            # self.tv.move(str(i + 2), str(i), "end")
            i += 2

        self.updateError()

    def updateError(self):
        self.errorText.delete("1.0", "end")
        if self.geomError:
            self.errorLst = [
                "Specified Geometry of Propellant is Invalid"
            ] + self.errorLst

        for line in self.errorLst:
            self.errorText.insert("end", line + "\n")
        self.errorLst = []

    def addSpecFrm(self, parent):
        specFrm = ttk.LabelFrame(parent, text="Design Summary")
        specFrm.grid(row=0, column=0, rowspan=2, sticky="nsew")
        specFrm.columnconfigure(0, weight=1)
        i = 0
        vinfText = " ".join(
            (
                "Velocity the shot would achieve if",
                "the barrel is extended to infinite length.",
            )
        )
        self.va, _, i = self.add12Disp(
            specFrm,
            i,
            "Asymptotic Vel.",
            "m/s",
            justify="right",
            infotext=vinfText,
        )
        teffText = " ".join(
            (
                "Thermal efficiency of the gun system, i.e.",
                "the amount of work done to both gas and",
                "projectile over chemical potential energy",
                "of propellant.",
            )
        )
        self.te, _, i = self.add12Disp(
            specFrm, i, "Thermal Eff.", "%", infotext=teffText
        )
        beffText = " ".join(
            (
                "Ballistic efficiency of the gun system, i.e.",
                "the amount of work done the projectile over",
                "chemical potential energy of propellant.",
            )
        )
        self.be, _, i = self.add12Disp(
            specFrm, i, "Ballistic Eff.", "%", infotext=beffText
        )
        self.cv, _, i = self.add12Disp(
            specFrm, i, "Chamber Volume", "m³", justify="right"
        )

    def addParFrm(self, parent):
        parFrm = ttk.LabelFrame(parent, text="Parameters")
        parFrm.grid(row=0, column=2, sticky="nsew")
        parFrm.columnconfigure(0, weight=3)
        parFrm.columnconfigure(1, weight=3)
        parFrm.columnconfigure(2, weight=1)

        # validation
        validationNN = parent.register(validateNN)
        validationGZ = parent.register(validateGZ)

        i = 0
        self.calmm, _, i = self.add3Input(
            parFrm, i, "Caliber", "mm", "50.0", validationGZ
        )
        self.tblmm, _, i = self.add3Input(
            parFrm, i, "Tube Length", "mm", "3500.0", validationGZ
        )
        self.shtkg, _, i = self.add3Input(
            parFrm, i, "Shot Mass", "kg", "1.0", validationGZ
        )

        chgText = " ".join(
            (
                "Mass of propellant charge to be used. Chamber",
                "volume is determined from this and load fraction.",
            )
        )
        self.chgkg, _, i = self.add3Input(
            parFrm,
            i,
            "Charge Mass",
            "kg",
            "1.0",
            validationGZ,
            infotext=chgText,
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

        specsText = " ".join(
            (
                "Nitration used for Nitrocellulose is between",
                "13.15%-13.25% for all preset propellants.",
                "Up to 14.14% is possible but has not been",
                "adopted by the industry due to cost.\n",
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
        """
        DNT coating can act to retard initial burning and 
        reduce the sensitivity of propellant behaviour to 
        initial manufacturing deviations"""

        CreateToolTip(propFrm, specsText)

        propFrm.rowconfigure(1, weight=1)
        propFrm.columnconfigure(0, weight=1)

        i += 1

        ttk.Label(parFrm, text="Grain Shape").grid(
            row=i,
            column=0,
            sticky="nsew",
            padx=2,
            pady=2,
        )

        i += 1

        geomVal = parent.register(self.updateGeom)

        # Create Dropdown menu
        self.dropGeom = ttk.Combobox(
            parFrm,
            values=self.geoOptions,
            state="readonly",
            justify="center",
            validate="focusin",
            validatecommand=(geomVal, "%P"),
        )

        self.dropGeom.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropGeom.current(0)
        self.dropGeom.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        self.dropGeom.configure(width=10)

        i += 1

        self.lengthPrimaryAs = StringVar()
        self.lengthPrimaryTip = StringVar()

        self.arcmm, _, i = self.add3Input(
            parFrm,
            i,
            self.lengthPrimaryAs,
            # "Arc Thickness",
            "mm",
            "1.0",
            validationGZ,
            infotext=self.lengthPrimaryTip,
        )

        self.lengthSecondaryAs = StringVar()
        self.lengthSecondaryTip = StringVar()

        self.permm, self.perw, i = self.add3Input(
            parFrm,
            i,
            self.lengthSecondaryAs,
            # "Perf Diameter",
            "mm",
            "0.0",
            validationNN,  # can be 0
            infotext=self.lengthSecondaryTip,
        )

        self.lengthRatioAs = StringVar()
        self.lengthRatioTip = StringVar()

        self.grlR, self.grlRw, i = self.add3Input(
            parFrm,
            i,
            self.lengthRatioAs,
            "",
            "2.5",
            validationNN,
            infotext=self.lengthRatioTip,
        )

        ldftext = " ".join(
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
                "the breech.\n",
                'This is different from the usual "loading',
                'density" quoted as the geometrical limit',
                "imposed by grain shape is already taken",
                "into account.",
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

        clrtext = " ".join(
            (
                "Chamber length ratio is the ratio between",
                "the length of reduced chamber (dividing",
                "the chamber volume with barrel cross",
                "section) to the actual chamber.",
            )
        )

        self.clr, _, i = self.add3Input(
            parFrm, i, "Chamber L.R.", "", "1.1", validationNN, infotext=clrtext
        )

        stpText = " ".join(
            (
                "Peak pressure that initially resists",
                "the movement of shot. This is made up of",
                "the rifling resisting the drive band or",
                "shell body, the static friction of the",
                "barrel, and the cartridge case holding onto",
                "the shot (greatly varying between 0.25-15MPa",
                "depending on caliber and desired RoF)\n",
                "For large caliber guns 25-30MPa is common.",
                "For rifles, 35-45MPa is common.",
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

    def updateSpec(self, *inp):
        self.specs.config(state="normal")
        compo = self.compositions[self.dropProp.get()]
        self.specs.delete("1.0", "end")
        for line in compo.desc.split(","):
            self.specs.insert("end", line + "\n")
        self.specs.config(state="disabled")
        # this updates the specification description

        return True

    def updateGeom(self, *inp):
        self.configEntry()
        geom = self.geometries[self.dropGeom.get()]
        arcText = " ".join(
            (
                "Specify the geometry of propellant grain",
                "using dimensions.\n",
                "Burning radially outward along the arc is",
                "the shortest path, and therefore is the",
                "primary geometrical determinant of burn rapidity.\n",
                "In theory micrometer level precision",
                "is possible. In practice, tolerance for indu-",
                "strial bulk production is in the range of",
                "0.15mm to 1mm depending on caliber arc",
                "thickness is generally found close to 1mm",
                "for small to intermediate calibers.\n",
                "Greater web thickness will lead to a longer burn",
                "time for the same load of propellant, giving",
                "a more gradual rise and lower peak pressure.",
            )
        )

        pDiaText = " ".join(
            (
                "Specify the geometry of propellant grain",
                "using dimensions.\n",
                "Perforations are formed by protrusions in the copper",
                "casting die, the misalignment of which contributes",
                "to deviations. Larger perforation diameter increase",
                "progressiveness of burn but contributes to a lower",
                "achievable loading density, with standard multi-perf grains",
                "coming with a hole diameter half that of web thickness.",
            )
        )

        diaText = " ".join(("Specify the diameter of the propellant grain.",))

        grlRtext = " ".join(
            (
                "Length to diameter ratio of the grain, usually",
                "around 2-2.5. A higher value facilitate progressive",
                "burning. If this is too low its possible for",
                "propellant to exhibit regressive burning behaviour,",
                "as the primary mode will be axial burning instead.",
                "Length is usually limited by propellant strength",
                "during the casting process.",
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

        heightText = " ".join(
            ("Specify the height of propellant rod or flake.",)
        )

        lengthRtext = " ".join(("Height to width ratio of the rod or flake.",))

        ratioEntryText = " ".join(
            (
                "Specify the geometry of propellant using",
                "ratios, scaled by the top entry. These",
                "entries are in use when Lock Geometry is",
                "enabled.",
            )
        )

        if geom == SimpleGeometry.SPHERE:
            self.lengthPrimaryAs.set("Diameter")
            self.lengthSecondaryAs.set("")
            self.lengthRatioAs.set("")
            self.primaryRatio.set("")
            self.secondaryRatio.set("")
            self.lengthPrimaryTip.set(diaText)
            self.lengthSecondaryTip.set("")
            self.lengthRatioTip.set("")
            self.ratioTip.set("")

        elif geom == SimpleGeometry.ROD:
            self.lengthPrimaryAs.set("Width")
            self.lengthSecondaryAs.set("Height")
            self.lengthRatioAs.set("Length/Width")
            self.primaryRatio.set("Width")
            self.secondaryRatio.set("Height")
            self.lengthPrimaryTip.set(widthText)
            self.lengthSecondaryTip.set(heightText)
            self.lengthRatioTip.set(lengthRtext)
            self.ratioTip.set(ratioEntryText)

        else:
            self.lengthPrimaryAs.set("Arc Thickness")
            self.lengthSecondaryAs.set("Perf Diameter")
            self.lengthRatioAs.set("Grain L/D")
            self.primaryRatio.set("Arc")
            self.secondaryRatio.set("P.D.")
            self.lengthPrimaryTip.set(arcText)
            self.lengthSecondaryTip.set(pDiaText)
            self.lengthRatioTip.set(grlRtext)
            self.ratioTip.set(ratioEntryText)

        return True

    def addOpsFrm(self, parent):
        opFrm = ttk.LabelFrame(parent, text="Options")
        opFrm.grid(row=1, column=2, sticky="nsew")

        opFrm.columnconfigure(1, weight=1)

        validationNN = parent.register(validateNN)

        i = 0

        self.primaryRatio = StringVar()
        self.ratioTip = StringVar()
        self.arcR, self.arcRw, i = self.add2Input(
            opFrm,
            i,
            0,
            self.primaryRatio,
            "1.0",
            validation=validationNN,
            infotext=self.ratioTip,
        )

        self.secondaryRatio = StringVar()
        self.perR, self.perRw, i = self.add2Input(
            opFrm,
            i,
            0,
            self.secondaryRatio,
            "0.5",
            validation=validationNN,
            infotext=self.ratioTip,
        )

        self.configEntry()

        self.arcmm.trace_add("write", self.arccallback)
        self.permm.trace_add("write", self.arccallback)
        self.arcR.trace_add("write", self.arccallback)
        self.perR.trace_add("write", self.arccallback)

        ttk.Checkbutton(
            opFrm,
            text="Lock\nGeom",
            variable=self.geoLocked,
            onvalue=1,
            offvalue=0,
            command=self.configEntry,
        ).grid(row=0, column=2, rowspan=2, sticky="nsew", padx=2, pady=2)

        samplbl = ttk.Label(opFrm, text="Sample")
        samplbl.grid(row=2, column=0, sticky="nsew", padx=2, pady=2)
        sampTxt = " ".join(
            ("The amount of sample points and domain to take them.", "")
        )
        CreateToolTip(samplbl, sampTxt)
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

        i = 3

        self.steps, _, i = self.add3Input(
            opFrm,
            i,
            "",
            "steps",
            "10",
            validation=validationNN,
            formatter=formatIntInput,
        )

        validationPI = parent.register(validatePI)

        tolText = " ".join(
            (
                "The maximum relative error, ε, that is allowed in the integrands",
                "for each component. Some components may have significantly less",
                "error than specified here, as shown under each entry.",
            )
        )
        self.accExp, _, i = self.add2Input(
            opFrm,
            i,
            0,
            "-log10(ε)",
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
        CreateToolTip(calButton, "Integrate system using RKF7(8) integrator")

    def configEntry(self):
        ratioEntrys = (self.arcRw, self.perRw)
        directEntrys = (self.perw,)
        geom = self.geometries[self.dropGeom.get()]
        if geom == SimpleGeometry.SPHERE:
            for en in ratioEntrys:
                en.config(state="disabled")
            for en in directEntrys:
                en.config(state="disabled")
            self.grlRw.config(state="disabled")
        else:
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
            self.grlRw.config(state="normal")

        self.arccallback(None, None, None)

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
            width=0,
        )
        self.errorText.grid(row=0, column=0, sticky="nsew")

    def wrapup(self):
        self.updateSpec()
        self.updateGeom()

    def addTblFrm(self, parent):
        columnList = [
            "Event",
            "Time",
            "Travel",
            "Burnup",
            "Velocity",
            "Pressure",
            "Temperature",
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
                stretch=1,
                width=width * 20,
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
    ):
        if isinstance(labelText, StringVar):
            lb = ttk.Label(parent, textvariable=labelText)
        else:
            lb = ttk.Label(parent, text=labelText)

        lb.grid(row=rowIndex, column=colIndex, sticky="nsew", padx=2, pady=2)
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

    def add11Disp(
        self,
        parent,
        rowIndex,
        labelText,
        default="0.0",
        entryWidth=5,
        justify="center",
        infotext=None,
    ):
        lb = ttk.Label(parent, text=labelText)
        lb.grid(
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
        en.grid(
            row=rowIndex + 1,
            column=0,
            columnspan=2,
            sticky="nsew",
            padx=2,
            pady=2,
        )
        if infotext is not None:
            CreateToolTip(lb, infotext)
        return e, en, rowIndex + 2

    def add12Disp(
        self,
        parent,
        rowIndex,
        labelText,
        unitText,
        default="0.0",
        entryWidth=5,
        justify="center",
        infotext=None,
    ):
        lb = ttk.Label(parent, text=labelText)
        lb.grid(
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
        if infotext is not None:
            CreateToolTip(lb, infotext)

        return e, en, rowIndex + 2

    def arccallback(self, var, index, mode):
        self.geomError = False

        if self.geoLocked.get() == 1:
            try:
                self.permm.set(
                    float(self.arcmm.get())
                    / float(self.arcR.get())
                    * float(self.perR.get())
                )
            except (ZeroDivisionError, ValueError):
                self.geomError = True
        else:
            try:
                self.perR.set(
                    float(self.permm.get())
                    / float(self.arcmm.get())
                    * float(self.arcR.get())
                )
            except (ZeroDivisionError, ValueError):
                self.geomError = True
        self.updateError()


if __name__ == "__main__":
    # high dpi scaling
    windll.shcore.SetProcessDpiAwareness(1)
    root = Tk()
    # one must supply the entire path
    loadfont(resolvepath("ui/Hack-Regular.ttf"))
    dpi = root.winfo_fpixels("1i")
    root.tk.call("tk", "scaling", "-displayof", ".", dpi / 72.0)
    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.12"))
    # root.tk.call("package", "require", "awdark")

    style = ttk.Style(root)
    style.theme_use("awdark")

    # ensure that the treeview rows are roughly the same height
    # regardless of dpi. on Windows, default is Segoe UI at 9 points
    # so the default row height should be around 12

    style.configure("Treeview", rowheight=round(12 * dpi / 72))
    style.configure("Treeview.Heading", font=("Hack", 8))
    style.configure("TButton", font=("Hack", 10, "bold"))
    style.configure("TLabelframe.Label", font=("Hack", 10, "bold"))
    style.configure("TNotebook.Tab", font=("Hack", 10))
    root.option_add("*Font", "Hack 8")

    # root.option_add("*TCombobox*Listbox*Font", "Hack 8")

    root.title("Phoenix's Internal Ballistics Solver v0.3")

    ibPanel = IB(root)
    root.mainloop()
