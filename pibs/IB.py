from tkinter import Frame, Menu, Text, filedialog, messagebox, StringVar, IntVar
from tkinter import Tk, ttk

import tkinter.font as tkFont
import traceback

from gun import Gun, DOMAIN_TIME, DOMAIN_LENG
from gun import POINT_START, POINT_BURNOUT, POINT_FRACTURE, POINT_EXIT
from gun import POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH
from gun import SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV
from recoiless import Recoiless

from prop import Propellant, GrainComp, GEOMETRIES, SimpleGeometry
from opt import Constrained
from optRecoiless import ConstrainedRecoiless
from tip import CreateToolTip
from lang import STRING

from misc import validateNN, validatePI, validateFLT
from misc import roundSig, formatMass, dot_aligned, toSI
from misc import formatIntInput
from misc import loadfont, resolvepath

from math import ceil, pi, log10

import matplotlib.pyplot as mpl
from matplotlib.figure import Figure
from labellines import labelLines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from multiprocessing import Process, Queue, freeze_support
from queue import Empty

import sys
import csv

import json

from ctypes import windll
import platform

from matplotlib import font_manager

from locWidget import Loc12Disp, Loc122Disp
from locWidget import Loc2Input, Loc3Input
from locWidget import LocLabelFrame, LocLabelCheck, LocDropdown

from material import MATERIALS


RECOILESS = "RECOILESS"
CONVENTIONAL = "CONVENTIONAL"
HIGHLOW = "HIGHLOW"

TYPES = {CONVENTIONAL: CONVENTIONAL, RECOILESS: RECOILESS}

CARTRIDGE = "CARTRIDGE"
TELESCOPED = "TELESCOPED"

AMMUNITIONS = {CARTRIDGE: CARTRIDGE, TELESCOPED: TELESCOPED}

SOLUTIONS = {
    SOL_LAGRANGE: SOL_LAGRANGE,
    SOL_PIDDUCK: SOL_PIDDUCK,
    SOL_MAMONTOV: SOL_MAMONTOV,
}
DOMAINS = {
    DOMAIN_TIME: DOMAIN_TIME,
    DOMAIN_LENG: DOMAIN_LENG,
}

CONTROLS = {
    POINT_PEAK_AVG: POINT_PEAK_AVG,
    POINT_PEAK_BREECH: POINT_PEAK_BREECH,
    POINT_PEAK_SHOT: POINT_PEAK_SHOT,
}

USE_LF = "USE_LF"
USE_CV = "USE_CV"
USE = {USE_LF: USE_LF, USE_CV: USE_CV}

FONTNAME = "Sarasa Fixed SC"
FONTSIZE = 9

GEOM_CONTEXT = {
    "font.size": FONTSIZE,
    "axes.titlesize": FONTSIZE,
    "axes.labelsize": FONTSIZE,
    "xtick.labelsize": FONTSIZE,
    "ytick.labelsize": FONTSIZE,
    "legend.fontsize": FONTSIZE,
    "figure.titlesize": FONTSIZE + 2,
    "lines.markersize": FONTSIZE / 4,
    "axes.axisbelow": True,
    "font.family": "Sarasa Fixed SC",
}

FIG_CONTEXT = {
    "font.size": FONTSIZE,
    "axes.titlesize": FONTSIZE,
    "axes.labelsize": FONTSIZE,
    "xtick.labelsize": FONTSIZE,
    "ytick.labelsize": FONTSIZE,
    "legend.fontsize": FONTSIZE,
    "figure.titlesize": FONTSIZE + 2,
    "lines.linewidth": 1,
    "font.weight": "bold",
    "lines.markersize": FONTSIZE / 2,
    "axes.axisbelow": False,
    "axes.labelweight": "bold",
    "yaxis.labellocation": "top",
    "font.family": "Sarasa Fixed SC",
}


class InteriorBallisticsFrame(Frame):
    def __init__(self, parent, menubar, dpi):
        ttk.Frame.__init__(self, parent)
        self.grid(row=0, column=0, sticky="nsew")
        self.LANG = StringVar(value=list(STRING.keys())[0])

        self.locs = []

        default_font = tkFont.Font(family=FONTNAME, size=FONTSIZE)
        self.option_add("*Font", default_font)

        self.queue = Queue()
        self.process = None

        self.dpi = dpi
        self.parent = parent
        self.forceUpdOnThemeWidget = []

        self.menubar = menubar

        fileMenu = Menu(menubar)
        menubar.add_cascade(
            label=self.getLocStr("fileLabel"), menu=fileMenu, underline=0
        )
        themeMenu = Menu(menubar)
        menubar.add_cascade(
            label=self.getLocStr("themeLabel"), menu=themeMenu, underline=0
        )
        debugMenu = Menu(menubar)
        menubar.add_cascade(
            label=self.getLocStr("debugLabel"), menu=debugMenu, underline=0
        )
        langMenu = Menu(menubar)
        menubar.add_cascade(label="Lang 语言", menu=langMenu, underline=0)

        self.fileMenu = fileMenu
        self.themeMenu = themeMenu
        self.debugMenu = debugMenu

        self.themeRadio = IntVar(value=0)
        self.useTheme()

        self.DEBUG = IntVar(value=0)

        fileMenu.add_command(
            label=self.getLocStr("saveLabel"), command=self.save, underline=0
        )
        fileMenu.add_command(
            label=self.getLocStr("loadLabel"), command=self.load, underline=0
        )
        fileMenu.add_command(
            label=self.getLocStr("exportLabel"),
            command=self.export,
            underline=0,
        )

        themeMenu.add_radiobutton(
            label=self.getLocStr("darkLabel"),
            variable=self.themeRadio,
            value=0,
            command=self.useTheme,
            underline=0,
        )
        themeMenu.add_radiobutton(
            label=self.getLocStr("lightLabel"),
            variable=self.themeRadio,
            value=1,
            command=self.useTheme,
            underline=0,
        )

        debugMenu.add_checkbutton(
            label=self.getLocStr("enableLabel"),
            variable=self.DEBUG,
            onvalue=1,
            offvalue=0,
            underline=0,
        )

        for lang in STRING.keys():
            langMenu.add_radiobutton(
                label=lang,
                variable=self.LANG,
                value=lang,
                command=self.changeLang,
                underline=0,
            )

        self.prop = None
        self.gun = None
        self.errorLst = []

        self.columnconfigure(1, weight=1)
        self.rowconfigure(2, weight=3)
        self.rowconfigure(3, weight=2)

        self.addNamePlate()
        # self.addTopFrm()
        self.addLeftFrm()
        self.addRightFrm()
        self.addErrFrm()

        self.addPlotFrm()
        self.addAuxFrm()
        self.addspecFrm()
        self.addTblFrm()

        # self.update_idletasks()
        self.addGeomPlot()
        # self.update_idletasks()
        self.addFigPlot()
        # self.update_idletasks()
        self.addAuxPlot()

        self.ambCallback()
        self.cvlfCallback()
        self.typeCallback()
        self.insetCallback()
        self.ctrlCallback()

        self.updateTable()
        self.updateSpec()
        self.updateGeom()

        self.forceUpdOnThemeWidget.append(self.errorText)
        self.forceUpdOnThemeWidget.append(self.specs)

        # self.bind("<Configure>", self.resizePlot)

        parent.bind("<Return>", lambda *_: self.onCalculate())
        parent.protocol("WM_DELETE_WINDOW", self.quit)

        self.tLid = None
        self.timedLoop()

    def getDescriptive(self):
        if self.gun is None:
            return "Unknown Design"
        else:
            kwargs = self.kwargs
            typ = kwargs["typ"]
            cal = kwargs["caliber"]
            blr = kwargs["lengthGun"] / cal
            car_len = kwargs["chamberVolume"] / (
                pi * (0.5 * cal) ** 2 * kwargs["chambrage"]
            )
            w = kwargs["shotMass"]
            return (
                "{:} {:.3g}x{:.4g}mm L{:.0f} ".format(
                    "{:.4g} g".format(w * 1e3)
                    if w < 1
                    else "{:.4g} kg".format(w),
                    cal * 1e3,
                    car_len * 1e3,
                    blr,
                )
                + typ
            )

    def save(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(
                self.getLocStr("excTitle"), self.getLocStr("noDataMsg")
            )
            return

        fileName = filedialog.asksaveasfilename(
            title=self.getLocStr("saveLabel"),
            filetypes=(("JSON file", "*.json"),),
            defaultextension=".json",
            initialfile=self.getDescriptive(),
        )

        if fileName == "":
            messagebox.showinfo(
                self.getLocStr("excTitle"), self.getLocStr("cancelMsg")
            )
            return

        try:
            with open(fileName, "r", encoding="utf-8") as file:
                comment = json.load(file)["Comment"]
        except Exception:
            comment = None  # either file DNE or some exception occured during reading.

        try:
            locValDict = {
                loc.getDescriptive(): str(loc.get())
                for loc in self.locs
                if hasattr(loc, "getDescriptive")
            }

            if comment is None:
                locValDict.update(
                    {
                        "Comment": "Comment written here will be preserved even if file is overwritten!"
                    }
                )

            else:
                locValDict.update({"Comment": comment})

            with open(fileName, "w", encoding="utf-8") as file:
                json.dump(locValDict, file, indent="\t", ensure_ascii=False)

            messagebox.showinfo(
                self.getLocStr("sucTitle"),
                self.getLocStr("savedLocMsg").format(fileName),
            )

        except Exception as e:
            messagebox.showinfo(self.getLocStr("excTitle"), str(e))

    def load(self):
        fileName = filedialog.askopenfilename(
            title=self.getLocStr("loadLabel"),
            filetypes=(("JSON File", "*.json"),),
            defaultextension=".json",
            initialfile=self.getDescriptive(),
        )
        if fileName == "":
            messagebox.showinfo("Exception Loading Design", "No File Selected")
            return

        try:
            locDict = {
                loc.getDescriptive(): loc
                for loc in self.locs
                if hasattr(loc, "getDescriptive")
            }

            with open(fileName, "r", encoding="utf-8") as file:
                fileDict = json.load(file)

            for key, value in fileDict.items():
                try:
                    locDict[key].set(value)
                except KeyError:
                    pass

        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()

            if self.DEBUG.get() == 1:
                errMsg = "".join(
                    traceback.format_exception(
                        exc_type, exc_value, exc_traceback
                    )
                )
            else:
                errMsg = str(e)
            messagebox.showinfo(self.getLocStr("excTitle"), errMsg)

    def export(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(
                self.getLocStr("excTitle"), self.getLocStr("noDataMsg")
            )
            return

        fileName = filedialog.asksaveasfilename(
            title=self.getLocStr("exportLabel"),
            filetypes=(("Comma Separated File", "*.csv"),),
            defaultextension=".csv",
            initialfile=self.getDescriptive(),
        )
        gunType = self.kwargs["typ"]
        columnList = self.getLocStr("columnList")[gunType]

        if fileName == "":
            messagebox.showinfo(
                self.getLocStr("excTitle"), self.getLocStr("cancelMsg")
            )
            return
        try:
            with open(fileName, "w", encoding="utf-8", newline="") as csvFile:
                csvWriter = csv.writer(
                    csvFile, delimiter=",", quoting=csv.QUOTE_MINIMAL
                )

                csvWriter.writerow(columnList)

                for line in self.tableData:
                    csvWriter.writerow(line)

            messagebox.showinfo(
                self.getLocStr("sucTitle"),
                self.getLocStr("savedLocMsg").format(fileName),
            )

        except Exception as e:
            messagebox.showinfo(self.getLocStr("excTitle"), str(e))

    def changeLang(self):
        self.menubar.entryconfig(0, label=self.getLocStr("fileLabel"))
        self.menubar.entryconfig(1, label=self.getLocStr("themeLabel"))
        self.menubar.entryconfig(2, label=self.getLocStr("debugLabel"))

        self.fileMenu.entryconfig(0, label=self.getLocStr("saveLabel"))
        self.fileMenu.entryconfig(1, label=self.getLocStr("loadLabel"))
        self.fileMenu.entryconfig(2, label=self.getLocStr("exportLabel"))

        self.themeMenu.entryconfig(0, label=self.getLocStr("darkLabel"))
        self.themeMenu.entryconfig(1, label=self.getLocStr("lightLabel"))
        self.debugMenu.entryconfig(0, label=self.getLocStr("enableLabel"))

        self.calcButtonTip.set(self.getLocStr("calcButtonText"))

        for locWidget in self.locs:
            locWidget.reLocalize()

        self.calButton.config(text=self.getLocStr("calcLabel"))

        self.updateTable()
        self.updateGeom()
        self.updateSpec()
        self.updateFigPlot()
        self.updateAuxPlot()

    def readTable(self, tag):
        i = [line[0] for line in self.tableData].index(tag)
        return self.tableData[i]

    def getLocStr(self, name, forceDefault=None):
        try:
            if forceDefault:
                raise KeyError
            return STRING[self.LANG.get()][name]
        except KeyError:
            try:
                return STRING["English"][name]
            except KeyError:
                return name

    def addNamePlate(self):
        nameFrm = LocLabelFrame(
            self, locKey="nameFrm", locFunc=self.getLocStr, allLLF=self.locs
        )
        nameFrm.grid(row=0, column=0, sticky="nsew", columnspan=2)

        nameFrm.columnconfigure(0, weight=1)
        nameFrm.rowconfigure(0, weight=1)

        name = StringVar(self)

        namePlate = ttk.Entry(
            nameFrm,
            textvariable=name,
            state="disabled",
            justify="left",
            font=(FONTNAME, FONTSIZE + 2),
        )
        namePlate.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        self.name = name

    """
    def addTopFrm(self):
        topFrm = ttk.Frame(self)
        topFrm.grid(row=1, column=1, sticky="nsew")

        topFrm.columnconfigure(0, weight=1)
        topFrm.rowconfigure(0, weight=1)

        pltOptnFrm = LocLabelFrame(
            topFrm,
            locKey="pltOptnFrm",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        pltOptnFrm.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
    """

    def addLeftFrm(self):
        leftFrm = ttk.Frame(self)
        leftFrm.grid(row=1, column=0, rowspan=5, sticky="nsew")
        leftFrm.columnconfigure(0, weight=1)
        leftFrm.rowconfigure(0, weight=1)

        parFrm = LocLabelFrame(
            leftFrm,
            locKey="parFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        parFrm.grid(row=0, column=0, sticky="nsew")
        parFrm.columnconfigure(0, weight=1)

        i = 0

        self.lx = Loc122Disp(
            parent=parFrm,
            row=i,
            labelLocKey="lxLabel",
            tooltipLocKey="calLxText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 3
        self.ammo = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="ammoLabel",
            unitText="",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2
        self.va = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="vaLabel",
            tooltipLocKey="vinfText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2
        self.te = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="teffLabel",
            tooltipLocKey="teffText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2
        self.be = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="beffLabel",
            tooltipLocKey="beffText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2
        self.ld = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="ldLabel",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2

        self.pa = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="paLabel",
            # unitText="m/s²",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )

        i += 2
        self.gm = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="gmLabel",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )

        i += 2
        self.bnm = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="bmLabel",
            tooltipLocKey="bmText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )

    def addRightFrm(self):
        """
        rightFrm
        |-mecFrm
        |-solFrm
        |-envFrm
        |-opFrm
          |-consFrm
          |-sampleFrm

        """
        rightFrm = ttk.Frame(self)
        rightFrm.grid(row=0, column=3, rowspan=6, sticky="nsew")
        rightFrm.columnconfigure(0, weight=1)
        rightFrm.rowconfigure(0, weight=1)

        validationNN = self.register(validateNN)
        validationPI = self.register(validatePI)

        mecFrm = LocLabelFrame(
            rightFrm,
            locKey="matFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        mecFrm.grid(row=1, column=0, sticky="nsew")
        mecFrm.columnconfigure(0, weight=1)

        self.dropMat = LocDropdown(
            parent=mecFrm,
            strObjDict=MATERIALS,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="matFrmLabel",
        )
        self.dropMat.grid(
            row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )

        self.dropMatTemp = LocDropdown(
            parent=mecFrm,
            strObjDict=self.dropMat.getObj().getTdict(),
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="matFrmTempLabel",
        )
        self.dropMatTemp.grid(
            row=1, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )

        self.dropMat.trace_add(
            "write",
            lambda *args: self.dropMatTemp.reset(
                self.dropMat.getObj().getTdict()
            ),
        )

        self.ssf = Loc2Input(
            parent=mecFrm,
            row=2,
            col=0,
            labelLocKey="sffLabel",
            default="1.35",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        self.isAf = LocLabelCheck(
            parent=mecFrm,
            labelLocKey="afLabel",
            row=3,
            col=0,
            locFunc=self.getLocStr,
            allLC=self.locs,
            columnspan=2,
        )

        solFrm = LocLabelFrame(
            rightFrm,
            locKey="solFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        solFrm.grid(row=2, column=0, sticky="nsew")
        solFrm.columnconfigure(0, weight=1)

        self.dropSoln = LocDropdown(
            parent=solFrm,
            strObjDict=SOLUTIONS,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="solFrmLabel",
        )
        self.dropSoln.grid(
            row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )

        envFrm = LocLabelFrame(
            rightFrm,
            locKey="envFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        envFrm.grid(row=3, column=0, sticky="nsew")
        envFrm.columnconfigure(0, weight=1)

        i = 0
        self.inAtmos = LocLabelCheck(
            parent=envFrm,
            labelLocKey="atmosLabel",
            row=i,
            col=0,
            locFunc=self.getLocStr,
            allLC=self.locs,
            columnspan=3,
        )
        self.inAtmos.trace_add("write", self.ambCallback)

        i += 1
        self.ambP = Loc3Input(
            parent=envFrm,
            row=i,
            col=0,
            labelLocKey="ambPresLabel",
            unitText="kPa",
            default="101.325",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.ambRho = Loc3Input(
            parent=envFrm,
            row=i,
            col=0,
            labelLocKey="ambRhoLabel",
            unitText="kg/m³",
            default="1.204",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.ambGam = Loc3Input(
            parent=envFrm,
            row=i,
            col=0,
            labelLocKey="ambGamLabel",
            default="1.400",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        opFrm = LocLabelFrame(
            rightFrm,
            locKey="opFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        opFrm.grid(row=4, column=0, sticky="nsew")
        opFrm.columnconfigure(0, weight=1)
        i = 0

        consFrm = LocLabelFrame(
            opFrm,
            locKey="consFrmLabel",
            style="SubLabelFrame.TLabelframe",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        consFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        consFrm.columnconfigure(0, weight=1)

        j = 0
        self.solve_W_Lg = LocLabelCheck(
            parent=consFrm,
            row=j,
            col=0,
            columnspan=3,
            default=0,
            labelLocKey="consButton",
            tooltipLocKey="useConsText",
            locFunc=self.getLocStr,
            allLC=self.locs,
        )

        self.solve_W_Lg.trace_add("write", self.ctrlCallback)

        j += 1

        self.lock_Lg = LocLabelCheck(
            parent=consFrm,
            row=j,
            col=0,
            columnspan=3,
            default=0,
            labelLocKey="lockButton",
            tooltipLocKey="lockText",
            locFunc=self.getLocStr,
            allLC=self.locs,
        )

        self.lock_Lg.trace_add("write", self.ctrlCallback)

        j += 1

        self.opt_lf = LocLabelCheck(
            parent=consFrm,
            row=j,
            col=0,
            columnspan=3,
            default=0,
            labelLocKey="minTVButton",
            tooltipLocKey="optLFText",
            locFunc=self.getLocStr,
            allLC=self.locs,
        )

        self.opt_lf.trace_add("write", self.ctrlCallback)

        j += 1
        self.vTgt = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="vTgtLabel",
            unitText="m/s",
            default="1200.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        j += 1
        self.pTgt = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="pTgtLabel",
            unitText="MPa",
            default="350.0",
            validation=validationNN,
            tooltipLocKey="pTgtText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        j += 1
        self.pControl = LocDropdown(
            parent=consFrm,
            strObjDict=CONTROLS,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="pTgtLabel",
        )

        self.pControl.grid(
            row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        j += 1
        self.minWeb = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="minWebLabel",
            unitText="μm",
            default="1.0",
            validation=validationNN,
            color="red",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.lgmax = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="maxLgLabel",
            unitText="m",
            default="100.0",
            validation=validationNN,
            color="red",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        i += 1

        sampleFrm = LocLabelFrame(
            opFrm,
            locKey="sampleFrmLabel",
            style="SubLabelFrame.TLabelframe",
            tooltipLocKey="sampText",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        sampleFrm.grid(
            row=i, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )
        sampleFrm.columnconfigure(1, weight=1)

        self.dropDomain = LocDropdown(
            parent=sampleFrm,
            strObjDict=DOMAINS,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="sampleFrmLabel",
        )

        j = 0
        self.dropDomain.grid(
            row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )

        j += 1
        self.step = Loc2Input(
            parent=sampleFrm,
            row=j,
            col=0,
            labelLocKey="stepLabel",
            default="100",
            validation=validationNN,
            formatter=formatIntInput,
            reverse=True,
            anchor="center",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.accExp = Loc2Input(
            parent=opFrm,
            row=i,
            col=0,
            labelLocKey="-log10(ε)",
            default="3",
            validation=validationPI,
            formatter=formatIntInput,
            color="red",
            tooltipLocKey="tolText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.pbar = ttk.Progressbar(opFrm, mode="indeterminate", maximum=100)
        self.pbar.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        i += 1
        errorFrm = LocLabelFrame(
            opFrm,
            locKey="errFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        errorFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        errorFrm.columnconfigure(0, weight=1)
        errorFrm.rowconfigure(0, weight=1)

        errScroll = ttk.Scrollbar(errorFrm, orient="vertical")
        errScroll.grid(row=0, column=1, sticky="nsew")
        self.errorText = Text(
            errorFrm,
            yscrollcommand=errScroll.set,
            wrap="word",
            height=6,
            width=0,
            font=(FONTNAME, FONTSIZE),
        )

        self.errorText.grid(row=0, column=0, sticky="nsew")

        i += 1

        self.calButton = ttk.Button(
            opFrm,
            text=self.getLocStr("calcLabel"),
            command=self.onCalculate,
        )
        self.calButton.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        opFrm.rowconfigure(i, weight=1)
        self.calcButtonTip = StringVar(value=self.getLocStr("calcButtonText"))
        CreateToolTip(self.calButton, self.calcButtonTip)

    def onCalculate(self):
        for loc in self.locs:
            try:
                loc.inhibit()
            except AttributeError:
                pass

        self.process = None

        try:
            constrain = self.solve_W_Lg.get() == 1
            lock = self.lock_Lg.get() == 1
            optimize = self.opt_lf.get() == 1
            debug = self.DEBUG.get() == 1
            atmosphere = self.inAtmos.get() == 1
            autofrettage = self.isAf.get() == 1

            gunType = self.typeOptn.getObj()
            telescoped = self.ammoOptn.getObj() == TELESCOPED

            if self.prop is None:
                raise ValueError("Invalid propellant.")

            if self.useCv.getObj() == USE_CV:
                chamberVolume = float(self.cvL.get()) * 1e-3
            else:
                chamberVolume = (
                    float(self.chgkg.get())
                    / self.prop.rho_p
                    / float(self.ldf.get())
                    * 100
                )

            chambrage = float(self.clr.get())
            chargeMass = float(self.chgkg.get())
            caliber = float(self.calmm.get()) * 1e-3

            gunLength = float(self.tblmm.get()) * 1e-3
            loadFraction = float(self.ldf.get()) * 1e-2

            if telescoped:  # converting to equivalent cartridged gun.
                maxInset = float(self.insetmm.get()) * 1e-3
                boreS = pi * 0.25 * caliber**2
                breechS = chambrage * boreS

                exactlyTelescopedV = breechS * maxInset

                if chamberVolume > exactlyTelescopedV:
                    inset = maxInset
                else:
                    inset = chamberVolume / breechS
                    # chamberVolume *= (chambrage - 1) / chambrage

                # gunLength += inset
                insetV = inset * boreS
                chamberVolume -= insetV
                loadFraction *= 1 + insetV / chamberVolume
            else:
                maxInset = 0

            self.kwargs = {
                "opt": optimize,
                "con": constrain,
                "deb": debug,
                "lock": lock,
                "typ": gunType,
                "dom": self.dropDomain.getObj(),
                "sol": self.dropSoln.getObj(),
                "control": self.pControl.getObj(),
                "structuralMaterial": self.dropMat.getObj().createMaterialAtTemp(
                    self.dropMatTemp.getObj()
                ),
                "structuralSafetyFactor": float(self.ssf.get()),
                "caliber": caliber,
                "shotMass": float(self.shtkg.get()),
                "propellant": self.prop,
                "grainSize": float(self.arcmm.get()) * 1e-3,
                "chargeMass": chargeMass,
                "chargeMassRatio": (
                    float(self.chgkg.get()) / float(self.shtkg.get())
                ),
                "chamberVolume": chamberVolume,
                "startPressure": float(self.stpMPa.get()) * 1e6,
                "lengthGun": gunLength,
                "chambrage": chambrage,  # chamber expansion
                "nozzleExpansion": float(
                    self.nozzExp.get()
                ),  # nozzle expansion
                "nozzleEfficiency": float(self.nozzEff.get())
                * 1e-2,  # nozzle efficiency
                "dragCoefficient": float(self.dgc.get())
                * 1e-2,  # drag coefficient
                "designPressure": float(self.pTgt.get())
                * 1e6,  # design pressure
                "designVelocity": float(self.vTgt.get()),  # design velocity
                "tol": 10 ** -int(self.accExp.get()),
                "minWeb": 1e-6 * float(self.minWeb.get()),
                "maxLength": float(self.lgmax.get()),
                "loadFraction": loadFraction,
                "step": int(self.step.get()),
                "maxInset": maxInset,
                "autofrettage": autofrettage,
            }

            if atmosphere:
                self.kwargs.update(
                    {
                        "ambientP": float(self.ambP.get()) * 1e3,
                        "ambientRho": float(self.ambRho.get()),
                        "ambientGamma": float(self.ambGam.get()),
                    }
                )
            else:
                self.kwargs.update(
                    {"ambientP": 0, "ambientRho": 0, "ambientGamma": 1}
                )

            self.process = Process(
                target=calculate, args=(self.queue, self.kwargs)
            )
            # self.pos = 0
            self.calButton.config(state="disabled")
            self.pbar.start(interval=10)
            self.process.start()

        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            self.gun = None
            self.errorLst.append("Exception when dispatching calculation:")
            if self.DEBUG.get():
                self.errorLst.append(
                    "".join(
                        traceback.format_exception(
                            exc_type, exc_value, exc_traceback
                        )
                    )
                )
            else:
                self.errorLst.append(str(e))

            self.updateTable()
            self.updateError()
            self.updateFigPlot()
            self.updateAuxPlot()

    def getValue(self):
        queue = self.queue

        try:
            (
                self.kwargs,
                self.gun,
                self.tableData,
                self.errorData,
                self.pressureTrace,
                self.structure,
                errorLst,
            ) = queue.get_nowait()

            self.errorLst.extend(errorLst)

        except Empty:
            return

        self.process = None
        kwargs = self.kwargs

        constrain = kwargs["con"]
        lock = kwargs["lock"]
        optimize = kwargs["opt"]
        gunType = kwargs["typ"]
        maxInset = kwargs["maxInset"]
        caliber = kwargs["caliber"]
        chambrage = kwargs["chambrage"]
        trueLF = kwargs["loadFraction"]

        boreS = pi * 0.25 * caliber**2
        breechS = boreS * chambrage
        exactlyTelescopedVreal = (breechS - boreS) * maxInset

        sigfig = int(-log10(kwargs["tol"])) + 1
        gun = self.gun
        # kwargs["lengthGun"] -= inset

        if kwargs["chamberVolume"] > exactlyTelescopedVreal:
            inset = maxInset
        else:
            inset = kwargs["chamberVolume"] / (breechS - boreS)

        insetV = boreS * inset

        kwargs["loadFraction"] /= (
            1 + insetV / kwargs["chamberVolume"]
        )  # convert the "true load fraction" to "apparent"
        kwargs[
            "chamberVolume"
        ] += insetV  # convert the "real volume" returned to "apparent volume"

        if gun is not None:
            if constrain:
                webmm = roundSig(kwargs["grainSize"] * 1e3, n=sigfig)
                self.arcmm.set(webmm)

                if not lock:
                    lgmm = roundSig(kwargs["lengthGun"] * 1e3, n=sigfig)
                    self.tblmm.set(lgmm)

                if optimize:
                    if self.useCv.getObj() == USE_CV:
                        self.cvL.set(
                            roundSig(kwargs["chamberVolume"] * 1e3, n=sigfig)
                        )
                    else:
                        self.ldf.set(
                            roundSig(kwargs["loadFraction"] * 100, n=sigfig)
                        )  # corrected "bulk" load fraction

            self.ld.set(
                toSI(trueLF * self.prop.rho_p, useSN=False).strip() + " kg/m³"
            )
            # true load density

            _, _, _, _, vg, *_ = self.readTable(POINT_EXIT)

            eta_t, eta_b = gun.getEff(vg)

            self.te.set(str(round(eta_t * 100, 2)) + " %")
            self.be.set(str(round(eta_b * 100, 2)) + " %")

            self.va.set(toSI(gun.v_j, unit="m/s"))

            cartridge_len = (
                kwargs["chamberVolume"] / breechS
            )  # is equivalent to chamber length

            self.lx.set(
                toSI(kwargs["lengthGun"] / caliber, unit="Cal"),
                toSI(
                    (kwargs["lengthGun"] + cartridge_len / kwargs["chambrage"])
                    / kwargs["caliber"],
                    unit="Cal",
                ),
            )

            self.ammo.set(
                "{:.3g} x {:.4g} mm".format(caliber * 1e3, cartridge_len * 1e3)
            )

            if gunType == CONVENTIONAL:
                ps = self.readTable(POINT_PEAK_SHOT)[7]
            elif gunType == RECOILESS:
                ps = self.readTable(POINT_PEAK_SHOT)[9]

            self.pa.set(toSI(ps * gun.S / gun.m, unit="m/s²"))

            self.name.set(self.getDescriptive())

            bore_mass = self.structure[0]

            self.gm.set(
                formatMass(bore_mass) if bore_mass is not None else "N/A"
            )

            breech_nozzle_mass = self.structure[2]

            self.bnm.set(
                "N/A"
                if breech_nozzle_mass is None
                else formatMass(breech_nozzle_mass)
            )

        self.updateTable()
        self.updateError()
        self.updateFigPlot()
        self.updateAuxPlot()

        self.pbar.stop()
        self.calButton.config(state="normal")

        for loc in self.locs:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    def addErrFrm(self):
        pass

    def addspecFrm(self):
        specFrm = LocLabelFrame(
            self,
            locKey="specFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        specFrm.grid(row=0, column=2, rowspan=6, sticky="nsew")
        specFrm.columnconfigure(0, weight=1)

        # validation
        validationNN = self.register(validateNN)
        validationFLT = self.register(validateFLT)

        i = 0

        self.typeOptn = LocDropdown(
            parent=specFrm,
            strObjDict=TYPES,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="typeLabel",
        )

        self.typeOptn.grid(
            row=i, column=0, sticky="nsew", padx=2, pady=2, columnspan=3
        )
        i += 1

        self.calmm = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="calLabel",
            unitText="mm",
            default="50.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.tblmm = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="tblLabel",
            unitText="mm",
            default="3500.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.shtkg = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="shtLabel",
            unitText="kg",
            default="1.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        self.ammoOptn = LocDropdown(
            parent=specFrm,
            strObjDict=AMMUNITIONS,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="ammoLabel",
        )

        self.ammoOptn.grid(
            row=i, column=0, sticky="nsew", padx=2, pady=2, columnspan=3
        )
        i += 1

        self.insetmm = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="insetLabel",
            unitText="mm",
            default="100.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        self.chgkg = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="chgLabel",
            unitText="kg",
            default="1.0",
            validation=validationNN,
            tooltipLocKey="chgText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        propFrm = LocLabelFrame(
            specFrm,
            locKey="propFrmLabel",
            style="SubLabelFrame.TLabelframe",
            locFunc=self.getLocStr,
            tooltipLocKey="specsText",
            allLLF=self.locs,
        )
        propFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        specFrm.rowconfigure(i, weight=1)

        propFrm.rowconfigure(1, weight=1)
        propFrm.columnconfigure(0, weight=1)

        self.dropProp = LocDropdown(
            parent=propFrm,
            strObjDict=GrainComp.readFile(
                resolvepath("data/propellants.csv")
            ),  # dict of composition.name (string) -> composition (object),
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="propFrmLabel",
        )
        self.dropProp.grid(
            row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )

        specScroll = ttk.Scrollbar(propFrm, orient="vertical")
        specScroll.grid(
            row=1,
            column=1,
            sticky="nsew",
            pady=2,
        )
        specHScroll = ttk.Scrollbar(propFrm, orient="horizontal")
        specHScroll.grid(
            row=2,
            column=0,
            sticky="nsew",
        )

        self.specs = Text(
            propFrm,
            # wrap=WORD,
            wrap="none",
            height=5,
            width=35,
            yscrollcommand=specScroll.set,
            xscrollcommand=specHScroll.set,
            font=(FONTNAME, FONTSIZE),
        )
        self.specs.grid(row=1, column=0, sticky="nsew")
        specScroll.config(command=self.specs.yview)
        specHScroll.config(command=self.specs.xview)

        i += 1
        specFrm.rowconfigure(i, weight=1)

        grainFrm = LocLabelFrame(
            specFrm,
            locKey="grainFrmLabel",
            style="SubLabelFrame.TLabelframe",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        grainFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        grainFrm.columnconfigure(0, weight=1)
        grainFrm.rowconfigure(0, weight=1)

        geomPlotFrm = LocLabelFrame(
            grainFrm,
            locKey="σ(Z)",
            style="SubLabelFrame.TLabelframe",
            locFunc=self.getLocStr,
            tooltipLocKey="geomPlotText",
            allLLF=self.locs,
        )

        j = 0
        geomPlotFrm.grid(
            row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        self.geomParentFrm = grainFrm
        self.geomPlotFrm = geomPlotFrm
        self.dropGeom = LocDropdown(
            parent=grainFrm,
            strObjDict=GEOMETRIES,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="grainFrmLabel",
        )

        j += 1
        self.dropGeom.grid(
            row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        j += 1
        self.arcmm = Loc3Input(
            parent=grainFrm,
            row=j,
            labelLocKey="",
            descLabelKey="Web",
            unitText="mm",
            default="1.0",
            validation=validationNN,
            tooltipLocKey="",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.burnRateFudge = Loc3Input(
            parent=grainFrm,
            color="red",
            row=j,
            labelLocKey="fudgeLabel",
            unitText="%",
            default="0.0",
            validation=validationFLT,
            tooltipLocKey="fudgeText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.grdR = Loc3Input(
            parent=grainFrm,
            row=j,
            labelLocKey="",
            descLabelKey="1/α",
            unitText="x",
            default="1.0",
            validation=validationNN,
            tooltipLocKey="",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.grlR = Loc3Input(
            parent=grainFrm,
            row=j,
            labelLocKey="",
            descLabelKey="1/β",
            unitText="x",
            default="2.5",
            validation=validationNN,
            tooltipLocKey="",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        i += 1
        self.useCv = LocDropdown(
            parent=specFrm,
            strObjDict=USE,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="cvlfLabel",
        )
        self.useCv.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        i += 1
        self.ldf = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="ldfLabel",
            unitText="%",
            default="50.0",
            validation=validationNN,
            tooltipLocKey="ldfText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.cvL = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="cvLabel",
            unitText="L",
            default="20",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        self.useCv.trace_add("write", self.cvlfCallback)
        self.cvL.trace_add("write", self.cvlfConsisCallback)
        self.ldf.trace_add("write", self.cvlfConsisCallback)
        self.chgkg.trace_add("write", self.cvlfConsisCallback)
        self.accExp.trace_add("write", self.cvlfConsisCallback)

        i += 1
        self.clr = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="clrLabel",
            unitText="x",
            default="1.5",
            validation=validationNN,
            tooltipLocKey="clrText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.dgc = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="dgcLabel",
            unitText="%",
            default="5.0",
            validation=validationNN,
            tooltipLocKey="dgcText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.stpMPa = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="stpLabel",
            unitText="MPa",
            default="10",
            validation=validationNN,
            tooltipLocKey="stpText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.nozzExp = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="nozzExpLabel",
            unitText="x",
            default="4",
            validation=validationNN,
            tooltipLocKey="nozzExpText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.nozzEff = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="nozzEffLabel",
            unitText="%",
            default="92.0",
            validation=validationNN,
            tooltipLocKey="nozzEffText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        self.dropProp.trace_add("write", self.updateSpec)
        self.dropGeom.trace_add("write", self.updateGeom)

        self.burnRateFudge.trace_add("write", self.callback)
        self.grdR.trace_add("write", self.callback)
        self.grlR.trace_add("write", self.callback)
        self.arcmm.trace_add("write", self.callback)

        self.typeOptn.trace_add("write", self.typeCallback)
        self.ammoOptn.trace_add("write", self.insetCallback)

    def addGeomPlot(self):
        geomPlotFrm = self.geomPlotFrm
        geomPlotFrm.columnconfigure(0, weight=1)
        geomPlotFrm.rowconfigure(0, weight=1)
        geomPlotFrm.grid_propagate(False)

        with mpl.rc_context(GEOM_CONTEXT):
            fig = Figure(dpi=96)
            self.geomFig = fig
            self.geomAx = fig.add_subplot(111)

            self.geomCanvas = FigureCanvasTkAgg(fig, master=geomPlotFrm)
            self.geomCanvas.get_tk_widget().grid(
                row=0, column=0, padx=0, pady=0, sticky="nsew"
            )
            fig.set_layout_engine(None)
            self.geomCanvas.draw_idle()

    def addPlotFrm(self):
        plotFrm = LocLabelFrame(
            self,
            locKey="plotFrmLabel",
            width=640,
            height=300,
            locFunc=self.getLocStr,
            tooltipLocKey="plotText",
            allLLF=self.locs,
        )
        plotFrm.grid(row=2, column=1, padx=2, pady=2, sticky="nsew")

        self.plotFrm = plotFrm

        for i in range(5):
            plotFrm.columnconfigure(i, weight=1)

        checks = []

        j = 1
        k = 0
        self.plotAvgP = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotAvgP",
            locFunc=self.getLocStr,
            allLC=checks,
        )
        k += 1
        self.plotBaseP = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotBaseP",
            locFunc=self.getLocStr,
            allLC=checks,
        )

        k += 1
        self.plotBreechNozzleP = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            descLabelKey="Breech/Nozzle Pressure",
            locFunc=self.getLocStr,
            allLC=checks,
        )
        k += 1
        self.plotStagP = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotStagP",
            locFunc=self.getLocStr,
            allLC=checks,
        )

        j += 1
        k = 0

        self.plotVel = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotVel",
            locFunc=self.getLocStr,
            allLC=checks,
        )
        k += 1
        self.plotNozzleV = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotNozzleV",
            locFunc=self.getLocStr,
            allLC=checks,
        )
        k += 1
        self.plotBurnup = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotBurnup",
            locFunc=self.getLocStr,
            allLC=checks,
        )
        k += 1
        self.plotEta = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotEta",
            locFunc=self.getLocStr,
            allLC=checks,
        )
        k += 1
        self.plotRecoil = LocLabelCheck(
            parent=plotFrm,
            row=j,
            col=k,
            labelLocKey="plotRecoil",
            locFunc=self.getLocStr,
            allLC=checks,
        )

        for check in checks:
            check.trace_add("write", self.updateFigPlot)

        self.locs.extend(checks)

    def addFigPlot(self):
        plotFrm = self.plotFrm
        plotFrm.columnconfigure(0, weight=1)
        plotFrm.rowconfigure(0, weight=1)
        plotFrm.grid_propagate(False)

        with mpl.rc_context(FIG_CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            axes = fig.add_subplot(111)

            ax = axes
            axP = ax.twinx()
            axv = ax.twinx()
            axF = ax.twinx()

            ax.yaxis.tick_right()
            axF.yaxis.tick_left()
            axv.yaxis.tick_left()

            ax.set_xlabel(" ")
            axP.spines.right.set_position(("data", 0.5))
            axF.spines.left.set_position(("data", 0.5))

            axP.yaxis.set_ticks(axP.get_yticks()[1:-1:])
            axF.yaxis.set_ticks(axF.get_yticks()[1:-1:])

            self.ax = ax
            self.axP = axP
            self.axv = axv
            self.axF = axF
            self.fig = fig

            self.pltCanvas = FigureCanvasTkAgg(fig, master=plotFrm)
            self.pltCanvas.get_tk_widget().grid(
                row=0, column=0, padx=2, pady=2, sticky="nsew", columnspan=5
            )

            self.pltCanvas.draw_idle()

    def addAuxFrm(self):
        auxFrm = LocLabelFrame(
            self,
            locKey="auxFrmLabel",
            width=640,
            height=200,
            locFunc=self.getLocStr,
            tooltipLocKey="auxText",
            allLLF=self.locs,
        )
        auxFrm.grid(row=3, column=1, padx=2, pady=2, sticky="nsew")
        self.auxFrm = auxFrm

        for i in range(2):
            auxFrm.columnconfigure(i, weight=1)

        auxChecks = []
        j = 1
        k = 0
        self.traceHull = LocLabelCheck(
            parent=auxFrm,
            row=j,
            col=k,
            labelLocKey="traceHull",
            locFunc=self.getLocStr,
            default=0,
            allLC=auxChecks,
        )

        k += 1
        self.tracePress = LocLabelCheck(
            parent=auxFrm,
            row=j,
            col=k,
            labelLocKey="tracePress",
            locFunc=self.getLocStr,
            allLC=auxChecks,
        )

        for check in auxChecks:
            check.trace_add("write", self.updateAuxPlot)

        self.locs.extend(auxChecks)

    def addAuxPlot(self):
        auxFrm = self.auxFrm
        auxFrm.columnconfigure(0, weight=1)
        auxFrm.rowconfigure(0, weight=1)

        auxFrm.grid_propagate(False)

        with mpl.rc_context(FIG_CONTEXT):
            fig = Figure(dpi=96, layout="constrained")

            axes = fig.add_subplot(111)

            ax = axes
            axH = ax.twinx()

            self.auxAx = ax  # auxiliary axes
            self.auxAxH = axH
            self.auxFig = fig
            self.auxAxH.yaxis.tick_right()
            self.auxCanvas = FigureCanvasTkAgg(fig, master=auxFrm)
            self.auxCanvas.get_tk_widget().grid(
                row=0, column=0, padx=2, pady=2, sticky="nsew", columnspan=2
            )

            self.auxCanvas.draw_idle()

    def timedLoop(self):
        if self.process is not None:
            self.getValue()

        self.tLid = self.after(100, self.timedLoop)

    def quit(self):
        root = self.parent
        if self.tLid is not None:
            root.after_cancel(self.tLid)
        root.quit()

    def updateFigPlot(self, *args):
        gun = self.gun
        if gun is None:
            with mpl.rc_context(FIG_CONTEXT):
                self.ax.cla()
                self.axP.cla()
                self.axv.cla()
                self.axF.cla()
                self.pltCanvas.draw_idle()
            return

        with mpl.rc_context(FIG_CONTEXT):
            self.ax.cla()
            self.axP.cla()
            self.axv.cla()
            self.axF.cla()

            vTgt = self.kwargs["designVelocity"]
            gunType = self.kwargs["typ"]
            dom = self.kwargs["dom"]

            xs, vs = [], []

            Pa_m = []
            Pas, Pss, Pbs, P0s = [], [], [], []
            Frs = []
            psis, etas = [], []
            vxs = []

            if gunType == CONVENTIONAL:
                for tag, t, l, psi, v, Pb, P, Ps, T in self.tableData:
                    if tag == POINT_PEAK_AVG:
                        if dom == DOMAIN_TIME:
                            xPeak = t * 1e3
                        elif dom == DOMAIN_LENG:
                            xPeak = l

                    if dom == DOMAIN_TIME:
                        xs.append(t * 1000)
                    elif dom == DOMAIN_LENG:
                        xs.append(l)

                    Fr = P * gun.S
                    vs.append(v)
                    Pas.append(P / 1e6)
                    Pss.append(Ps / 1e6)
                    Pbs.append(Pb / 1e6)
                    Frs.append(Fr / 1e6)
                    psis.append(psi)

            elif gunType == RECOILESS:
                for (
                    tag,
                    t,
                    l,
                    psi,
                    v,
                    vx,
                    Px,
                    P0,
                    P,
                    Ps,
                    T,
                    eta,
                ) in self.tableData:
                    if tag == POINT_PEAK_AVG:
                        if dom == DOMAIN_TIME:
                            xPeak = t * 1e3
                        elif dom == DOMAIN_LENG:
                            xPeak = l
                        Pa_m.append("^")
                    elif tag == "*":
                        Pa_m.append("x")
                    else:
                        Pa_m.append(",")

                    if dom == DOMAIN_TIME:
                        xs.append(t * 1000)
                    elif dom == DOMAIN_LENG:
                        xs.append(l)

                    Fr = P * gun.S * (1 - gun.C_f * gun.S_j_bar)
                    vs.append(v)
                    vxs.append(vx)
                    Pas.append(P / 1e6)
                    Pss.append(Ps / 1e6)
                    Pbs.append(Px / 1e6)
                    P0s.append(P0 / 1e6)
                    Frs.append(Fr / 1e6)
                    psis.append(psi)
                    etas.append(eta)

            self.axP.spines.right.set_position(("data", xPeak))
            self.axF.spines.left.set_position(("data", xPeak))

            if self.plotBreechNozzleP.get():
                self.axP.plot(
                    xs,
                    Pbs,
                    c="xkcd:goldenrod",
                    label=self.getLocStr("figBreech")
                    if gunType == CONVENTIONAL
                    else self.getLocStr("figNozzleP")
                    if gunType == RECOILESS
                    else "",
                )

            if gunType == RECOILESS:
                if self.plotStagP.get():
                    self.axP.plot(
                        xs,
                        P0s,
                        "seagreen",
                        label=self.getLocStr("figStagnation"),
                    )

                if self.plotNozzleV.get():
                    self.axv.plot(
                        xs, vxs, "royalblue", label=self.getLocStr("figNozzleV")
                    )

                if self.plotEta.get():
                    self.ax.plot(
                        xs, etas, "crimson", label=self.getLocStr("figOutflow")
                    )

            if self.plotAvgP.get():
                self.axP.plot(
                    xs, Pas, "tab:green", label=self.getLocStr("figAvgP")
                )

            if self.plotBaseP.get():
                self.axP.plot(
                    xs, Pss, "yellowgreen", label=self.getLocStr("figShotBase")
                )

            Pd = float(self.pTgt.get())
            self.axP.axhline(
                Pd,
                c="tab:green",
                linestyle=":",
                label=self.getLocStr("figTgtP"),
            )

            if self.plotVel.get():
                self.axv.plot(
                    xs,
                    vs,
                    "tab:blue",
                    label=self.getLocStr("figShotVel"),
                )
            self.axv.axhline(
                vTgt,
                c="tab:blue",
                linestyle=":",
                label=self.getLocStr("figTgtV"),
            )

            if self.plotBurnup.get():
                self.ax.plot(
                    xs, psis, c="tab:red", label=self.getLocStr("figPsi")
                )

            if self.plotRecoil.get():
                if gunType == CONVENTIONAL:
                    self.axF.plot(
                        xs,
                        Frs,
                        c="tab:green",
                        label=self.getLocStr("figRecoil"),
                        linestyle="dotted",
                    )
                elif gunType == RECOILESS:
                    self.axF.plot(
                        xs,
                        tuple(-v for v in Frs),
                        c="tab:green",
                        label="-" + self.getLocStr("figRecoil"),
                        linestyle="dotted",
                    )

            linesLabeled = []
            for lines, xvals in zip(
                (
                    self.axP.get_lines(),
                    self.ax.get_lines(),
                    self.axv.get_lines(),
                    self.axF.get_lines(),
                ),
                (
                    (0.2 * xs[-1] + 0.8 * xPeak, xs[-1]),
                    (0, xs[-1]),
                    (xPeak, 0.2 * xs[-1] + 0.8 * xPeak),
                    (0, xPeak),
                ),
            ):
                labelLines(lines, align=True, xvals=xvals, outline_width=4)
                linesLabeled.append(lines)

            self.ax.set_xlim(left=0, right=xs[-1])
            pmax = max(Pas + Pbs + Pss + P0s)
            self.axP.set(ylim=(0, pmax * 1.1))
            self.axF.set(ylim=(0, pmax * gun.S * 1.1))
            self.axv.set(ylim=(0, max(vs + vxs) * 1.15))
            self.ax.set_ylim(bottom=0, top=1.05)

            fmax = pmax * gun.S

            self.axP.yaxis.set_ticks(
                [v for v in self.axP.get_yticks() if v <= pmax][1:]
            )
            self.axF.yaxis.set_ticks(
                [v for v in self.axF.get_yticks() if v <= fmax][1:]
            )

            self.ax.yaxis.tick_right()
            self.axF.yaxis.tick_left()
            self.axv.yaxis.tick_left()

            tkw = dict(size=4, width=1.5)

            self.ax.tick_params(axis="y", colors="tab:red", **tkw)
            self.axv.tick_params(axis="y", colors="tab:blue", **tkw)
            self.axP.tick_params(axis="y", colors="tab:green", **tkw)
            self.axF.tick_params(axis="y", colors="tab:green", **tkw)
            self.ax.tick_params(axis="x", **tkw)

            if dom == DOMAIN_TIME:
                self.ax.set_xlabel(self.getLocStr("figTimeDomain"))
            elif dom == DOMAIN_LENG:
                self.ax.set_xlabel(self.getLocStr("figLengDomain"))

            self.fig.set_layout_engine("constrained")

            self.pltCanvas.draw_idle()

    def updateAuxPlot(self, *args):
        gun = self.gun
        if gun is None:
            with mpl.rc_context(FIG_CONTEXT):
                self.auxAx.cla()
                self.auxAxH.cla()
                self.auxCanvas.draw_idle()
            return

        with mpl.rc_context(FIG_CONTEXT):
            self.auxAx.cla()
            self.auxAxH.cla()
            gunType = self.kwargs["typ"]

            pTrace = self.pressureTrace

            cmap = mpl.colormaps["afmhot"]

            x_max = 0
            y_max = 0
            T_max = max(trace[2] for trace in pTrace)
            T_min = min(trace[2] for trace in pTrace)

            for i, (tag, psi, T, trace) in enumerate(pTrace[::-1]):
                if tag != "":
                    continue
                if self.themeRadio.get():
                    color = cmap(1 - (T - T_min) / (T_max - T_min))
                else:
                    color = cmap((T - T_min) / (T_max - T_min))

                alpha = None
                linestyle = None

                x, y = zip(*trace)
                y = [v * 1e-6 for v in y]
                x_max = max(x_max, max(x))
                y_max = max(y_max, max(y))

                if self.tracePress.get():
                    self.auxAx.plot(x, y, c=color, alpha=alpha, ls=linestyle)

            self.auxAx.set_xlim(left=0, right=x_max)
            self.auxAx.set_ylim(bottom=0, top=y_max * 1.15)

            l_c = gun.l_c
            chi_k = gun.chi_k
            r = 0.5 * gun.caliber * 1e3
            self.auxAx.plot(
                (l_c, l_c), (0, y_max * 1.15), c="grey", ls="dotted", zorder=1.5
            )

            tkw = dict(size=4, width=1.5)
            self.auxAx.tick_params(axis="y", colors="tab:green", **tkw)
            self.auxAx.tick_params(axis="x", **tkw)

            self.auxAx.set_xlabel(self.getLocStr("figAuxDomain"))

            HTrace = self.structure[1]

            if HTrace is not None and self.traceHull.get():
                xHull, rHull = zip(*HTrace)
                rHull = [r * 1e3 for r in rHull]
                inline = [chi_k**0.5 * r if x < l_c else r for x in xHull]
                rRim = rHull[0]

                BTrace = self.structure[3]
                if gunType == CONVENTIONAL and BTrace is not None:
                    xBreech, rBreech = zip(*BTrace)
                    rBreech = [r * 1e3 for r in rBreech]
                    rFace = rBreech[-1]

                    self.auxAxH.plot(xBreech, rBreech, c="tab:green")
                    self.auxAxH.fill_between(
                        xBreech,
                        rBreech,
                        [0 for _ in xBreech],
                        alpha=0.5 if self.tracePress.get() else 0.8,
                        color="tab:green",
                    )
                    self.auxAx.set_xlim(left=min(xBreech))

                    xHull = list(xBreech) + list(xHull)
                    rHull = [rRim for _ in xBreech] + rHull
                    inline = [rFace for _ in xBreech] + inline

                elif gunType == RECOILESS and BTrace is not None:
                    xNozzle, rIn, rOut = zip(*BTrace)
                    rIn = [r * 1e3 for r in rIn]
                    rOut = [r * 1e3 for r in rOut]

                    self.auxAxH.plot(xNozzle, rIn, c="tab:green")
                    self.auxAxH.plot(xNozzle, rOut, c="tab:green")

                    self.auxAxH.fill_between(
                        xNozzle,
                        rIn,
                        rOut,
                        alpha=0.5 if self.tracePress.get() else 0.8,
                        color="tab:green",
                    )

                    self.auxAx.set_xlim(left=min(xNozzle))

                self.auxAxH.plot(xHull, rHull, c="tab:blue")
                self.auxAxH.plot(xHull, inline, c="tab:blue")

                self.auxAxH.fill_between(
                    xHull,
                    rHull,
                    inline,
                    alpha=0.5 if self.tracePress.get() else 0.8,
                    color="tab:blue",
                )

            self.auxAxH.set_ylim(bottom=0)

            self.auxFig.set_layout_engine("constrained")
            self.auxCanvas.draw_idle()

    def addTblFrm(self):
        tblFrm = LocLabelFrame(
            self, locKey="tblFrmLabel", locFunc=self.getLocStr, allLLF=self.locs
        )
        tblFrm.grid(row=4, column=1, sticky="nsew")

        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(0, weight=1)
        # configure the numerical
        self.tv = ttk.Treeview(
            tblFrm, selectmode="browse", height=8
        )  # this set the nbr. of values
        self.tv.grid(row=0, column=0, sticky="nsew")

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

    def updateSpec(self, *args):
        self.specs.config(state="normal")
        compo = self.dropProp.getObj()
        self.specs.delete("1.0", "end")
        t_Font = tkFont.Font(font=self.specs.cget("font"))
        width = self.specs.winfo_width() // t_Font.measure("m")
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} K {:}\n".format(
                self.getLocStr("TvDesc"),
                compo.T_v,
                self.getLocStr("isochorDesc"),
            ),
        ),

        self.specs.insert(
            "end",
            "{:}: {:>4.0f} kg/m³\n".format(
                self.getLocStr("densityDesc"), compo.rho_p
            ),
        )
        isp = compo.getIsp()
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} m/s {:>3.0f} s\n".format(
                self.getLocStr("vacISPDesc"), isp, isp / 9.805
            ),
        )
        isp = compo.getIsp(50)
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} m/s {:>3.0f} s\n{:}\n".format(
                self.getLocStr("atmISPDesc"),
                isp,
                isp / 9.805,
                self.getLocStr("pRatioDesc"),
            ),
        )
        self.specs.insert("end", "{:}:\n".format(self.getLocStr("brDesc")))
        for p in (100e6, 200e6, 300e6):
            self.specs.insert(
                "end",
                "{:>12}".format(toSI(compo.getLBR(p), unit="m/s", dec=3))
                + " @ {:>12}\n".format(toSI(p, unit="Pa", dec=3)),
            )

        self.specs.insert(
            "end",
            "-" * width + "\n",
        )
        for line in compo.desc.split(","):
            self.specs.insert("end", line + "\n")
        self.specs.config(state="disabled")
        # this updates the specification description

        self.callback()
        self.cvlfConsisCallback()  # update the chamber volume / load fraction with current data

        return True

    def updateGeom(self, *args):
        geom = self.dropGeom.getObj()
        if geom == SimpleGeometry.SPHERE:
            self.grlR.remove()
            self.grdR.remove()

        elif geom == SimpleGeometry.CYLINDER:
            self.grlR.restore()
            self.grdR.remove()

        else:
            self.grlR.restore()
            self.grdR.restore()

        if geom == SimpleGeometry.SPHERE:
            self.arcmm.reLocalize("diamLabel", "diaText")
            self.grlR.reLocalize("", "")
            self.grdR.reLocalize("", "")

        elif geom == SimpleGeometry.ROD:
            self.arcmm.reLocalize("widthLabel", "widthText")
            self.grlR.reLocalize("ltwLabel", "rodRText")
            self.grdR.reLocalize("htwLabel", "heightRText")

        elif geom == SimpleGeometry.CYLINDER:
            self.arcmm.reLocalize("diamLabel", "diaText")
            self.grlR.reLocalize("ltdLabel", "cylLRText")
            self.grdR.reLocalize("", "")

        else:
            self.arcmm.reLocalize("athLabel", "arcText")
            self.grlR.reLocalize("ltdLabel", "perfLRText")
            self.grdR.reLocalize("pdtalLabel", "pDiaRText")

        self.callback()

        return True

    def updateGeomPlot(self):
        with mpl.rc_context(GEOM_CONTEXT):
            N = 10
            prop = self.prop
            Zb = prop.Z_b
            self.geomAx.cla()
            if prop is not None:
                xs = [i / N for i in range(N)]
                ys = [prop.f_sigma_Z(x) for x in xs]

                if Zb > 1:
                    xs.extend((1, 1))
                    ys.extend(prop.f_ullim())

                xs.append(Zb)
                ys.append(prop.f_sigma_Z(Zb))

                xs.append(xs[-1])
                ys.append(0)

                self.geomAx.plot(xs, ys)
                self.geomAx.grid(
                    which="major", color="grey", linestyle="dotted"
                )
                self.geomAx.minorticks_on()
                self.geomAx.set_xlim(left=0, right=min(prop.Z_b, 2))
                self.geomAx.xaxis.set_ticks(
                    [i * 0.5 for i in range(ceil(min(prop.Z_b, 2) / 0.5) + 1)]
                )
                self.geomAx.set_ylim(bottom=0, top=max(ys))
                self.geomAx.yaxis.set_ticks(
                    [i * 0.5 for i in range(ceil(max(ys) / 0.5) + 1)]
                )

            self.geomFig.set_layout_engine("constrained")
            self.geomCanvas.draw_idle()

    def updateError(self):
        self.errorText.delete("1.0", "end")
        for line in self.errorLst:
            self.errorText.insert("end", line + "\n")
        self.errorLst = []

    def updateTable(self):
        self.tv.delete(*self.tv.get_children())

        try:
            gunType = self.kwargs["typ"]
            tableData, errorData = self.tableData, self.errorData
        except AttributeError:
            gunType = self.typeOptn.getObj()
            tableData, errorData = [], []

        if gunType == CONVENTIONAL:
            # fmt: off
            useSN = (
                False, False, False, True, False, False, False, False, True
            )
            # fmt: on

            units = (None, "s", "m", None, "m/s", "Pa", "Pa", "Pa", "K")

        elif gunType == RECOILESS:
            # fmt: off
            useSN = (
                False, False, False, True, False, False, False, False, False,
                False, True, True
            )
            units = (
                None, "s", "m", None, "m/s", "m/s", "Pa", "Pa", "Pa", "Pa", "K",
                None
            )
            # fmt: on

        tableData = dot_aligned(
            tableData,
            units=units,
            useSN=useSN,
        )

        errorData = dot_aligned(errorData, units=units, useSN=useSN)

        columnList = self.getLocStr("columnList")[gunType]
        self.tv["columns"] = columnList
        self.tv["show"] = "headings"
        self.tv.tag_configure(POINT_PEAK_AVG, foreground="#2ca02c")
        self.tv.tag_configure(POINT_PEAK_BREECH, foreground="orange")
        self.tv.tag_configure(POINT_PEAK_SHOT, foreground="yellow green")
        self.tv.tag_configure(POINT_BURNOUT, foreground="red")
        self.tv.tag_configure(POINT_FRACTURE, foreground="brown")
        self.tv.tag_configure(POINT_EXIT, foreground="steel blue")
        self.tv.tag_configure(POINT_START, foreground="steel blue")
        self.tv.tag_configure("*", foreground="tan")

        t_Font = tkFont.Font(family=FONTNAME, size=FONTSIZE)

        self.tv.tag_configure("monospace", font=t_Font)
        self.tv.tag_configure("error", font=t_Font, foreground="dim gray")

        # we use a fixed width font so any char will do
        fontWidth, _ = t_Font.measure("m"), t_Font.metrics("linespace")

        winWidth = self.tv.winfo_width()
        width = int(winWidth / len(self.tv["columns"]))

        for i, column in enumerate(columnList):  # foreach column
            self.tv.heading(
                i, text=column, anchor="e"
            )  # let the column heading = column name
            self.tv.column(
                column,
                stretch=True,  # will adjust to window resizing
                width=width,
                minwidth=fontWidth * 14,
                anchor="e",
            )
        for i, (row, erow) in enumerate(zip(tableData, errorData)):
            self.tv.insert(
                "", "end", str(i + 1), values=row, tags=(row[0], "monospace")
            )
            self.tv.insert(
                str(i + 1),
                "end",
                str(-i - 1),
                values=tuple("±" + e if e != erow[0] else e for e in erow),
                tags="error",
            )
            self.tv.move(str(-i - 1), str(i + 1), "end")

    def callback(self, *args):
        """
        updates the propellant object on write to the ratio entry fields
        and, on changing the propellant or geometrical specification.

        Double calling is due to value validation, no workaround has been
        found at this time!
        """

        geom = self.dropGeom.getObj()
        compo = self.dropProp.getObj()

        try:
            self.prop = Propellant(
                compo,
                geom,
                float(self.grdR.get()),
                float(self.grlR.get()),
                float(self.burnRateFudge.get())
                * 1e-2,  # convert percentage to float
            )
            self.updateGeomPlot()

        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            self.prop = None
            if self.DEBUG.get():
                self.errorLst.append(
                    "".join(
                        traceback.format_exception(
                            exc_type, exc_value, exc_traceback
                        )
                    )
                )
            else:
                self.errorLst.append(str(e))

        self.updateError()

    def typeCallback(self, *args):
        gunType = self.typeOptn.getObj()

        if gunType == CONVENTIONAL:
            self.nozzExp.remove()
            self.nozzEff.remove()

            self.plotBreechNozzleP.reLocalize("plotBreechP")
            self.bnm.reLocalize("bmLabel", "bmText")

            self.dropSoln.enable()

            for widget in (self.plotEta, self.plotNozzleV):
                widget.disable()

        elif gunType == RECOILESS:
            self.nozzExp.restore()
            self.nozzEff.restore()

            self.plotBreechNozzleP.reLocalize("plotNozzleP")
            self.bnm.reLocalize("nmLabel", "")

            self.dropSoln.setByObj(SOL_LAGRANGE)
            self.dropSoln.disable()

            for widget in (self.plotEta, self.plotNozzleV):
                widget.enable()

    def ctrlCallback(self, *args):
        if self.solve_W_Lg.get() == 0:
            self.vTgt.disable()
            self.pTgt.disable()
            self.opt_lf.disable()
            self.lock_Lg.disable()
            self.minWeb.disable()
            self.lgmax.disable()
            self.pControl.disable()
            self.tblmm.enable()

        else:
            if self.lock_Lg.get() == 1:
                self.vTgt.disable()
                self.opt_lf.disable()
                self.tblmm.enable()
            else:
                self.vTgt.enable()
                self.opt_lf.enable()
                self.tblmm.disable()

            if self.opt_lf.get() == 1:
                self.lock_Lg.disable()
            else:
                self.lock_Lg.enable()

            self.pTgt.enable()
            self.minWeb.enable()
            self.lgmax.enable()
            self.pControl.enable()

    def cvlfConsisCallback(self, *args):
        try:
            sigfig = int(self.accExp.get()) + 1
            if self.useCv.getObj() == USE_CV:  # use Cv
                cv = float(self.cvL.get()) / 1e3
                w = float(self.chgkg.get())
                rho = self.prop.rho_p
                self.ldf.set(roundSig(w / cv / rho * 100, n=sigfig))

            else:  # using load fraction
                w = float(self.chgkg.get())
                lf = float(self.ldf.get()) / 100
                rho = self.prop.rho_p
                self.cvL.set(roundSig((w / rho / lf) * 1e3, n=sigfig))

        except (ZeroDivisionError, ValueError):
            return

    def ambCallback(self, *args):
        self.ambP.enable() if self.inAtmos.get() else self.ambP.disable()
        self.ambRho.enable() if self.inAtmos.get() else self.ambRho.disable()
        self.ambGam.enable() if self.inAtmos.get() else self.ambGam.disable()

    def cvlfCallback(self, *args):
        useCv = self.useCv.getObj() == USE_CV

        self.ldf.disable() if useCv else self.ldf.enable()
        self.cvL.enable() if useCv else self.cvL.disable()

    def insetCallback(self, *args):
        isCartridge = self.ammoOptn.getObj() == CARTRIDGE
        self.insetmm.remove() if isCartridge else self.insetmm.restore()

    def useTheme(self):
        style = ttk.Style(self)
        choice = self.themeRadio.get()
        if choice == 0:
            style.theme_use("awdark")
        elif choice == 1:
            style.theme_use("awlight")

        self.setTheme()

    def setTheme(self):
        dpi = self.dpi

        style = ttk.Style(self)
        # ensure that the treeview rows are roughly the same height
        # regardless of dpi. on Windows, default is Segoe UI at 9 points
        # so the default row height should be around 12

        style.configure(
            "Treeview", rowheight=round(12 * (FONTSIZE / 8) * dpi / 72.0)
        )
        style.configure("Treeview.Heading", font=(FONTNAME, FONTSIZE))
        style.configure("TButton", font=(FONTNAME, FONTSIZE + 2, "bold"))
        style.configure(
            "TLabelframe.Label", font=(FONTNAME, FONTSIZE + 2, "bold")
        )
        style.configure(
            "SubLabelFrame.TLabelframe.Label",
            font=(FONTNAME, FONTSIZE + 2),
        )
        style.configure("TCheckbutton", font=(FONTNAME, FONTSIZE))

        style.configure(
            "TNotebook.Tab",
            padding=[10, 2],
            font=(FONTNAME, FONTSIZE + 2, "bold"),
        )

        bgc = str(style.lookup("TFrame", "background"))
        fgc = str(style.lookup("TFrame", "foreground"))
        ebgc = str(style.lookup("TCombobox", "fieldbackground"))

        GEOM_CONTEXT.update(
            {
                "xtick.color": fgc,
                "ytick.color": fgc,
                "axes.edgecolor": fgc,
                "axes.facecolor": ebgc,
                "figure.facecolor": bgc,
                "figure.edgecolor": fgc,
            }
        )

        FIG_CONTEXT.update(
            {
                "figure.facecolor": bgc,
                "figure.edgecolor": fgc,
                "axes.edgecolor": fgc,
                "axes.facecolor": ebgc,
                "axes.labelcolor": fgc,
                "text.color": fgc,
                "xtick.color": fgc,
                "ytick.color": fgc,
            }
        )

        # some widgets also needs to be manually updated
        for w in self.forceUpdOnThemeWidget:
            w.config(background=ebgc, foreground=fgc)

        try:
            self.fig.set_facecolor(bgc)
            self.geomFig.set_facecolor(bgc)
            self.auxFig.set_facecolor(bgc)

            # this is necessary because twinned axis does not necessarily
            # follow the rcParams
            for ax in (
                self.ax,
                self.axv,
                self.axP,
                self.geomAx,
                self.axF,
                self.auxAx,
                self.auxAxH,
            ):
                ax.set_facecolor(ebgc)
                ax.spines["top"].set_color(fgc)
                ax.spines["bottom"].set_color(fgc)
                ax.spines["left"].set_color(fgc)
                ax.spines["right"].set_color(fgc)

            self.updateGeomPlot()
            self.updateFigPlot()
            self.updateAuxPlot()

        except AttributeError:
            pass
        try:
            self.tv.tag_configure("oddrow", background=bgc)
            self.tv.tag_configure("evenrow", background=ebgc)
        except AttributeError:
            pass

        self.update_idletasks()


def calculate(
    queue,
    kwargs,
):
    gunType = kwargs["typ"]
    constrain = kwargs["con"]
    optimize = kwargs["opt"]
    lock = kwargs["lock"]
    debug = kwargs["deb"]

    try:
        if constrain:
            if gunType == CONVENTIONAL:
                constrained = Constrained(**kwargs)
            elif gunType == RECOILESS:
                constrained = ConstrainedRecoiless(**kwargs)

            if optimize:
                l_f, e_1, l_g = constrained.findMinV(**kwargs)
                kwargs.update({"loadFraction": l_f})
            else:
                e_1, l_g = constrained.solve(**kwargs, known_bore=lock)

            kwargs.update({"grainSize": 2 * e_1})

            if not lock:
                kwargs.update({"lengthGun": l_g})

            chamberVolume = (
                kwargs["chargeMass"]
                / kwargs["propellant"].rho_p
                / kwargs["loadFraction"]
            )
            kwargs.update({"chamberVolume": chamberVolume})

        if gunType == CONVENTIONAL:
            gun = Gun(**kwargs)
        elif gunType == RECOILESS:
            gun = Recoiless(**kwargs)

        tableData, errorData, pressureTrace, structure = gun.integrate(**kwargs)
        errorReport = []

    except Exception as e:
        gun = None
        tableData = []
        errorData = []
        pressureTrace = []
        structure = []

        errorReport = ["Exception while calculating:"]
        exc_type, exc_value, exc_traceback = sys.exc_info()
        if debug:
            errorReport.append(
                "".join(
                    traceback.format_exception(
                        exc_type, exc_value, exc_traceback
                    )
                )
            )
        else:
            errorReport.append(str(e))

    queue.put(
        (
            kwargs,
            gun,
            tableData,
            errorData,
            pressureTrace,
            structure,
            errorReport,
        )
    )


def main():
    freeze_support()

    # if check avoids hackery when not profiling
    # Optional; hackery *seems* to work fine even when not profiling, it's just wasteful
    import cProfile

    if sys.modules["__main__"].__file__ == cProfile.__file__:
        import IB  # Imports you again (does *not* use cache or execute as __main__)

        globals().update(
            vars(IB)
        )  # Replaces current contents with newly imported stuff
        sys.modules[
            "__main__"
        ] = IB  # Ensures pickle lookups on __main__ find matching version

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
    myappid = "Phoenix.Interior.Ballistics.Solver"  # arbitrary string
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    root = Tk()

    root.iconbitmap(resolvepath("ui/logo.ico"))

    loadfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"), False, True)
    font_manager.fontManager.addfont(
        resolvepath("ui/sarasa-fixed-sc-regular.ttf")
    )
    # from tkinter import font
    # print(font.families())

    dpi = root.winfo_fpixels("1i")

    # Tk was originally developed for a dpi of 72
    # root.tk.call("tk", "scaling", "-displayof", ".", dpi / 72.0)
    scale = 1.0 * dpi / 72.0
    root.tk.call("tk", "scaling", scale)

    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.12"))

    root.option_add("*tearOff", False)
    root.title("PIBS v0.4.8")
    menubar = Menu(root)
    root.config(menu=menubar)
    """
    tabControl = ttk.Notebook(root)
    tabControl.pack(expand=1, fill="both", side="left")
    ibFrame = IB(tabControl, menubar, dpi, scale)
    tabControl.add(ibFrame, text="INTERIOR")
    """
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    InteriorBallisticsFrame(root, menubar, dpi)

    # center(root)
    root.update_idletasks()
    root.minsize(root.winfo_width(), root.winfo_height())  # set minimum size
    root.state("zoomed")  # maximize window
    root.mainloop()


if __name__ == "__main__":
    # print(__name__)
    main()
