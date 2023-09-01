from tkinter import (
    Frame,
    Menu,
    Text,
    filedialog,
    messagebox,
    StringVar,
    IntVar,
    Tk,
    ttk,
)

import tkinter.font as tkFont
import traceback

from gun import Gun, DOMAIN_TIME, DOMAIN_LENG
from gun import (
    POINT_PEAK,
    POINT_BURNOUT,
    POINT_FRACTURE,
    POINT_EXIT,
    POINT_PEAK_SHOT,
    POINT_PEAK_BREECH,
)
from gun import SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV
from recoiless import Recoiless

from prop import Propellant, GrainComp, GEOMETRIES, SimpleGeometry
from opt import Constrained
from optRecoiless import ConstrainedRecoiless
from tip import CreateToolTip
from lang import STRING

from misc import (
    toSI,
    validateNN,
    validatePI,
    formatFloatInput,
    formatIntInput,
    dot_aligned,
    resolvepath,
    roundSig,
    center,
    loadfont,
)
from math import ceil, pi, log10

import matplotlib.pyplot as mpl
from matplotlib.figure import Figure
from labellines import labelLines, labelLine
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from multiprocessing import Process, Queue
from queue import Empty

import sys
import csv
import datetime

from ctypes import windll
import platform
import multiprocessing

from matplotlib import font_manager


RECOILESS = "Recoiless Gun"
CONVENTIONAL = "Conventional Gun"


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
    # "path.simplify_threshold": 1,
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
    # "path.simplify_threshold": 1,
    "font.family": "Sarasa Fixed SC",
}


class IB(Frame):
    def __init__(self, parent, dpi, scale):
        Frame.__init__(self, parent)
        self.LANG = StringVar(value=list(STRING.keys())[0])

        default_font = tkFont.Font(family=FONTNAME, size=FONTSIZE)
        self.option_add("*Font", default_font)

        self.queue = Queue()
        self.process = None

        # self.tableData = None
        self.outflowData = None

        self.pos = -1
        self.dpi = dpi
        self.parent = parent
        self.forceUpdOnThemeWidget = []
        self.scale = scale

        menubar = Menu(parent)

        self.menubar = menubar

        fileMenu = Menu(menubar)
        menubar.add_cascade(
            label=self.getString("fileLabel"), menu=fileMenu, underline=0
        )
        themeMenu = Menu(menubar)
        menubar.add_cascade(
            label=self.getString("themeLabel"), menu=themeMenu, underline=0
        )
        debugMenu = Menu(menubar)
        menubar.add_cascade(
            label=self.getString("debugLabel"), menu=debugMenu, underline=0
        )
        langMenu = Menu(menubar)
        menubar.add_cascade(label="Lang 语言", menu=langMenu, underline=0)
        solMenu = Menu(menubar)
        menubar.add_cascade(
            label=self.getString("solLabel"), menu=solMenu, underline=0
        )

        self.fileMenu = fileMenu
        self.themeMenu = themeMenu
        self.debugMenu = debugMenu
        self.solMenu = solMenu

        self.themeRadio = IntVar(value=0)
        self.useTheme()

        self.DEBUG = IntVar(value=0)

        self.soln = IntVar(value=0)

        self.useCv = IntVar(value=0)

        fileMenu.add_command(
            label=self.getString("saveLabel"), command=self.save, underline=0
        )
        fileMenu.add_command(
            label=self.getString("loadLabel"), command=self.load, underline=0
        )
        fileMenu.add_command(
            label=self.getString("exportLabel"),
            command=self.export,
            underline=0,
        )

        themeMenu.add_radiobutton(
            label=self.getString("darkLabel"),
            variable=self.themeRadio,
            value=0,
            command=self.useTheme,
            underline=0,
        )
        themeMenu.add_radiobutton(
            label=self.getString("lightLabel"),
            variable=self.themeRadio,
            value=1,
            command=self.useTheme,
            underline=0,
        )

        debugMenu.add_checkbutton(
            label=self.getString("enableLabel"),
            variable=self.DEBUG,
            onvalue=1,
            offvalue=0,
            underline=0,
        )

        solMenu.add_radiobutton(
            label=self.getString("SOL_LAGRANGE"),
            variable=self.soln,
            value=0,
            underline=0,
        )
        solMenu.add_radiobutton(
            label=self.getString("SOL_PIDDUCK"),
            variable=self.soln,
            value=1,
            underline=0,
        )
        solMenu.add_radiobutton(
            label=self.getString("SOL_MAMONTOV"),
            variable=self.soln,
            value=2,
            underline=0,
        )
        solMenu.add_separator()
        """
        solMenu.add_checkbutton(
            label=self.getString("atmosLabel"),
            variable=self.inAtmos,
            onvalue=1,
            offvalue=0,
            underline=0,
        )
        solMenu.add_separator()
        """
        solMenu.add_checkbutton(
            label=self.getString("useLFLabel"), variable=self.useCv, onvalue=0
        )
        solMenu.add_checkbutton(
            label=self.getString("useCVLabel"), variable=self.useCv, onvalue=1
        )

        self.useCv.trace_add("write", self.cvlfCallback)

        for lang in STRING.keys():
            langMenu.add_radiobutton(
                label=lang,
                variable=self.LANG,
                value=lang,
                command=self.changeLang,
                underline=0,
            )

        parent.config(menu=menubar)

        self.compositions = GrainComp.readFile(
            resolvepath("data/propellants.csv")
        )  # dict of composition.name (string) -> composition (object)
        self.geometries = {self.getString(k): v for k, v in GEOMETRIES.items()}
        # dict of localized name (string) -> geometry (object)

        self.prop = None
        self.gun = None
        self.errorLst = []

        self.propOptions = tuple(self.compositions.keys())

        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(1, weight=1)
        self.addTopFrm(parent)
        self.addLeftFrm(parent)
        self.addRightFrm(parent)
        self.addErrFrm(parent)

        self.addPlotFrm(parent)
        self.addParFrm(parent)
        self.addTblFrm(parent)
        parent.update_idletasks()
        self.addGeomPlot()
        parent.update_idletasks()
        self.addFigPlot()

        self.ambCallback()
        self.cvlfCallback()
        self.typeCallback()
        self.ctrlCallback()

        self.updateSpec()
        self.updateGeom()

        self.forceUpdOnThemeWidget.append(self.errorText)
        self.forceUpdOnThemeWidget.append(self.specs)

        parent.bind("<Configure>", self.resizePlot)
        root.protocol("WM_DELETE_WINDOW", self.quit)
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
                self.getString("excTitle"), self.getString("noDataMsg")
            )
            return

        fileName = filedialog.asksaveasfilename(
            title=self.getString("saveLabel"),
            filetypes=(("Gun Design File", "*.gun"),),
            defaultextension=".gun",
            initialfile=self.getDescriptive(),
        )
        if fileName == "":
            messagebox.showinfo(
                self.getString("excTitle"), self.getString("cancelMsg")
            )
            return

        commentLines = []
        try:
            with open(fileName, "r", encoding="utf-8", newline="\n") as file:
                lines = file.readlines()

            isComment = False
            for line in lines:
                if line[:11] == "COMMENT END":
                    isComment = False
                if isComment:
                    commentLines.append(line)
                if line[:11] == "COMMENT BEG":
                    isComment = True

        except Exception:
            pass  # either file DNE or some exception occured during reading.

        try:
            kwargs = self.kwargs
            sigfig = int(-log10(kwargs["tol"]))

            with open(fileName, "w", encoding="utf-8", newline="\n") as file:
                file.write("{:>45}\n".format(self.getDescriptive()))
                file.write(
                    "LAST MODIFIED  {:>30}\n".format(
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
                    )
                )
                file.write("COMMENT BEG\n")
                for line in commentLines:
                    file.write(line)
                file.write("COMMENT END\n")
                file.write("DESIGN BEG\n")
                for desc, value, unit in zip(
                    (
                        "SIGNIFICANT FIGURE",
                        "TYPE",
                        "CALIBER",
                        "TUBE LENGTH",
                        "SHOT MASS",
                        "CHAMBER VOLUME",
                        "CHAMBRAGE",
                        "ENGRAVING PRESSURE",
                        "NOZZLE EFFICIENCY",
                        "NOZZLE EXPANSION R.",
                        "CHARGE MASS",
                        "PROPELLANT",
                        "GEOMETRY",
                        "2xWEB",
                        "1/ALPHA",
                        "1/BETA",
                        "DRAG COEFFICIENT",
                    ),
                    (
                        sigfig,
                        kwargs["typ"],
                        kwargs["caliber"] * 1e3,
                        kwargs["lengthGun"] * 1e3,
                        kwargs["shotMass"],
                        kwargs["chamberVolume"] * 1e3,
                        kwargs["chambrage"],
                        kwargs["startPressure"] * 1e-6,
                        kwargs["nozzleEfficiency"] * 1e2,
                        kwargs["nozzleExpansion"],
                        kwargs["chargeMass"],
                        kwargs["propellant"].composition.name,
                        kwargs["propellant"].geometry.name,
                        kwargs["grainSize"] * 1e3,
                        kwargs["propellant"].R1,
                        kwargs["propellant"].R2,
                        kwargs["dragCoefficient"] * 1e2,
                    ),
                    (
                        "",
                        "",
                        "mm",
                        "mm",
                        "kg",
                        "L",
                        "",
                        "MPa",
                        "%",
                        "",
                        "kg",
                        "",
                        "",
                        "mm",
                        "",
                        "",
                        "%",
                    ),
                ):
                    file.write(
                        "{:<20}{:>20} {:<4}\n".format(
                            desc,
                            roundSig(value, sigfig)
                            if isinstance(value, float)
                            else value,
                            unit,
                        )
                    )
                file.write("DESIGN END\n")

            messagebox.showinfo(
                self.getString("sucTitle"),
                self.getString("savedLocMsg").format(fileName),
            )

        except Exception as e:
            messagebox.showinfo(self.getString("excTitle"), str(e))

    def load(self):
        fileName = filedialog.askopenfilename(
            title=self.getString("loadLabel"),
            filetypes=(("Gun Design File", "*.gun"),),
            defaultextension=".gun",
            initialfile=self.getDescriptive(),
        )
        if fileName == "":
            messagebox.showinfo("Exception Loading Design", "No File Selected")
            return

        try:
            fileKwargs = {}
            with open(fileName, "r", encoding="utf-8", newline="\n") as file:
                lines = file.readlines()
            isDesign = False
            for line in lines:
                if line[:10] == "DESIGN END":
                    isDesign = False
                if isDesign:
                    desc = line[:20].strip()
                    try:
                        if desc == "SIGNIFICANT FIGURE":
                            val = int(line[20:40].strip())
                        else:
                            val = float(line[20:40].strip())
                    except ValueError:
                        val = line[20:40].strip()
                    fileKwargs.update({desc: val})
                if line[:10] == "DESIGN BEG":
                    isDesign = True

            self.typeOptn.set(
                self.typeOptn["values"][
                    (CONVENTIONAL, RECOILESS).index(fileKwargs["TYPE"])
                ]
            )
            self.calmm.set(fileKwargs["CALIBER"])
            self.tblmm.set(fileKwargs["TUBE LENGTH"])

            self.shtkg.set(fileKwargs["SHOT MASS"])
            self.chgkg.set(fileKwargs["CHARGE MASS"])
            self.useCv.set(1)
            self.cvL.set(fileKwargs["CHAMBER VOLUME"])
            self.clr.set(fileKwargs["CHAMBRAGE"])
            self.stpMPa.set(fileKwargs["ENGRAVING PRESSURE"])
            self.nozzEff.set(fileKwargs["NOZZLE EFFICIENCY"])
            self.nozzExp.set(fileKwargs["NOZZLE EXPANSION R."])

            self.dropGeom.set(
                self.dropGeom["values"][
                    list(GEOMETRIES.keys()).index(fileKwargs["GEOMETRY"])
                ]
            )  # since self.geomtries maps the localized name to object, while
            # we are saving the non-localized .name.
            self.dropProp.set(fileKwargs["PROPELLANT"])
            self.arcmm.set(fileKwargs["2xWEB"])
            self.grdR.set(fileKwargs["1/ALPHA"])
            self.grlR.set(fileKwargs["1/BETA"])

            self.dgc.set(fileKwargs["DRAG COEFFICIENT"])
            self.accExp.set(fileKwargs["SIGNIFICANT FIGURE"])

        except Exception as e:
            messagebox.showinfo(self.getString("excTitle"), str(e))

    def export(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(
                self.getString("excTitle"), self.getString("noDataMsg")
            )
            return

        fileName = filedialog.asksaveasfilename(
            title=self.getString("exportLabel"),
            filetypes=(("Comma Separated File", "*.csv"),),
            defaultextension=".csv",
            initialfile=self.getDescriptive(),
        )

        if fileName == "":
            messagebox.showinfo(
                self.getString("excTitle"), self.getString("cancelMsg")
            )
            return
        try:
            gunType = self.kwargs["typ"]
            with open(fileName, "w", encoding="utf-8", newline="") as csvFile:
                csvWriter = csv.writer(
                    csvFile, delimiter=",", quoting=csv.QUOTE_MINIMAL
                )

                if gunType == CONVENTIONAL:
                    headers = (
                        "",
                        "T - s",
                        "L - m",
                        "V - m/s",
                        "P_l=L - Pa",
                        "1/L ∫{0->L} P_l dl - Pa",
                        "P_l=0 - Pa",
                        "ψ - 1",
                        "T - K",
                    )
                elif gunType == RECOILESS:
                    headers = (
                        "",
                        "T - s",
                        "L - m",
                        "V - m/s",
                        "V_l=0 - m/s",
                        "P_l=L - Pa",
                        "1/L ∫{0->L} P_l dl - Pa",
                        "P_v=0 - Pa",
                        "P_l=0 - Pa",
                        "ψ - 1",
                        "η - 1",
                        "T - K",
                    )

                csvWriter.writerow(headers)

                for line in self.tableData:
                    tag, t, l, psi, v, P, T, eta = line

                    if gunType == CONVENTIONAL:
                        Pb, Pt = gun.toPsPb(l, P)
                        csvWriter.writerow((tag, t, l, v, Pb, P, Pt, psi, T))
                    elif gunType == RECOILESS:
                        Ps, P0, Px, vx = gun.toPsP0PxVx(l, v, P, T, eta)
                        csvWriter.writerow(
                            (tag, t, l, v, vx, Ps, P, P0, Px, psi, eta, T)
                        )

            messagebox.showinfo(
                self.getString("sucTitle"),
                self.getString("savedLocMsg").format(fileName),
            )

        except Exception as e:
            messagebox.showinfo(self.getString("excTitle"), str(e))

    def ambCallback(self, *args):
        self.ambPw.config(state="normal" if self.inAtmos.get() else "disabled")
        self.ambRhow.config(
            state="normal" if self.inAtmos.get() else "disabled"
        )
        self.ambGamw.config(
            state="normal" if self.inAtmos.get() else "disabled"
        )

    def cvlfCallback(self, *args):
        useCv = self.useCv.get()

        self.ldfw.config(state="disabled" if useCv else "normal")
        self.cvLw.config(state="normal" if useCv else "disabled")

    def changeLang(self):
        self.menubar.entryconfig(0, label=self.getString("fileLabel"))
        self.menubar.entryconfig(1, label=self.getString("themeLabel"))
        self.menubar.entryconfig(2, label=self.getString("debugLabel"))
        self.menubar.entryconfig(4, label=self.getString("solLabel"))

        self.fileMenu.entryconfig(0, label=self.getString("saveLabel"))
        self.fileMenu.entryconfig(1, label=self.getString("loadLabel"))
        self.fileMenu.entryconfig(2, label=self.getString("exportLabel"))

        self.themeMenu.entryconfig(0, label=self.getString("darkLabel"))
        self.themeMenu.entryconfig(1, label=self.getString("lightLabel"))
        self.debugMenu.entryconfig(0, label=self.getString("enableLabel"))

        self.solMenu.entryconfig(0, label=self.getString("SOL_LAGRANGE"))
        self.solMenu.entryconfig(1, label=self.getString("SOL_PIDDUCK"))
        self.solMenu.entryconfig(2, label=self.getString("SOL_MAMONTOV"))
        """
        # separator @ 3
        self.solMenu.entryconfig(4, label=self.getString("atmosLabel"))
        """
        # seaparator @ 3
        self.solMenu.entryconfig(4, label=self.getString("useCVLabel"))
        self.solMenu.entryconfig(5, label=self.getString("useLFLabel"))

        self.calLb.config(text=self.getString("calLabel"))
        self.tblLb.config(text=self.getString("tblLabel"))
        self.shtLb.config(text=self.getString("shtLabel"))
        self.chgLb.config(text=self.getString("chgLabel"))
        self.ldfLb.config(text=self.getString("ldfLabel"))
        self.clrLb.config(text=self.getString("clrLabel"))
        self.dgcLb.config(text=self.getString("dgcLabel"))
        self.stpLb.config(text=self.getString("stpLabel"))
        self.vTgtLb.config(text=self.getString("vTgtLabel"))
        self.pTgtLb.config(text=self.getString("pTgtLabel"))
        self.minWebLb.config(text=self.getString("minWebLabel"))
        self.lgmaxLb.config(text=self.getString("maxLgLabel"))
        self.ldLb.config(text=self.getString("ldLabel"))
        self.lxLb.config(text=self.getString("lxLabel"))
        self.vaLb.config(text=self.getString("vaLabel"))
        self.teLb.config(text=self.getString("teffLabel"))
        self.beLb.config(text=self.getString("beffLabel"))
        self.cvLLb.config(text=self.getString("cvLabel"))
        self.ambPLb.config(text=self.getString("ambPresLabel"))
        self.ambRhoLb.config(text=self.getString("ambRhoLabel"))
        self.ambGamLb.config(text=self.getString("ambGamLabel"))
        self.ammoLb.config(text=self.getString("ammoLabel"))
        self.outflowToLb.config(text=self.getString("outflowToLabel"))

        self.nozzExpLb.config(text=self.getString("nozzExpLabel"))
        self.nozzEffLb.config(text=self.getString("nozzEffLabel"))

        # self.psmaxLb.config(text=self.getString("psmaxLabel"))
        # self.pamaxLb.config(text=self.getString("pamaxLabel"))
        # self.pbmaxLb.config(text=self.getString("pbmaxLabel"))

        self.useConstraintTip.set(self.getString("useConsText"))
        self.optimizeLFTip.set(self.getString("optLFText"))
        self.chgTip.set(self.getString("chgText"))
        self.vinfTip.set(self.getString("vinfText"))
        self.lxTip.set(self.getString("calLxText"))
        self.geomPlotTip.set(self.getString("geomPlotText"))
        self.teffTip.set(self.getString("teffText"))
        self.beffTip.set(self.getString("beffText"))
        self.specsTip.set(self.getString("specsText"))
        self.ldfTip.set(self.getString("ldfText"))
        self.clrTip.set(self.getString("clrText"))
        self.dgcTip.set(self.getString("dgcText"))
        self.stpTip.set(self.getString("stpText"))
        self.nozzExpTip.set(self.getString("nozzExpText"))
        self.nozzEffTip.set(self.getString("nozzEffText"))
        self.tolTip.set(self.getString("tolText"))
        self.sampleTip.set(self.getString("sampText"))
        self.calcButtonTip.set(self.getString("calcButtonText"))
        self.pTgtTip.set(self.getString("pTgtText"))
        self.plotTip.set(self.getString("plotText"))

        self.tblFrm.config(text=self.getString("tblFrmLabel"))
        self.plotFrm.config(text=self.getString("plotFrmLabel"))
        self.errorFrm.config(text=self.getString("errFrmLabel"))
        self.parFrm.config(text=self.getString("parFrmLabel"))
        self.specFrm.config(text=self.getString("specFrmLabel"))
        self.opFrm.config(text=self.getString("opFrmLabel"))
        self.consFrm.config(text=self.getString("consFrmLabel"))
        self.pltOptnFrm.config(text=self.getString("pltOptnFrm"))
        self.sampleFrm.config(text=self.getString("sampleFrmLabel"))
        self.propFrm.config(text=self.getString("propFrmLabel"))
        self.grainFrm.config(text=self.getString("grainFrmLabel"))
        self.envFrm.config(text=self.getString("envFrmLabel"))
        self.outflowFrm.config(text=self.getString("outflowFrmLabel"))

        self.useConstraint.config(text=self.getString("consButton"))
        self.optimizeLF.config(text=self.getString("minTVButton"))

        for i, columnName in enumerate(self.getString("columnList")):
            self.tv.heading(i, text=columnName)

        self.plotAvgPCheck.config(text=self.getString("plotAvgP"))
        self.plotBasePCheck.config(text=self.getString("plotBaseP"))
        self.plotBreechNozzlePCheck.config(
            text=self.getString("plotBreechNozzleP")
        )
        self.plotStagPCheck.config(text=self.getString("plotStagP"))
        self.plotVelCheck.config(text=self.getString("plotVel"))
        self.plotNozzleVCheck.config(text=self.getString("plotNozzleV"))
        self.plotBurnupCheck.config(text=self.getString("plotBurnup"))
        self.plotEtaCheck.config(text=self.getString("plotEta"))
        self.plotOutflowCheck.config(text=self.getString("plotOutflow"))

        self.plotRecoilCheck.config(text=self.getString("plotRecoil"))

        self.inAtmosCheck.config(text=self.getString("atmosLabel"))

        self.stepLb.config(text=self.getString("stepLabel"))
        self.calButton.config(text=self.getString("calcLabel"))

        gunTypeIndex = self.typeOptn["values"].index(self.gunType.get())
        self.typeOptn.config(
            values=(
                self.getString("CONVENTIONAL"),
                self.getString("RECOILESS"),
            )
        )
        self.typeOptn.current(gunTypeIndex)

        domainIndex = self.dropOptn["values"].index(self.dropOptn.get())
        self.dropOptn.config(
            values=(
                self.getString("DOMAIN_TIME"),
                self.getString("DOMAIN_LENG"),
            )
        )
        self.dropOptn.current(domainIndex)

        geomIndex = self.dropGeom["values"].index(self.dropGeom.get())
        self.geometries = {self.getString(k): v for k, v in GEOMETRIES.items()}
        # mapping of the localized string to the object.
        self.dropGeom.config(values=tuple(self.geometries.keys()))
        self.dropGeom.current(geomIndex)

        self.updateGeom()
        self.updateSpec()
        self.updateFigPlot()

    def readTable(self, tag):
        i = [line[0] for line in self.tableData].index(tag)
        return self.tableData[i]

    def getString(self, name):
        try:
            return STRING[self.LANG.get()][name]
        except KeyError:
            return STRING["English"][name]

    def addTopFrm(self, parent):
        topFrm = ttk.Frame(parent)
        topFrm.grid(row=0, column=0, sticky="nsew")

        topFrm.columnconfigure(0, weight=1)
        topFrm.rowconfigure(0, weight=1)

        pltOptnFrm = ttk.LabelFrame(topFrm, text=self.getString("pltOptnFrm"))
        pltOptnFrm.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        self.pltOptnFrm = pltOptnFrm

        for i in range(10):
            pltOptnFrm.columnconfigure(i, weight=1)

        j = 0

        self.pbar = ttk.Progressbar(
            pltOptnFrm, mode="indeterminate", maximum=100
        )
        self.pbar.grid(
            row=j, column=0, columnspan=10, sticky="nsew", padx=2, pady=2
        )

        j += 1

        self.plotAvgP = IntVar(value=1)
        self.plotAvgPCheck = ttk.Checkbutton(
            pltOptnFrm, text=self.getString("plotAvgP"), variable=self.plotAvgP
        )
        self.plotAvgPCheck.grid(row=j, column=0, sticky="nsew")

        self.plotBaseP = IntVar(value=1)
        self.plotBasePCheck = ttk.Checkbutton(
            pltOptnFrm,
            text=self.getString("plotBaseP"),
            variable=self.plotBaseP,
        )
        self.plotBasePCheck.grid(row=j, column=1, sticky="nsew")

        self.plotBreechNozzleP = IntVar(value=1)
        self.plotBreechNozzlePCheck = ttk.Checkbutton(
            pltOptnFrm,
            text=self.getString("plotBreechNozzleP"),
            variable=self.plotBreechNozzleP,
        )
        self.plotBreechNozzlePCheck.grid(row=j, column=2, sticky="nsew")

        self.plotStagP = IntVar(value=1)
        self.plotStagPCheck = ttk.Checkbutton(
            pltOptnFrm,
            text=self.getString("plotStagP"),
            variable=self.plotStagP,
        )
        self.plotStagPCheck.grid(row=j, column=3, sticky="nsew")

        self.plotVel = IntVar(value=1)
        self.plotVelCheck = ttk.Checkbutton(
            pltOptnFrm, text=self.getString("plotVel"), variable=self.plotVel
        )
        self.plotVelCheck.grid(row=j, column=4, sticky="nsew")

        self.plotNozzleV = IntVar(value=1)
        self.plotNozzleVCheck = ttk.Checkbutton(
            pltOptnFrm,
            text=self.getString("plotNozzleV"),
            variable=self.plotNozzleV,
        )
        self.plotNozzleVCheck.grid(row=j, column=5, sticky="nsew")

        self.plotBurnup = IntVar(value=1)
        self.plotBurnupCheck = ttk.Checkbutton(
            pltOptnFrm,
            text=self.getString("plotBurnup"),
            variable=self.plotBurnup,
        )
        self.plotBurnupCheck.grid(row=j, column=6, sticky="nsew")

        self.plotEta = IntVar(value=1)
        self.plotEtaCheck = ttk.Checkbutton(
            pltOptnFrm, text=self.getString("plotEta"), variable=self.plotEta
        )
        self.plotEtaCheck.grid(row=j, column=7, sticky="nsew")

        self.plotOutflow = IntVar(value=0)
        self.plotOutflowCheck = ttk.Checkbutton(
            pltOptnFrm,
            text=self.getString("plotOutflow"),
            variable=self.plotOutflow,
        )
        self.plotOutflowCheck.grid(row=j, column=8, sticky="nsew")

        self.plotRecoil = IntVar(value=0)
        self.plotRecoilCheck = ttk.Checkbutton(
            pltOptnFrm,
            text=self.getString("plotRecoil"),
            variable=self.plotRecoil,
        )
        self.plotRecoilCheck.grid(row=j, column=9, sticky="nsew")

        self.plotAvgP.trace_add("write", self.updateFigPlot)
        self.plotBaseP.trace_add("write", self.updateFigPlot)
        self.plotBreechNozzleP.trace_add("write", self.updateFigPlot)
        self.plotStagP.trace_add("write", self.updateFigPlot)
        self.plotEta.trace_add("write", self.updateFigPlot)
        self.plotNozzleV.trace_add("write", self.updateFigPlot)
        self.plotBurnup.trace_add("write", self.updateFigPlot)
        self.plotVel.trace_add("write", self.updateFigPlot)
        self.plotOutflow.trace_add("write", self.updateFigPlot)
        self.plotRecoil.trace_add("write", self.updateFigPlot)

    def addLeftFrm(self, parent):
        leftFrm = ttk.Frame(parent)
        leftFrm.grid(row=0, column=1, rowspan=4, sticky="nsew")
        leftFrm.columnconfigure(0, weight=1)
        leftFrm.rowconfigure(0, weight=1)

    def addRightFrm(self, parent):
        rightFrm = ttk.Frame(parent)
        rightFrm.grid(row=0, column=3, rowspan=4, sticky="nsew")
        rightFrm.columnconfigure(0, weight=1)
        rightFrm.rowconfigure(0, weight=1)

        specFrm = ttk.LabelFrame(rightFrm, text=self.getString("specFrmLabel"))
        specFrm.grid(row=0, column=0, sticky="nsew")
        specFrm.columnconfigure(0, weight=1)
        self.specFrm = specFrm

        i = 0

        self.lxTip = StringVar(value=self.getString("calLxText"))
        self.lxLb, self.lx, self.tlx, _, _, i = self.add122Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("lxLabel"),
            unitText_up="Cal",
            unitText_dn="Cal",
            infotext=self.lxTip,
        )

        self.ammoLb, self.ammo, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("ammoLabel"),
            unitText="",
            default="0.00 x 0.000 mm",
        )

        self.vinfTip = StringVar(value=self.getString("vinfText"))
        self.vaLb, self.va, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("vaLabel"),
            infotext=self.vinfTip,
        )

        self.teffTip = StringVar(value=self.getString("teffText"))
        self.teLb, self.te, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("teffLabel"),
            unitText="%",
            infotext=self.teffTip,
        )

        self.beffTip = StringVar(value=self.getString("beffText"))
        self.beLb, self.be, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("beffLabel"),
            unitText="%",
            infotext=self.beffTip,
        )

        self.ldLb, self.ld, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("ldLabel"),
            unitText="kg/m³",
            # justify="right",
        )
        """
        self.pamaxLb, self.pamax, _, i = self.add12Disp(
            parent=specFrm, rowIndex=i, labelText=self.getString("pamaxLabel")
        )
        self.pbmaxLb, self.pbmax, _, i = self.add12Disp(
            parent=specFrm, rowIndex=i, labelText=self.getString("pbmaxLabel")
        )
        self.psmaxLb, self.psmax, _, i = self.add12Disp(
            parent=specFrm, rowIndex=i, labelText=self.getString("psmaxLabel")
        )
        """

        validationNN = parent.register(validateNN)
        validationPI = parent.register(validatePI)

        outflowFrm = ttk.LabelFrame(
            rightFrm, text=self.getString("outflowFrmLabel")
        )
        outflowFrm.grid(row=1, column=0, sticky="nsew")
        outflowFrm.columnconfigure(0, weight=1)
        self.outflowFrm = outflowFrm

        i = 0
        self.outflowToLb, self.outflowTo, _, _, i = self.add3Input(
            parent=outflowFrm,
            rowIndex=i,
            labelText=self.getString("outflowToLabel"),
            validation=validationNN,
            unitText="ms  ",
            default="1.0",
        )

        envFrm = ttk.LabelFrame(rightFrm, text=self.getString("envFrmLabel"))
        envFrm.grid(row=2, column=0, sticky="nsew")
        envFrm.columnconfigure(0, weight=1)
        self.envFrm = envFrm

        i = 0
        self.inAtmos = IntVar(value=1)
        self.inAtmos.trace_add("write", self.ambCallback)
        self.inAtmosCheck = ttk.Checkbutton(
            envFrm, text=self.getString("atmosLabel"), variable=self.inAtmos
        )
        self.inAtmosCheck.grid(row=i, column=0, columnspan=3, sticky="nsew")
        i += 1

        self.ambPLb, self.ambP, self.ambPw, _, i = self.add3Input(
            parent=envFrm,
            rowIndex=i,
            colIndex=0,
            labelText=self.getString("ambPresLabel"),
            unitText="kPa",
            default="101.325",
            validation=validationNN,
        )

        self.ambRhoLb, self.ambRho, self.ambRhow, _, i = self.add3Input(
            parent=envFrm,
            rowIndex=i,
            colIndex=0,
            labelText=self.getString("ambRhoLabel"),
            unitText="kg/m³",
            default="1.204",
            validation=validationNN,
        )

        self.ambGamLb, self.ambGam, self.ambGamw, _, i = self.add3Input(
            parent=envFrm,
            rowIndex=i,
            colIndex=0,
            labelText=self.getString("ambGamLabel"),
            default="1.400",
            validation=validationNN,
        )

        opFrm = ttk.LabelFrame(rightFrm, text=self.getString("opFrmLabel"))
        opFrm.grid(row=3, column=0, sticky="nsew")
        opFrm.columnconfigure(1, weight=1)
        self.opFrm = opFrm

        i = 0

        consFrm = ttk.LabelFrame(
            opFrm,
            text=self.getString("consFrmLabel"),
            style="SubLabelFrame.TLabelframe",
        )
        consFrm.grid(
            row=i, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )
        self.consFrm = consFrm
        j = 0

        self.vTgtLb, self.vTgt, _, _, j = self.add3Input(
            parent=consFrm,
            rowIndex=j,
            colIndex=0,
            labelText=self.getString("vTgtLabel"),
            unitText="m/s",
            default="1500.0",
            validation=validationNN,
        )

        self.pTgtTip = StringVar(value=self.getString("pTgtText"))
        self.pTgtLb, self.pTgt, _, _, j = self.add3Input(
            parent=consFrm,
            rowIndex=j,
            colIndex=0,
            labelText=self.getString("pTgtLabel"),
            unitText="MPa",
            default="350.0",
            validation=validationNN,
            infotext=self.pTgtTip,
        )

        self.minWebLb, self.minWeb, self.minWebw, _, j = self.add3Input(
            parent=consFrm,
            rowIndex=j,
            colIndex=0,
            labelText=self.getString("minWebLabel"),
            unitText="μm",
            default="1.0",
            validation=validationNN,
            color="red",
        )
        self.lgmaxLb, self.lgmax, self.lgmaxw, _, j = self.add3Input(
            parent=consFrm,
            rowIndex=j,
            colIndex=0,
            labelText=self.getString("maxLgLabel"),
            unitText="m",
            default="1000.0",
            validation=validationNN,
            color="red",
        )

        j += 1
        self.solve_W_Lg = IntVar()
        self.solve_W_Lg.set(0)
        self.useConstraint = ttk.Checkbutton(
            consFrm, text=self.getString("consButton"), variable=self.solve_W_Lg
        )
        self.useConstraint.grid(row=j, column=0, columnspan=3, sticky="nsew")
        self.solve_W_Lg.trace_add("write", self.ctrlCallback)

        self.useConstraintTip = StringVar(value=self.getString("useConsText"))
        CreateToolTip(self.useConstraint, self.useConstraintTip)

        j += 1
        self.opt_lf = IntVar()

        self.optimizeLF = ttk.Checkbutton(
            consFrm, text=self.getString("minTVButton"), variable=self.opt_lf
        )
        self.optimizeLF.grid(row=j, column=0, columnspan=3, sticky="nsew")
        self.optimizeLFTip = StringVar(value=self.getString("optLFText"))
        CreateToolTip(self.optimizeLF, self.optimizeLFTip)
        i += 1

        sampleFrm = ttk.LabelFrame(
            opFrm,
            text=self.getString("sampleFrmLabel"),
            style="SubLabelFrame.TLabelframe",
        )
        sampleFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        sampleFrm.columnconfigure(0, weight=1)
        sampleFrm.columnconfigure(1, weight=1)
        self.sampleFrm = sampleFrm

        self.dropOptn = ttk.Combobox(
            sampleFrm,
            values=(
                self.getString("DOMAIN_TIME"),
                self.getString("DOMAIN_LENG"),
            ),
            state="readonly",
            justify="center",
        )

        j = 0
        self.dropOptn.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropOptn.current(0)
        self.dropOptn.grid(
            row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2
        )
        # self.dropOptn.configure(width=0)

        j += 1

        self.stepLb, self.step, _, j = self.add2Input(
            parent=sampleFrm,
            rowIndex=j,
            colIndex=0,
            labelText=self.getString("stepLabel"),
            default="100",
            validation=validationNN,
            formatter=formatIntInput,
            reverse=True,
            anchor="center",
        )

        self.sampleTip = StringVar(value=self.getString("sampText"))
        CreateToolTip(sampleFrm, self.sampleTip)

        i += 1

        self.tolTip = StringVar(value=self.getString("tolText"))
        self.accExpLb, self.accExp, _, i = self.add2Input(
            parent=opFrm,
            rowIndex=i,
            colIndex=0,
            labelText="-log10(ε)",
            default="4",
            validation=validationPI,
            formatter=formatIntInput,
            color="red",
            infotext=self.tolTip,
        )

        self.calButton = ttk.Button(
            opFrm,
            text=self.getString("calcLabel"),
            # command=self.calculate,  # underline=0
            command=self.onCalculate,
        )
        self.calButton.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        opFrm.rowconfigure(i, weight=1)
        self.calcButtonTip = StringVar(value=self.getString("calcButtonText"))
        CreateToolTip(self.calButton, self.calcButtonTip)

    def onCalculate(self):
        constrain = self.solve_W_Lg.get() == 1
        optimize = self.opt_lf.get() == 1
        debug = self.DEBUG.get() == 1
        atmosphere = self.inAtmos.get() == 1

        invGunTypeLookup = {
            self.getString("CONVENTIONAL"): CONVENTIONAL,
            self.getString("RECOILESS"): RECOILESS,
        }
        gunType = invGunTypeLookup[self.typeOptn.get()]

        invDomainLookup = {
            self.getString("DOMAIN_TIME"): DOMAIN_TIME,
            self.getString("DOMAIN_LENG"): DOMAIN_LENG,
        }

        sols = (SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV)

        self.kwargs = {
            "opt": optimize,
            "con": constrain,
            "deb": debug,
            "typ": gunType,
            "dom": invDomainLookup[self.dropOptn.get()],
            "sol": sols[self.soln.get()],
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

        self.kwargs.update({"outflowTo": float(self.outflowTo.get()) * 1e-3})

        self.process = None

        try:
            if self.prop is None:
                raise ValueError("Invalid propellant.")

            if self.useCv.get():
                chamberVolume = float(self.cvL.get()) * 1e-3
            else:
                chamberVolume = (
                    float(self.chgkg.get())
                    / self.prop.rho_p
                    / float(self.ldf.get())
                    * 100
                )

            self.kwargs.update(
                {
                    "caliber": float(self.calmm.get()) * 1e-3,
                    "shotMass": float(self.shtkg.get()),
                    "propellant": self.prop,
                    "grainSize": float(self.arcmm.get()) * 1e-3,
                    "chargeMass": float(self.chgkg.get()),
                    "chargeMassRatio": (
                        float(self.chgkg.get()) / float(self.shtkg.get())
                    ),
                    "chamberVolume": chamberVolume,
                    "startPressure": float(self.stpMPa.get()) * 1e6,
                    "lengthGun": float(self.tblmm.get()) * 1e-3,
                    "chambrage": float(self.clr.get()),  # chamber expansion
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
                    "loadFraction": 1e-2 * float(self.ldf.get()),
                    "step": int(self.step.get()),
                }
            )

            self.process = Process(
                target=calculate,
                args=(
                    self.queue,
                    self.kwargs,
                ),
            )
            self.pos = 0
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

            self.tv.delete(*self.tv.get_children())
            self.updateError()
            self.updateFigPlot()

    def getValue(self):
        queue = self.queue

        try:
            if self.pos == 0:
                self.kwargs = queue.get_nowait()
                self.pos += 1

            if self.pos == 1:
                self.gun = queue.get_nowait()
                self.pos += 1

            if self.pos == 2:
                self.intgRecord = queue.get_nowait()
                self.pos += 1

            if self.pos == 3:
                self.tableData = queue.get_nowait()
                self.pos += 1

            if self.pos == 4:
                self.errorData = queue.get_nowait()
                self.pos += 1

            if self.pos == 5:
                self.errorLst.extend(queue.get_nowait())
                self.pos += 1

        except Empty:
            return

        constrain = self.kwargs["con"]
        optimize = self.kwargs["opt"]
        gunType = self.kwargs["typ"]

        self.pos = -1
        kwargs = self.kwargs
        sigfig = int(-log10(kwargs["tol"]))

        if self.gun is not None:
            if constrain:
                webmm = roundSig(1e3 * kwargs["grainSize"], n=sigfig)
                self.arcmm.set(webmm)

                lgmm = roundSig(kwargs["lengthGun"] * 1e3, n=sigfig)
                self.tblmm.set(lgmm)

                if optimize:
                    if self.useCv.get():
                        self.cvL.set(
                            roundSig(kwargs["chamberVolume"] * 1e3, n=sigfig)
                        )
                    else:
                        self.ldf.set(
                            roundSig(kwargs["loadFraction"] * 100, n=sigfig)
                        )

            self.ld.set(
                toSI(kwargs["loadFraction"] * self.prop.rho_p, useSN=True)
            )

            _, te, lg, _, vg, pe, Te, _ = self.readTable(POINT_EXIT)

            if gunType == CONVENTIONAL:
                _, pb = self.gun.toPsPb(lg, pe)
                self.outflowData = self.gun.hugoniot(
                    Te,
                    n=kwargs["step"],
                    t_offset=te,  # tol=kwargs["tol"],
                    t_to=kwargs["outflowTo"],
                )
            elif gunType == RECOILESS:
                self.outflowData = None

            eta_t, eta_b = self.gun.getEff(vg)

            self.te.set(round(eta_t * 100, 1))
            self.be.set(round(eta_b * 100, 1))

            """
            _, tp, lp, _, vp, pp, Tp, etap = self.readTable(POINT_PEAK)

            self.pamax.set(toSI(pp, unit="Pa"))

            _, tp, lp, _, vp, pp, Tp, etap = self.readTable(POINT_PEAK_BREECH)

            if gunType == CONVENTIONAL:
                _, pb = self.gun.toPsPb(lp, pp)

            elif gunType == RECOILESS:
                _, _, pb, _ = self.gun.toPsP0PxVx(lp, vp, pp, Tp, etap)

            self.pbmax.set(toSI(pb, unit="Pa"))

            _, tp, lp, _, vp, pp, Tp, etap = self.readTable(POINT_PEAK_SHOT)

            if gunType == CONVENTIONAL:
                ps, _ = self.gun.toPsPb(lp, pp)

            elif gunType == RECOILESS:
                ps, _, _, _ = self.gun.toPsP0PxVx(lp, vp, pp, Tp, etap)

            self.psmax.set(toSI(ps, unit="Pa"))
            """

            self.lx.set(toSI(kwargs["lengthGun"] / kwargs["caliber"]))
            self.tlx.set(
                toSI(
                    (kwargs["lengthGun"] + self.gun.l_0 / kwargs["chambrage"])
                    / kwargs["caliber"]
                )
            )
            self.va.set(toSI(self.gun.v_j, unit="m/s"))

            caliber = kwargs["caliber"]
            cartridge_len = kwargs["chamberVolume"] / (
                pi * (0.5 * caliber) ** 2 * kwargs["chambrage"]
            )
            self.ammo.set(
                "{:.3g} x {:.4g} mm".format(caliber * 1e3, cartridge_len * 1e3)
            )

        # populate the table with new values

        self.tv.delete(*self.tv.get_children())
        useSN = (False, False, False, True, False, False, True, True)
        units = (None, "s", "m", None, "m/s", "Pa", "K", None)
        tableData = dot_aligned(
            self.tableData,
            units=units,
            useSN=useSN,
        )

        errorData = dot_aligned(self.errorData, units=units, useSN=useSN)

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
        self.updateError()
        self.updateFigPlot()

        self.pbar.stop()
        self.calButton.config(state="normal")
        self.process = None

    def addErrFrm(self, parent):
        errorFrm = ttk.LabelFrame(parent, text=self.getString("errFrmLabel"))
        errorFrm.grid(row=3, column=0, sticky="nsew")
        errorFrm.columnconfigure(0, weight=1)
        errorFrm.rowconfigure(0, weight=1)
        self.errorFrm = errorFrm

        errScroll = ttk.Scrollbar(errorFrm, orient="vertical")
        errScroll.grid(row=0, column=1, sticky="nsew")
        self.errorText = Text(
            errorFrm,
            yscrollcommand=errScroll.set,
            wrap="word",
            height=5,
            width=0,
            font=(FONTNAME, FONTSIZE),
        )

        self.errorText.grid(row=0, column=0, sticky="nsew")

    def addParFrm(self, parent):
        parFrm = ttk.LabelFrame(parent, text=self.getString("parFrmLabel"))
        parFrm.grid(row=0, column=2, rowspan=4, sticky="nsew")
        parFrm.columnconfigure(0, weight=1)
        self.parFrm = parFrm
        # validation
        validationNN = parent.register(validateNN)

        i = 0
        self.gunType = StringVar()
        self.typeOptn = ttk.Combobox(
            parFrm,
            textvariable=self.gunType,
            values=(
                self.getString("CONVENTIONAL"),
                self.getString("RECOILESS"),
            ),
            state="readonly",
            justify="center",
        )
        self.typeOptn.option_add("*TCombobox*Listbox.Justify", "center")
        self.typeOptn.current(0)
        self.typeOptn.grid(
            row=i, column=0, sticky="nsew", padx=2, pady=2, columnspan=3
        )

        i += 1
        self.calLb, self.calmm, _, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("calLabel"),
            unitText="mm",
            default="50.0",
            validation=validationNN,
        )
        self.tblLb, self.tblmm, _, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("tblLabel"),
            unitText="mm",
            default="3500.0",
            validation=validationNN,
        )
        self.shtLb, self.shtkg, _, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("shtLabel"),
            unitText="kg",
            default="1.0",
            validation=validationNN,
        )

        self.chgTip = StringVar(value=self.getString("chgText"))
        self.chgLb, self.chgkg, _, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("chgLabel"),
            unitText="kg",
            default="1.0",
            validation=validationNN,
            infotext=self.chgTip,
        )

        # allow propellant specification to grow
        parFrm.rowconfigure(i, weight=4)

        propFrm = ttk.LabelFrame(
            parFrm,
            text=self.getString("propFrmLabel"),
            style="SubLabelFrame.TLabelframe",
        )
        propFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        propFrm.rowconfigure(1, weight=1)
        propFrm.columnconfigure(0, weight=1)
        # propFrm.columnconfigure(1, weight=1)
        self.propFrm = propFrm

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

        self.specsTip = StringVar(value=self.getString("specsText"))
        CreateToolTip(propFrm, self.specsTip)

        i += 1

        # allow grain frame to grow
        parFrm.rowconfigure(i, weight=1)

        grainFrm = ttk.LabelFrame(
            parFrm,
            text=self.getString("grainFrmLabel"),
            style="SubLabelFrame.TLabelframe",
        )
        grainFrm.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        grainFrm.columnconfigure(0, weight=1)
        grainFrm.rowconfigure(0, weight=1)
        self.grainFrm = grainFrm

        geomPlotFrm = ttk.LabelFrame(
            grainFrm,
            text="σ(Z)",
            style="SubLabelFrame.TLabelframe",
        )

        geomPlotFrm.grid(
            row=0,
            column=0,
            columnspan=3,
            sticky="nsew",
            padx=2,
            pady=2,
        )

        self.geomPlotTip = StringVar(value=self.getString("geomPlotText"))
        CreateToolTip(geomPlotFrm, self.geomPlotTip)

        self.geomParentFrm = grainFrm
        self.geomPlotFrm = geomPlotFrm

        j = 1

        # Create Dropdown menu
        self.currGeom = StringVar()

        self.dropGeom = ttk.Combobox(
            grainFrm,
            textvariable=self.currGeom,
            values=tuple(self.geometries.keys()),
            state="readonly",
            justify="center",
        )

        self.dropGeom.option_add("*TCombobox*Listbox.Justify", "center")
        self.dropGeom.current(0)

        self.dropGeom.grid(
            row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        self.lengthPrimaryAs = StringVar()
        self.lengthPrimaryTip = StringVar()

        j += 1
        self.arcLb, self.arcmm, _, _, j = self.add3Input(
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
        self.grdRLb, self.grdR, self.grdRw, _, j = self.add3Input(
            parent=grainFrm,
            rowIndex=j,
            labelText=self.ratioAs,
            unitText="x",
            default="1.0",
            validation=validationNN,
            infotext=self.lengthSecondaryTip,
        )

        self.lengthRatioAs = StringVar()
        self.lengthRatioTip = StringVar()

        self.grlRLb, self.grlR, self.grlRw, _, j = self.add3Input(
            parent=grainFrm,
            rowIndex=j,
            labelText=self.lengthRatioAs,
            unitText="x",
            default="2.5",
            validation=validationNN,
            infotext=self.lengthRatioTip,
        )

        i += 1
        self.ldfTip = StringVar(value=self.getString("ldfText"))
        self.ldfLb, self.ldf, self.ldfw, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("ldfLabel"),
            unitText="%",
            default="50.0",
            validation=validationNN,
            infotext=self.ldfTip,
        )

        self.cvLLb, self.cvL, self.cvLw, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("cvLabel"),
            unitText="L",
            default="20",
            validation=validationNN,
        )

        self.cvL.trace_add("write", self.cvlfConsisCallback)
        self.ldf.trace_add("write", self.cvlfConsisCallback)
        self.chgkg.trace_add("write", self.cvlfConsisCallback)
        self.accExp.trace_add("write", self.cvlfConsisCallback)

        self.clrTip = StringVar(value=self.getString("clrText"))
        self.clrLb, self.clr, _, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("clrLabel"),
            unitText="x",
            default="1.5",
            validation=validationNN,
            infotext=self.clrTip,
        )

        self.dgcTip = StringVar(value=self.getString("dgcText"))
        self.dgcLb, self.dgc, _, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("dgcLabel"),
            unitText="%",
            default="5.0",
            validation=validationNN,
            infotext=self.dgcTip,
        )

        self.stpTip = StringVar(value=self.getString("stpText"))
        self.stpLb, self.stpMPa, _, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("stpLabel"),
            unitText="MPa",
            default="10",
            validation=validationNN,
            infotext=self.stpTip,
        )

        self.nozzExpTip = StringVar(value=self.getString("nozzExpText"))
        (
            self.nozzExpLb,
            self.nozzExp,
            self.nozzExpw,
            self.nozzExpU,
            i,
        ) = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("nozzExpLabel"),
            unitText="x",
            default="4",
            validation=validationNN,
            infotext=self.nozzExpTip,
        )

        self.nozzEffTip = StringVar(value=self.getString("nozzEffText"))
        (
            self.nozzEffLb,
            self.nozzEff,
            self.nozzEffw,
            self.nozzEffU,
            i,
        ) = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("nozzEffLabel"),
            unitText="%",
            default="92.0",
            validation=validationNN,
            infotext=self.nozzEffTip,
        )

        self.currProp.trace_add("write", self.updateSpec)
        self.currGeom.trace_add("write", self.updateGeom)

        self.grdR.trace_add("write", self.callback)
        self.grlR.trace_add("write", self.callback)
        self.arcmm.trace_add("write", self.callback)

        self.gunType.trace_add("write", self.typeCallback)

    def addGeomPlot(self):
        geomPlotFrm = self.geomPlotFrm

        # _, _, width, height = self.geomPlotFrm.bbox("insert")

        geomPlotFrm.config(height=0.5 * geomPlotFrm.winfo_width())

        width = geomPlotFrm.winfo_width() - 2
        height = geomPlotFrm.winfo_height() - 2

        geomPlotFrm.columnconfigure(0, weight=1)
        geomPlotFrm.rowconfigure(0, weight=1)

        """
        geomPlotFrm.config(width=width, height=width)
        # we lock the frame the plot is put in
        geomPlotFrm.grid_propagate(False)
        # and lock it there
        If we set the global dpi awareness for this window
        to false (0), this will result in a correctly sized
        window, and everything would be fine, ironically.
        (in this case, dpi = dpi)
        """

        dpi = self.dpi
        with mpl.rc_context(GEOM_CONTEXT):
            fig = Figure(
                figsize=(width / dpi, max(height / dpi, 0.5 * width / dpi)),
                dpi=96,
                layout="constrained",
            )
            self.geomFig = fig
            self.geomAx = fig.add_subplot(111)

            self.geomCanvas = FigureCanvasTkAgg(fig, master=geomPlotFrm)
            self.geomCanvas.get_tk_widget().grid(
                row=0, column=0, padx=0, pady=0, sticky="nsew"
            )
            fig.set_layout_engine(None)
            self.geomCanvas.draw_idle()

        geomPlotFrm.grid_propagate(False)
        # last ditch effort to prevent blowing up the frame

    def addPlotFrm(self, parent):
        plotFrm = ttk.LabelFrame(
            parent, text=self.getString("plotFrmLabel"), width=640, height=480
        )
        plotFrm.grid(row=1, column=0, sticky="nsew")
        plotFrm.columnconfigure(0, weight=1)
        plotFrm.rowconfigure(0, weight=1)

        self.plotFrm = plotFrm
        self.plotTip = StringVar(value=self.getString("plotText"))
        CreateToolTip(self.plotFrm, self.plotTip)

    def addFigPlot(self):
        plotFrm = self.plotFrm

        # this will force a resize event to ensure the inserted
        # graph is of the correct size.

        # this is necessary here because winfo_width() will return
        # valid values with simply update_idletask()
        width = plotFrm.winfo_width() - 6
        # in pixels, -2 to account for label frame border thickness
        # additional -4 for padding
        height = plotFrm.winfo_height() - 6
        # technically we also need to account for the height of the
        # label frame text, but since its screen dependent its not
        # really possible.

        dpi = self.dpi
        with mpl.rc_context(FIG_CONTEXT):
            fig = Figure(
                figsize=(width / dpi, height / dpi),
                dpi=96,
                layout="constrained",
            )
            # fig.subplots_adjust(bottom=0.1)

            axes = fig.add_subplot(111)

            ax = axes
            axP = ax.twinx()
            axv = ax.twinx()
            axF = ax.twinx()

            ax.yaxis.tick_right()
            axF.yaxis.tick_left()
            ax.set_xlabel(" ")

            axv.spines.right.set_position(("axes", 1.0 + 45 * dpi / 96 / width))
            axP.spines.right.set_position(("data", 0.5))

            axP.yaxis.set_ticks(axP.get_yticks()[1:-1:])

            self.ax = ax
            self.axP = axP
            self.axv = axv
            self.axF = axF
            self.fig = fig

            self.pltCanvas = FigureCanvasTkAgg(fig, master=plotFrm)
            self.pltCanvas.get_tk_widget().grid(
                row=0, column=0, padx=2, pady=2, sticky="nsew"
            )

            self.pltCanvas.draw_idle()

        """
        For some reason the above specification is not strictly adhered to
        in practice, this might be a bug on tkAgg or matplotlib backend.
        """
        plotFrm.grid_propagate(False)

    def resizePlot(self, event):
        # we use the bbox method here as it has already accounted for padding
        # so no adjustment here is necessary
        _, _, width, height = self.plotFrm.bbox("insert")

        dpi = self.dpi
        with mpl.rc_context(FIG_CONTEXT):
            self.axv.spines.right.set_position(
                ("axes", 1 + 45 * dpi / 96 / width)
            )

    def timedLoop(self):
        if self.pos >= 0:  # and not self.process.is_alive():
            self.getValue()

        self.tLid = self.parent.after(100, self.timedLoop)

    def quit(self):
        root = self.parent
        if self.tLid is not None:
            root.after_cancel(self.tLid)
        root.quit()
        # root.destroy()

    def updateFigPlot(self, *args):
        with mpl.rc_context(FIG_CONTEXT):
            self.ax.cla()
            self.axP.cla()
            self.axv.cla()
            self.axF.cla()

            dpi = self.dpi
            size = self.fig.get_size_inches() * self.fig.dpi

            self.axv.spines.right.set_position(
                ("axes", 1 + 45 * dpi / 96 / size[0])
            )
            gun = self.gun
            if gun is not None:
                vTgt = self.kwargs["designVelocity"]
                gunType = self.kwargs["typ"]
                dom = self.kwargs["dom"]
                tol = self.kwargs["tol"]

                xs, vs = [], []
                Pas, Pss, Pbs, P0s = [], [], [], []
                Frs = []  # recoil forces
                psis, etas = [], []
                vxs = []

                for i, (t, (l, psi, v, p, T, eta)) in enumerate(
                    self.intgRecord
                ):
                    if dom == DOMAIN_TIME:
                        xs.append(t * 1000)
                    elif dom == DOMAIN_LENG:
                        xs.append(l)

                    vs.append(v)
                    Pas.append(p / 1e6)

                    if gunType == CONVENTIONAL:
                        Ps, Pb = gun.toPsPb(l, p)
                        P0 = 0
                        vx = 0
                        Fr = p * gun.S

                    elif gunType == RECOILESS:
                        Ps, P0, Pb, vx = gun.toPsP0PxVx(l, v, p, T, eta)
                        Fr = p * gun.S - gun.C_f * gun.S_j * p

                    Pss.append(Ps / 1e6)
                    Pbs.append(Pb / 1e6)
                    P0s.append(P0 / 1e6)
                    vxs.append(vx)

                    Frs.append(Fr / 1e6)

                    psis.append(psi)
                    etas.append(eta)

                if self.plotVel.get():
                    self.axv.scatter(
                        xs, vs, color="tab:blue", marker="s", s=FONTSIZE
                    )
                if self.plotAvgP.get():
                    self.axP.scatter(
                        xs, Pas, color="tab:green", marker="s", s=FONTSIZE
                    )
                if self.plotBurnup.get():
                    self.ax.scatter(
                        xs, psis, color="tab:red", marker="s", s=FONTSIZE
                    )

                xPeak = 0
                for i, (tag, t, l, psi, v, p, T, eta) in enumerate(
                    self.tableData
                ):
                    if dom == DOMAIN_TIME:
                        x = t * 1000
                    elif dom == DOMAIN_LENG:
                        x = l
                    xs.append(x)
                    if tag == POINT_PEAK:
                        xPeak = x
                    vs.append(v)
                    Pas.append(p / 1e6)

                    if gunType == CONVENTIONAL:
                        Ps, Pb = gun.toPsPb(l, p)
                        P0 = 0
                        vx = 0
                        Fr = p * gun.S

                    elif gunType == RECOILESS:
                        Ps, P0, Pb, vx = gun.toPsP0PxVx(l, v, p, T, eta)
                        Fr = p * gun.S - gun.C_f * gun.S_j * p

                    Pss.append(Ps / 1e6)
                    Pbs.append(Pb / 1e6)
                    P0s.append(P0 / 1e6)
                    vxs.append(vx)
                    psis.append(psi)
                    etas.append(eta)

                    Frs.append(Fr / 1e6)

                self.axP.spines.right.set_position(("data", xPeak))

                (xs, vs, vxs, Pas, Pss, Pbs, P0s, psis, etas, Frs) = zip(
                    *sorted(
                        zip(xs, vs, vxs, Pas, Pss, Pbs, P0s, psis, etas, Frs),
                        key=lambda line: line[0],
                    )
                )

                if self.plotVel.get():
                    self.axv.plot(
                        xs,
                        vs,
                        "tab:blue",
                        label=self.getString("figShotVel"),
                        marker=".",
                        alpha=1,
                    )
                self.axv.scatter(
                    xs[-1],
                    vTgt,
                    s=FONTSIZE**2,
                    c="tab:blue",
                    marker=5,
                    alpha=1,
                )

                if self.plotBreechNozzleP.get():
                    self.axP.plot(
                        xs,
                        Pbs,
                        c="xkcd:goldenrod",
                        label=self.getString("figBreech")
                        if gunType == CONVENTIONAL
                        else self.getString("figNozzleP")
                        if gunType == RECOILESS
                        else "",
                        linestyle="dashed",
                        alpha=0.75,
                    )

                if gunType == RECOILESS:
                    if self.plotStagP.get():
                        self.axP.plot(
                            xs,
                            P0s,
                            "seagreen",
                            label=self.getString("figStagnation"),
                            linestyle="dashed",
                            alpha=0.75,
                        )

                    if self.plotNozzleV.get():
                        self.axv.plot(
                            xs,
                            vxs,
                            "royalblue",
                            label=self.getString("figNozzleV"),
                            alpha=0.75,
                            linestyle="dashed",
                        )

                    if self.plotEta.get():
                        self.ax.plot(
                            xs,
                            etas,
                            "crimson",
                            label=self.getString("figOutflow"),
                            alpha=0.75,
                            linestyle="dashed",
                        )

                if self.plotAvgP.get():
                    self.axP.plot(
                        xs,
                        Pas,
                        "tab:green",
                        label=self.getString("figAvgP"),
                        marker=".",
                        alpha=1,
                    )

                if self.plotBaseP.get():
                    self.axP.plot(
                        xs,
                        Pss,
                        "yellowgreen",
                        label=self.getString("figShotBase"),
                        linestyle="dashed",
                        alpha=0.75,
                    )

                Pd = float(self.pTgt.get())
                self.axP.scatter(
                    xPeak,
                    Pd,
                    s=FONTSIZE**2,
                    c="tab:green",
                    alpha=1,
                    marker=5,  # caret right:5
                    label=self.getString("figTgtP"),
                )

                if self.plotBurnup.get():
                    self.ax.plot(
                        xs,
                        psis,
                        c="tab:red",
                        label=self.getString("figPsi"),
                        marker=".",
                        alpha=1,
                    )

                if self.plotRecoil.get():
                    self.axF.plot(
                        xs,
                        Frs,
                        c="tab:green",
                        label=self.getString("figRecoil"),
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
                        (0, xPeak),
                        (xPeak, 0.2 * xs[-1] + 0.8 * xPeak),
                        (0, xPeak),
                    ),
                ):
                    labelLines(lines, align=True, xvals=xvals, outline_width=4)
                    linesLabeled.append(lines)

                xmax = xs[-1]

                Fo = None

                if (
                    self.outflowData is not None
                    and dom == DOMAIN_TIME
                    and self.plotOutflow.get()
                ):
                    xo = []
                    Fo = []

                    Cf = Recoiless.getCf(gun.theta + 1, 1, tol)
                    # outflow force is multiplied by thrust factor

                    for t, Q, F, M in self.outflowData:
                        xo.append(t * 1e3)
                        Fo.append(F * Cf / 1e6)

                    xmax = xo[-1]

                    if self.plotRecoil.get():
                        (l,) = self.axF.plot(
                            xo,
                            Fo,
                            c="tab:green",
                            label=self.getString("figRecoil"),
                            linestyle="dotted",
                        )
                        labelLine(
                            l,
                            x=0.5 * (xs[-1] + xmax),
                            align=True,
                            outline_width=4,
                        )

                    self.ax.axvline(x=xs[-1], color="grey", ls=":")

                _, ti, li, _, vi, pi, Ti, etai = self.readTable(POINT_PEAK_SHOT)

                if gunType == CONVENTIONAL:
                    pi, _ = gun.toPsPb(li, pi)
                elif gunType == RECOILESS:
                    pi, _, _, _ = gun.toPsP0PxVx(li, vi, pi, Ti, etai)

                if self.plotBaseP.get():
                    self.axP.scatter(
                        ti * 1e3 if dom == DOMAIN_TIME else li,
                        pi / 1e6,
                        marker="+",
                        s=FONTSIZE**2,
                        c="yellowgreen",
                    )

                _, ti, li, _, vi, pi, Ti, etai = self.readTable(
                    POINT_PEAK_BREECH
                )

                if gunType == CONVENTIONAL:
                    _, pi = gun.toPsPb(li, pi)
                elif gunType == RECOILESS:
                    _, _, pi, _ = gun.toPsP0PxVx(li, vi, pi, Ti, etai)

                if self.plotBreechNozzleP.get():
                    self.axP.scatter(
                        ti * 1e3 if dom == DOMAIN_TIME else li,
                        pi / 1e6,
                        marker="+",
                        s=FONTSIZE**2,
                        c="xkcd:goldenrod",
                    )

                self.ax.set_xlim(left=0, right=xmax)
                self.ax.set_ylim(bottom=0, top=1.05)
                pmax = max(Pas + Pbs + Pss + P0s)
                self.axP.set(ylim=(0, pmax * 1.05))
                self.axv.set(ylim=(0, max(vs) * 1.05))
                self.axF.set(ylim=(0, pmax * gun.S * 1.05))

                tkw = dict(size=4, width=1.5)
                self.ax.yaxis.tick_right()
                self.axF.yaxis.tick_left()
                self.ax.tick_params(axis="y", colors="tab:red", **tkw)
                self.axv.tick_params(axis="y", colors="tab:blue", **tkw)
                self.axP.tick_params(axis="y", colors="tab:green", **tkw)
                self.axF.tick_params(axis="y", colors="tab:green", **tkw)
                self.ax.tick_params(axis="x", **tkw)

                if dom == DOMAIN_TIME:
                    self.ax.set_xlabel(self.getString("figTimeDomain"))
                elif dom == DOMAIN_LENG:
                    self.ax.set_xlabel(self.getString("figLengDomain"))

            else:
                self.axP.spines.right.set_position(("axes", 0.5))
                self.ax.set_xlabel(" ")

            self.axP.yaxis.set_ticks(self.axP.get_yticks()[1:-1:])

            self.fig.set_layout_engine("constrained")
            self.pltCanvas.draw_idle()

    def addTblFrm(self, parent):
        columnList = self.getString("columnList")
        tblFrm = ttk.LabelFrame(parent, text=self.getString("tblFrmLabel"))
        tblFrm.grid(row=2, column=0, sticky="nsew")

        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(0, weight=1)
        # configure the numerical
        self.tblFrm = tblFrm
        self.tv = ttk.Treeview(
            tblFrm, selectmode="browse", height=10
        )  # this set the nbr. of values
        self.tv.grid(row=0, column=0, sticky="nsew")

        self.tv["columns"] = columnList
        self.tv["show"] = "headings"
        self.tv.tag_configure(POINT_PEAK, foreground="#2ca02c")
        self.tv.tag_configure(POINT_PEAK_BREECH, foreground="orange")
        self.tv.tag_configure(POINT_PEAK_SHOT, foreground="yellow green")
        self.tv.tag_configure(POINT_BURNOUT, foreground="red")
        self.tv.tag_configure(POINT_FRACTURE, foreground="brown")

        t_Font = tkFont.Font(family=FONTNAME, size=FONTSIZE)

        self.tv.tag_configure("monospace", font=t_Font)
        self.tv.tag_configure("error", font=t_Font, foreground="grey")

        # we use a fixed width font so any char will do
        fontWidth, _ = t_Font.measure("m"), t_Font.metrics("linespace")

        for i, column in enumerate(columnList):  # foreach column
            self.tv.heading(
                i, text=column, anchor="e"
            )  # let the column heading = column name
            self.tv.column(
                column,
                stretch=True,  # will adjust to window resizing
                width=fontWidth * 16,
                minwidth=fontWidth * 16,
                anchor="e",
            )

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

        """
        self.tv.configure(
            yscrollcommand=vertscroll.set
        )  # assign the scrollbar to the Treeview Widget
        """

    def updateSpec(self, *args):
        self.specs.config(state="normal")
        compo = self.compositions[self.dropProp.get()]
        self.specs.delete("1.0", "end")
        t_Font = tkFont.Font(font=self.specs.cget("font"))
        width = self.specs.winfo_width() // t_Font.measure("m")
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} K {:}\n".format(
                self.getString("TvDesc"),
                compo.T_v,
                self.getString("isochorDesc"),
            ),
        ),

        self.specs.insert(
            "end",
            "{:}: {:>4.0f} kg/m^3\n".format(
                self.getString("densityDesc"), compo.rho_p
            ),
        )
        isp = compo.getIsp()
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} m/s {:>3.0f} s\n".format(
                self.getString("vacISPDesc"), isp, isp / 9.805
            ),
        )
        isp = compo.getIsp(50)
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} m/s {:>3.0f} s\n{:}\n".format(
                self.getString("atmISPDesc"),
                isp,
                isp / 9.805,
                self.getString("pRatioDesc"),
            ),
        )
        self.specs.insert("end", "{:}:\n".format(self.getString("brDesc")))
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
            self.lengthPrimaryAs.set(self.getString("diamLabel"))

            self.lengthRatioAs.set("")
            self.ratioAs.set("")

            self.lengthPrimaryTip.set(self.getString("diaText"))
            self.lengthRatioTip.set("")
            self.lengthSecondaryTip.set("")

        elif geom == SimpleGeometry.ROD:
            self.lengthPrimaryAs.set(self.getString("widtLabel"))
            self.lengthRatioAs.set(self.getString("ltwLabel"))
            self.ratioAs.set(self.getString("htwLabel"))

            self.lengthPrimaryTip.set(self.getString("widthText"))
            self.lengthRatioTip.set(self.getString("rodRText"))
            self.lengthSecondaryTip.set(self.getString("heightRText"))

        elif geom == SimpleGeometry.CYLINDER:
            self.lengthPrimaryAs.set(self.getString("diamLabel"))
            self.lengthRatioAs.set(self.getString("ltdLabel"))
            self.ratioAs.set("")

            self.lengthPrimaryTip.set(self.getString("diaText"))
            self.lengthRatioTip.set(self.getString("cylLRText"))
            self.lengthSecondaryTip.set("")

        else:
            self.lengthPrimaryAs.set(self.getString("athLabel"))
            self.lengthRatioAs.set(self.getString("ltdLabel"))
            self.ratioAs.set(self.getString("pdtalLabel"))

            self.lengthPrimaryTip.set(self.getString("arcText"))
            self.lengthRatioTip.set(self.getString("perfLRText"))
            self.lengthSecondaryTip.set(self.getString("pDiaRText"))

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

    def callback(self, *args):
        """
        updates the propellant object on write to the ratio entry fields
        and, on changing the propellant or geometrical specification.

        Double calling is due to value validation, no workaround has been
        found at this time!
        """

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
        invGunTypeLookup = {
            self.getString("CONVENTIONAL"): CONVENTIONAL,
            self.getString("RECOILESS"): RECOILESS,
        }
        gunType = invGunTypeLookup[self.typeOptn.get()]

        if gunType == CONVENTIONAL:
            self.nozzExpw.config(state="disabled")
            self.nozzEffw.config(state="disabled")

            self.nozzExpLb.grid_remove()
            self.nozzExpw.grid_remove()
            self.nozzExpU.grid_remove()
            self.nozzEffLb.grid_remove()
            self.nozzEffw.grid_remove()
            self.nozzEffU.grid_remove()

        elif gunType == RECOILESS:
            self.nozzExpw.config(state="normal")
            self.nozzEffw.config(state="normal")

            self.nozzExpLb.grid()
            self.nozzExpw.grid()
            self.nozzExpU.grid()
            self.nozzEffLb.grid()
            self.nozzEffw.grid()
            self.nozzEffU.grid()

    def ctrlCallback(self, *args):
        if self.solve_W_Lg.get() == 0:
            self.optimizeLF.config(state="disabled")
            self.minWebw.config(state="disabled")
            self.lgmaxw.config(state="disabled")
        else:
            self.optimizeLF.config(state="normal")
            self.minWebw.config(state="normal")
            self.lgmaxw.config(state="normal")

    def cvlfConsisCallback(self, *args):
        try:
            sigfig = int(self.accExp.get())
            if self.useCv.get():  # use Cv
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

    def add2Input(
        self,
        parent,
        rowIndex=0,
        colIndex=0,
        labelText="",
        default="1.0",
        validation=None,
        entryWidth=15,
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
        return lb, e, en, rowIndex + 1

    def add3Input(
        self,
        parent,
        rowIndex=0,
        colIndex=0,
        labelText="",
        unitText="",
        default="0.0",
        validation=None,
        entryWidth=15,
        formatter=formatFloatInput,
        color=None,
        infotext=None,
    ):
        lb, e, en, _ = self.add2Input(
            parent=parent,
            rowIndex=rowIndex,
            colIndex=colIndex,
            labelText=labelText,
            default=default,
            validation=validation,
            entryWidth=entryWidth,
            formatter=formatter,
            color=color,
            infotext=infotext,
        )
        ulb = ttk.Label(parent, text=unitText)
        ulb.grid(
            row=rowIndex, column=colIndex + 2, sticky="nsew", padx=2, pady=2
        )

        return lb, e, en, ulb, rowIndex + 1

    def add12Disp(
        self,
        parent,
        rowIndex=0,
        colIndex=0,
        labelText="",
        unitText="",
        default="0.0",
        entryWidth=15,
        justify="center",
        infotext=None,
        reverse=False,
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
            row=rowIndex + 1,
            column=colIndex + (1 if reverse else 0),
            sticky="nsew",
            padx=2,
            pady=2,
        )
        ttk.Label(parent, text=unitText).grid(
            row=rowIndex + 1,
            column=colIndex + (0 if reverse else 1),
            sticky="nsew",
            padx=2,
            pady=2,
        )
        if infotext is not None:
            CreateToolTip(lb, infotext)

        return lb, e, en, rowIndex + 2

    def add122Disp(
        self,
        parent,
        rowIndex=0,
        colIndex=0,
        labelText="",
        unitText_up="",
        unitText_dn="",
        default_up="0.0",
        default_dn="0.0",
        entryWidth=15,
        justify_up="center",
        justify_dn="center",
        infotext=None,
        reverse=False,
    ):
        lb, e, en, rowIndex = self.add12Disp(
            parent=parent,
            rowIndex=rowIndex,
            colIndex=colIndex,
            labelText=labelText,
            unitText=unitText_up,
            default=default_up,
            entryWidth=entryWidth,
            justify=justify_up,
            infotext=infotext,
            reverse=reverse,
        )
        e2 = StringVar(parent)
        e2.default = default_dn
        e2.set(default_dn)
        parent.rowconfigure(rowIndex, weight=0)
        en2 = ttk.Entry(
            parent,
            textvariable=e2,
            width=entryWidth,
            state="disabled",
            justify=justify_dn,
        )
        en2.grid(
            row=rowIndex + 1,
            column=colIndex + (1 if reverse else 0),
            sticky="nsew",
            padx=2,
            pady=2,
        )
        ttk.Label(parent, text=unitText_dn).grid(
            row=rowIndex + 1,
            column=colIndex + (0 if reverse else 1),
            sticky="nsew",
            padx=2,
            pady=2,
        )

        return lb, e, e2, en, en2, rowIndex + 2

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

            # this is necessary because twinned axis does not necessarily
            # follow the rcParams
            for ax in (self.ax, self.axv, self.axP, self.geomAx, self.axF):
                ax.set_facecolor(ebgc)
                ax.spines["top"].set_color(fgc)
                ax.spines["bottom"].set_color(fgc)
                ax.spines["left"].set_color(fgc)
                ax.spines["right"].set_color(fgc)

            self.updateGeomPlot()
            self.updateFigPlot()

        except AttributeError:
            pass

        self.update_idletasks()


def calculate(
    queue,
    kwargs,
):
    tableData = []
    errorData = []
    errorReport = []
    intgRecord = []

    gunType = kwargs["typ"]
    constrain = kwargs["con"]
    optimize = kwargs["opt"]
    debug = kwargs["deb"]

    sigfig = int(-log10(kwargs["tol"]))

    try:
        if constrain:
            if gunType == CONVENTIONAL:
                constrained = Constrained(**kwargs)
            elif gunType == RECOILESS:
                constrained = ConstrainedRecoiless(**kwargs)

            if optimize:
                l_f, e_1, l_g = constrained.findMinV(**kwargs)
                kwargs.update({"loadFraction": roundSig(l_f, n=sigfig)})
            else:
                e_1, l_g = constrained.solve(**kwargs)

            kwargs.update({"grainSize": roundSig(2 * e_1, n=sigfig)})
            kwargs.update({"lengthGun": roundSig(l_g, n=sigfig)})

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

        tableData, errorData = gun.integrate(
            **kwargs,
            record=intgRecord,
        )

    except Exception as e:
        gun = None
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

    queue.put(kwargs)
    queue.put(gun)
    queue.put(intgRecord)
    queue.put(tableData)
    queue.put(errorData)
    queue.put(errorReport)


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
    myappid = "Phoenix.Internal Ballistics.Solver.046"  # arbitrary string
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

    root.title("PIBS v0.4.6")

    ibPanel = IB(root, dpi, scale)

    center(root)

    root.minsize(root.winfo_width(), root.winfo_height())  # set minimum size
    root.state("zoomed")  # maximize window

    root.mainloop()
