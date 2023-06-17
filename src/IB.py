GEOM_CONTEXT = {
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.titlesize": 10,
    # "figure.autolayout": True,
    "lines.markersize": 2,
    "axes.axisbelow": True,
    # "path.simplify_threshold": 1,
}

FIG_CONTEXT = {
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 12,
    # "figure.autolayout": True,
    "lines.linewidth": 1,
    "font.weight": "bold",
    "lines.markersize": 4,
    "axes.axisbelow": False,
    "axes.labelweight": "bold",
    "yaxis.labellocation": "top",
    # "path.simplify_threshold": 1,
}

from tkinter import *
from tkinter import ttk
import tkinter.font as tkFont
import traceback

from gun import Gun, DOMAIN_TIME, DOMAIN_LENG
from gun import POINT_PEAK, POINT_BURNOUT, POINT_FRACTURE
from recoiless import Recoiless

from prop import Propellant, GrainComp, GEOMETRIES, SimpleGeometry
from opt import Constrained
from optRecoiless import ConstrainedRecoiless
from tip import ToolTip, CreateToolTip
from lang import STRING

from misc import (
    toSI,
    validateNN,
    validatePI,
    formatFloatInput,
    formatIntInput,
    dot_aligned,
    resolvepath,
)
from math import ceil, floor, log10

import matplotlib.pyplot as mpl
from matplotlib.figure import Figure
from labellines import labelLine, labelLines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Cursor


from multiprocessing import Process, Queue
from queue import Empty

import sys

RECOILESS = "Recoiless Gun"
CONVENTIONAL = "Conventional Gun"


class IB(Frame):
    def __init__(self, parent, dpi):
        Frame.__init__(self, parent)
        self.LANG = StringVar(value=list(STRING.keys())[0])
        self.queue = Queue()
        self.process = None
        self.pos = -1
        self.dpi = dpi
        self.parent = parent
        self.forceUpdOnThemeWidget = []

        menubar = Menu(parent)

        self.menubar = menubar

        themeMenu = Menu(menubar)
        menubar.add_cascade(label=self.getString("themeLabel"), menu=themeMenu)
        debugMenu = Menu(menubar)
        menubar.add_cascade(label=self.getString("debugLabel"), menu=debugMenu)
        langMenu = Menu(menubar)
        menubar.add_cascade(label="Lang 语言", menu=langMenu)

        self.themeMenu = themeMenu
        self.debugMenu = debugMenu

        self.themeRadio = IntVar()
        self.themeRadio.set(1)
        self.useTheme()

        self.DEBUG = IntVar()
        self.DEBUG.set(1)

        themeMenu.add_radiobutton(
            label=self.getString("darkLabel"),
            variable=self.themeRadio,
            value=1,
            command=self.useTheme,
        )
        themeMenu.add_radiobutton(
            label=self.getString("lightLabel"),
            variable=self.themeRadio,
            value=2,
            command=self.useTheme,
        )

        debugMenu.add_checkbutton(
            label=self.getString("enableLabel"),
            variable=self.DEBUG,
            onvalue=1,
            offvalue=0,
        )

        for lang in STRING.keys():
            langMenu.add_radiobutton(
                label=lang,
                variable=self.LANG,
                value=lang,
                command=self.changeLang,
            )

        parent.config(menu=menubar)

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
        parent.rowconfigure(1, weight=1)

        self.addTopBar(parent)

        self.addRightFrm(parent)
        self.addErrFrm(parent)
        self.addPlotFrm(parent)
        self.addParFrm(parent)
        self.addTblFrm(parent)

        parent.update_idletasks()
        self.addGeomPlot()

        parent.update_idletasks()
        self.addFigPlot()

        self.updateSpec()
        self.updateGeom()

        self.forceUpdOnThemeWidget.append(self.errorText)
        self.forceUpdOnThemeWidget.append(self.specs)

        parent.bind("<Configure>", self.resizePlot)

        # self.resized = False
        self.timedLoop()

    def changeLang(self):
        self.menubar.entryconfig(0, label=self.getString("themeLabel"))
        self.menubar.entryconfig(1, label=self.getString("debugLabel"))
        self.themeMenu.entryconfig(0, label=self.getString("darkLabel"))
        self.themeMenu.entryconfig(1, label=self.getString("lightLabel"))
        self.debugMenu.entryconfig(0, label=self.getString("enableLabel"))

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

        self.nozzExpLb.config(text=self.getString("nozzExpLabel"))
        self.nozzEffLb.config(text=self.getString("nozzEffLabel"))

        self.chgTip.set(self.getString("chgText"))
        self.vinfTip.set(self.getString("vinfText"))
        self.lxTip.set(self.getString("calLxTxt"))
        self.geomPlotTip.set(self.getString("geomPlotTxt"))
        self.teffTip.set(self.getString("teffText"))
        self.beffTip.set(self.getString("beffText"))

        self.tblFrm.config(text=self.getString("tblFrmLabel"))
        self.plotFrm.config(text=self.getString("plotFrmLabel"))
        self.errorFrm.config(text=self.getString("errFrmLabel"))
        self.parFrm.config(text=self.getString("parFrmLabel"))
        self.specFrm.config(text=self.getString("specFrmLabel"))
        self.opFrm.config(text=self.getString("opFrmLabel"))
        self.consFrm.config(text=self.getString("consFrmLabel"))
        self.topFrm.config(text=self.getString("pltOptnFrm"))

        self.useConstraint.config(text=self.getString("consButton"))
        self.optimizeLF.config(text=self.getString("minTVButton"))
        self.sampleFrm.config(text=self.getString("sampleFrmLabel"))

        self.lxLb.config(text=self.getString("lxLabel"))
        self.vaLb.config(text=self.getString("vaLabel"))
        self.pPLb.config(text=self.getString("pPLabel"))
        self.teLb.config(text=self.getString("teffLabel"))
        self.beLb.config(text=self.getString("beffLabel"))
        self.cvLb.config(text=self.getString("cvLabel"))
        self.ldpLb.config(text=self.getString("ldLabel"))

        self.propFrm.config(text=self.getString("propFrmLabel"))
        self.grainFrm.config(text=self.getString("grainFrmLabel"))

        for i, columnName in enumerate(self.getString("columnList")):
            self.tv.heading(i, text=columnName)

    def getString(self, name):
        try:
            return STRING[self.LANG.get()][name]
        except KeyError:
            return STRING["English"][name]

    def addTopBar(self, parent):
        topFrm = ttk.LabelFrame(parent, text=self.getString("pltOptnFrm"))
        topFrm.grid(row=0, column=0, sticky="nsew")
        self.topFrm = topFrm

        for i in range(8):
            topFrm.columnconfigure(i, weight=1)

        self.pbar = ttk.Progressbar(topFrm, mode="indeterminate", maximum=100)
        self.pbar.grid(
            row=0, column=0, columnspan=8, sticky="nsew", padx=2, pady=2
        )

        self.plotAvgP = IntVar(value=1)

        ttk.Checkbutton(
            topFrm, text="Length Avg. P.", variable=self.plotAvgP
        ).grid(row=1, column=0, sticky="nsew")
        self.plotBaseP = IntVar(value=1)
        ttk.Checkbutton(
            topFrm, text="Shot Base P.", variable=self.plotBaseP
        ).grid(row=1, column=1, sticky="nsew")

        self.plotBreechNozzleP = IntVar(value=1)
        ttk.Checkbutton(
            topFrm, text="Breech/Nozzle P.", variable=self.plotBreechNozzleP
        ).grid(row=1, column=2, sticky="nsew")

        self.plotStagP = IntVar(value=1)
        ttk.Checkbutton(
            topFrm, text="Stagnation P.", variable=self.plotStagP
        ).grid(row=1, column=3, sticky="nsew")

        self.plotVel = IntVar(value=1)
        ttk.Checkbutton(topFrm, text="Shot Vel.", variable=self.plotVel).grid(
            row=1, column=4, sticky="nsew"
        )

        self.plotNozzleV = IntVar(value=1)
        ttk.Checkbutton(
            topFrm, text="Nozzle Throat Vel.", variable=self.plotNozzleV
        ).grid(row=1, column=5, sticky="nsew")

        self.plotBurnup = IntVar(value=1)
        ttk.Checkbutton(topFrm, text="Burnup", variable=self.plotBurnup).grid(
            row=1, column=6, sticky="nsew"
        )
        self.plotEta = IntVar(value=1)
        ttk.Checkbutton(topFrm, text="Escape", variable=self.plotEta).grid(
            row=1, column=7, sticky="nsew"
        )

        self.plotAvgP.trace_add("write", self.updateFigPlot)
        self.plotBaseP.trace_add("write", self.updateFigPlot)
        self.plotBreechNozzleP.trace_add("write", self.updateFigPlot)
        self.plotStagP.trace_add("write", self.updateFigPlot)
        self.plotEta.trace_add("write", self.updateFigPlot)
        self.plotNozzleV.trace_add("write", self.updateFigPlot)
        self.plotBurnup.trace_add("write", self.updateFigPlot)
        self.plotVel.trace_add("write", self.updateFigPlot)
        """
        self.summary = ttk.Entry(topFrm, justify="right", state="disabled")
        self.summary.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        """

    def addRightFrm(self, parent):
        rightFrm = ttk.Frame(parent)
        rightFrm.grid(row=0, column=2, rowspan=4, sticky="nsew")
        rightFrm.columnconfigure(0, weight=1)
        rightFrm.rowconfigure(0, weight=1)

        specFrm = ttk.LabelFrame(rightFrm, text=self.getString("specFrmLabel"))
        specFrm.grid(row=0, column=0, sticky="nsew")
        specFrm.columnconfigure(0, weight=1)
        self.specFrm = specFrm

        i = 0

        self.lxTip = StringVar(value=self.getString("calLxTxt"))
        self.lxLb, self.lx, self.tlx, _, _, i = self.add122Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("lxLabel"),
            unitText_up="Cal",
            unitText_dn="Cal",
            justify_up="right",
            justify_dn="right",
            infotext=self.lxTip,
        )

        self.vinfTip = StringVar(value=self.getString("vinfText"))
        self.vaLb, self.va, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("vaLabel"),
            unitText="m/s",
            justify="right",
            infotext=self.vinfTip,
        )

        self.pPLb, self.ptm, self.pbm, _, _, i = self.add122Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("pPLabel"),
            unitText_up="Pa",
            unitText_dn="Pa",
            justify_up="right",
            justify_dn="right",
            infotext=self.getString("pMaxTxt"),
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
        self.cvLb, self.cv, _, i = self.add12Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("cvLabel"),
            unitText="m³",
            justify="right",
        )
        self.ldpLb, self.ldp, self.ld, _, _, i = self.add122Disp(
            parent=specFrm,
            rowIndex=i,
            labelText=self.getString("ldLabel"),
            unitText_up="%",
            unitText_dn="kg/m³",
            justify_dn="right",
        )

        opFrm = ttk.LabelFrame(rightFrm, text=self.getString("opFrmLabel"))
        opFrm.grid(row=1, column=0, sticky="nsew")
        opFrm.columnconfigure(1, weight=1)
        self.opFrm = opFrm

        validationNN = parent.register(validateNN)
        validationPI = parent.register(validatePI)

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

        self.vTgtLb, self.vTgt, _, j = self.add3Input(
            parent=consFrm,
            rowIndex=0,
            colIndex=0,
            labelText=self.getString("vTgtLabel"),
            unitText="m/s",
            default="1500.0",
            validation=validationNN,
        )

        self.pTgtLb, self.pTgt, _, j = self.add3Input(
            parent=consFrm,
            rowIndex=j,
            colIndex=0,
            labelText=self.getString("pTgtLabel"),
            unitText="MPa",
            default="350.0",
            validation=validationNN,
            infotext=self.getString("pTgtTxt"),
        )

        self.minWebLb, self.minWeb, _, j = self.add3Input(
            parent=consFrm,
            rowIndex=j,
            colIndex=0,
            labelText=self.getString("minWebLabel"),
            unitText="μm",
            default="1.0",
            validation=validationNN,
            color="red",
        )
        self.lgmaxLb, self.lgmax, _, j = self.add3Input(
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
        self.solve_W_Lg.trace_add("write", self.setCD)

        CreateToolTip(self.useConstraint, self.getString("useConsTxt"))

        j += 1
        self.opt_lf = IntVar()
        self.optimizeLF = ttk.Checkbutton(
            consFrm, text=self.getString("minTVButton"), variable=self.opt_lf
        )
        self.optimizeLF.grid(row=j, column=0, columnspan=3, sticky="nsew")
        self.setCD()

        CreateToolTip(self.optimizeLF, self.getString("optLFTxt"))
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
            values=self.domainOptions,
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

        self.stepsLb, self.steps, _, j = self.add2Input(
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

        CreateToolTip(sampleFrm, self.getString("sampTxt"))

        i += 1

        self.accExpLb, self.accExp, _, i = self.add2Input(
            parent=opFrm,
            rowIndex=i,
            colIndex=0,
            labelText="-log10(ε)",
            default="3",
            validation=validationPI,
            formatter=formatIntInput,
            color="red",
            infotext=self.getString("tolText"),
        )

        self.calButton = ttk.Button(
            opFrm,
            text="CALCULATE",
            # command=self.calculate,  # underline=0
            command=self.onCalculate,
        )
        self.calButton.grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )

        opFrm.rowconfigure(i, weight=1)
        CreateToolTip(
            self.calButton, "Integrate system using RKF7(8) integrator"
        )

    def onCalculate(self):
        constrain = self.solve_W_Lg.get() == 1
        optimize = self.opt_lf.get() == 1
        debug = self.DEBUG.get() == 1
        gunType = self.typeOptn.get()

        self.tableData = []
        self.errorData = []
        self.intgRecord = []
        self.kwargs = {
            "opt": optimize,
            "con": constrain,
            "deb": debug,
            "typ": gunType,
            "dom": self.dropOptn.get(),
        }
        self.process = None

        try:
            chamberVolume = (
                float(self.chgkg.get())
                / self.prop.rho_p
                / self.prop.maxLF
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
                    "chargeMassRatio": float(self.chgkg.get())
                    / float(self.shtkg.get()),
                    "chamberVolume": chamberVolume,
                    "startPressure": float(self.stpMPa.get()) * 1e6,
                    "lengthGun": float(self.tblmm.get()) * 1e-3,
                    "chamberExpansion": float(
                        self.clr.get()
                    ),  # chamber expansion
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
                    "step": int(self.steps.get()),
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
        constrain = self.kwargs["con"]
        optimize = self.kwargs["con"]
        gunType = self.kwargs["typ"]

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

        self.pos = -1
        kwargs = self.kwargs

        if self.gun is not None:
            chamberVolume = kwargs["chamberVolume"]

            self.cv.set(toSI(chamberVolume, useSN=True))

            if constrain:
                webmm = round(
                    1e3 * kwargs["grainSize"],
                    3 - int(floor(log10(abs(1e3 * kwargs["grainSize"])))),
                )
                self.arcmm.set(webmm)

                lgmm = round(
                    kwargs["lengthGun"] * 1e3,
                    3 - int(floor(log10(abs(kwargs["lengthGun"] * 1000)))),
                )

                self.tblmm.set(lgmm)
                if optimize:
                    lfpercent = round(
                        kwargs["loadFraction"] * 100,
                        3
                        - int(floor(log10(abs(kwargs["loadFraction"] * 100)))),
                    )
                    self.ldf.set(lfpercent)

            self.ldp.set(
                round(self.prop.maxLF * kwargs["loadFraction"] * 100, 1)
            )
            self.ld.set(
                toSI(
                    self.prop.maxLF * kwargs["loadFraction"] * self.prop.rho_p,
                    useSN=True,
                )
            )

            i = [i[0] for i in self.tableData].index("SHOT EXIT")
            vg = self.tableData[i][4]
            te, be = self.gun.getEff(vg)
            self.te.set(round(te * 100, 1))
            self.be.set(round(te / self.gun.phi * 100, 1))

            i = [i[0] for i in self.tableData].index("PEAK PRESSURE")
            _, tp, lp, _, vp, pp, Tp, etap = self.tableData[i]

            if gunType == CONVENTIONAL:
                Pb, Pt = self.gun.toPbPt(lp, pp)
                self.ptm.set(toSI(Pt))
                self.pbm.set(toSI(Pb))
            else:
                Pb, P0, Px, _ = self.gun.toPbP0PxVx(lp, vp, pp, Tp, etap)
                self.ptm.set(toSI(P0))
                self.pbm.set(toSI(Pb))

            self.lx.set(toSI(kwargs["lengthGun"] / kwargs["caliber"]))
            self.tlx.set(
                toSI(
                    (
                        kwargs["lengthGun"]
                        + self.gun.l_0 / kwargs["chamberExpansion"]
                    )
                    / kwargs["caliber"]
                )
            )
            self.va.set(toSI(self.gun.v_j))

        self.tv.delete(*self.tv.get_children())
        useSN = (False, False, False, True, False, False, True, True)
        units = (None, "s", "m", None, "m/s", "Pa", "K", None)
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
            wrap=WORD,
            height=5,
            width=0,
            font=("Hack", 8),
        )

        self.errorText.grid(row=0, column=0, sticky="nsew")

    def addParFrm(self, parent):
        parFrm = ttk.LabelFrame(parent, text=self.getString("parFrmLabel"))
        parFrm.grid(row=0, column=1, rowspan=4, sticky="nsew")
        parFrm.columnconfigure(0, weight=1)
        self.parFrm = parFrm
        # validation
        validationNN = parent.register(validateNN)

        i = 0
        self.gunType = StringVar()
        self.typeOptn = ttk.Combobox(
            parFrm,
            textvariable=self.gunType,
            values=(CONVENTIONAL, RECOILESS),
            state="readonly",
            justify="center",
        )
        self.typeOptn.option_add("*TCombobox*Listbox.Justify", "center")
        self.typeOptn.current(0)
        self.typeOptn.grid(
            row=i, column=0, sticky="nsew", padx=2, pady=2, columnspan=3
        )

        i += 1
        self.calLb, self.calmm, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("calLabel"),
            unitText="mm",
            default="50.0",
            validation=validationNN,
        )
        self.tblLb, self.tblmm, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("tblLabel"),
            unitText="mm",
            default="3500.0",
            validation=validationNN,
        )
        self.shtLb, self.shtkg, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("shtLabel"),
            unitText="kg",
            default="1.0",
            validation=validationNN,
        )

        self.chgTip = StringVar(value=self.getString("chgText"))
        self.chgLb, self.chgkg, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("chgLabel"),
            unitText="kg",
            default="1.0",
            validation=validationNN,
            infotext=self.chgTip,
        )

        # allow propellant specification to grow
        parFrm.rowconfigure(i, weight=3)

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
        # self.dropProp.configure(width=10)
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
            height=10,
            width=36,
            yscrollcommand=specScroll.set,
            xscrollcommand=specHScroll.set,
            font=("Hack", 8),
        )
        self.specs.grid(row=1, column=0, sticky="nsew")
        specScroll.config(command=self.specs.yview)
        specHScroll.config(command=self.specs.xview)

        CreateToolTip(propFrm, self.getString("specsText"))

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
            width=100,
            height=100,
        )  # set an arbitrary initial size here
        # to prevent matplotlib giving a fit when adding geomFig.

        geomPlotFrm.grid(
            row=0,
            column=0,
            columnspan=3,
            sticky="nsew",
            padx=2,
            pady=2,
        )

        self.geomPlotTip = StringVar(value=self.getString("geomPlotTxt"))
        CreateToolTip(geomPlotFrm, self.geomPlotTip)

        self.geomParentFrm = grainFrm
        self.geomPlotFrm = geomPlotFrm

        j = 1

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
        self.arcLb, self.arcmm, _, j = self.add3Input(
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
        self.grdRLb, self.grdR, self.grdRw, j = self.add3Input(
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

        self.grlRLb, self.grlR, self.grlRw, j = self.add3Input(
            parent=grainFrm,
            rowIndex=j,
            labelText=self.lengthRatioAs,
            unitText="x",
            default="2.5",
            validation=validationNN,
            infotext=self.lengthRatioTip,
        )

        i += 1

        self.ldfLb, self.ldf, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("ldfLabel"),
            unitText="%",
            default="50.0",
            validation=validationNN,
            infotext=self.getString("ldftext"),
        )

        self.clrLb, self.clr, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("clrLabel"),
            unitText="x",
            default="1.5",
            validation=validationNN,
            infotext=self.getString("clrtext"),
        )

        self.dgcLb, self.dgc, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("dgcLabel"),
            unitText="%",
            default="5.0",
            validation=validationNN,
            infotext=self.getString("dgctext"),
        )

        self.stpLb, self.stpMPa, _, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("stpLabel"),
            unitText="MPa",
            default="10",
            validation=validationNN,
            infotext=self.getString("stpText"),
        )

        self.nozzExpLb, self.nozzExp, self.nozzExpw, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("nozzExpLabel"),
            unitText="x",
            default="4",
            validation=validationNN,
            infotext=self.getString("nozzExpTxt"),
        )

        self.nozzEffLb, self.nozzEff, self.nozzEffw, i = self.add3Input(
            parent=parFrm,
            rowIndex=i,
            labelText=self.getString("nozzEffLabel"),
            unitText="%",
            default="92.0",
            validation=validationNN,
            infotext=self.getString("nozzEffTxt"),
        )

        self.currProp.trace_add("write", self.updateSpec)
        self.currGeom.trace_add("write", self.updateGeom)

        self.grdR.trace_add("write", self.callback)
        self.grlR.trace_add("write", self.callback)
        self.arcmm.trace_add("write", self.callback)

        self.gunType.trace_add("write", self.typeCallback)
        self.typeCallback()

    def addGeomPlot(self):
        geomParentFrm = self.geomParentFrm
        geomPlotFrm = self.geomPlotFrm
        geomPlotFrm.columnconfigure(0, weight=1)
        geomPlotFrm.rowconfigure(0, weight=1)

        # _, _, width, height = self.geomPlotFrm.bbox("insert")
        width = geomPlotFrm.winfo_width() - 2
        height = geomPlotFrm.winfo_height() - 2

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

            self.geomCanvas = FigureCanvasTkAgg(
                fig,
                master=geomPlotFrm,
            )
            self.geomCanvas.get_tk_widget().grid(
                row=0, column=0, padx=0, pady=0, sticky="nsew"
            )

    def addPlotFrm(self, parent):
        plotFrm = ttk.LabelFrame(parent, text=self.getString("plotFrmLabel"))
        plotFrm.grid(row=1, column=0, sticky="nsew")
        plotFrm.columnconfigure(0, weight=1)
        plotFrm.rowconfigure(0, weight=1)

        self.plotFrm = plotFrm

    def addFigPlot(self):
        plotFrm = self.plotFrm
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

            ax.yaxis.tick_right()
            ax.set_xlabel("Domain")

            axv.spines.right.set_position(("axes", 1.0 + 40 * dpi / 96 / width))
            axP.spines.right.set_position(("data", 0.5))

            axP.yaxis.set_ticks(axP.get_yticks()[1:-1:])

            self.ax = ax
            self.axP = axP
            self.axv = axv
            self.fig = fig

            self.pltCanvas = FigureCanvasTkAgg(fig, master=plotFrm)
            self.pltCanvas.get_tk_widget().grid(
                row=0, column=0, padx=2, pady=2, sticky="nsew"
            )

    def resizePlot(self, event):
        # we use the bbox method here as it has already accounted for padding
        # so no adjustment here is necessary

        _, _, width, height = self.plotFrm.bbox("insert")

        dpi = self.dpi

        with mpl.rc_context(FIG_CONTEXT):
            self.axv.spines.right.set_position(
                ("axes", 1 + 40 * dpi / 96 / width)
            )

        # self.resized = True

    def timedLoop(self):
        """
        if self.resized:
            self.updateSpec(None, None, None)
            self.resized = False
        """

        if self.pos >= 0:  # and not self.process.is_alive():
            self.getValue()

        self.parent.after(100, self.timedLoop)

    def updateFigPlot(self, *args):
        with mpl.rc_context(FIG_CONTEXT):
            gun = self.gun
            try:
                gunType = self.kwargs["typ"]
            except AttributeError:
                gunType = CONVENTIONAL

            self.ax.cla()
            self.axP.cla()
            self.axv.cla()

            self.axv.set_axisbelow(False)
            self.axP.set_axisbelow(False)
            dpi = self.dpi
            size = self.fig.get_size_inches() * self.fig.dpi

            self.axv.spines.right.set_position(
                ("axes", 1 + 40 * dpi / 96 / size[0])
            )

            if gun is not None:
                xs = []
                vs = []
                Ps = []
                Pbs = []
                Pts = []
                P0s = []
                Pxs = []
                psis = []
                etas = []
                vxs = []
                dom = self.dropOptn.get()

                for i, (t, (l, psi, v, p, T, eta)) in enumerate(
                    self.intgRecord
                ):
                    if dom == DOMAIN_TIME:
                        xs.append(t * 1000)
                    elif dom == DOMAIN_LENG:
                        xs.append(l)
                    vs.append(v)
                    Ps.append(p / 1e6)
                    if gunType == CONVENTIONAL:
                        Pb, Pt = gun.toPbPt(l, p)
                        P0, Px = 0, 0
                        vx = 0
                    else:
                        Pb, P0, Px, vx = gun.toPbP0PxVx(l, v, p, T, eta)
                        Pt = 0
                    Pbs.append(Pb / 1e6)
                    Pts.append(Pt / 1e6)
                    P0s.append(P0 / 1e6)
                    Pxs.append(Px / 1e6)
                    vxs.append(vx)
                    psis.append(psi)
                    etas.append(eta)

                if self.plotVel.get():
                    self.axv.scatter(xs, vs, color="tab:blue", marker="s", s=8)
                if self.plotAvgP.get():
                    self.axP.scatter(xs, Ps, color="tab:green", marker="s", s=8)
                if self.plotBurnup.get():
                    self.ax.scatter(xs, psis, color="tab:red", marker="s", s=8)

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
                    Ps.append(p / 1e6)

                    if gunType == CONVENTIONAL:
                        Pb, Pt = gun.toPbPt(l, p)
                        P0, Px = 0, 0
                        vx = 0
                    else:
                        Pb, P0, Px, vx = gun.toPbP0PxVx(l, v, p, T, eta)
                        Pt = 0

                    Pbs.append(Pb / 1e6)
                    Pts.append(Pt / 1e6)
                    P0s.append(P0 / 1e6)
                    Pxs.append(Px / 1e6)
                    vxs.append(vx)
                    psis.append(psi)
                    etas.append(eta)

                self.axP.spines.right.set_position(("data", xPeak))

                self.ax.set_xlim(left=0, right=xs[-1])
                self.ax.set_ylim(bottom=0, top=1.05)

                self.axP.set(ylim=(0, max(Ps + Pbs + Pts + P0s + Pxs) * 1.05))
                self.axv.set(ylim=(0, max(vs) * 1.05))

                (xs, vs, vxs, Ps, Pbs, Pts, P0s, Pxs, psis, etas) = zip(
                    *sorted(
                        zip(xs, vs, vxs, Ps, Pbs, Pts, P0s, Pxs, psis, etas),
                        key=lambda line: line[0],
                    )
                )
                if self.plotVel.get():
                    self.axv.plot(
                        xs,
                        vs,
                        "tab:blue",
                        label="Shot Velocity\nm/s",
                        marker=".",
                        alpha=1,
                    )
                vd = float(self.vTgt.get())
                self.axv.scatter(
                    xs[-1],
                    vd,
                    8**2,
                    "tab:blue",
                    marker=5,
                    alpha=1,
                )

                if self.typeOptn.get() == CONVENTIONAL:
                    if self.plotBreechNozzleP.get():
                        self.axP.plot(
                            xs,
                            Pts,
                            "seagreen",
                            label="Breech Face",
                            linestyle="dashed",
                            alpha=0.75,
                        )
                else:
                    if self.plotStagP.get():
                        self.axP.plot(
                            xs,
                            P0s,
                            "yellowgreen",
                            label="Stagnation",
                            linestyle="dashed",
                            alpha=0.75,
                        )

                    if self.plotBreechNozzleP.get():
                        self.axP.plot(
                            xs,
                            Pxs,
                            "seagreen",
                            label="Nozz. Throat P.",
                            linestyle="dashdot",
                            alpha=0.75,
                        )
                    if self.plotNozzleV.get():
                        self.axv.plot(
                            xs,
                            vxs,
                            "royalblue",
                            label="Nozz. Throat Vel.",
                            alpha=0.75,
                            linestyle="dotted",
                        )

                    if self.plotEta.get():
                        self.ax.plot(
                            xs,
                            etas,
                            "tab:orange",
                            label="Outflow Frac.",
                            alpha=0.75,
                            linestyle="dotted",
                        )

                if self.plotAvgP.get():
                    self.axP.plot(
                        xs,
                        Ps,
                        "tab:green",
                        label="Avg. Pressure\nMPa",
                        marker=".",
                        alpha=1,
                    )

                if self.plotBaseP.get():
                    self.axP.plot(
                        xs,
                        Pbs,
                        "olive",
                        label="Shot Base",
                        linestyle="dotted",
                        alpha=0.75,
                    )

                Pd = float(self.pTgt.get())
                self.axP.scatter(
                    xPeak,
                    Pd,
                    s=8**2,
                    c="tab:green",
                    alpha=1,
                    marker=5,  # caret right:5
                    label="P. Target",
                )
                if self.plotBurnup.get():
                    self.ax.plot(
                        xs,
                        psis,
                        "tab:red",
                        label="Volume Burnup",
                        marker=".",
                        alpha=1,
                    )

                tkw = dict(size=4, width=1.5)
                self.ax.yaxis.tick_right()
                self.ax.tick_params(axis="y", colors="tab:red", **tkw)
                self.axv.tick_params(axis="y", colors="tab:blue", **tkw)
                self.axP.tick_params(axis="y", colors="tab:green", **tkw)
                self.ax.tick_params(axis="x", **tkw)

                if dom == DOMAIN_TIME:
                    self.ax.set_xlabel("Time - ms")
                elif dom == DOMAIN_LENG:
                    self.ax.set_xlabel("Length - m")

                for lines, xvals in zip(
                    (
                        self.axP.get_lines(),
                        self.ax.get_lines(),
                        self.axv.get_lines(),
                    ),
                    (
                        (0.2 * xs[-1] + 0.8 * xPeak, xs[-1]),
                        (0, xPeak),
                        (xPeak, 0.2 * xs[-1] + 0.8 * xPeak),
                    ),
                ):
                    labelLines(lines, align=True, xvals=xvals)

            else:
                self.axP.spines.right.set_position(("axes", 0.5))
                self.ax.set_xlabel("Domain")

            self.axP.yaxis.set_ticks(self.axP.get_yticks()[1:-1:])
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
        self.tv.tag_configure(POINT_PEAK, foreground="orange")
        self.tv.tag_configure(POINT_BURNOUT, foreground="red")
        self.tv.tag_configure(POINT_FRACTURE, foreground="brown")

        t_Font = tkFont.Font(family="hack", size=8)

        self.tv.tag_configure("monospace", font=t_Font)
        self.tv.tag_configure("error", font=("hack", 8), foreground="grey")

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
        vertscroll.grid(row=0, column=1, sticky="nsew")

        """
        horzscroll = ttk.Scrollbar(tblFrm, orient="horizontal")
        horzscroll.configure(command=self.tv.xview)
        horzscroll.grid(row=1, column=0, sticky="nsew")
        self.tv.configure(
            yscrollcommand=vertscroll.set, xscrollcommand=horzscroll.set
        )  # assign the scrollbar to the Treeview Widget
        """

    def updateSpec(self, *args):
        self.specs.config(state="normal")
        compo = self.compositions[self.dropProp.get()]
        self.specs.delete("1.0", "end")
        t_Font = tkFont.Font(family="hack", size=8)
        width = self.specs.winfo_width() // t_Font.measure("m")
        self.specs.insert(
            "end", " Adb.Temp: {:>4.0f} K (Isochoric)\n".format(compo.T_v)
        ),

        self.specs.insert(
            "end", "  Density: {:>4.0f} kg/m^3\n".format(compo.rho_p)
        )
        isp = compo.getIsp()
        self.specs.insert(
            "end",
            " Isp(Vac): {:>4.0f} m/s {:>3.0f} s\n".format(isp, isp / 9.805),
        )
        isp = compo.getIsp(50)
        self.specs.insert(
            "end",
            " Isp(Atm): {:>4.0f} m/s {:>3.0f} s (Pc:Pa=50)\n".format(
                isp, isp / 9.805
            ),
        )
        self.specs.insert("end", "Burn rate:\n")
        for p in (100e6, 200e6, 300e6):
            self.specs.insert(
                "end",
                "{:>13}".format(toSI(compo.getLBR(p), unit="m/s", dec=3))
                + " @ {:>13}\n".format(toSI(p, unit="Pa", dec=3)),
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
            self.lengthPrimaryAs.set("Diameter")

            self.lengthRatioAs.set("")
            self.ratioAs.set("")

            self.lengthPrimaryTip.set(self.getString("diaText"))
            self.lengthRatioTip.set("")
            self.lengthSecondaryTip.set("")

        elif geom == SimpleGeometry.ROD:
            self.lengthPrimaryAs.set("Width")
            self.lengthRatioAs.set("Length / Width")
            self.ratioAs.set("Height / Width")

            self.lengthPrimaryTip.set(self.getString("widthText"))
            self.lengthRatioTip.set(self.getString("rodRtext"))
            self.lengthSecondaryTip.set(self.getString("heightRtext"))

        elif geom == SimpleGeometry.CYLINDER:
            self.lengthPrimaryAs.set("Diameter")

            self.lengthRatioAs.set("Length / Diameter")
            self.ratioAs.set("")

            self.lengthPrimaryTip.set(self.getString("diaText"))
            self.lengthRatioTip.set(self.getString("cylLRtext"))
            self.lengthSecondaryTip.set("")

        else:
            self.lengthPrimaryAs.set("Arc Thickness")
            self.lengthRatioAs.set("Length / Diameter")
            self.ratioAs.set("Perf.Dia. / A.Th.")

            self.lengthPrimaryTip.set(self.getString("arcText"))
            self.lengthRatioTip.set(self.getString("perfLRtext"))
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

                # xs, ys = zip(*sorted(zip(xs, ys), key=lambda x: x[0]))
                self.geomAx.plot(xs, ys)
                self.geomAx.grid(
                    which="major", color="grey", linestyle="dotted"
                )
                self.geomAx.minorticks_on()
                self.geomAx.set_xlim(left=0, right=prop.Z_b)
                self.geomAx.xaxis.set_ticks(
                    [i * 0.5 for i in range(ceil(prop.Z_b / 0.5) + 1)]
                )
                self.geomAx.set_ylim(bottom=0, top=max(ys))
                self.geomAx.yaxis.set_ticks(
                    [i * 0.5 for i in range(ceil(max(ys) / 0.5) + 1)]
                )

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

        Double calling is due value validation, no workaround has been found
        at this time!
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
        gunType = self.typeOptn.get()
        if gunType == CONVENTIONAL:
            self.nozzExpw.config(state="disabled")
            self.nozzEffw.config(state="disabled")
        else:
            self.nozzExpw.config(state="normal")
            self.nozzEffw.config(state="normal")

    def setCD(self, *args):
        if self.solve_W_Lg.get() == 0:
            self.optimizeLF.config(state="disabled")
        else:
            self.optimizeLF.config(state="normal")

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
        return lb, e, en, rowIndex + 1

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
        lb, e, en, _ = self.add2Input(
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
        return lb, e, en, rowIndex + 1

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
        rowIndex,
        colIndex=0,
        labelText="",
        unitText_up="",
        unitText_dn="",
        default_up="0.0",
        default_dn="0.0",
        entryWidth=5,
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
        if self.themeRadio.get() == 1:
            style.theme_use("awdark")
        else:
            style.theme_use("awlight")
        self.setTheme()

    def setTheme(self):
        dpi = self.dpi
        style = ttk.Style(self)
        # ensure that the treeview rows are roughly the same height
        # regardless of dpi. on Windows, default is Segoe UI at 9 points
        # so the default row height should be around 12

        style.configure("Treeview", rowheight=round(12 * dpi / 72.0))
        style.configure("Treeview.Heading", font=("Hack", 9))
        style.configure("TButton", font=("Hack", 10, "bold"))
        style.configure("TLabelframe.Label", font=("Hack", 11, "bold"))
        style.configure("TCheckbutton", font=("Hack", 9))
        style.configure("SubLabelFrame.TLabelframe.Label", font=("Hack", 10))
        # style.configure("TNotebook.Tab", font=("Hack", 10))

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
            }
        )

        FIG_CONTEXT.update(
            {
                "figure.facecolor": bgc,
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

            for ax in (self.ax, self.axv, self.axP, self.geomAx):
                ax.set_facecolor(ebgc)
                ax.spines["top"].set_color(fgc)
                ax.spines["bottom"].set_color(fgc)
                ax.spines["left"].set_color(fgc)
                ax.spines["right"].set_color(fgc)

            self.updateGeomPlot()
            self.updateFigPlot()

        except AttributeError as e:
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
    try:
        if constrain:
            if gunType == CONVENTIONAL:
                constrained = Constrained(**kwargs)
            else:
                constrained = ConstrainedRecoiless(**kwargs)

            if optimize:
                l_f, e_1, l_g = constrained.findMinV(**kwargs)

                kwargs.update(
                    {
                        "loadFraction": round(
                            l_f, 3 - int(floor(log10(abs(l_f))))
                        )
                    }
                )

            else:
                e_1, l_g = constrained.solve(**kwargs)

            kwargs.update(
                {
                    "grainSize": round(
                        2 * e_1, 3 - int(floor(log10(abs(2 * e_1))))
                    )
                }
            )
            kwargs.update(
                {"lengthGun": round(l_g, 3 - int(floor(log10(abs(l_g)))))}
            )

            chamberVolume = (
                kwargs["chargeMass"]
                / kwargs["propellant"].rho_p
                / kwargs["propellant"].maxLF
                / kwargs["loadFraction"]
            )

            kwargs.update({"chamberVolume": chamberVolume})

        else:
            pass

        if gunType == CONVENTIONAL:
            gun = Gun(**kwargs)
        else:
            gun = Recoiless(**kwargs)

        tableData, errorData = gun.integrate(
            steps=kwargs["step"],
            dom=kwargs["dom"],
            tol=kwargs["tol"],
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
