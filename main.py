from tkinter import *
from tkinter import ttk
from idlelib.tooltip import Hovertip
from gun import *


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


def formatInput(event):
    v = event.widget.get()
    if v == "":
        event.widget.insert(0, event.widget.default)
    else:
        event.widget.delete(0, END)
        event.widget.insert(0, float(v))


def dot_aligned(matrix):
    transposed = []
    # print(matrix)

    for seq in zip(*matrix):
        # print(seq)
        snums = [str(n) for n in seq]
        dots = [s.find(".") for s in snums]
        m = max(dots)
        transposed.append(tuple(" " * (m - d) + s for s, d in zip(snums, dots)))

    return tuple(zip(*transposed))


class IB(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        parent.title("Phoenix's Internal Ballistics Solver")
        # parent.geometry("600x400")

        self.compositions = GrainComp.readFile("propellants.csv")
        self.geometries = {i.desc: i for i in Geometry}

        self.propOptions = tuple(self.compositions.keys())
        self.geoOptions = tuple(self.geometries.keys())

        self.prop = None
        self.gun = None
        self.tableData = []
        self.errorLst = []
        self.geomError = False

        # datatype of menu text
        self.clickedProp = StringVar()
        # initial menu text
        self.clickedProp.set(self.propOptions[0])

        self.clickedGeom = StringVar()
        self.clickedGeom.set(self.geoOptions[0])
        self.propBanner = StringVar()
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
        compo = self.compositions[self.clickedProp.get()]
        geom = self.geometries[self.clickedGeom.get()]

        self.tableData = []

        try:
            self.prop = Propellant(
                compo,
                geom,
                float(self.webmm.get()) / 1000,
                float(self.permm.get()) / 1000,
                float(self.grlR.get()),
            )
        except Exception as e:
            print(e)
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

        except Exception as e:
            self.gun = None
            self.errorLst.append("Exception when defining guns:")
            self.errorLst.append(str(e))

        if self.gun is not None:
            try:
                self.tableData = self.gun.integrate()
            except Exception as e:
                self.errorLst.append("Exception integrating:")
                self.errorLst.append(str(e))

        self.tv.delete(*self.tv.get_children())
        self.tableData = dot_aligned(self.tableData)

        for row in self.tableData:
            self.tv.insert("", "end", values=row, tags=(row[0],))

        self.propBanner.set("")

        self.updateError()

    def updateError(self):
        self.errorText.delete("1.0", "end")
        if self.geomError:
            self.errorLst.append("Invalid geometry")

        for line in self.errorLst:
            self.errorText.insert("end", line + "\n")
        self.errorLst = []

    def addSpecFrm(self, parent):
        specFrm = LabelFrame(parent, text="Design Summary")
        specFrm.grid(row=0, column=0, rowspan=2, sticky="nsew")

        propSpec = Label(specFrm, textvariable=self.propBanner).grid(
            row=0, column=0
        )

    def addParFrm(self, parent):
        parFrm = LabelFrame(parent, text="Parameters")
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

        Label(parFrm, text="Propellant").grid(row=4, column=0)
        # Create Dropdown menu
        self.dropProp = OptionMenu(parFrm, self.clickedProp, *self.propOptions)
        self.dropProp.grid(row=4, column=1, columnspan=2, sticky="nsew")
        self.dropProp.configure(width=25)
        Label(parFrm, text="Grain Shape").grid(row=5, column=0)
        # Create Dropdown menu
        self.dropGeom = OptionMenu(parFrm, self.clickedGeom, *self.geometries)
        self.dropGeom.grid(row=5, column=1, columnspan=2, sticky="nsew")
        self.dropGeom.configure(width=25)
        i += 2

        parFrm.rowconfigure(4, weight=0)
        parFrm.rowconfigure(5, weight=0)

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
            parFrm, i, "Shot Start Pressure", "MPa", "30", validationNN
        )

    def addOpsFrm(self, parent):
        opFrm = LabelFrame(parent, text="Operation")
        opFrm.grid(row=1, column=2, sticky="nsew")

        opFrm.columnconfigure(0, weight=1)
        opFrm.columnconfigure(1, weight=1)
        opFrm.columnconfigure(1, weight=1)
        opFrm.rowconfigure(2, weight=1)

        validationNN = parent.register(validateNN)
        i = 0
        self.webR, webRw, i = self.add2Input(
            opFrm, i, 1, "W.Th. ", "2.0", validationNN, 10
        )
        self.perR, perRw, i = self.add2Input(
            opFrm, i, 1, "P.Dia. ", "1.0", validationNN, 10
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

        Checkbutton(
            opFrm,
            text="Lock Geometry",
            variable=self.geoLocked,
            onvalue=1,
            offvalue=0,
            command=configEntry,
        ).grid(row=0, column=0, rowspan=2, sticky="nsew")

        Button(opFrm, text="Calculate", command=self.calculate).grid(
            row=2, column=0, columnspan=3, sticky="nsew"
        )

    def addErrFrm(self, parent):
        errorFrm = LabelFrame(parent, text="Errors")
        errorFrm.grid(row=1, column=1, sticky="nsew")
        errorFrm.columnconfigure(0, weight=1)
        errorFrm.rowconfigure(0, weight=1)

        errScroll = Scrollbar(errorFrm, orient="vertical")
        errScroll.grid(row=0, column=1, sticky="nsew")
        self.errorText = Text(
            errorFrm,
            yscrollcommand=errScroll.set,
            wrap=WORD,
            height=0,
            width=60,
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
        tblFrm = LabelFrame(parent, text="Result Table")
        tblFrm.grid(row=0, column=1, sticky="nsew")
        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(0, weight=1)

        self.tv = ttk.Treeview(tblFrm, selectmode="browse")
        self.tv.grid(row=0, column=0, sticky="nsew")

        self.tv["columns"] = columnList
        self.tv["show"] = "headings"
        self.tv.tag_configure("PEAK PRESSURE", foreground="orange")

        for column in columnList:  # foreach column
            self.tv.heading(
                column, text=column
            )  # let the column heading = column name
            self.tv.column(column, stretch=1, width=0)

        vertscroll = Scrollbar(tblFrm, orient="vertical")  # create a scrollbar
        vertscroll.configure(command=self.tv.yview)  # make it vertical
        vertscroll.grid(row=0, column=1, sticky="nsew")

        horzscroll = Scrollbar(tblFrm, orient="horizontal")
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
        entryWidth=0,
    ):
        Label(parent, text=labelText).grid(
            row=rowIndex, column=colIndex, sticky="nsew"
        )
        parent.rowconfigure(rowIndex, weight=0)
        e = StringVar(parent)
        e.set(default)
        en = Entry(
            parent,
            textvariable=e,
            validate="key",
            validatecommand=(validation, "%P"),
            width=entryWidth,
        )
        en.default = default
        en.grid(row=rowIndex, column=colIndex + 1, sticky="nsew")
        en.bind("<FocusOut>", formatInput)
        return e, en, rowIndex + 1

    def add3Input(
        self,
        parent,
        rowIndex,
        labelText,
        unitText,
        default="0.0",
        validation=None,
        entryWidth=0,
    ):
        Label(parent, text=labelText).grid(
            row=rowIndex, column=0, sticky="nsew"
        )
        parent.rowconfigure(rowIndex, weight=0)
        e = StringVar(parent)
        e.set(default)
        en = Entry(
            parent,
            textvariable=e,
            validate="key",
            validatecommand=(validation, "%P"),
            width=entryWidth,
        )
        en.default = default
        en.grid(row=rowIndex, column=1, sticky="nsew")
        en.bind("<FocusOut>", formatInput)
        Label(parent, text=unitText).grid(row=rowIndex, column=2, sticky="nsew")
        return e, en, rowIndex + 1

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
    root = Tk()
    ibPanel = IB(root)
    # allow parent window resizing

    # root.rowconfigure(0, weight=1)
    # root.columnconfigure(1, weight=1)

    root.mainloop()
