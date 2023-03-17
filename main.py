from tkinter import *
from tkinter import ttk

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


class IB(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        parent.title("Phoenix's Internal Ballistics Solver")
        # ------------------------ requisite file IO ---------------------------------

        compositions = GrainComp.readFile("propellants.csv")
        gemoetries = {i.desc: i for i in Geometry}

        propOptions = tuple(compositions.keys())
        geoOptions = tuple(gemoetries.keys())

        self.prop = None
        self.gun = None
        self.tableData = []
        self.errorLst = []
        self.geomError = False

        columnList = [
            "Time/s",
            "Travel/m",
            "Burnup/1",
            "Velocity/ms^-1",
            "Pressure/Pa",
        ]

        # ------------------------ initialize the GUI -------------------------------

        # datatype of menu text
        clickedProp = StringVar()
        # initial menu text
        clickedProp.set(propOptions[0])

        clickedGeom = StringVar()
        clickedGeom.set(geoOptions[0])
        propBanner = StringVar()
        geoLocked = IntVar()
        errorBanner = StringVar()

        def updateError():
            if self.geomError:
                self.errorLst.append("Invalid geometry")
            errorBanner.set("\n".join(self.errorLst))
            self.errorLst = []

        # ----------------------- specs panel ----------------------------------

        specFrm = LabelFrame(parent, text="Design Summary")
        specFrm.grid(row=0, column=0, sticky="nsew")

        propSpec = Label(specFrm, textvariable=propBanner).grid(row=0, column=0)

        errorFrm = LabelFrame(parent, text="Errors")
        errorFrm.grid(row=1, column=0, sticky="nsew")

        errorLbl = Label(
            errorFrm, textvariable=errorBanner, anchor="w", justify="left"
        ).grid(row=0, column=0)

        # ------------------------ initialize the table frame --------------------

        tblFrm = LabelFrame(parent, text="Result Table")
        tblFrm.grid(row=0, rowspan=2, column=1, sticky="nsew")
        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(0, weight=1)

        tv = ttk.Treeview(tblFrm, selectmode="browse")
        tv.grid(row=0, column=0, sticky="nsew")

        tv["columns"] = columnList
        tv["show"] = "headings"

        for column in columnList:  # foreach column
            tv.heading(
                column, text=column
            )  # let the column heading = column name
            tv.column(
                column, width=100, stretch=1
            )  # set the columns size to 50px

        # tv.place(relheight=1, relwidth=1)
        treescroll = Scrollbar(tblFrm, orient="vertical")  # create a scrollbar
        treescroll.configure(command=tv.yview)  # make it vertical
        tv.configure(
            yscrollcommand=treescroll.set
        )  # assign the scrollbar to the Treeview Widget

        # ----------------------- input panel ------------------------------
        treescroll.grid(row=0, column=1, sticky="nsew")
        parFrm = LabelFrame(
            parent,
            text="Parameters",
        )
        parFrm.grid(row=0, column=2, sticky="nsew")

        # validation
        validationNN = parent.register(validateNN)

        def addParInput(parent, rowIndex, labelText, unitText, default="0.0"):
            Label(parent, text=labelText).grid(row=rowIndex, column=0)
            e = StringVar(parent)
            e.set(default)
            en = Entry(
                parent,
                textvariable=e,
                validate="key",
                validatecommand=(validationNN, "%P"),
            )
            en.default = default
            en.grid(row=rowIndex, column=1)
            en.bind("<FocusOut>", formatInput)
            Label(parent, text=unitText).grid(row=rowIndex, column=2)
            return e, en, rowIndex + 1

        i = 0
        calmm, _, i = addParInput(parFrm, i, "Caliber", "mm")
        tblmm, _, i = addParInput(parFrm, i, "Tube Length", "mm")
        shtkg, _, i = addParInput(parFrm, i, "Shot Mass", "kg")
        chgkg, _, i = addParInput(parFrm, i, "Charge Mass", "kg")

        Label(parFrm, text="Propellant").grid(row=4, column=0)
        # Create Dropdown menu
        dropProp = OptionMenu(parFrm, clickedProp, *propOptions).grid(
            row=4, column=1, columnspan=2, sticky="nsew"
        )
        Label(parFrm, text="Grain Shape").grid(row=5, column=0)
        # Create Dropdown menu
        dropGeom = OptionMenu(parFrm, clickedGeom, *gemoetries).grid(
            row=5, column=1, columnspan=2, sticky="nsew"
        )
        i += 2

        webmm, webw, i = addParInput(parFrm, i, "Web Thickness", "mm")
        permm, perw, i = addParInput(parFrm, i, "Perf Diameter", "mm")
        grlR, grlRw, i = addParInput(parFrm, i, "Grain L/D", "", "2.25")

        ldf, _, i = addParInput(parFrm, i, "Load Factor", "%", "50.0")
        stpMPa, _, i = addParInput(parFrm, i, "Shot Start Pressure", "MPa")

        def calculate(event=None):
            # force an immediate redraw after calculation
            compo = compositions[clickedProp.get()]
            geom = gemoetries[clickedGeom.get()]

            self.tableData = []

            try:
                self.prop = Propellant(
                    compo,
                    geom,
                    float(webmm.get()) / 1000,
                    float(permm.get()) / 1000,
                    float(grlR.get()),
                )
            except Exception as e:
                print(e)
                self.prop = None
                self.errorLst.append(
                    "Exception when defining propellant:\n{:}".format(e)
                )

            try:
                chamberVolume = (
                    float(chgkg.get())
                    / self.prop.rho_p
                    / self.prop.maxLF
                    / float(ldf.get())
                    * 100
                )
                self.gun = Gun(
                    float(calmm.get()) / 1000,
                    float(shtkg.get()),
                    self.prop,
                    float(chgkg.get()),
                    chamberVolume,
                    float(stpMPa.get()) * 1e6,
                    float(tblmm.get()) / 1000,
                )
                self.tableData = self.gun.integrate()
            except Exception as e:
                print(e)
                self.gun = None
                self.errorLst.append(
                    "Exception when defining guns:\n{:}".format(e)
                )

            tv.delete(*tv.get_children())

            def dot_align(v):
                decimalLeft = 12
                decimalRight = 6
                dotIndex = "{:,}".format(v).find(".")
                spacePad = max(decimalLeft - dotIndex, 0)
                return (
                    " " * spacePad
                    + "{:,}".format(v)[: dotIndex + 4]
                    + " "
                    + format(v)[dotIndex + 4 : dotIndex + 7]
                )

            # foreach row of pokemon data insert the row into the treeview.
            for row in self.tableData:
                tv.insert("", "end", values=tuple(dot_align(v) for v in row))

            propBanner.set("")

            updateError()

        # ----------------------- operation panel -----------------------

        opFrm = LabelFrame(parent, text="Operation")
        opFrm.grid(row=1, column=2, sticky="nsew")

        opFrm.columnconfigure(0, weight=1)
        opFrm.rowconfigure(4, weight=1)

        def addRatInput(parent, rowIndex, labelText, default="1.0"):
            Label(parent, text=labelText).grid(row=rowIndex, column=0)
            e = StringVar(parent)
            e.set(default)
            en = Entry(
                parent,
                textvariable=e,
                validate="key",
                validatecommand=(validationNN, "%P"),
            )
            en.default = default
            en.grid(row=rowIndex, column=1)
            en.bind("<FocusOut>", formatInput)
            return e, en, rowIndex + 1

        i = 1
        webR, webRw, i = addRatInput(opFrm, i, "W.Th. ", "2.0")
        perR, perRw, i = addRatInput(opFrm, i, "P.Dia. ", "1.0")

        ratioEntrys = (webRw, perRw)
        directEntrys = (perw,)

        def configEntry():
            if geoLocked.get() == 0:
                for en in ratioEntrys:
                    en.config(state="disabled")
                for en in directEntrys:
                    en.config(state="normal")
            else:
                for en in ratioEntrys:
                    en.config(state="normal")
                for en in directEntrys:
                    en.config(state="disabled")

        configEntry()

        def webcallback(var, index, mode):
            if any(
                (
                    webmm.get() == "",
                    webR.get() == "",
                    perR.get() == "",
                )
            ):
                pass
            else:
                if geoLocked.get() == 1:
                    try:
                        self.geomError = False
                        permm.set(
                            float(webmm.get())
                            / float(webR.get())
                            * float(perR.get())
                        )
                    except ZeroDivisionError:
                        self.geomError = True

        webmm.trace_add("write", webcallback)
        webR.trace_add("write", webcallback)
        perR.trace_add("write", webcallback)
        grlR.trace_add("write", webcallback)

        Checkbutton(
            opFrm,
            text="Lock Geometry",
            variable=geoLocked,
            onvalue=1,
            offvalue=0,
            command=configEntry,
        ).grid(row=0, column=0, columnspan=2, sticky="nsew")

        Button(opFrm, text="Calculate", command=calculate).grid(
            row=4, column=0, columnspan=2, sticky="nsew"
        )

        parent.bind("<space>", calculate)


if __name__ == "__main__":
    root = Tk()
    ibPanel = IB(root)
    # allow parent window resizing
    root.columnconfigure(1, weight=1)
    root.rowconfigure(1, weight=1)

    root.mainloop()
