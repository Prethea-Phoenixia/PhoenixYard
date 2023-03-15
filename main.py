from tkinter import *
from tkinter import ttk

from gun import *


def validate(inp):
    if inp == "":
        return True
    try:
        float(inp)
    except ValueError:
        return False
    return True


def formatInput(event):
    v = event.widget.get()
    if v == "":
        event.widget.insert(0, "0.0")
    else:
        event.widget.delete(0, END)
        event.widget.insert(0, float(v))


def popup(issue):
    popupWindow = Toplevel()
    popupWindow.wm_title("Issue")
    popupWindow.overrideredirect(True)
    Label(popupWindow, text=issue).grid(row=0, column=0)
    # Button(popupWindow, text="Okay", command=popupWindow.destroy).grid(row=1, column=0)
    popupWindow.after(2000, lambda *args: popupWindow.destroy())


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

        columnList = ["Time/s", "Travel/m", "Burnup/1", "Velocity/ms^-1", "Pressure/Pa"]

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
            tv.heading(column, text=column)  # let the column heading = column name
            tv.column(column, width=100)  # set the columns size to 50px

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
        validation = parent.register(validate)

        def addParInput(parent, rowIndex, labelText, unitText):
            Label(parent, text=labelText).grid(row=rowIndex, column=0)
            e = StringVar(parent)
            e.set("0.0")
            en = Entry(
                parent,
                textvariable=e,
                validate="key",
                validatecommand=(validation, "%P"),
            )
            en.grid(row=rowIndex, column=1)
            en.bind("<FocusOut>", formatInput)
            Label(parent, text=unitText).grid(row=rowIndex, column=2)
            return e, rowIndex + 1

        i = 0
        calmm, i = addParInput(parFrm, i, "Caliber", "mm")
        tblmm, i = addParInput(parFrm, i, "Tube Length", "mm")
        shtkg, i = addParInput(parFrm, i, "Shot Mass", "kg")
        chgkg, i = addParInput(parFrm, i, "Charge Mass", "kg")

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

        permm, i = addParInput(parFrm, i, "Perf Diameter", "mm")
        webmm, i = addParInput(parFrm, i, "Web Thickness", "mm")
        grlmm, i = addParInput(parFrm, i, "Grain Length", "mm")

        ldf, i = addParInput(parFrm, i, "Load Factor", "%")
        stpMPa, i = addParInput(parFrm, i, "Shot Start Pressure", "MPa")
        """
        def Load_data():
            
        """

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
                    float(grlmm.get()) / 1000,
                )
            except:
                self.prop = None

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
                self.tableData = gun.integrate()
            except:
                self.gun = None

            errors = []
            if self.prop is None:
                errors.append("Unable to define: Propellant")

            if self.gun is None:
                errors.append("Unable to define: Gun")

            errorBanner.set("\n".join(errors))

            tv.delete(*tv.get_children())

            # foreach row of pokemon data insert the row into the treeview.
            for row in self.tableData:
                tv.insert("", "end", values=row)

            propBanner.set("")

        # ----------------------- operation panel -----------------------

        opFrm = LabelFrame(parent, text="Operation")
        opFrm.grid(row=1, column=2, sticky="nsew")

        opFrm.columnconfigure(0, weight=1)
        opFrm.rowconfigure(1, weight=1)

        Checkbutton(
            opFrm, text="Lock Geometry", variable=geoLocked, onvalue=1, offvalue=0
        ).grid(row=0, column=0, sticky="nsew")

        Button(opFrm, text="Calculate", command=calculate).grid(
            row=1, column=0, sticky="nsew"
        )

        parent.bind("<space>", calculate)


if __name__ == "__main__":
    root = Tk()
    ibPanel = IB(root)
    # allow parent window resizing
    root.columnconfigure(1, weight=1)
    root.rowconfigure(1, weight=1)

    root.mainloop()
