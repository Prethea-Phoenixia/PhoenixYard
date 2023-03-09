import tkinter

from tkinter import *
from tkinter import ttk


def main():
    root = Tk()

    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    root.title("Phoenix Yard")

    pltCvs = Canvas(root, width=640, height=480, bg="white", highlightthickness=0)
    pltCvs.grid(row=0, rowspan=2, column=1, sticky=N + S + E + W)

    # --------------------- initialize variable -----------------------------

    calmm = 100  # 100mm
    tOrx = False  # false -> t plot, true -> x plot

    specFrm = ttk.Frame(root, padding=10)
    specFrm.grid(row=0, column=0, rowspan=2)

    ttk.Label(specFrm, text="test").grid(row=0, column=0)

    # ----------------------- input panel ------------------------------

    parFrm = ttk.Frame(root, padding=10)
    parFrm.grid(row=0, column=2)

    ttk.Label(parFrm, text="Caliber").grid(row=0, column=0)
    ttk.Entry(parFrm, textvariable=calmm).grid(row=0, column=1)
    ttk.Label(parFrm, text="mm").grid(row=0, column=2)

    def redrawGraph(event):

        pltCvs.delete("all")

        pxWidth = pltCvs.winfo_width()
        pxHeight = pltCvs.winfo_height()
        print(pxWidth, pxHeight)
        # x1,y1,x2,y2
        # horizontal axis
        pltCvs.create_line(0.1 * pxWidth, 0.9 * pxHeight, 0.9 * pxWidth, 0.9 * pxHeight)
        # vertical axis
        pltCvs.create_line(0.1 * pxWidth, 0.1 * pxHeight, 0.1 * pxWidth, 0.9 * pxHeight)

        hzLabel = pltCvs.create_text((0.91 * pxWidth, 0.9 * pxHeight), anchor="w")
        pltCvs.itemconfig(hzLabel, text="t" if tOrx else "x")

    def calculate(event=None):
        # force an immediate redraw after calculation
        redrawGraph(None)

    def changePlotType():
        nonlocal tOrx
        tOrx = not tOrx
        redrawGraph(None)

    # ----------------------- operation panel -----------------------

    opFrm = ttk.Frame(root, padding=10)
    opFrm.grid(row=1, column=2)

    ttk.Button(opFrm, text="x<>t", command=changePlotType).grid(row=0, column=0)
    ttk.Button(opFrm, text="Calculate", command=calculate).grid(row=0, column=1)

    root.bind("<Configure>", redrawGraph)
    root.bind("<space>", calculate)

    root.mainloop()


if __name__ == "__main__":
    main()
