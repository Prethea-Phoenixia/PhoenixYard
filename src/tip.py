from tkinter import *
import tkinter.font as tkFont
from math import ceil


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

useConsTxt = " ".join(
    (
        "Constrain the design to specified muzzle velocity and peak pressure",
        "by varying the web thickness and tube length. Currently solved",
        "solution is correct for chamber length ratio of 1.0x.",
    )
)

optLFTxt = " ".join(
    (
        "Optimize load fraction such that constrained design will result in",
        "the solution that minimize tube volume (including chamber).",
    )
)

calLxTxt = " ".join(
    (
        "Tube length, commonly expressed in literature as L/xx.",
        "Top is measured from shot start, bottom is measured from breechface.",
    )
)


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
