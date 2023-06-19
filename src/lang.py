ENGLISH = {
    "themeLabel": "Theme",
    "lightLabel": "Light",
    "darkLabel": "Dark",
    "debugLabel": "Debug Info",
    "enableLabel": "Enabled",
    "vTgtLabel": "V. Tgt.",
    "pTgtLabel": "P. Tgt.",
    "minWebLabel": "Min. W.",
    "maxLgLabel": "Max. L.",
    "calLabel": "Caliber",
    "tblLabel": "Tube Length",
    "shtLabel": "Shot Mass",
    "chgLabel": "Charge Mass",
    "ldfLabel": "Load Factor",
    "clrLabel": "Chamber L.R.",
    "dgcLabel": "Drag coe.",
    "stpLabel": "Start P.",
    "nozzExpLabel": "Nozz. Exp.",
    "nozzEffLabel": "Nozz. Eff.",
    "tblFrmLabel": "Result Table",
    "plotFrmLabel": "Plot",
    "errFrmLabel": "Exceptions",
    "parFrmLabel": "Parameters",
    "specFrmLabel": "Design Summary",
    "opFrmLabel": "Operations",
    "consFrmLabel": "Constraints",
    "consButton": "Constrain Design",
    "minTVButton": "Minimize Tube Volume",
    "sampleFrmLabel": "Sampling",
    "pltOptnFrm": "Plot Option",
    "lxLabel": "Length Ratio",
    "vaLabel": "Asymptotic Vel.",
    "pPLabel": "Peak Pressure",
    "teffLabel": "Thermal Eff.",
    "beffLabel": "Ballistic Eff.",
    "cvLabel": "Chamber Volume",
    "ldLabel": "Loading Density",
    "propFrmLabel": "Propellant",
    "grainFrmLabel": "Grain Geometry",
    "diamLabel": "Dia.",
    "widtLabel": "Width",
    "ltwLabel": "Len. / W.",
    "htwLabel": "H. / W.",
    "ltdLabel": "Len. / Dia.",
    "athLabel": "Arc Thickness",
    "pdtalLabel": "P.D. / A.Th.",
    "plotAvgP": "Length Avg. P.",
    "plotBaseP": "Shot Base P.",
    "plotBreechNozzleP": "Breech / Nozzle P.",
    "plotStagP": "Stagnation P.",
    "plotVel": "Shot Vel.",
    "plotNozzleV": "Nozzle Throat Vel.",
    "plotBurnup": "Burnup",
    "plotEta": "Escape",
    "stepsLabel": "Steps",
    "calcLabel": "CALCULATE",
    "CONVENTIONAL": "Conventional Gun",
    "RECOILESS": "Recoiless Gun",
    "DOMAIN_TIME": "Time",
    "DOMAIN_LENG": "Length",
    "TvDesc": "Adb.Temp",
    "isochorDesc": "(Isochoric)",
    "densityDesc": " Density",
    "vacISPDesc": "Isp(Vac)",
    "atmISPDesc": "Isp(Atm)",
    "pRatioDesc": "(Pc:Pa=50)",
    "brDesc": "Burnrate",
    "chgText": " ".join(
        (
            "Mass of propellant charge to be used. Lagrangian approximation",
            "is decently held up until a charge/shot mass ratio of 1. Higher",
            "ratios tends to excite pressure waves inducing large pressure",
            "spikes, requiring 1-D models to resolve.",
        )
    ),
    "vinfText": " ".join(
        (
            "Velocity the shot would achieve if",
            "the barrel is extended to infinite length.",
        )
    ),
    "teffText": " ".join(
        (
            "Thermal efficiency of the gun system, i.e.",
            "the amount of work done to both gas and",
            "projectile over potential from propellant",
        )
    ),
    "beffText": " ".join(
        (
            "Ballistic efficiency of the gun system, i.e.",
            "the amount of work done the projectile over",
            "chemical potential energy of propellant.",
        )
    ),
    "specsText": " ".join(
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
    ),
    "geomPlotTxt": " ".join(
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
    ),
    "ldftext": " ".join(
        (
            "Percentage of chamber volume filled by",
            "the outlines of the grain, also known",
            "as 'packing density'. Value of 0-100 % ",
            "(not inclusive) are supported, although",
            "it is highly unlikely that values greater",
            "than 60% could be achieved in practice,",
            "due to packing behaviour when loading the",
            "cartridge.\n",
            "A high value is also undesirable for",
            "causing excessive peak pressure as well",
            "as moving the pressure spike closer to",
            "the breech.\n",
            "Only internal voids are accounted for,",
            "thus this is not simply the usually quoted",
            "'load fraction' which is by weight.",
        )
    ),
    "arcText": " ".join(
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
    ),
    "pDiaRText": " ".join(
        (
            "Specify diameter of perforation over arc width.\n",
            "Perforations are formed by protrusions in the",
            "copper casting die. Standard multi-perf grain tends",
            "to come in the range of 0.5-1, whereas longer, single",
            "perf, or tubular grains come to around 1.33 for the same.",
        )
    ),
    "diaText": " ".join(("Specify the diameter of the propellant grain.",)),
    "perfLRtext": " ".join(
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
    ),
    "cylLRtext": " ".join(
        (
            "Specify length to diameter ratio of the grain.",
            "Cylindrical or tubular propellant can be made rather long",
        )
    ),
    "rodRtext": " ".join(
        (
            "Specify the length to width ratio of propellant rod or flake.",
            "Can be quite long for tape like propellant.",
        )
    ),
    "widthText": " ".join(
        (
            "Specify the width of propellant rod or flake.\n",
            "This value is used to scale all other dimension",
            "for rod or flake like propellants. Does not have",
            "to be the smallest dimension.",
        )
    ),
    "heightRtext": " ".join(
        ("Specify the height to width ratio of propellant rod or flake.",)
    ),
    "tolText": " ".join(
        (
            "The maximum relative error, ε, that is allowed in the integrands",
            "for each component. Some components may have significantly less",
            "error than specified here, shown under each entry in the table.",
        )
    ),
    "stpText": " ".join(
        (
            "Peak pressure that initially resists",
            "the movement of shot. This is made up of",
            "the rifling engraving the driving band or",
            "shell body, the static friction of the",
            "barrel, and the cartridge case crimping onto",
            "the shot (greatly varying between 0.25 - 100MPa",
            "depending on caliber and desired RoF)\n",
            "For recoiless weapon, the nozzle open pressure",
            "is also taken to be this value, for simplicity,",
            "and such that the balance of momentum is conserved",
            "throughout.",
        )
    ),
    "clrtext": " ".join(
        (
            "Chamber length ratio is defined as the ratio of the actual chamber cross section",
            "to that of the barrel, or put another way, the length a chamber would have if its",
            "in line with the barrel to that of the actual chamber length. This parameters is",
            "used to correct for chamberage effect, which is most significant at the start of",
            "IB cycle. Currently the chamberage correction for conventional guns employ a length",
            "averaged approach. Recoiless guns are (currenlty) not corrected for chamberage effect.",
        )
    ),
    "dgctext": " ".join(
        (
            "Drag coefficient, or the pressure induced by barrel",
            "friction divided by the shot base pressure. Currently",
            "2%-7% is reported for rifled weapons, with smaller",
            "calibers on the higher end and larger caliber shots",
            "with driving band on the lower end.",
        )
    ),
    "sampTxt": " ".join(
        (
            "Samples are taken equidistantly along specified domain.",
            "Sampling is done after the system has been solved and thus",
            "does not influence the accuracy of characteristc points.",
            "This can be used to sanity-check and validate calculations",
            "made to an accuracy specification.",
        )
    ),
    "pTgtTxt": " ".join(
        (
            "Pressure design target: length averaged, Lagrangian",
            "chamber pressure.",
        )
    ),
    "useConsTxt": " ".join(
        (
            "Constrain the design to specified muzzle velocity and peak pressure",
            "by controlling the web thickness and tube length. Currently solved",
            "solution is correct for chamber length ratio of 1.0x.",
        )
    ),
    "optLFTxt": " ".join(
        (
            "Find the optimum load fraction for this charge mass, minimizing",
            "tube volume (including chamber), in addition to the above constraints.",
            "The minimum volume solution is not necessarily appropriate for indirect",
            "fire guns, as the proximity of burnout point with the muzzle end tends",
            "to exacerbate velocity dispersion due to uneven burning, and the effect",
            "of ambient temperature on propellant force.",
        )
    ),
    "calLxTxt": " ".join(
        (
            "Tube length, commonly expressed in literature as L/xx.",
            "Top is measured from shot start, bottom is measured from breechface.",
        )
    ),
    "pMaxTxt": " ".join(
        (
            "Peak pressure, when measured from the breech, and when measured at",
            "the base of projectile. Breech pressure is around 1.0-1.2x of the",
            "result of copper crusher test, with the factor usually taken to be",
            "1.12x. These values are tabulated at the peak of space-averaged",
            "pressure, which does not co-occur when the effect of chamberage is",
            "considered (Chamber Expansion not equal to 1.0x).",
        )
    ),
    "nozzExpTxt": " ".join(
        (
            "Area expansion ratio of the recoiless gun's rear nozzle, or the ratio",
            "between the cross section area of the nozzle end, and the area of the",
            "throat. This is used to size the throat, or opening, with larger nozzle",
            "expansion ratio",
            "resulting in more efficient nozzle, and less gas leakage, however with",
            "quickly diminishing returns. Usually chosen to be around 4.",
        )
    ),
    "nozzEffTxt": " ".join(
        (
            "Efficiency of the nozzle, accouting for the less efficient geometry used",
            "in real nozzles to simplify production, where a short nozzle of the straight",
            "walled type can be reasonably estimated to be within a few percent of unity,",
            "but for the effect of high pressure gas that is usually not corrected for in",
            "rocket theories. Therefore a more conservative estimate of 92% is adopted",
        )
    ),
    "columnList": [
        "Event",
        "Time",
        "Travel",
        "Burnup",
        "Velocity",
        "Avg. Pressure",
        "Avg. Temp.",
        "Outflow Fract.",
    ],
    "SEVEN_PERF_CYLINDER": "7 Perf. Cylinder",
    "SEVEN_PERF_ROSETTE": "7 Perf. Rosette Prism",
    "FOURTEEN_PERF_ROSETTE": "14 Perf. Rosette Prism",
    "NINETEEN_PERF_ROSETTE": "19 Perf. Rosette Prism",
    "NINETEEN_PERF_CYLINDER": "19 Perf. Cylinder",
    "NINETEEN_PERF_HEXAGON": "19 Perf. Hexagonal Prism",
    "NINETEEN_PERF_ROUNDED_HEXAGON": "19 Perf. Rounded Hex. Prism",
    "SPHERE": "Sphere",
    "ROD": "Strip / Flake",
    "CYLINDER": "Cylinder",
    "TUBE": "Tube",
    "figTimeDomain": "Time - ms",
    "figLengDomain": "Travel - m",
}
CHINESE = {
    "themeLabel": "主题",
    "lightLabel": "浅色",
    "darkLabel": "深色",
    "debugLabel": "调试信息",
    "enableLabel": "启用",
    "vTgtLabel": "设计弹速",
    "pTgtLabel": "最大膛压",
    "minWebLabel": "最小弧厚",
    "maxLgLabel": "最长身管",
    "calLabel": "口径",
    "tblLabel": "身管长",
    "shtLabel": "弹重",
    "chgLabel": "药量",
    "ldfLabel": "装药系数",
    "clrLabel": "药室扩大系数",
    "dgcLabel": "阻力系数",
    "stpLabel": "挤进压强",
    "nozzExpLabel": "喷管面积比",
    "nozzEffLabel": "喷管效率",
    "tblFrmLabel": "计算结果打表",
    "plotFrmLabel": "计算结果作图",
    "errFrmLabel": "程序运行例外",
    "parFrmLabel": "内弹道参量",
    "specFrmLabel": "系统性能指标",
    "opFrmLabel": "程序运算控制",
    "consFrmLabel": "反算设置",
    "consButton": "反算弧厚身管长",
    "minTVButton": "求最小膛容解",
    "sampleFrmLabel": "采样控制",
    "pltOptnFrm": "作图控制",
    "lxLabel": "倍径",
    "vaLabel": "极限弹速",
    "pPLabel": "最大压强",
    "teffLabel": "热机效率",
    "beffLabel": "弹道效率",
    "cvLabel": "药室容积",
    "ldLabel": "装填密度",
    "propFrmLabel": "火药",
    "grainFrmLabel": "药粒形状参数",
    "diamLabel": "直径",
    "widtLabel": "宽度",
    "ltwLabel": "长宽比",
    "htwLabel": "高宽比",
    "ltdLabel": "药粒长径比",
    "athLabel": "弧厚",
    "pdtalLabel": "孔径弧厚比",
    "plotAvgP": "空间平均压强",
    "plotBaseP": "弹底压强",
    "plotBreechNozzleP": "膛底压强",
    "plotStagP": "滞止点压强",
    "plotVel": "弹速",
    "plotNozzleV": "喉口流出速度",
    "plotBurnup": "火药燃尽比例",
    "plotEta": "燃气流出比例",
    "stepsLabel": "采样点",
    "calcLabel": "计算",
    "CONVENTIONAL": "常规火炮",
    "RECOILESS": "无后座炮",
    "DOMAIN_TIME": "时间",
    "DOMAIN_LENG": "空间",
    "chgText": "".join(
        (
            "拉格朗日假设是简化弹道学问题的经典假设。该假设下膛内密度始终沿炮膛均匀分布，火药燃速",
            "以膛内空间平均压强计算。在拉格朗日假设下计算的结果，可以较为准确地描述装药量小于",
            "弹重的情况。装药量较大时，一方面膛内压力波传递与反射对于膛内气体影响较大，破坏了",
            "压强分布假设；另一方面火药膛内运动较为明显，不能简单地认为在平均压力点燃烧。",
            "研究这些问题涉及到对膛内气体两相流动问题的精确解。",
        )
    ),
    "vinfText": "".join(("该装药条件下的极限速度（vj），即炮管长度趋于无穷时弹速的极限。")),
    "calLxTxt": "".join(
        (
            "炮管倍径，即炮管长度与口径的比例。西方文献中常用 L/xx表示。\n",
            "上：身管倍径，至弹底起始位置；下：身管（含药室）倍径，至膛底",
        )
    ),
    "geomPlotTxt": "".join(
        (
            "相对燃烧表面σ=S/S1，关于相对厚度Z=e/e1的函数。曲线递增代表增面燃烧，曲线递减代表减面燃烧。",
            "简单形状火药燃烧全过程Z从0到1，为减面燃烧；多孔火药的燃烧分裂前Z从0到1，有条件地增面燃烧，",
            "分裂后Z从1到Zb，由于剩余体积较小，分裂出棒状体几何形状复杂，以等体积圆柱处理，燃烧体现减面性。",
            "增面燃烧有利于减缓压强上升，降低最大压强，提高压强充满率，对于武器设计十分有利。",
        )
    ),
    "teffText": "".join(
        (
            "对于火药在膛内完全燃烧的情况，火炮热机效率(γg')是出膛瞬间火药气体温度（Tg），",
            "与火药恒容绝热火焰温度（T1）的函数，γg'=1-T/T1。热机功率联系了出膛速度（vg）与极限速度（vj），",
            "即vg=sqrt(γg')*vj。当Tg趋近于0时，出膛速度趋近于极限速度（vj）。",
        )
    ),
    "beffText": "".join(
        (
            "火炮弹道效率（γg）定义为热机效率除以次要功系数（φ），即γg=γg'/φ，",
            "也是炮弹出膛动能占火药做功潜能的比例，是衡量火炮系统火药利用效率的主要指标之一。",
        )
    ),
    "columnList": [
        "特征点",
        "历时",
        "行程",
        "火药燃去比例",
        "弹速",
        "空间平均压强",
        "平均燃气温度",
        "燃气流出比例",
    ],
    "SEVEN_PERF_CYLINDER": "圆柱形七孔柱",
    "SEVEN_PERF_ROSETTE": "花边形七孔柱",
    "FOURTEEN_PERF_ROSETTE": "花边十四孔柱",
    "NINETEEN_PERF_ROSETTE": "花边十九孔柱",
    "NINETEEN_PERF_CYLINDER": "圆柱形十九孔柱",
    "NINETEEN_PERF_HEXAGON": "正六边形十九孔柱",
    "NINETEEN_PERF_ROUNDED_HEXAGON": "圆角六边形十九孔柱",
    "SPHERE": "球体",
    "ROD": "立方体状",
    "CYLINDER": "圆柱",
    "TUBE": "管状",
    "TvDesc": "绝热燃烧温度",
    "isochorDesc": "（定容）",
    "densityDesc": "　　　　密度",
    "vacISPDesc": "　　真空比冲",
    "atmISPDesc": "　　大气比冲",
    "pRatioDesc": "　（压强比50）",
    "brDesc": "燃速",
    "figTimeDomain": "时间 / ms",
    "figLengDomain": "行程 / m",
}
STRING = {"English": ENGLISH, "中文": CHINESE}
