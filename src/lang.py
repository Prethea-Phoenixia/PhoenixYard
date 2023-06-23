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
    "cylLRtext": " ".join(("Specify length to diameter ratio of the grain.",)),
    "rodRtext": " ".join(
        ("Specify the length to width ratio of propellant rod or flake.",)
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
            "Chamber length ratio is defined as the ratio of the actual chamber",
            "cross section to that of the barrel, or put another way, the length",
            "a chamber would have if its in line with the barrel to that of the",
            "actual chamber length. This parameters is used to correct for",
            "chamberage effect, which is most significant at the start of IB cycle.",
            "Currently the chamberage correction for conventional guns employ a length"
            "averaged approach. Recoiless guns are (currenlty) not corrected for",
            "chamberage effect.",
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
    "calcButtonTxt": "Integrate system using RKF7(8) integrator",
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
    "figShotVel": "Shot Velocity\nm/s",
    "figBreech": "Breech Face",
    "figStagnation": "Stagnation",
    "figNozzleP": "Nozz. Throat P.",
    "figNozzleV": "Nozz. Throat Vel.",
    "figOutflow": "Outflow Frac.",
    "figAvgP": "Avg. Pressure\nMPa",
    "figShotBase": "Shot Base",
    "figTgtP": "P. Target",
    "figPsi": "Volume Burnup",
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
    "plotBurnup": "火药燃去比例",
    "plotEta": "燃气流出比例",
    "stepsLabel": "采样点",
    "calcLabel": "计算",
    "CONVENTIONAL": "常规火炮",
    "RECOILESS": "无后坐炮",
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
    "specsText": "".join(
        (
            "   单基火药一般由硝化纤维作为主要成分，加入少许稳定剂、塑化剂，提高化学稳定性，降低燃速爆温。",
            "其生产工艺成熟，成本可控，且具有良好的机械性能。\n",
            "   双基火药在硝化纤维为主的基础上，溶入相当比例的硝化甘油作为胶化剂，同时积极参与火药燃烧。",
            "在提高火药力以外还有改善燃烧产物化学平衡，减少发烟的作用。除硝化甘油外，亦有二甘醇二硝酸酯",
            "作为胶化剂的技术路线，可减少对硝化甘油生产所需的酒精的消耗。少量添加其他成分作稳定剂。",
            "双基火药更高的爆温也同时提高了燃气的侵蚀性，提高性能的同时，对系统寿命有消极影响。",
            "基于二甘醇二硝酸酯的双基火药挥发性更强，热稳定性更弱。\n",
            "   三基火药成分以硝基胍为主，硝化纤维与硝化甘油以一定比例溶解而成。硝基胍在燃烧中起到改善膛内",
            "欠氧状态，降低产物可燃物比例，起到在不提高爆温的前提下提高火药力的作用。与双基火药一样，",
            "同样可以采用其他胶化剂。三基火药通常机械性能欠佳，生产时需要低温缓慢挤出成型，膛内较高压强或",
            "冲击波条件下容易发生药粒碎裂的情况。此外，硝基胍合成流程较长，原材料成本较高，因此产量有限。",
            "因此，也有使用混合硝酸酯替代硝基胍，以降低造假，减少依赖的一些探索。\n",
            "   硝胺火药以环三亚甲基三硝胺（亦称黑索金）为主，混溶一定比例的硝基火药制成。常规硝胺炸药药力",
            "高，研磨成微米微粒后燃烧速率与火药敏感性可控。基于环四亚甲基四硝胺（亦称奥托金）的硝铵火药药力",
            "更高，但燃烧表现不稳定。",
        )
    ),
    "geomPlotTxt": "".join(
        (
            "相对燃烧表面σ=S/S1，关于相对厚度Z=e/e1的函数。曲线递增代表增面燃烧，曲线递减代表减面燃烧。",
            "简单形状火药燃烧全过程Z从0到1，为减面燃烧；多孔火药的燃烧分裂前Z从0到1，有条件地增面燃烧，",
            "分裂后Z从1到Zb，由于剩余体积较小，分裂出棒状体几何形状复杂，以等体积圆柱处理，燃烧体现减面性。",
            "增面燃烧有利于减缓压强上升，降低最大压强，提高压强充满率，改善膛内受力情况。",
        )
    ),
    "ldftext": "".join(
        (
            "装填系数，即药粒外轮廓对于药膛空间充满比例。对于无内孔的药粒，等于装填密度(Δ)与",
            "火药密度(ρp)的比例。对于有内孔的药粒，再乘上药粒内空间填充率。该系数反应了药室内",
            "药包，药筒，组合模块等装药结构、药粒真实堆叠对于空间的利用效率。过低的装填系数可能造",
            "成火药出现低压异常燃烧现象，严重情况下无法克服挤进压力将弹头推出炮管。过高的装填系数",
            "装填时可能需要用力挤压，使药粒变形碎裂；在击发时也难以保证传火效率，无法均匀点火，",
            "在大装药量下还容易诱发燃烧压力波。经验上，小口径武器由于药粒形状简单，装药量小，可以",
            "采用较大的装填系数。装药量大，采用多孔火药的大口径武器，装填系数普遍较低。",
        )
    ),
    "arcText": "".join(
        (
            "   设置火药药粒的弧厚。弧厚又称肉厚，对于简单火药，为平行表面间最短的距离；对于多孔火药，",
            "是相邻两孔内表面间最短距离、也是最外层孔的内表面与药粒外表面相平行的最短距离。在火药",
            "所有表面同时点火，以各向同速的平行层燃烧的模型下，肉厚是最有意义的关于药粒的几何参量。\n",
            "   计算结果表明，最大压强和燃尽点距离对于肉厚最敏感；对于管长于燃尽点的",
            "内弹道系统，弹速与肉厚的相关性较小。\n",
            "   多孔火药的加工精度一般在0.15mm-1um之间，因此通常弧厚不小于1mm，否则批次间差距过大，不利于保证",
            "一致性。简单形状火药相对而言可以加工，筛选成更小的几何尺寸。因此，在小口径武器或无坐力炮等需要较小尺",
            "寸药粒的情况下通常使用简单形状药粒。",
        )
    ),
    "pDiaRText": "".join(
        (
            "   设置火药药粒孔径对于弧厚的比例。\n",
            "   生产过程中，多孔火药一般是以粘稠悬浊液、胶体的形式，在压强下挤入黄铜模具翻模制造。",
            "火药中孔洞尺寸，形状由模具中的柱状突出物半径，截面形状，偏心率等决定。此外，火药具有一定的",
            "可塑性，因此拆模，包装，运输过程中的挤压形变，也会对于孔径产生一定影响。\n",
            "   标准的多孔火药中，该值一般取值在0.5到1之间。对于较长的柱状火药，该值可以取到1.33。",
        )
    ),
    "diaText": "".join(("设置火药药粒的直径。",)),
    "perfLRtext": "".join(
        (
            "设置火药药粒的长径比。多孔火药长径比一般在1.82-3.57之间。该比例越大，火药燃烧增面性越强。",
            "当取值过小，以至于药粒柱高小于弧厚时，药粒可能出现减面燃烧的现象。",
            "一般而言较长的药粒难以同时着火，需要每隔一定间隔开槽。",
            "运动性越差，发射时停留在在膛底，作用压强越高。",
        )
    ),
    "cylLRtext": "".join(("设置、管柱状火药的长径比",)),
    "rodRtext": "".join(("设置长方体火药长宽比。",)),
    "widthText": "".join(("设置长方体火药宽度。不需要一定是最小跨度边。",)),
    "heightRtext": "".join(("设置火药的高宽比例。",)),
    "tolText": "".join(
        (
            "积分器最大允许相对绝对误差，ε。龙格-库塔-菲尔伯格自适应积分器根据每一步各项的取值，导数，",
            "对常微分方程生成七阶，八阶两个估测。这两个估测的差值作为局部误差，再根据步长与积分域的比例",
            "外推，除以该点取值，估测全局相对误差。将该值与用户设定的误差大小作比较，按一定算法调整步长",
            "直到该条件得到满足。重复以上直到积分器达到停止条件，或遇到数值奇点。",
        )
    ),
    "stpText": "".join(
        (
            "   设置挤进压强。内弹道计算中常见的化简是认为在第一阶段，火药在弹后空间定容燃烧，直到克服挤进",
            "压强。此时，弹头瞬间完成挤进炮膛，与膛线结合，开始运动，火药燃烧进入由常微方程描述的第二阶段。",
            "因为该参量实际上反映的是弹道简化，因此物理意义并不十分明确。\n",
            "   历史上常用静压法测定挤进压强，即用压力机缓慢将弹头挤入膛线所需的最大静压。一般而言，对于大口径",
            "线膛炮10-20MPa，对小口径武器取30-40MPa。虽然较大的挤进压强造成了最大压强增加，不利于弹道系统；",
            "但更高的挤进压强有利于改善点火条件，提高点火一致性，因此即便是滑膛炮也需要一定的挤进压强。\n",
            "   实际挤进是一个运动的过程，弹头在挤进过程中，膛压可达到最大膛压的六、七成，挤进完成时火药燃烧，",
            "弹头速度都达到出膛的一成左右，耗时大致是火药燃烧时间的四分之一到三分之一。\n",
            "   对无坐力火炮，为满足严格无坐力条件，喷口打开压强与挤进取相同数值。",
        )
    ),
    "clrtext": "".join(
        (
            "药室扩大系数指的是药室截面积与炮管截面积的比例。膛压较大时，定装弹药过长的药室存在抽壳困难",
            "的情况，后膛过长也不利于火炮装填作业。\n",
            "药室缩进产生流速突变对于弹道内循环的起始阶段影响最大。这里按照坡膛面积、流速突变考虑，对于常",
            "规火炮，以修正次要功系数的方法体现这样的效应。对于无坐力炮，由于涉及压强弹速较低，目前并未",
            "进行缩进修正。",
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
    "calcButtonTxt": "调用龙格-库塔-菲尔伯格7-8阶积分器求解系统。",
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
    "figShotVel": "弹速\nm/s",
    "figBreech": "膛底",
    "figStagnation": "滞止点",
    "figNozzleP": "喷管喉压",
    "figNozzleV": "喷管喉速",
    "figOutflow": "流出比例",
    "figAvgP": "平均压强\nMPa",
    "figShotBase": "弹底",
    "figTgtP": "设计压强",
    "figPsi": "燃去比例",
}
STRING = {"English": ENGLISH, "中文": CHINESE}
