ENGLISH = {
    "fontName": "Sarasa Mono SC",
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
    "athLabel": "Arc Width",
    "pdtalLabel": "P.D. / A.W.",
    "plotAvgP": "Avg. P.",
    "plotBaseP": "Base P.",
    "plotBreechNozzleP": "Breech P.",
    "plotStagP": "Stag. P.",
    "plotVel": "Shot Vel.",
    "plotNozzleV": "Throat Vel.",
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
            "Mass of propellant charge to be used. This program employ",
            "calssical Lagrangian approximation for calculation,"
            "implying a space-invariant density for propellant gasses,",
            "and combustion occuring at space-average pressure.\n"
            "This assumption is almost correct for charge to shot mass",
            "ratios -> 0. As charge to shot mass ratio increase, pressure",
            "wave phenomenon becomes significant. More sophisticated models",
            "are required to take account of such phenomenon, although",
            "empirically this program is rather accurate for cases with",
            "shot to mass ratio of close to 1 (e.g. M829A1).",
        )
    ),
    "vinfText": " ".join(
        (
            "Velocity the shot would achieve if",
            "the barrel is extended to infinite length,",
            "or equivalently, if the gas is allowed to",
            "expand to temperature of absolute zero.",
        )
    ),
    "teffText": " ".join(
        (
            "Thermal efficiency of the gun system, i.e.",
            "the amount of work done to both gas and",
            "projectile over potential from propellant.",
            "This is equivalent to 1 substracting the",
            "average temperature of exhaust gas with",
            "that of its adiabatic flame temperature.",
        )
    ),
    "beffText": " ".join(
        (
            "Ballistic efficiency of the gun system, i.e.",
            "the amount of work done the projectile over",
            "chemical potential energy of propellant,",
            "or the thermal efficiency divided with the",
            "secondary work factor.",
        )
    ),
    "specsText": " ".join(
        (
            "Specify the propellant composition used.\n",
            "Pure Nitrocellulose is mixed with inert plasticizer",
            "and stabilizer to reduce burn rate and lower",
            "flame temperature, forming single based",
            "propellant. These have excellent mechanical strength,",
            "low cost and high availability.\n",
            "Double based propellant is formed when energetic",
            "Nitroglycerin is used as gelatinizer instead.",
            "While more energetic and less sooty,",
            "it also burns hotter and erodes barrel more.",
            "DEGDN can be substituted at the cost of thermal",
            "stability if supply of NG is unavailable.\n",
            "Triple base propellants primarily consist of",
            "Nitroguanidine, improving combustion chemistry,",
            "when mixed with conventional propellant, generating",
            "higher energy content while keeping flame",
            "temperature low. However, it is mechanically",
            "the weakest, and Nitroguanidine availability is",
            "limited. Mixed Nitrate Ether propellants are",
            "available that delivers similar level of performance",
            "with higher availability.\n"
            "Nitramine propellants mix nitramines like RDX or HMX",
            "with small fraction of organic nitro propellant.",
            "Nitramines are normally explosives, but their",
            "burn characteristc can be modified for use as",
            "propellant by grounding to micrometer sized",
            "particles. HMX is preferable to RDX in terms",
            "of propellant force, but burn less controllably.",
            "Nitramine propellants tends to be considerably",
            "erosive.",
        )
    ),
    "geomPlotText": " ".join(
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
    "ldfText": " ".join(
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
            "bulk production is in the range of 1μm - 0.15mm,",
            "depending on sources.\n",
            "Arc Width is generally found close to 1mm",
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
    "perfLRText": " ".join(
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
    "cylLRText": " ".join(("Specify length to diameter ratio of the grain.",)),
    "rodRText": " ".join(
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
    "heightRText": "Specify the height to width ratio of propellant rod or flake.",
    "tolText": " ".join(
        (
            "The maximum tolerated error, ε, for each component of the",
            "integrands\n"
            "The absolute difference between the 8th order estimate and the",
            "7th order estimate, as constructed from the Runge-Kutta-",
            "Fehlberg 7(8) method, is taken to be the local truncation error.",
            "This is scaled by the ratio between integration range and step-",
            "size to estimate global truncation error. Finally, this is",
            "divided by the current value for each integrand to to derive the",
            "relative error.\n",
            "The exact error critera is constructed such that for values",
            "between (0,1) the maximum allowed error is exactly ε, and for",
            "higher values x, ε*x. The underlying SoE is expressed in scaled",
            "unitless formulation, i.e. with scaled velocity from [0,1) and",
            "scaled length & time (0,+inf).\n",
            "Defualt value corresponds to ε=1e-3, if value jumps are observed",
            "on plot, adjust this value appropriately.",
        )
    ),
    "stpText": " ".join(
        (
            "Starting Pressure, parameter used to model projectile engraving",
            "to simplfy internal ballistic solutions.\n",
            "The engraving process is often simplified to be instantaneouse."
            "That is, propellant burns in an isochoric fashion until the",
            "pressure at base of projectile is enough to overcome starting",
            "pressure, at which point starting pressure vanish, projectile",
            "starts moving as descirbed by the system of ordinary",
            "differential equation.\n",
            "Historically start pressure was determined as the highest",
            "pressure recorded when a projectile is pressed down into",
            "the barrel. Typical values for rifled small arms are in the",
            "30-40MPa range while for larger caliber rifled guns, values",
            "are reported between 10-20MPa. Although a higher start pressure",
            "somewhat increase the peak pressure, an initial resistance is",
            "necessary to improve the consistency of propellant ignition.",
            "Thus, even for smoothbore guns it is necessary provide some",
            "starting resistance.\n"
            "It is generally accepted now that projectile engraving is",
            "actually a dynamic process, with the projectile already moving",
            "at some speed. The resistance is also correspondingly higher,",
            "with detailed calculation showing it to be rather close to the",
            "peak pressure. ",
        )
    ),
    "clrText": " ".join(
        (
            "Chamber length ratio, or the average cross section of chamber",
            "over that of the barrel.\n",
            "Given the density of propellant is much less than that of the",
            "shot, having chamber that is wider than the barrel improves",
            "compactness with higher propellant loading. However, this",
            "introduce the chamberage effect, or the increase in speed of",
            "the gasses as it is compressed into the smaller barrel, most",
            "prominent at the beginning of the internal ballistic cycle.\n",
            "For conventional gun, the chamberage effect is corrected",
            "by averaging the influence of chamberage over the length of",
            "barrel and increasing the secondary work factor by a",
            "corresponding amount. For recoiless gun, due to the slow",
            "velocity of shot and low pressure, currently the chamberage",
            "effect is not corrected for.",
        )
    ),
    "dgcText": " ".join(
        (
            "Drag coefficient is used to model all resistive forces",
            "that are proportional to shot base pressure.\n"
            "2%-7% is commonly reported, with small arms that engages with",
            "the rifling with the entire projectile on the high end,",
            "and large caliber projectiles with driving band on the lower",
            "end.",
        )
    ),
    "sampText": " ".join(
        (
            "The domain that is sampled equidistantly form shot start until",
            "projectile exit.\n",
            "Sampling is done after the calculation for characteristc points,",
            "and in the case of time-domain calculations, independent of the",
            "former, enabling its use as a sanity check for the adaptive",
            "results. Length-domain sampling is not independent of",
            "characteristc point calculations, though.",
        )
    ),
    "pTgtText": "Maximum length-averaged chamber pressure.",
    "useConsText": " ".join(
        (
            "Solve the required arc thickness and barrel length to achieve",
            "peak pressure and muzzle velocity constraint, to the specified",
            "accuracy.\n",
            "Note while numerical precision is arbitrary,",
            "calculated result should be interpreted within the limites of",
            "reasonable manufacturability.",
        )
    ),
    "optLFText": " ".join(
        (
            "Find the optimum load fraction for this charge mass, minimizing",
            "tube volume (including chamber), under the constraint set out",
            "above. The all-burnt point is usually solved to be close, or",
            "outside of the muzzle, tending to cause large muzzle velcoity",
            "disperson, as propellant burn rates can be heavily influenced",
            "by factors such as ambient temperature.\n",
            "The minimum volume solution is most appropriate for direct fire,",
            "high velocity guns, such as tank guns and anti-tank guns, where",
            "compactness and mobility are cruicial, and accuracy is not too",
            "much impacted by velocity inconsistency between shots.\n",
            "For guns intended for indirect fire or long range fire, for",
            "example artillery and naval guns, an extended barrel that keeps",
            "the all-burnt point within 80% of the barrel length with top",
            "charge is more appropriate.",
        )
    ),
    "calLxText": " ".join(
        (
            "Tube length, commonly expressed in literature as L/xx.",
            "Top is measured from shot start, bottom is measured from",
            "breechface.",
        )
    ),
    "nozzExpText": " ".join(
        (
            "Area expansion ratio of the recoiless gun's rear nozzle, or the",
            "cross section of the nozzle end over that of the throat. This is",
            "used to size the throat to achieve recoiless condiiton while",
            "the projectile is still in bore.\n",
            "Usually taken to be 2-4. Higher expansion ratios are only",
            "marginally more efficient while incurring additional penalty to",
            "weight.",
        )
    ),
    "nozzEffText": " ".join(
        (
            "Efficiency of nozzle for recoiless gun. From rocket theory,",
            "it is possible to construct nozzle of decent efficiency",
            "(as high as 95%) using simple gemometries, for example",
            "a straight walled truncated cone with a 15 deg half angle.",
            "Covolume correction is probably unnecessary due to the low",
            "pressure involved in its internal ballistics, which is on",
            "the order of rocket engines.\n",
            "Nozzle Efficiency is introduced to account for geometrical",
            "inefficiencies, wear as would be expected from a service",
            "weapon, and finally to average out the effect of thrust factor",
            "(for muzzle) inreacsing as the projectile leaves the tube.\n",
            "Typically taken to be around 92%, this implies that an actual",
            "recoiless gun will tend to slightly recoil forward while",
            "projectile is in bore, and backward after it has left.",
        )
    ),
    "calcButtonText": "Integrate system using RKF7(8) integrator",
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
    "SEVEN_PERF_ROSETTE": "7 Perf. Rosette",
    "FOURTEEN_PERF_ROSETTE": "14 Perf. Rosette",
    "NINETEEN_PERF_ROSETTE": "19 Perf. Rosette",
    "NINETEEN_PERF_CYLINDER": "19 Perf. Cylinder",
    "NINETEEN_PERF_HEXAGON": "19 Perf. Hexagon",
    "NINETEEN_PERF_ROUNDED_HEXAGON": "19 Perf. Rounded Hex.",
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
    "fontName": "Sarasa Mono SC",
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
    "geomPlotText": "".join(
        (
            "相对燃烧表面σ=S/S1，关于相对厚度Z=e/e1的函数。曲线递增代表增面燃烧，曲线递减代表减面燃烧。",
            "简单形状火药燃烧全过程Z从0到1，为减面燃烧；多孔火药的燃烧分裂前Z从0到1，有条件地增面燃烧，",
            "分裂后Z从1到Zb，由于剩余体积较小，分裂出棒状体几何形状复杂，以等体积圆柱处理，燃烧体现减面性。",
            "增面燃烧有利于减缓压强上升，降低最大压强，提高压强充满率，改善膛内受力情况。",
        )
    ),
    "ldfText": "".join(
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
    "perfLRText": "".join(
        (
            "设置火药药粒的长径比。多孔火药长径比一般在1.82-3.57之间。该比例越大，火药燃烧增面性越强。",
            "当取值过小，以至于药粒柱高小于弧厚时，药粒可能出现减面燃烧的现象。",
            "一般而言较长的药粒难以同时着火，需要每隔一定间隔开槽。",
            "运动性越差，发射时停留在在膛底，作用压强越高。",
        )
    ),
    "cylLRText": "".join(("设置、管柱状火药的长径比",)),
    "rodRText": "".join(("设置长方体火药长宽比。",)),
    "widthText": "".join(("设置长方体火药宽度。不需要一定是最小跨度边。",)),
    "heightRText": "".join(("设置火药的高宽比例。",)),
    "tolText": "".join(
        (
            "   积分器最大允许误差，ε。龙格-库塔-菲尔伯格自适应积分器根据每一步各项的取值，导数，",
            "对常微分方程生成七阶，八阶两个估测。这两个估测的差值作为局部误差，再根据步长与积分域的比例",
            "外推，除以该点取值，估测全局相对误差。将该值与用户设定的误差大小作比较，按一定算法调整步长",
            "直到该条件得到满足。重复以上直到积分器达到停止条件，或遇到数值奇点。\n",
            "   鉴于实际积分的常微分系统采用归一化无量纲形式，各项取值为：速度［0，1），行程、时间",
            "［0，正无穷）。\n",
            "   默认取值为对应ε=1E-3，若由图可见采样点与特征点差距较大，提示需要提高精度计算。",
        )
    ),
    "stpText": "".join(
        (
            "   设置挤进压强。内弹道计算中常见的化简是认为在第一阶段，火药在弹后空间定容燃烧，直到克服挤进",
            "压强。此时，弹头瞬间完成挤进炮膛，与膛线结合，开始运动，火药燃烧进入由常微方程描述的第二阶段。",
            "因为该参量实际上反映的是弹道简化，因此物理意义并不十分明确。\n",
            "   历史上常用静压法测定挤进压强，即用压力机缓慢将弹头挤入膛线所需的最大静压。一般而言，对于大口径",
            "线膛炮10-20MPa，对小口径武器取30-40MPa。虽然较大的挤进压强一定程度上造成最大压强增加，",
            "但更高的挤进压强有利于改善点火条件，提高点火一致性，因此即便是滑膛炮也需要一定的挤进压强。\n",
            "   实际挤进是一个运动的过程，弹头在挤进过程中，膛压可达到最大膛压的六、七成，挤进完成时火药燃烧，",
            "弹头速度都达到出膛的一成左右，耗时大致是火药燃烧时间的四分之一到三分之一。\n",
            "   对无坐力火炮，为满足严格无坐力条件，喷口打开压强与挤进取相同数值。",
        )
    ),
    "clrText": "".join(
        (
            "药室扩大系数指的是药室截面积与炮管截面积的比例。膛压较大时，定装弹药过长的药室存在抽壳困难",
            "的情况，后膛过长也不利于火炮装填作业，并造成点火困难。\n",
            "药室缩进产生流速突变对于弹道内循环的起始阶段影响最大。这里按照坡膛面积、流速突变考虑，对于常",
            "规火炮，以修正次要功系数的方法体现这样的效应。对于无坐力炮，由于涉及压强弹速较低，未对缩进",
            "效应进行修正。",
        )
    ),
    "dgcText": "".join(
        (
            "阻力系数，即作用于弹头的一切与弹底压强相关的阻力，与弹底压强的比例。对于线膛武器，",
            "取值从2%到7%不等，其中小口径武器相对较高。",
        )
    ),
    "sampText": "".join(
        (
            "   采样时将给定域等间隔均分的数目。程序中，特征点的计算先行进行，与采样设置无关。采样点的计",
            "算在时间域上独立于特征点计算，可以用于特征点计算结果的验证。\n",
            "   需要注意的是，由于基于长度的常微分系统在零点上有数值奇点，因此在长度域上的采样点计算与",
            "特征点计算并不独立。",
        )
    ),
    "pTgtText": "设计最大空间平均压强",
    "useConsText": "".join(
        (
            "根据设计最大压强与设计弹速，按照精度要求，反算满足条件的火药弧厚，身管长度。",
            "需要注意的是，反解虽然可以按任意数值精度进行，但计算的结果不一定具有实际应用价值。",
            "实际情况下要求内弹道系统能够在火药，炮管加工偏差的情况下仍然满足技战术要求。",
        )
    ),
    "optLFText": "".join(
        (
            "   求该装药量下，满足最大设计压强和出躺弹速的的最小膛容解（含药室）。\n",
            "   最小膛容解下火药燃尽点常常接近或位于炮口，实际使用中，火药燃速容易受环境温度，批次生产药力差异，",
            "弹药装填条件影响，使得弹速的可重复性较差。这样的解一般适用于射程较短的高速直射火炮，例如坦克炮，反坦克炮；",
            "对于这些武器弹速的差异对于精度影响有限。而对于榴弹炮等射程较远，精度要求较高的武器，应适当延长炮管至燃尽",
            "点，使其在身管长度80%或以内。",
        )
    ),
    "calLxText": "".join(
        (
            "身管倍径，即将火炮身管长度表述为口径的倍率\n",
            "上为从弹底运动起点至炮口，下为从炮膛膛底计。",
        )
    ),
    "nozzExpText": "".join(
        (
            "   无坐力炮喷管的扩张系数，或喷管尾与喉部的面积比。结合无坐力条件可以确定喉口的截面积。\n",
            "   无坐力炮的喷管类似于火箭，其推力系数（喷管喷出口的动压相比燃气平均压强的放大比例）随",
            "扩张系数的提升而增大，但总体上收益是边缘递减的。更高的推力系数使满足无坐力条件所需的喉口",
            "面积减小，对于提高无坐力炮的内弹道效率有积极作用，但是降低了系统的机动性以及紧凑性，从总",
            "体设计上必须仔细考量得失。一般无坐力炮取2-4。",
        )
    ),
    "nozzEffText": "".join(
        (
            "喷管效率。参照火箭理论的结果，选用合适的扩张角下，即使简单台体的喷管也可以达到95%的",
            "推进效率。由于无坐力炮一般涉及的膛压较低，没有太大的必要针对余容进行修正。考虑到真实",
            "喷管的磨损漏气情况，同时为达到快速再发射的目的，使得弹在膛内时少量前坐，出膛后少量",
            "后坐，引入喷管效率进行修正。一般无坐力炮该值取92%左右。",
        )
    ),
    "calcButtonText": "调用龙格-库塔-菲尔伯格7-8阶积分器求解系统。",
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
