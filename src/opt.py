class Constrained:
    def __init__(
        self,
        caliber,
        shotMass,
        startPressure,
        chamberExpansion,
        dragCoe,
        designPressure,
        designVelocity,
    ):
        # constants for constrained designs
        self.caliber = caliber
        self.shotMass = shotMass
        self.startPressure = startPressure
        self.chamberExpansion = chamberExpansion
        self.dragCoe = dragCoe

        # design limits
        self.designPressure = designPressure
        self.designVelocity = designVelocity
