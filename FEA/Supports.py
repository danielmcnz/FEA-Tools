import numpy as np

from .Vector import Vec2


class Support:
    def __init__(self, pos : Vec2):
        self.pos : Vec2 = pos

        self.x_dof : int = 0
        self.y_dof : int = 0
        self.moment_dof : int = 0

class PinSupport(Support):

    def __init__(self, pos : Vec2):
        Support.__init__(self, pos)
        
        self.moment_dof = 1

class FixedSupport(Support):

    def __init__(self, pos : Vec2):
        Support.__init__(self, pos)


class RollerSupportX(Support):

    def __init__(self, pos : Vec2):
        Support.__init__(self, pos)

        self.x_dof = 1
        self.moment_dof = 1


class RollerSupportY(Support):

    def __init__(self, pos : Vec2):
        Support.__init__(self, pos)

        self.y_dof = 1
        self.moment_dof = 1


class BarElementSupport(Support):

    def __init__(self, pos : Vec2):
        Support.__init__(self, pos)

        self.x_dof = 1
        self.y_dof = 1