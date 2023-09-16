import numpy as np

from .Vector import Vec2

class Direction:
    HORIZONTAL : Vec2 = Vec2(1, 0)
    VERTICAL : Vec2 = Vec2(0, 1)


class Support:
    def __init__(self, pos : Vec2, direction : Vec2):
        self.pos : Vec2 = pos

        self.x_dof : int = 0
        self.y_dof : int = 0
        self.moment_dof : int = 0

        self.direction : Vec2 = direction

class PinSupport(Support):

    def __init__(self, pos : Vec2, direction : Vec2):
        Support.__init__(self, pos, direction)
        
        self.moment_dof = 1

class FixedSupport(Support):

    def __init__(self, pos : Vec2, direction : Vec2):
        Support.__init__(self, pos, direction)


class RollerSupport(Support):

    def __init__(self, pos : Vec2, direction : Vec2):
        Support.__init__(self, pos, direction)

        if direction.x == 1:
            self.x_dof = 1
        elif direction.y == 1:
            self.y_dof = 1
        else:
            raise ValueError("Direction must be either horizontal or vertical")

        self.moment_dof = 1


class BarElementSupport(Support):

    def __init__(self, pos : Vec2, direction : Vec2):
        Support.__init__(self, pos, direction)

        self.x_dof = 1
        self.y_dof = 1