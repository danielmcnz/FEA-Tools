from ..FEA import *

def example1():
    E = 200e9
    I = 6e-6
    A = 2e-4
    L1 = 3
    L2 = 5
    L3 = 5
    L4 = 6
    angle1 = 0
    angle2 = 53.1
    angle3 = -36.9
    angle4 = 0

    supports = [
        PinSupport(Vec2(0, 0)),
        RollerSupportY(Vec2(6, 0)),
        RollerSupportY(Vec2(6, 4))
    ]

    elements = [
        Element([Vec2(3, 4), Vec2(6, 4)], E, I, L1, A, angle1),
        Element([Vec2(0, 0), Vec2(3, 4)], E, I, L2, A, angle2),
        Element([Vec2(3, 4), Vec2(6, 0)], E, I, L3, A, angle3),
        Element([Vec2(0, 0), Vec2(6, 0)], E, I, L4, A, angle4, UDL=(0, -29e3, 1, 1))
    ]

    structure = Structure(elements, supports)

    structure.plot_structure(5, 20, True)



def example2():
    E = 200e9
    I = 1e-4
    A = 2e-4
    L1 = 4
    L2 = 3
    L3 = 4
    angle1 = -30
    angle2 = 90
    angle3 = 0

    supports = [
        FixedSupport(Vec2(0, 0)),
        PinSupport(Vec2(4, 0)),
        FixedSupport(Vec2(4-3.464, 3+2))
    ]

    elements = [
        Element([Vec2(4-3.464, 3+2), Vec2(4, 3)], E, I, L1, A, angle1),
        Element([Vec2(4, 0), Vec2(4, 3)], E, I, L2, A, angle2, LVL=(0, 1000*9.81*4*3, 1, 1)),
        Element([Vec2(0, 0), Vec2(4, 0)], E, I, L3, A, angle3, UDL=(0, -1000*9.81*4, 1, 1))
    ]

    structure = Structure(elements, supports)

    structure.plot_structure(20, 20, True)