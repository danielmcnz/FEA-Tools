from ..FEA import *

def two_elements(display=True):
    supports = [
        FixedSupport(Vec2(0, 10)),
        RollerSupportX(Vec2(0, 0)),
        PinSupport(Vec2(10, 0))
    ]

    elements = [
        Element([Vec2(0, 10), Vec2(0, 0)], 200e9, 5e-4, 10, 1e-5, -90),
        Element([Vec2(0, 0), Vec2(10, 0)], 200e9, 5e-4, 10, 1e-5, 0)
    ]

    Q = np.array([[0], [140e3], [0]])

    structure = Structure(elements, supports, Q)

    if display:
        structure.plot_structure(100, 20, deflections=True, annotations=False)