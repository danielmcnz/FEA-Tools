from ..FEA import *

def bar_element(display=True):
    E = 200e9
    L1 = 10
    L2 = 10 * 1.41
    A = np.pi * (0.1 ** 2) / 4
    angle1 = 0
    angle2 = -135

    supports = [
        FixedSupport(Vec2(0, 10)),
        BarElementSupport(Vec2(10, 10)),
        FixedSupport(Vec2(0, 0)),
    ]

    elements = [
        Element([Vec2(0, 10), Vec2(10, 10)], E, 0, L1, A, angle1),
        Element([Vec2(0, 0), Vec2(10, 10)], E, 0, L2, A, angle2)
    ]

    # external force vector
    Q = np.array([
        [0], 
        [100e3]
    ]) 

    structure = Structure(elements, supports, Q)

    if display:
        structure.plot_structure(1000, 20, True)