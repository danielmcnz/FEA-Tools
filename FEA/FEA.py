import numpy as np
import matplotlib.pyplot as plt

from FEA.Vector import Vec2
from FEA.Supports import (Support, PinSupport, FixedSupport, RollerSupport)
from FEA.Element import (Element, Node, GLOB_DOF)
from FEA.Structure import Structure