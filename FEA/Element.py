import numpy as np
import matplotlib.pyplot as plt

# from FEA.FrameElement import FrameElement
# from FEA.BarElement import BarElement

class Element:
    
    def __init__(self):

        # if(self._element_type == FrameElement):
        
        self.structural_deflections : np.ndarray
        self.element_deflections : np.ndarray

        self.UDL : int                  = None
        self.LVL : int                  = None
        self.point_load : int           = None

        self.UDL_forces : np.ndarray    = None
        self.UDL_f_eq : np.ndarray      = None
        self.UDL_F_eq : np.ndarray      = None
        self.LVL_forces : np.ndarray    = None
        self.LVL_f_eq : np.ndarray      = None
        self.LVL_F_eq : np.ndarray      = None
        self.point_load_forces : np.ndarray = None
        self.PL_f_shear : np.ndarray    = None
        self.PL_f_axial : np.ndarray    = None
        self.PL_F_shear : np.ndarray    = None
        self.PL_F_axial : np.ndarray    = None

        # elif(self._element_type == BarElement):

        self.nodal_displacements : np.ndarray

        self.E : int
        self.L : int
        self.A : int
        self.angle : int

        self.assembly_mat : np.ndarray

        self.lambda_mat : np.ndarray

        self.local_stiffness : np.ndarray
        self.local_stiffness_hat : np.ndarray
        self.global_stiffness : np.ndarray

        self.global_force : np.ndarray
        self.local_force : np.ndarray

        self.strain : np.ndarray


    def calculate_global_stiffness(self) -> None:
        pass


    def calculate_force_vector(self, q) -> None:
        pass

    def calculate_deflections(self, x) -> None:
        pass


    def plot_element(self, nodes : np.ndarray, displacement_magnitude : int, n_points : int, q : np.ndarray = None) -> np.ndarray:
        pass