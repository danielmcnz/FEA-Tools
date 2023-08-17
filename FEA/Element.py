import numpy as np

# from FEA.FrameElement import FrameElement
# from FEA.BarElement import BarElement

class Element:
    
    def __init__(self):

        # self._element_type = None

        # if(self._element_type == FrameElement):
        self.structural_deflections : np.ndarray
        self.element_deflections : np.ndarray
        # elif(self._element_type == BarElement):
        self.nodal_displacements : np.ndarray
        self.strain : np.ndarray

        self.E : int
        self.L : int
        self.A : int
        self.angle : int

        self.UDL : int
        self.point_load : int

        self.UDL_forces : np.ndarray
        self.point_load_forces : np.ndarray

        self.assembly_mat : np.ndarray

        self.lambda_mat : np.ndarray

        self.local_stiffness : np.ndarray
        self.local_stiffness_hat : np.ndarray
        self.global_stiffness : np.ndarray

        self.global_force : np.ndarray
        self.local_force : np.ndarray

    
    # def set_type(self, element_type):
    #     self._element_type = type(element_type)


    def calculate_global_stiffness(self) -> None:
        pass


    def calculate_force_vector(self, q) -> None:
        pass


    def plot_element(self, nodes : np.ndarray, displacement_magnitude : int, n_points : int, q : np.ndarray = None) -> np.ndarray:
        pass