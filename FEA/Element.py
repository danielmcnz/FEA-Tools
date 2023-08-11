import numpy as np

class Element:
    
    def __init__(self):
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


    def calculate_global_stiffness(self) -> None:
        pass


    def calculate_force_vector(self, q) -> None:
        pass