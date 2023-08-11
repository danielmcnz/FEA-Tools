import numpy as np
from typing import List

from FEA.Element import Element


class BarElement(Element):
    """
    A bar element for FEA analysis.

    Attributes
    ----------
    E : int
        The Young's modulus of the element.
    L : int
        The length of the element.
    A : int
        The cross-sectional area of the element.
    angle : int
        The angle of the element.
    assembly_mat : np.ndarray
        The assembly matrix for the element.
    lambda_mat : np.ndarray
        The lambda matrix for the element.
    local_stiffness : np.ndarray
        The local stiffness matrix for the element.
    local_stiffness_hat : np.ndarray
        The local hat stiffness matrix for the element.
    global_stiffness : np.ndarray
        The global stiffness matrix for the element.
    global_force_vector : np.ndarray
        The global force vector for the element.
    local_force_vector : np.ndarray
        The local force vector for the element.
    strain : np.ndarray
        The strain in the element.
    
    Methods
    -------
    __to_local()
        Returns the local stiffness matrix for a bar element.
    __to_global()
        Returns the global stiffness matrix for a bar element.
    __to_structure()
        Returns the global stiffness matrix for a bar element.
    calculate_global_stiffness()
        Calculates the global stiffness matrix for a bar element.
    calculate_force_vector(q)
        Calculates the force vector for a bar element.
    get_lambda_mat(angle)
        calculates the lambda matrix
    local_to_global_deflections(element_angle, local_deflections)
        calculates the global deflections from local deflections
    global_to_local_deflections(element_angle, global_deflections)
        calculates the local deflections from global deflections
    local_to_global_forces(lambda_mat, local_forces)
        calculates the global forces from local forces
    """

    def __init__(self, assembly_mat : np.ndarray, E : int, L : int, A : int, angle : int) -> None:
        """
        Parameters
        ----------
        assembly_mat : np.ndarray
            The assembly matrix for the element.
        E : int
            The Young's modulus of the element.
        L : int
            The length of the element.
        A : int
            The cross-sectional area of the element.
        angle : int
            The angle of the element.
        
        Returns
        -------
        None
        """

        super().__init__()

        # assert(assembly_mat.shape == (2, 4))

        self.E : int = E
        self.L : int = L
        self.A : int = A
        self.angle : int = angle

        self.assembly_mat : np.ndarray = assembly_mat

        self.lambda_mat : np.ndarray = None

        self.local_stiffness : np.ndarray = None
        self.local_stiffness_hat : np.ndarray = None
        self.global_stiffness : np.ndarray = None

        self.global_force : np.ndarray = None
        self.local_force : np.ndarray = None
        self.strain : np.ndarray = None


    def __to_local(self) -> np.ndarray:
        """
        Returns the local stiffness matrix for a bar element.

        Returns
        -------
        np.ndarray
            The local stiffness matrix.
        """

        mat = np.array([
            [1, -1],
            [-1, 1]]
        )

        self.local_stiffness = ((self.E * self.A) / self.L) * mat

        return self.local_stiffness
    

    def __to_global(self) -> List[np.ndarray]:
        """
        Returns the global stiffness matrix for a bar element.

        Returns
        -------
        List[np.ndarray]
            The local hat stiffness matrix and the lambda matrix.
        """

        if(self.local_stiffness is None):
            raise AssertionError("local  stiffness and local stiffness hat matrices cannot be None. First call __to_local().")

        if(self.local_stiffness.shape != (2, 2)):
            raise ValueError("local_stiffness must be a 2x2 matrix")
        
        self.lambda_mat = BarElement.get_lambda_mat(self.angle)

        self.local_stiffness_hat = self.lambda_mat.T @ self.local_stiffness @ self.lambda_mat

        return self.local_stiffness_hat, self.lambda_mat # global_stiffnes, lambda
    

    def __to_structure(self) -> np.ndarray:
        """
        Returns the global stiffness matrix for a bar element.
        
        Returns
        -------
        np.ndarray
            The global stiffness matrix for the element.
        """

        if(self.local_stiffness_hat is None and self.lambda_mat):
            raise AssertionError("global_stiffness matrix cannot be None. First call __to_global().")

        self.global_stiffness = self.assembly_mat @ self.local_stiffness_hat @ self.assembly_mat.T

        return self.global_stiffness

    
    def calculate_global_stiffness(self) -> None:
        """
        Calculates the global stiffness matrix for a bar element.
        
        Returns
        -------
        None
        """

        self.__to_local()
        self.__to_global()
        self.__to_structure()

    
    def calculate_force_vector(self, q) -> None:
        """
        Calculates the force vector for a bar element.
        
        Parameters
        ----------
        q : np.ndarray
            The solution to the equation of motion.
            
        Returns
        -------
        None
        """

        if(self.global_stiffness is None):
            raise AssertionError("global_stiffness matrix cannot be None. First call calculate_global_stiffness().")

        self.global_force = self.local_stiffness_hat @ (self.assembly_mat.T @ q) # in global co-ordinates

        nodal_displacement_vector = self.lambda_mat @ self.assembly_mat.T @ q # in local co-ordinates

        self.strain = (nodal_displacement_vector[:,0][1] - nodal_displacement_vector[:,0][0]) / self.L

        self.local_force = self.local_stiffness @ nodal_displacement_vector


    @staticmethod
    def get_lambda_mat(angle : int) -> np.ndarray:
        """
        calculates the lambda matrix

        Parameters
        ----------
        angle : int
            The angle of the element.

        Returns
        -------
        np.ndarray
            The lambda matrix.        
        """

        return np.array([
            [np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle)), 0, 0],
            [0, 0, np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))]]
        )

    
    @staticmethod
    def local_to_global_deflections(element_angle : int, local_deflections : np.ndarray) -> np.ndarray:
        """
        calculates the global deflections from local deflections
        
        Parameters
        ----------
        element_angle : int
            The angle of the element.
        local_deflections : np.ndarray
            The local deflections of the element.
        
        Returns
        -------
        np.ndarray
            The global deflections.
        """

        transformation_mat = BarElement.get_lambda_mat(element_angle)

        return transformation_mat @ local_deflections

    
    @staticmethod
    def global_to_local_deflections(element_angle : int, global_deflections : np.ndarray) -> np.ndarray:
        """
        calculates the local deflections from global deflections
        
        Parameters
        ----------
        element_angle : int
            The angle of the element.
        global_deflections : np.ndarray
            The global deflections of the element.
        
        Returns
        -------
        np.ndarray
            The local deflections.
        """

        transformation_mat = BarElement.get_lambda_mat(element_angle).T

        return transformation_mat @ global_deflections
    

    @staticmethod
    def local_to_global_forces(lambda_mat : np.ndarray, local_forces : np.ndarray) -> np.ndarray:
        """
        calculates the global forces from local forces
        
        Parameters
        ----------
        lambda_mat : np.ndarray
            The lambda matrix of the element.
        local_forces : np.ndarray
            The local forces of the element.
        
        Returns
        -------
        np.ndarray
            The global forces.
        """

        return lambda_mat.T @ local_forces