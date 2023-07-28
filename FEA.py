import numpy as np
from typing import List


class Element:
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
    bar_element_mat : np.ndarray
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

        self.E = E
        self.L = L
        self.A = A
        self.angle = angle

        self.assembly_mat : np.ndarray = assembly_mat

        self.bar_element_mat : np.ndarray = None # lambda

        self.local_stiffness : np.ndarray = None
        self.local_stiffness_hat : np.ndarray = None
        self.global_stiffness : np.ndarray = None

        self.global_force_vector : np.ndarray = None
        self.local_force_vector : np.ndarray = None
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
        
        self.bar_element_mat = np.array([
            [np.cos(np.deg2rad(self.angle)), np.sin(np.deg2rad(self.angle)), 0, 0],
            [0, 0, np.cos(np.deg2rad(self.angle)), np.sin(np.deg2rad(self.angle))]]
        )

        self.local_stiffness_hat = self.bar_element_mat.T @ self.local_stiffness @ self.bar_element_mat

        return self.local_stiffness_hat, self.bar_element_mat # global_stiffnes, lambda
    

    def __to_structure(self) -> np.ndarray:
        """
        Returns the global stiffness matrix for a bar element.
        
        Returns
        -------
        np.ndarray
            The global stiffness matrix for the element.
        """

        if(self.local_stiffness_hat is None and self.bar_element_mat):
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

        self.global_force_vector = self.local_stiffness_hat @ (self.assembly_mat.T @ q) # in global co-ordinates

        nodal_displacement_vector = self.bar_element_mat @ self.assembly_mat.T @ q # in local co-ordinates

        self.strain = (nodal_displacement_vector[:,0][1] - nodal_displacement_vector[:,0][0]) / self.L

        self.local_force_vector = self.local_stiffness @ nodal_displacement_vector


class Structure:
    """
    A structure for FEA analysis, made up on a collection of bar elements.

    Attributes
    ----------
    elements : List[Element]
        The list of elements in the structure.
    external_force_vector : np.ndarray
        The external force vector for the structure.
    q : np.ndarray
        The solution to the equation of motion.

    Methods
    -------
    __solve_EOM(K_g)
        Solves the equation of motion K_g * q = F for q.
    solve()
        Solves for the structure elements and finds the force vectors locally and globally.
    """

    def __init__(self, elements : List[Element], external_force_vector : np.ndarray) -> None:
        """
        Parameters
        ----------
        elements : List[Element]
            The list of elements in the structure.
        external_force_vector : np.ndarray
            The external force vector for the structure.
            
        Returns
        -------
        None
        """
        
        self.elements : List[Element] = elements
        self.external_force_vector : np.ndarray = external_force_vector

        # equation for motions from K_g * q = F
        self.q : np.ndarray = None


    def __solve_EOM(self, K_g : np.ndarray) -> np.ndarray:
        """
        Solves the equation of motion K_g * q = F for q.

        Parameters
        ----------
        K_g : np.ndarray
            The global stiffness matrix.

        Returns
        -------
        np.ndarray
            The solution to the equation of motion.
        """

        if(self.external_force_vector.shape != (K_g.shape[0], 1)):
            raise ValueError("external_force_vector must be a column vector of size K_g.shape[0]")
        
        self.q = np.linalg.solve(K_g, self.external_force_vector)

        return self.q
    

    def solve(self) -> None:
        """
        Solves for the structure elements and finds the force vectors locally and globally.

        Returns
        -------
        None
        """

        total_stiffness : np.ndarray = np.zeros((self.external_force_vector.shape[0], self.external_force_vector.shape[0]))

        for element in self.elements:
            element.calculate_global_stiffness()
            total_stiffness += element.global_stiffness

        self.__solve_EOM(total_stiffness)

        for element in self.elements:
            element.calculate_force_vector(self.q)