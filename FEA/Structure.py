import numpy as np
from typing import List
import matplotlib.pyplot as plt


from FEA.Element import Element

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

    def __init__(self, elements : List[Element], external_force_vector : np.ndarray, element_points : np.ndarray = None) -> None:
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
        self.element_points : np.ndarray = element_points
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


    def single_deflection_display(self, magnification : int, deflection_index : int) -> None:
        """
        Displays the structure, where the first element point is the deflected point

        Returns
        -------
        None
        """

        deflected_points = self.element_points.copy()

        print(deflected_points[0])
        deflected_points[deflection_index][0] += self.q[0] * magnification
        deflected_points[deflection_index][1] += self.q[1] * magnification

        plt.plot(self.element_points[:,0], self.element_points[:,1], 'b')
        plt.plot(deflected_points[:,0], deflected_points[:,1], '--', color='r')

        plt.show()