import numpy as np
from typing import List, Dict
import matplotlib.pyplot as plt

from FEA.Element import Element, Node, GLOB_DOF
from FEA.Supports import Support


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
    _solve_EOM(K_g)
        Solves the equation of motion K_g * q = F for q.
    solve()
        Solves for the structure elements and finds the force vectors locally and globally.
    """

    def __init__(self, elements : List[Element], supports : List[Support], external_force_vector : np.ndarray) -> None:
        """
        Parameters
        ----------
        elements : List[Element]
            The list of elements in the structure.
        supports : List[Support]
            The list of supports in the structure.
        external_force_vector : np.ndarray
            The external force vector for the structure.
            
        Returns
        -------
        None
        """

        GLOB_DOF.cur_index = 0
        
        self.elements : List[Element] = elements
        self.supports : List[Support] = supports
        self.external_force_vector : np.ndarray = external_force_vector.astype(np.float64)
        self.total_stiffness : np.ndarray = None

        # equation for motions from K_g * q = F
        self.q : np.ndarray = None

    
    def _calculate_asm_mats(self):
        dofs : List[Node] = []
        n_dofs : int = 0
        for element in self.elements:
            element._find_nodes(self.supports)

            for node in element.nodes:
                if node not in dofs:
                    dofs.append(node)
                    n_dofs += node.x.value + node.y.value + node.moment.value

        for node in dofs:
            for dof in node:
                if(dof.value):
                    dof.assign_dof_index()
                
        # this is very scuffed, need to optimise this
        for i in range(len(self.elements)):
            for j in range(len(self.elements[i].nodes)):
                if self.elements[i].nodes[j] not in dofs:
                    for n in dofs:
                        if n.pos == self.elements[i].nodes[j].pos:
                            self.elements[i].nodes[j] = n
        
        for element in self.elements:
            element._calculate_asm_mat(n_dofs)


    def _solve_EOM(self, K_g : np.ndarray) -> np.ndarray:
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

        self._calculate_asm_mats()

        self.total_stiffness : np.ndarray = np.zeros((self.external_force_vector.shape[0], self.external_force_vector.shape[0]))

        for element in self.elements:
            element.calculate_global_stiffness()
            self.total_stiffness += element.global_stiffness

        for element in self.elements:
            if(element.UDL_forces is not None):
                self.external_force_vector += element.UDL_forces
            if(element.LVL_forces is not None):
                self.external_force_vector += element.LVL_forces
            if(element.point_load_forces is not None):
                self.external_force_vector += element.point_load_forces

        self._solve_EOM(self.total_stiffness)

        for element in self.elements:
            element.calculate_force_vector(self.q)

        np.set_printoptions(formatter={'float': '{:0.5e}'.format})


    def plot_structure(self, displacement_magnitude : int, resolution : int) -> None:
        """
        Plots the structure.

        Parameters
        ----------
        nodes : np.ndarray
            The nodes of the structure, made up of the elements of all the nodes.
        displacement_magnitude : int
            The magnitude to increase the displacements visually by.
        resolution : int
            The number of points to plot between each node.

        Returns
        -------
        None
        """

        i = 0

        for element in self.elements:
            element.plot_element(displacement_magnitude, resolution)
            i += 1

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))

        plt.grid()
        plt.legend(by_label.values(), by_label.keys())
        plt.show()

        plt.show()