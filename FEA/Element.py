from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from typing import List

from .Supports import Support
from .Vector import Vec2


class GLOB_DOF:
    cur_index : int = 0

    def __init__(self, name : str = "", value : int = 1):
        self.value : int = value

        self.name : str = name

        self.index = -1

    def assign_dof_index(self):
        self.index = GLOB_DOF.cur_index
        GLOB_DOF.cur_index += 1


    def __add__(self, other : GLOB_DOF) -> int:
        return self.value + other.value

    
    def __add__(self, other : GLOB_DOF) -> int:
        return self.value - other.value
    

    def __eq__(self, other : GLOB_DOF):
        return self.value == other.value and self.index == other.index


    def __call__(self) -> int:
        return self.value
    
    
    def __str__(self) -> str:
        return "N/A" if self.value == 0 else str(self.index)
    

class Node():
    def __init__(self, pos : Vec2, x : GLOB_DOF, y : GLOB_DOF, moment : GLOB_DOF):
        self.pos : Vec2 = pos
        self.x : GLOB_DOF = x
        self.y : GLOB_DOF = y
        self.moment : GLOB_DOF = moment
        

    
    def __str__(self):
        return "Node -> pos: " + str(self.pos) + " dof: (" + str(self.x) + ", " + str(self.y) + ", " + str(self.moment) + ")"


    def __iter__(self):
        self._iter_list = [self.x, self.y, self.moment]
        return iter(self._iter_list)
    

    def __next__(self):
        return next(self._iter_list)


    def __eq__(self, other : Node):
        if isinstance(other, Node):
            return self.pos == other.pos and self.x == other.x and self.y == other.y and self.moment == other.moment
        return False
    

    def __hash__(self) -> int:
        return hash((self.pos, self.x, self.y, self.moment))


    def __bool__(self, other : Node):
        return self.__eq__(self, other)

    
    def to_global(self):
        pass


class Element:
    """
    An element for FEA analysis.

    Attributes
    ----------
    E : int
        The Young's modulus of the element.
    I : int
        The second moment of area of the element.
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
    
    Methods
    -------
    _to_local()
        Returns the local stiffness matrix for an element.
    _to_global()
        Returns the global stiffness matrix for an element.
    _to_structure()
        Returns the global stiffness matrix for an element.
    calculate_global_stiffness()
        Calculates the global stiffness matrix for an element.
    calculate_force_vector(q)
        Calculates the force vector for an element.
    get_lambda_mat(angle)
        calculates the lambda matrix
    """

    element_index : int = 0
     
    def __init__(self, node_pos : List[Vec2], E : int, I : int, L : int, A : int, angle : int, UDL : tuple = (0, 0, 1, 1), LVL : tuple = (0, 0, 1, 1), point_load : tuple = (0, 0, 0, 1, 1)):
        """
        Parameters
        ----------
        assembly_mat : np.ndarray
            The assembly matrix for the element.
        E : int
            The Young's modulus of the element.
        I : int
            The second moment of area of the element.
        L : int
            The length of the element.
        A : int
            The cross-sectional area of the element.
        angle : int
            The angle of the element.
        UDL : tuple
            The uniform distributed load if present on the element.
            
            format : (angle of horizontal of element to horizontal of UDL (CCW), magnitude of load, x sign, y sign)
        LVL : tuple
            The linearly distributed load if present on the element.

            format : (angle of horizontal of element to horizontal of LVL  (CCW), magnitude of load, x sign, y sign)
        point_load : tuple
            The point load if present on the element.
            
            format : (distance from node 1 (a), angle of horizontal of element to horizontal of PL (CCW), magnitude of load, x sign, y sign)

            More comments on angle: get the x and y with the point load directions, then draw triangle with point load and direction.
            This will give the +ve angle within this triangle
        Returns
        -------
        None
        """
        Element.element_index += 1

        self.id : int = Element.element_index
        

        self.node_pos : List[Vec2] = node_pos
        self.nodes : List[Node] = [None, None]
        self.n_element_dofs : int = 6
        
        self.E : int = E
        self.I : int = I
        self.L : int = L
        self.A : int = A
        self.angle : int = angle

        self.UDL : int = UDL
        self.LVL : int = LVL
        self.point_load : tuple = point_load

        self.UDL_forces : np.ndarray            = None
        self.UDL_f_shear : np.ndarray           = None
        self.UDL_f_axial : np.ndarray           = None
        self.UDL_F_shear : np.ndarray           = None
        self.UDL_F_axial : np.ndarray           = None

        self.LVL_forces : np.ndarray            = None
        self.LVL_f_shear : np.ndarray           = None
        self.LVL_f_axial : np.ndarray           = None
        self.LVL_F_axial : np.ndarray           = None
        self.LVL_F_shear : np.ndarray           = None

        self.point_load_forces : np.ndarray     = None
        self.PL_f_shear : np.ndarray            = None
        self.PL_f_axial : np.ndarray            = None
        self.PL_F_shear : np.ndarray            = None
        self.PL_F_axial : np.ndarray            = None

        self.assembly_mat : np.ndarray          = None

        self.lambda_mat : np.ndarray            = None

        self.local_stiffness : np.ndarray       = None
        self.local_stiffness_hat : np.ndarray   = None
        self.global_stiffness : np.ndarray      = None

        self.stress : np.ndarray                = None
        self.strain : np.ndarray                = None

        self.local_forces : np.ndarray          = None
        self.global_forces : np.ndarray         = None

        self.local_deflections : np.ndarray     = None
        self.global_deflections : np.ndarray    = None

    
    def _find_nodes(self, supports : List[Support]) -> np.ndarray:

        self.nodes = [None, None]

        for support in supports:
            if support.pos == self.node_pos[0]:
                self.nodes[0] = Node(self.node_pos[0], GLOB_DOF('x', support.x_dof), GLOB_DOF('y', support.y_dof), GLOB_DOF('moment', support.moment_dof))
            if support.pos == self.node_pos[1]:
                self.nodes[1] = Node(self.node_pos[1], GLOB_DOF('x', support.x_dof), GLOB_DOF('y', support.y_dof), GLOB_DOF('moment', support.moment_dof))

        for i in range(len(self.nodes)):
            if self.nodes[i] is None:
                self.nodes[i] = Node(self.node_pos[i], GLOB_DOF('x'), GLOB_DOF('y'), GLOB_DOF('moment'))
                


    def _calculate_asm_mat(self, n_structure_dofs : int):
        self.assembly_mat = np.zeros((n_structure_dofs, self.n_element_dofs), dtype=np.int32)

        for i in range(len(self.nodes)):
            index = i * 3

            if index <= self.n_element_dofs:
                if self.nodes[i].x.index >= 0:
                    self.assembly_mat[self.nodes[i].x.index][index] = self.nodes[i].x.value

                if self.nodes[i].y.index >= 0:
                    self.assembly_mat[self.nodes[i].y.index][index+1] = self.nodes[i].y.value
                    
                if self.nodes[i].moment.index >= 0:
                    self.assembly_mat[self.nodes[i].moment.index][index+2] = self.nodes[i].moment.value


    def _to_local(self) -> np.ndarray:
        """
        creates the local stiffness matrix for an element

        Returns
        -------
        np.ndarray
            The local stiffness matrix.
        """

        if(self.I == 0):
            self.I = 1e-12
        beta = (self.A * self.L ** 2) / self.I

        mat = np.array([
            [beta , 0       , 0          , -beta, 0        , 0          ],
            [0    , 12      , 6*self.L   , 0    , -12      , 6*self.L   ],
            [0    , 6*self.L, 4*self.L**2, 0    , -6*self.L, 2*self.L**2],
            [-beta, 0       , 0          , beta , 0        , 0          ],
            [0    , -12     , -6*self.L  , 0    , 12       , -6*self.L  ],
            [0    , 6*self.L, 2*self.L**2, 0    , -6*self.L, 4*self.L**2]
        ])

        self.local_stiffness = ((self.E * self.I) / (self.L**3)) * mat

        return self.local_stiffness
    

    def _to_global(self) -> np.ndarray:
        """
        creates the global stiffness matrix for an element

        Returns
        -------
        np.ndarray
            The global stiffness matrix.
        """

        self.lambda_mat = Element.get_lambda_mat(self.angle)

        self.local_stiffness_hat = self.lambda_mat.T @ self.local_stiffness @ self.lambda_mat

        return self.local_stiffness_hat


    def _to_structure(self) -> np.ndarray:
        """
        creates the structural stiffness matrix for an element

        Returns
        -------
        np.ndarray
            The structural stiffness matrix.
        """

        if(self.local_stiffness_hat is None):
            raise AssertionError("global_stiffness matrix cannot be None. First call _to_global().")

        self.global_stiffness = self.assembly_mat @ self.local_stiffness_hat @ self.assembly_mat.T

        return self.global_stiffness
    

    def _UDL_forces(self) -> None:
        """
        Calculates the uniformly distributed load forces if present on the element.

        Returns
        -------
        None
        """

        if(self.UDL == 0):
            return

        angle = self.UDL[0]
        udl = self.UDL[1]
        shear_sign = self.UDL[2]
        axial_sign = self.UDL[3]
        
        shear_mag = shear_sign * udl * np.cos(np.deg2rad(angle))
        axial_mag = axial_sign * udl * np.sin(np.deg2rad(angle))

        self.UDL_f_shear = np.array([
            [0],
            [self.L / 2],
            [self.L ** 2 / 12],
            [0],
            [self.L / 2],
            [-self.L ** 2 / 12]
        ])

        self.UDL_f_axial = np.array([
            [self.L / 2],
            [0],
            [0],
            [self.L / 2],
            [0],
            [0]
        ])

        self.UDL_f_shear *= shear_mag
        self.UDL_f_axial *= axial_mag

        self.UDL_F_shear = self.lambda_mat.T @ self.UDL_f_shear
        self.UDL_F_axial = self.lambda_mat.T @ self.UDL_f_axial

        self.UDL_forces = self.assembly_mat @ (self.UDL_F_shear + self.UDL_F_axial)


    def _LVL_forces(self) -> None:
        """
        Calculates the linearly distributed load forces if present on the element.

        Returns
        -------
        None
        """

        if(self.LVL == 0):
            return

        angle = self.LVL[0]
        lvl = self.LVL[1]
        shear_sign = self.LVL[2]
        axial_sign = self.LVL[3]        

        shear_mag = shear_sign * lvl * np.cos(np.deg2rad(angle))
        axial_mag = axial_sign * lvl * np.sin(np.deg2rad(angle))

        self.LVL_f_shear = np.array([
            [0],
            [3 * self.L / 20],
            [self.L ** 2 / 30],
            [0],
            [7 * self.L / 20],
            [-self.L ** 2 / 20]
        ])

        self.LVL_f_axial = np.array([
            [self.L / 2],
            [0],
            [0],
            [self.L / 2],
            [0],
            [0]
        ])

        self.LVL_f_shear *= shear_mag
        self.LVL_f_axial *= axial_mag

        self.LVL_F_shear = self.lambda_mat.T @ self.LVL_f_shear
        self.LVL_F_axial = self.lambda_mat.T @ self.LVL_f_axial

        self.LVL_forces = self.assembly_mat @ (self.LVL_F_shear + self.LVL_F_axial)


    def _point_load_forces(self) -> None:
        """
        Calculates the point load forces if present on the element.
            This includes the shear and axial equivalent forces (locally and globally), and then
            creates the force vectors.

        Returns
        -------
        None
        """

        if(self.point_load[2] == 0):
            return

        a = self.point_load[0]
        angle = self.point_load[1]
        point_load = self.point_load[2]
        x_sign = self.point_load[3]
        y_sign = self.point_load[4]

        shear_mag = y_sign * point_load * np.cos(np.deg2rad(angle))
        axial_mag = x_sign * point_load * np.sin(np.deg2rad(angle))

        self.PL_f_shear = np.array([
            [0],
            [1 - 3 * a ** 2 / self.L ** 2 + 2 * a ** 3 / self.L ** 3],
            [a ** 3 / self.L ** 2 - 2 * a ** 2 / self.L + a],
            [0],
            [3 * a ** 2 / self.L ** 2 - 2 * a ** 3 / self.L ** 3],
            [a ** 3 / self.L ** 2 - a ** 2 / self.L]
        ])

        self.PL_f_shear *= shear_mag

        self.PL_f_axial = np.array([
            [1 - a / self.L],
            [0],
            [0],
            [a / self.L],
            [0],
            [0]
        ])

        self.PL_f_axial *= axial_mag

        self.PL_F_shear = self.lambda_mat.T @ self.PL_f_shear
        self.PL_F_axial = self.lambda_mat.T @ self.PL_f_axial

        self.point_load_forces = self.assembly_mat @ (self.PL_F_shear + self.PL_F_axial)


    def calculate_global_stiffness(self) -> None:
        """
        Calculates the global stiffness matrix for an element.
        
        Returns
        -------
        None
        """

        self._to_local()
        self._to_global()
        self._to_structure()

        self._UDL_forces()
        self._LVL_forces()
        self._point_load_forces()


    @staticmethod
    def get_local_deflections(lmbda : np.ndarray, D_e : np.ndarray = None, A : np.ndarray = None, q : np.ndarray = None) -> np.ndarray:
        error_msg = """
Incorrect format. Below are possible equations:

[lmbda, D_e] -> d_e = lmbda @ D_e
[lmbda, A, q] -> d_e = lmbda @ A.T @ q
"""

        if D_e is not None:
            return lmbda @ D_e
        elif A is not None and q is not None:
            return lmbda @ A.T @ q
        else:
            raise ValueError(error_msg)

    
    @staticmethod
    def get_global_deflections(A : np.ndarray, q : np.ndarray) -> np.ndarray:
        return A.T @ q


    @staticmethod
    def get_local_forces(K_e : np.ndarray, d_e : np.ndarray = None, lmbda : np.ndarray = None, D_e : np.ndarray = None) -> np.ndarray:
        error_msg = """
Incorrect format. Below are possible equations:

[K_e, d_e] -> f_e = K_e @ d_e
[K_e, lmbda, D_e] -> f_e = K_e @ lmbda @ D_e
"""

        if d_e is not None:
            return K_e @ d_e
        elif lmbda is not None and D_e is not None:
            return K_e @ lmbda @ D_e
        else:
            raise ValueError(error_msg)


    @staticmethod
    def get_global_forces(K_e_hat : np.ndarray, D_e : np.ndarray) -> np.ndarray:
        return K_e_hat @ D_e


    def calculate_force_vector(self, q) -> None:
        """
        Calculates the force vector for an element.
        
        Parameters
        ----------
        q : np.ndarray
            The solution to the equation of motion.
            
        Returns
        -------
        None
        """

        self.global_deflections = Element.get_global_deflections(self.assembly_mat, q)

        self.local_deflections = Element.get_local_deflections(self.lambda_mat, self.global_deflections)

        self.local_forces = Element.get_local_forces(self.local_stiffness, self.local_deflections)

        self.global_forces = Element.get_global_forces(self.local_stiffness_hat, self.global_deflections)


    def calculate_stresses_and_strains(self):

        self.strain = (self.local_deflections[3:] - self.local_deflections[:3]) / self.L

        self.stress = (self.local_forces[3:] - self.local_forces[:3]) / self.A
        
    
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

        c = np.cos(np.deg2rad(angle))
        s = np.sin(np.deg2rad(angle))

        lambda_mat_3x3 = np.array([[c , s, 0],
                                   [-s, c, 0],
                                   [0 , 0, 1]
        ])

        lambda_mat = np.zeros((6, 6))

        lambda_mat[0:3, 0:3] = lambda_mat_3x3
        lambda_mat[3:6, 3:6] = lambda_mat_3x3

        return lambda_mat
    

    def calculate_deflections(self, x, q : np.ndarray = None) -> np.ndarray:
        """
        Calculates the Xg and Yg global deflections of an element.

        Parameters
        ----------
        x : np.ndarray
            The x values to calculate the deflections at (can be single value, or array).

        Returns
        -------
        np.ndarray
            The Xg and Yg global deflections.
        """

        if q is not None:
            d_e = Element.get_local_deflections(self.lambda_mat, A=self.assembly_mat, q=q)
        else:
            d_e = self.local_deflections

        phi_1 = (1 - x / self.L)
        phi_2 = x / self.L

        N_1 = 1 - 3 * (x ** 2) / (self.L ** 2) + 2 * (x ** 3) / (self.L ** 3)
        N_2 = (x ** 3) / self.L ** 2 - 2 * (x ** 2) / self.L + x
        N_3 = 3 * (x ** 2) / (self.L ** 2) - 2 * (x ** 3) / (self.L ** 3)
        N_4 = (x ** 3) / (self.L ** 2) - (x ** 2) / self.L

        element_axial_displacement = phi_1 * d_e[0] + phi_2 * d_e[3]
        element_transverse_displacement = N_1 * d_e[1] + N_2 * d_e[2] + N_3 * d_e[4] + N_4 * d_e[5]

        deflections_XG = element_axial_displacement * np.cos(np.deg2rad(self.angle)) - element_transverse_displacement * np.sin(np.deg2rad(self.angle))
        deflections_YG = element_axial_displacement * np.sin(np.deg2rad(self.angle)) + element_transverse_displacement * np.cos(np.deg2rad(self.angle))

        return np.array([deflections_XG, deflections_YG])
    

    def plot_DOF(self) -> None:
        """
        """

        font_size = 20

        width = 16
        height = 9

        fig, ax = plt.subplots(figsize=(width, height))

        ax.plot([self.nodes[0].pos.x, self.nodes[1].pos.x], [self.nodes[0].pos.y, self.nodes[1].pos.y], 'b')

        arrow_space_from_element = Vec2(0.15, 0.15)
        arrow_len = 0.5
        
        kw = dict(arrowstyle="Simple, tail_width=0.5, head_width=4, head_length=8", color='k')

        # ------------------------------- #
        #   plot element no. and angle    #
        # ------------------------------- #

        mid_point = (self.nodes[1].pos + self.nodes[0].pos) / 2

        ax.annotate(f"element {self.id}", xy=(mid_point.x, mid_point.y), xytext=(-110, 170), textcoords='offset points', color='black', fontsize=font_size)
        ax.annotate(f"$\u03B1^{self.id}$={str(self.angle)}", xy=(mid_point.x, mid_point.y), xytext=(-110, 150), textcoords='offset points', color='black', fontsize=font_size)

        # ------------------------------- #

        # ------------------------------- #
        #  Plot local coordinate system   #
        # ------------------------------- #

        start_point = self.nodes[0].pos
        end_point = self.nodes[1].pos

        # Calculate the direction vector of the element
        direction = end_point - start_point

        midpoint = (end_point - start_point).abs() / 2
        midpoint.x += min(start_point.x, end_point.x)
        midpoint.y += min(start_point.y, end_point.y)

        # Normalize the direction vector
        normalized_direction = direction.normalize()

        # Calculate the perpendicular vector for the local y-axis
        perpendicular = normalized_direction.perpendicular()

        midpoint += Vec2(0.15, 0.15) * perpendicular

        # Calculate the coordinates for the local axes
        start_x_axis = (midpoint.x, midpoint.y)
        end_x_axis = (midpoint.x + arrow_len * normalized_direction.x,
                        midpoint.y + arrow_len * normalized_direction.y)
        start_y_axis = (midpoint.x, midpoint.y)
        end_y_axis = (midpoint.x + arrow_len * perpendicular.x,
                        midpoint.y + arrow_len * perpendicular.y)

        arrow_x = patches.FancyArrowPatch(
            start_x_axis, 
            end_x_axis, 
            **kw
        )

        arrow_y = patches.FancyArrowPatch(
            start_y_axis, 
            end_y_axis,
            **kw
        )

        ax.annotate("$x^e$", xy=end_x_axis, xytext=(5, 0), textcoords='offset points', color='black', fontsize=font_size)
        ax.annotate("$y^e$", xy=end_y_axis, xytext=(0, 5), textcoords='offset points', color='black', fontsize=font_size)

        ax.add_patch(arrow_x)
        ax.add_patch(arrow_y)
        
        # ------------------------------- #


        # ------------------------------- #
        # Plot degrees of freedom arrows  #
        # ------------------------------- #

        for i in range(len(self.nodes)):
            p_dist = (self.nodes[i].pos - self.nodes[1-i].pos).normalize()

            p : Vec2 = self.nodes[i].pos

            x_init : float = p.x
            y_init : float = p.y

            if self.nodes[i].pos.x > self.nodes[1-i].pos.x:
                x_init += arrow_space_from_element.x
            else:
                x_init -= arrow_len + arrow_space_from_element.x

            if self.nodes[i].pos.y > self.nodes[1-i].pos.y:
                y_init = p.y + arrow_space_from_element.y
            else:
                y_init = p.y

            if self.nodes[i].x:
                a = patches.FancyArrowPatch(
                    (x_init, y_init), 
                    (x_init + arrow_len, y_init),
                    **kw
                )

                ax.annotate("D"+str(i*2+(i+1)), xy=(x_init + arrow_len, y_init), xytext=(0, 0), textcoords='offset points', color='black', fontsize=font_size)

                plt.gca().add_patch(a)

            if self.nodes[i].y:

                a = patches.FancyArrowPatch(
                    (x_init, y_init), 
                    (x_init, y_init + arrow_len), 
                    **kw
                )

                ax.annotate("D"+str(i*2+(i+2)), xy=(x_init, y_init + arrow_len), xytext=(0, 0), textcoords='offset points', color='black', fontsize=font_size)

                plt.gca().add_patch(a)

            if self.nodes[i].moment:
                center = (x_init, y_init)
                radius = 0.2

                start_angle = 0
                end_angle = 270  # 3/4 of a full circle

                arc = patches.Arc(center, radius*2, radius*2, angle=0, theta1=start_angle, theta2=end_angle,
                                linewidth=2, color='black')

                # Calculate the coordinates of the end point of the arc
                end_angle_rad = np.radians(end_angle)
                end_x = center[0] + radius * np.cos(end_angle_rad)
                end_y = center[1] + radius * np.sin(end_angle_rad)

                # Create an arrow patch at the end of the arc
                arrow = patches.FancyArrowPatch(
                    (end_x - 0.01, end_y), 
                    (end_x + 0.05, end_y),
                    **kw
                )

                ax.annotate("D"+str(i*2+(i+3)), xy=(end_x - 0.05, end_y), xytext=(0, -20), textcoords='offset points', color='black', fontsize=font_size)

                plt.gca().add_patch(arc)
                plt.gca().add_patch(arrow)
            
        # ------------------------------- #

        ax.axis('off')
        ax.axis('equal')

        plt.show()



    def plot_element(self, displacement_magnitude : int, resolution : int, axes, deflections : int, external_deflections : np.ndarray) -> None:
        """
        Plots the element in both deflected and undeflected form.

        Parameters
        ----------
        nodes : np.ndarray
            The nodes of the element.
        displacement_magnitude : int
            The magnitude to increase the displacements visually by.
        resolution : int
            The number of points to plot the element with.

        Returns
        -------
        None
        """
        
        x = np.linspace(0, self.L, resolution)

        deflections_XG, deflections_YG = self.calculate_deflections(x)

        x_undeflected = np.linspace(self.nodes[0].pos.x, self.nodes[1].pos.x, resolution)
        y_undeflected = np.linspace(self.nodes[0].pos.y, self.nodes[1].pos.y, resolution)

        x_deflected = x_undeflected + deflections_XG * displacement_magnitude
        y_deflected = y_undeflected + deflections_YG * displacement_magnitude

        axes.plot(x_undeflected, y_undeflected, 'b', label="Undeflected")

        if deflections:
            axes.plot(x_deflected, y_deflected, 'r', label="Deflected")

        if external_deflections is not None:
            deflections_XG, deflections_YG = self.calculate_deflections(x, external_deflections)
            
            x_deflected = x_undeflected + deflections_XG * displacement_magnitude
            y_deflected = y_undeflected + deflections_YG * displacement_magnitude
        
            axes.plot(x_deflected, y_deflected, 'g', label="External Deflections")
