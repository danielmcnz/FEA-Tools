import numpy as np
from typing import List, Dict
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from .Element import Element, Node, GLOB_DOF
from .Supports import Support, RollerSupport, PinSupport, FixedSupport, Direction
from .Vector import Vec2


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

    def __init__(self, elements : List[Element], supports : List[Support], Q : np.ndarray = None) -> None:
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

        self.global_nodes : List[Node] = []
        self.n_dofs : int = 0

        self._original_Q = Q
        
        self.elements : List[Element] = elements
        self.supports : List[Support] = supports
        self.Q : np.ndarray = None
        self.total_stiffness : np.ndarray = None

        self.q : np.ndarray = None

        self.solve()

    
    def _calculate_asm_mats(self):
        self.n_dofs = 0
        for element in self.elements:
            element._find_nodes(self.supports)

            for node in element.nodes:
                if node not in self.global_nodes:
                    self.global_nodes.append(node)
                    self.n_dofs += node.x.value + node.y.value + node.moment.value

        for node in self.global_nodes:
            for dof in node:
                if(dof.value):
                    dof.assign_dof_index()
                
        # this is very scuffed, need to optimise this
        for i in range(len(self.elements)):
            for j in range(len(self.elements[i].nodes)):
                if self.elements[i].nodes[j] not in self.global_nodes:
                    for n in self.global_nodes:
                        if n.pos == self.elements[i].nodes[j].pos:
                            self.elements[i].nodes[j] = n
        
        for element in self.elements:
            element._calculate_asm_mat(self.n_dofs)


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

        if(self.Q.shape != (K_g.shape[0], 1)):
            raise ValueError("external_force_vector must be a column vector of size K_g.shape[0]")

        self.q = np.linalg.solve(K_g, self.Q)

        return self.q
    

    def solve(self) -> None:
        """
        Solves for the structure elements and finds the force vectors locally and globally.

        Returns
        -------
        None
        """
        
        GLOB_DOF.cur_index = 0
        Element.element_index = 0

        self.global_nodes.clear()

        self._calculate_asm_mats()

        if self._original_Q is None:
            self._original_Q = np.zeros((self.n_dofs, 1))
        self.Q = self._original_Q.astype(np.float64)

        self.total_stiffness : np.ndarray = np.zeros((self.Q.shape[0], self.Q.shape[0]))

        for element in self.elements:
            element.calculate_global_stiffness()
            self.total_stiffness += element.global_stiffness

        for element in self.elements:
            if(element.UDL_forces is not None):
                self.Q += element.UDL_forces
            if(element.LVL_forces is not None):
                self.Q += element.LVL_forces
            if(element.point_load_forces is not None):
                self.Q += element.point_load_forces

        self._solve_EOM(self.total_stiffness)

        for element in self.elements:
            element.calculate_force_vector(self.q)
            element.calculate_stresses_and_strains()


    def plot_structure(self, displacement_magnitude : int, resolution : int, deflections : int = False, annotations : int = False, deflection_annotations : int = False, width : int = 15, height : int = 5, external_deflections : np.ndarray = None) -> None:
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

        fig, axes = plt.subplots(figsize=(width, height))

        for element in self.elements:
            element.plot_element(displacement_magnitude, resolution, axes, deflections, external_deflections)
            i += 1

        if deflections:
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))

            axes.legend(by_label.values(), by_label.keys())

        # ------------------------------- #
        #        plot the supports        #
        # ------------------------------- #
        for support in self.supports:
                pos = support.pos

                support_height = 0.8
                support_width = 0.6
                
                roller_radius = 0.1

                ground_height = Vec2(pos.x, pos.y)
                ground_dashes = 7
                
                if support.direction == Direction.HORIZONTAL.abs():
                    ground_height.x -= support.direction.x * support_height

                    # plot the triangle
                    axes.plot([pos.x, pos.x - support.direction.x * support_height], [pos.y, pos.y - support_width / 2], color='black')
                    axes.plot([pos.x, pos.x - support.direction.x * support_height], [pos.y, pos.y + support_width / 2], color='black')
                    axes.plot([pos.x - support.direction.x * support_height, pos.x - support.direction.x * support_height], [pos.y - support_width / 2, pos.y + support_width / 2], color='black')
                    
                    # add roller specific stuff
                    if type(support) == RollerSupport:
                        ground_height.x -= 2 * support.direction.x * roller_radius

                        r1 = patches.Circle((pos.x - support.direction.x * support_height - roller_radius, pos.y), roller_radius, edgecolor='black', facecolor='none')
                        r2 = patches.Circle((pos.x - support.direction.x * support_height - roller_radius, pos.y - 2 * roller_radius), roller_radius, edgecolor='black', facecolor='none')
                        r3 = patches.Circle((pos.x - support.direction.x * support_height - roller_radius, pos.y + 2 * roller_radius), roller_radius, edgecolor='black', facecolor='none')

                        plt.gca().add_patch(r1)
                        plt.gca().add_patch(r2)
                        plt.gca().add_patch(r3)

                    # plot ground plane
                    axes.plot([ground_height.x, ground_height.x], [pos.y - support_width, pos.y + support_width], color='black')


                elif support.direction == Direction.VERTICAL.abs():
                    ground_height.y -= support.direction.y * support_height

                    # plot the triangle
                    axes.plot([pos.x, pos.x - support_width / 2], [pos.y, pos.y - support.direction.y * support_height], color='black')
                    axes.plot([pos.x, pos.x + support_width / 2], [pos.y, pos.y - support.direction.y * support_height], color='black')
                    axes.plot([pos.x - support_width / 2, pos.x + support_width / 2], [pos.y - support.direction.y * support_height, pos.y - support.direction.y * support_height], color='black')

                    # add roller specific stuff
                    if type(support) == RollerSupport:
                        ground_height.y -= 2 * support.direction.y * roller_radius

                        r1 = patches.Circle((pos.x, pos.y - support.direction.y * support_height - roller_radius), roller_radius, edgecolor='black', facecolor='none')
                        r2 = patches.Circle((pos.x - 2 * roller_radius, pos.y - support.direction.y * support_height - roller_radius), roller_radius, edgecolor='black', facecolor='none')
                        r3 = patches.Circle((pos.x + 2 * roller_radius, pos.y - support.direction.y * support_height - roller_radius), roller_radius, edgecolor='black', facecolor='none')

                        plt.gca().add_patch(r1)
                        plt.gca().add_patch(r2)
                        plt.gca().add_patch(r3)

                    # plot ground plane
                    axes.plot([pos.x - support_width, pos.x + support_width], [ground_height.y, ground_height.y], color='black')

                    # plot ground "dashes"
                    for i in range(ground_dashes):
                        axes.plot([pos.x - support_width+0.1 + (2 * support_width / ground_dashes) * i, pos.x - support_width + (2 * support_width / ground_dashes) * i], [ground_height.y, ground_height.y - 0.1], color='black')

        # ------------------------------- #


        if annotations or deflection_annotations:
            font_size = 10

            arrow_space_from_element = Vec2(0.01, 0.05)
            arrow_len = 0.5
            
            kw = dict(arrowstyle="Simple, tail_width=0.5, head_width=4, head_length=8", color='k')

            

            for element in self.elements:
                mid_point = (element.nodes[1].pos + element.nodes[0].pos) / 2
                axes.annotate(f"E{element.id}", xy=(mid_point.x, mid_point.y), xytext=(-10, 10), textcoords='offset points', color='black', fontsize=font_size)
            
            if deflection_annotations:
                print("deflection annotations not currently supported")

            # ------------------------------- #
            # Plot degrees of freedom arrows  #
            # ------------------------------- #

            for i in range(len(self.global_nodes)):
                # p_dist = (self.global_nodes[i].pos - self.global_nodes[1-i].pos).normalize()

                p : Vec2 = self.global_nodes[i].pos

                x_init : float = p.x + arrow_space_from_element.x
                y_init : float = p.y + arrow_space_from_element.y

                if self.global_nodes[i].x.index >= 0:
                    if annotations:
                        x_dof = patches.FancyArrowPatch(
                            (x_init, y_init), 
                            (x_init + arrow_len, y_init),
                            **kw
                        )

                        axes.annotate("q"+str(self.global_nodes[i].x.index+1), xy=(x_init + arrow_len, y_init), xytext=(0, 0), textcoords='offset points', color='black', fontsize=font_size)

                        plt.gca().add_patch(x_dof)

                if self.global_nodes[i].y.index >= 0:
                    if annotations:
                        y_dof = patches.FancyArrowPatch(
                            (x_init, y_init), 
                            (x_init, y_init + arrow_len), 
                            **kw
                        )

                        axes.annotate("q"+str(self.global_nodes[i].y.index+1), xy=(x_init, y_init + arrow_len), xytext=(0, 0), textcoords='offset points', color='black', fontsize=font_size)

                        plt.gca().add_patch(y_dof)

                if self.global_nodes[i].moment.index >= 0:
                    if annotations:
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

                        
                        axes.annotate("q"+str(self.global_nodes[i].moment.index+1), xy=(end_x - 0.05, end_y), xytext=(-27, -2), textcoords='offset points', color='black', fontsize=font_size)

                        plt.gca().add_patch(arc)
                        plt.gca().add_patch(arrow)
                
            # ------------------------------- #

        axes.axis('off')
        axes.axis('equal')

        plt.show()