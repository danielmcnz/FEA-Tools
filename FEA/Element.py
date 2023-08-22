import numpy as np
import matplotlib.pyplot as plt


class Element():
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
     
    def __init__(self, assembly_mat : np.ndarray, E : int, I : int, L : int, A : int, angle : int, UDL : tuple = (0, 0, 1, 1), LVL : tuple = (0, 0, 1, 1), point_load : tuple = (0, 0, 0, 1, 1)):
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
        
        self.E : int = E
        self.I : int = I
        self.L : int = L
        self.A : int = A
        self.angle : int = angle

        self.UDL : int = UDL
        self.LVL : int = LVL
        self.point_load : tuple = point_load

        self.UDL_forces : np.ndarray        = None
        self.UDL_f_shear : np.ndarray       = None
        self.UDL_f_axial : np.ndarray       = None
        self.UDL_F_shear : np.ndarray       = None
        self.UDL_F_axial : np.ndarray       = None

        self.LVL_forces : np.ndarray        = None
        self.LVL_f_shear : np.ndarray       = None
        self.LVL_f_axial : np.ndarray       = None
        self.LVL_F_axial : np.ndarray       = None
        self.LVL_F_shear : np.ndarray       = None

        self.point_load_forces : np.ndarray = None
        self.PL_f_shear : np.ndarray        = None
        self.PL_f_axial : np.ndarray        = None
        self.PL_F_shear : np.ndarray        = None
        self.PL_F_axial : np.ndarray        = None

        self.assembly_mat : np.ndarray = assembly_mat

        self.lambda_mat : np.ndarray            = None

        self.local_stiffness : np.ndarray       = None
        self.local_stiffness_hat : np.ndarray   = None
        self.global_stiffness : np.ndarray      = None
        self.strain : np.ndarray                = None

        self.local_force : np.ndarray           = None
        self.global_force : np.ndarray          = None

        self.element_deflections : np.ndarray       = None
        self.structural_deflections : np.ndarray    = None

    
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

        if(self.global_stiffness is None):
            raise AssertionError("global_stiffness matrix cannot be None. First call calculate_global_stiffness().")

        self.structural_deflections = self.assembly_mat.T @ q

        self.element_deflections = self.lambda_mat @ self.structural_deflections

        self.local_force = self.local_stiffness @ self.element_deflections
        
        self.strain = (self.element_deflections[3] - self.element_deflections[0]) / self.L

        self.global_force = self.local_stiffness_hat @ self.structural_deflections

    
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
    

    def calculate_deflections(self, x) -> np.ndarray:
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

        phi_1 = (1 - x / self.L)
        phi_2 = x / self.L

        N_1 = 1 - 3 * (x ** 2) / (self.L ** 2) + 2 * (x ** 3) / (self.L ** 3)
        N_2 = (x ** 3) / self.L ** 2 - 2 * (x ** 2) / self.L + x
        N_3 = 3 * (x ** 2) / (self.L ** 2) - 2 * (x ** 3) / (self.L ** 3)
        N_4 = (x ** 3) / (self.L ** 2) - (x ** 2) / self.L

        element_axial_displacement = phi_1 * self.element_deflections[0] + phi_2 * self.element_deflections[3]
        element_transverse_displacement = N_1 * self.element_deflections[1] + N_2 * self.element_deflections[2] + N_3 * self.element_deflections[4] + N_4 * self.element_deflections[5]

        deflections_XG = element_axial_displacement * np.cos(np.deg2rad(self.angle)) - element_transverse_displacement * np.sin(np.deg2rad(self.angle))
        deflections_YG = element_axial_displacement * np.sin(np.deg2rad(self.angle)) + element_transverse_displacement * np.cos(np.deg2rad(self.angle))

        return np.array([deflections_XG, deflections_YG])


    def plot_element(self, nodes : np.ndarray, displacement_magnitude : int, resolution : int) -> None:
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

        x_undeflected = np.linspace(nodes[0][0], nodes[1][0], resolution)
        y_undeflected = np.linspace(nodes[0][1], nodes[1][1], resolution)

        x_deflected = x_undeflected + deflections_XG * displacement_magnitude
        y_deflected = y_undeflected + deflections_YG * displacement_magnitude

        plt.plot(x_undeflected, y_undeflected, 'b', label="Undeflected")
        plt.plot(x_deflected, y_deflected, 'r', label="Deflected")
