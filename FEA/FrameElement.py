import numpy as np

from FEA.Element import Element


class FrameElement(Element):
    """
    A frame element for FEA analysis.

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
    __to_local()
        Returns the local stiffness matrix for a frame element.
    __to_global()
        Returns the global stiffness matrix for a frame element.
    __to_structure()
        Returns the global stiffness matrix for a frame element.
    calculate_global_stiffness()
        Calculates the global stiffness matrix for a frame element.
    calculate_force_vector(q)
        Calculates the force vector for a frame element.
    get_lambda_mat(angle)
        calculates the lambda matrix
    """
     
    def __init__(self, assembly_mat : np.ndarray, E : int, I : int, L : int, A : int, angle : int, UDL : int = 0, point_load : int = 0):
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

        Returns
        -------
        None
        """
        
        super().__init__()
        
        self.E : int = E
        self.I : int = I
        self.L : int = L
        self.A : int = A
        self.angle : int = angle

        self.UDL : int = UDL
        self.point_load : int = point_load

        self.UDL_forces : np.ndarray = None
        self.point_load_forces : np.ndarray = None

        self.assembly_mat : np.ndarray = assembly_mat

        self.lambda_mat : np.ndarray = None

        self.local_stiffness : np.ndarray = None
        self.local_stiffness_hat : np.ndarray = None
        self.global_stiffness : np.ndarray = None

        self.element_deflections : np.ndarray = None
        self.structural_deflections : np.ndarray = None
        self.local_force : np.ndarray = None
        self.global_force : np.ndarray = None

    
    def __to_local(self) -> np.ndarray:
        """
        creates the local stiffness matrix for a frame element

        Returns
        -------
        np.ndarray
            The local stiffness matrix.
        """

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
    

    def __to_global(self) -> np.ndarray:
        """
        creates the global stiffness matrix for a frame element

        Returns
        -------
        np.ndarray
            The global stiffness matrix.
        """

        self.lambda_mat = FrameElement.get_lambda_mat(self.angle)

        self.local_stiffness_hat = self.lambda_mat.T @ self.local_stiffness @ self.lambda_mat

        return self.local_stiffness_hat


    def __to_structure(self) -> np.ndarray:
        """
        creates the structural stiffness matrix for a frame element

        Returns
        -------
        np.ndarray
            The structural stiffness matrix.
        """

        if(self.local_stiffness_hat is None):
            raise AssertionError("global_stiffness matrix cannot be None. First call __to_global().")

        self.global_stiffness = self.assembly_mat @ self.local_stiffness_hat @ self.assembly_mat.T

        return self.global_stiffness
    

    def _UDL_forces(self) -> None:
        """
        """

        f_eq = np.array([
            [0],
            [self.UDL * self.L / 2],
            [self.UDL * self.L ** 2 / 12],
            [0],
            [self.UDL * self.L / 2],
            [-self.UDL * self.L ** 2 / 12]
        ])

        F_eq = self.lambda_mat.T @ f_eq

        self.UDL_forces = self.assembly_mat.T @ F_eq


    def _point_load_forces(self) -> None:
        """
        """

        f_eq = np.array([
            [0],
            [self.point_load / 2],
            [self.point_load * self.L / 8],
            [0],
            [self.point_load / 2],
            [-self.point_load * self.L / 8]
        ])

        F_eq = self.lambda_mat.T @ f_eq

        self.point_load_forces = self.assembly_mat.T @ F_eq


    def calculate_global_stiffness(self) -> None:
        """
        Calculates the global stiffness matrix for a frame element.
        
        Returns
        -------
        None
        """

        self.__to_local()
        self.__to_global()
        self.__to_structure()

        self._UDL_forces()
        self._point_load_forces()


    def calculate_force_vector(self, q) -> None:
        """
        Calculates the force vector for a frame element.
        
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


    def plot_element(self, nodes : np.ndarray, displacement_magnitude : int, n_points : int, q : np.ndarray = None) -> np.ndarray:
        """
        """
        
        x = np.linspace(0, self.L, n_points)

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

        x_undeflected = np.linspace(nodes[0][0], nodes[1][0], n_points)
        y_undeflected = np.linspace(nodes[0][1], nodes[1][1], n_points)

        x_deflected = x_undeflected + deflections_XG * displacement_magnitude
        y_deflected = y_undeflected + deflections_YG * displacement_magnitude

        return np.array([[x_undeflected, y_undeflected], [x_deflected, y_deflected]])