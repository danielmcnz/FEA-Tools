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
     
    def __init__(self, assembly_mat : np.ndarray, E : int, I : int, L : int, A : int, angle : int):
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

        self.assembly_mat : np.ndarray = assembly_mat

        self.lambda_mat : np.ndarray = None

        self.local_stiffness : np.ndarray = None
        self.local_stiffness_hat : np.ndarray = None
        self.global_stiffness : np.ndarray = None

        self.element_deflection : np.ndarray = None
        self.structural_deflection : np.ndarray = None
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

        self.structural_deflection = self.assembly_mat.T @ q

        self.element_deflection = self.lambda_mat @ self.structural_deflection

        self.local_force = self.local_stiffness @ self.element_deflection

        self.global_force = self.local_stiffness_hat @ self.structural_deflection

    
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