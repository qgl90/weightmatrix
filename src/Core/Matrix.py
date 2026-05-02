import numpy as np 
#from typing import Self
from .Precision import array_d
import itertools
from typing import TypeVar

TMatrix = TypeVar("TMatrix", bound =  "Matrix")
#### Matrix 5x5 class representing covariance 
class Matrix :     
    """
    A class used to represent a Covariance Matrix 

    Attributes
    ----------
    matrix : np.array( (5,5) ) 
        a 5x5 matrix
    Methods
    -------
    says(sound=None)
        Prints the animals name and what sound it makes
    """ 
    #2D array 5x5 matrix double-precision-number    
    def __init__(self, matrix = np.zeros( (5,5), "double")):    
        """constructor passing a 5x5 matrix, by default a 0-matrix set"""
        self.matrix  = matrix
    def Transp( self)  : 
        """
        Transpose the matrix
        """
        return Matrix(self.matrix.transpose())        
    def Prod( col : np.array ) -> np.array :
        """
        Product matrix with array M * vector = vector
        """        
        return np.matmul(np.matrix, col)
    def Inverse( self )   :
        """Compute a new Matrix inverting it M^{-1}*M = Unity. Return M^{-1}"""
        try : 
            #if cannot invert, return a pseudo inverse
            rr =  Matrix( np.linalg.inv(self.matrix))
        except :
            rr =  Matrix( np.linalg.pinv(self.matrix))            
        return rr
    def AddMS_TX(self,  th2 : float) -> None :
        """
        Add Multiple-Scattering on xz plane effect, where th2 is the m-scattering term squared 
        scat2 = math.pow(14.*qop,2)* usePlanes[idx].Thickness  * math.sqrt( 1. + Tx_At_Plane**2 + Ty_At_Plane**2 )
        """
        for i, j in itertools.product( range(5), range(5)):
            if( i == 2 or j ==2 or i <j ): continue #skip the matrix[2][2] which is the tx , tx diagonal term      
            self.matrix[i][j]+= th2 * ( self.matrix[i][j] * self.matrix[2][2]  -  self.matrix[i][2] * self.matrix[j][2])                                        
        fact = 1./(1.+th2*self.matrix[2][2]);    
        for i, j in itertools.product( range(5), range(5)):
            if i < j : continue 
            self.matrix[i][j] = self.matrix[i][j]*fact        
        for i, j in itertools.product( range(5), range(5)):
            if i < j : continue
            self.matrix[j][i]=self.matrix[i][j]
        return
    def AddMS_TY( self, th2 : float) -> None :
        """
        Add Multiple-Scattering on yz plane effect, where th2 is the m-scattering term squared 
        scat2 = math.pow(14.*qop,2)* usePlanes[idx].Thickness  * math.sqrt( 1. + Tx_At_Plane**2 + Ty_At_Plane**2 )
        """
        for i, j in itertools.product( range(5), range(5)):            
            if i ==3 or j == 3 or i < j  : continue #skip the matrix[3][3] which is the ty , ty diagonal term 
            if i < j          : continue            
            self.matrix[i][j] = self.matrix[i][j]+ th2 * ( self.matrix[i][j] * self.matrix[3][3]  -  self.matrix[i][3] * self.matrix[j][3])
        fact = 1./(1.+th2*self.matrix[3][3]);
        for i, j in itertools.product( range(5), range(5)):
            if i < j : continue 
            self.matrix[i][j]*= fact        
        for i, j in itertools.product( range(5), range(5)):
            if i < j : continue
            self.matrix[j][i]=self.matrix[i][j] 
    def AddMS(self, th2 : float ) -> None :
        ## TODO : add the tx, ty argument 
        ## th2 = math.pow(14.*qop,2)* usePlanes[idx].Thickness  * math.sqrt( 1. + Tx_At_Plane**2 + Ty_At_Plane**2 ) ?
        """
        zx + zy plane addition of multiple scattering to *this 
        From Pierre Mail    
                I have the formulae for the mult. scatt. contribution for non-negligible tx and ty. 
                First, you have to rad length by sqrt(1+tx^2+ty^2) if the material plane is perp. to z (or, for any orientation, by 1/sin(a), 
                where a is the angle between the track and the plane). 
                Then you add to the error matrix the following terms:
                on (tx,tx) : (1+tx^2+ty^2)* (1+tx^2) * s2
                on (ty,ty) : (1+tx^2+ty^2)* (1+ty^2) * s2
                on (tx,ty) : (1+tx^2+ty^2)* tx * ty  * s2
                where s2 = (14/p)^2.Lrad, where Lrad is the mumber of rad. lengths, accounting for the factor above.
                If you neglect the correlation terms, you van use the same functions as defore. Else you can code the matrix operation (1+WC)^(-1).W with the functions in Matrix5.cc
        """   
        self.AddMS_TX(self, th2, tx_state, ty_state)
        self.AddMS_TY(self, th2, tx_state, ty_state)
    def Scaled( self , scale : float = 1.) : 
        return Matrix( self.matrix * scale)
    def Propagate( self , delta : float, curv : float) -> None :
        """
        Propagate This Matrix over a delta, curvature
        """
        d13 = -delta
        d24 = -delta
        d15 = curv * delta*delta/2.
        d35 = -curv *delta        
        der = Matrix( array_d( [ [ 1., 0. ,  d13, 0., d15],\
                                 [ 0., 1. ,  0., d24, 0. ],\
                                 [ 0., 0. ,  0.,  0., d35],\
                                 [ 0., 0. ,  0.,  1.,  0.],\
                                 [ 0., 0. ,  0.,  0.,  1.]]))
        self.matrix =  (der.Transp()  * (self.matrix)) * der 
        return    
    def __add__(self, other : TMatrix)   -> TMatrix:
        """add 2 matrix objects"""
        return Matrix( self.matrix + other.matrix)
    def __sub__(self, other: TMatrix)    -> TMatrix:
        """subtract 2 matrix objects"""
        return Matrix(self.matrix - other.matrix)
    def __mul__(self, other : TMatrix) :
        """multiply 2 matrix objects"""        
        return np.matmul(self.matrix, other.matrix)
    #multiply by a single scalar
    def __mul__(self, other : float) :
        """multiply matrix by float number"""
        return Matrix( self.matrix * other)        
    def __mul__(self, other : TMatrix) -> TMatrix:
        """multiply matrix with matrix"""
        return Matrix(np.matmul(self.matrix, other.matrix))

    def __str__(self):
        return f"{self.matrix}"