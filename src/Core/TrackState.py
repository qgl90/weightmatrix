import numpy as np
from math import cosh  
import numpy as np
from termcolor import colored
def PT(p, eta):
    return p / cosh(eta)
class TrackState :
    """
    Represents the state of a track at a given position in LHCb coordinates
    X,Y,Z 
    TX = dx/dz 
    TY = dy/dz 
    QOP = q/p of the track
    """    
    def __init__(self, x : float , y : float , z :float ,\
                 tx :float = None  , ty  : float = None , qop    : float = None ,\
                 eta:float = None  , phi : float = None , charge :float = None  , momentum : float = None 
                 ) -> None : 
        self.z        = z
        self.hasState = False
        if tx != None and ty != None and qop != None : 
            self.state = np.array( [ x,y,tx,ty,qop], dtype=float )
            self.hasState = True
        else :              
            self.hasState = False
            raise ValueError("Invalid state, fix me")        
    def X(self)  : return self.state[0]
    def Y(self)  : return self.state[1]
    def TX(self) : return self.state[2]
    def TY(self) : return self.state[3]
    def QOP(self): return self.state[4]
    def Z(self)  : return self.z    
    def arrState(self) : return self.state
    def __repr__(self):
        if not(self.hasState) : return colored(f"INVALID, NO STATE at Z = {self.z}", "red")
        else :                  return colored(f"<TrackState Z = {self.z},\n \t (x,y,tx,ty,qop)= {self.state[0]},{self.state[1]},{self.state[2]},{self.state[3]},{self.state[4]}>", "yellow")
    def __str__(self):
        if not(self.hasState) : return colored(f"INVALID STATE Z = {self.z}", "red")
        else :                  return colored(f"<TrackState Z = {self.z},\n \t (x,y,tx,ty,qop)= {self.state[0]},{self.state[1]},{self.state[2]},{self.state[3]},{self.state[4]}>", "yellow")

    def xyz(self) : return np.array([ self.state[0], self.state[1], self.z])