
from hepunits.units import centimeter2, gram
from Utils.Logger import Logger
"""
A simple python script with definitions of materials used in the code 
"""
rad_lengths = {    
    "air" : { 
        "radlen"    :  36.62     * gram /centimeter2 , 
        "radlen_mm" :  3.039E+05
    },
    "He" : {
        "radlen" : 94.32* gram /centimeter2 ,
        "radlen_mm" : 56.71E+05
    }
}
def radlen( material : str ) -> float :
    if material not in rad_lengths.keys():
        Logger.error(f"'{material}' not in available rad_lengths dictionary, check materials.py file and add it")
        raise Exception(f"'{material}' not in available rad_lengths dictionary, check materials.py file and add it")
    return rad_lengths[material]["radlen_mm"]