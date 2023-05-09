import math
import numpy as np
from numpy.linalg import norm
import vg
import os
 
# This code attempts to implement the flexible water model
# Potential is currently wrong. Need to fix again.

def calculate_molecule_properties(coordinates):
    O1 = coordinates[0]
    H11 = coordinates[1]
    H21 = coordinates[2]
    O2 = coordinates[3]
    H12 = coordinates[4]
    H22 = coordinates[5]

    r11 = np.sqrt(np.sum((O1-H11)**2, axis=0))
    r21 = np.sqrt(np.sum((O1-H21)**2, axis=0))
    r12 = np.sqrt(np.sum((O2-H12)**2, axis=0))
    r22 = np.sqrt(np.sum((O2-H22)**2, axis=0))

    OH11 = O1 - H11
    OH21 = O1 - H21
    OH12 = O2 - H12
    OH22 = O2 - H22
    angle1 = vg.angle(OH11, OH21)
    angle2 = vg.angle(OH12, OH22)

    rO = np.sqrt(np.sum((O1-O2)**2, axis=0))

    print("r1 (O1-H11):", np.round(r11, 8), "Å")
    print("r2 (O1-H12):", np.round(r21, 8), "Å")
    print("r3 (O2-H21):", np.round(r12, 8), "Å")
    print("r4 (O2-H22):", np.round(r22, 8), "Å")
    print("angle1:", np.round(angle1, 8), "°")
    print("angle2:", np.round(angle2, 8), "°")
    print("O-O:", np.round(rO, 8), "Å")