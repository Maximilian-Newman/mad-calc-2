import xfoil
from matplotlib import pyplot as plt
import math
import numpy as np

XFOIL_FAIL_TOLERANCE = 2
KINEMATIC_VISCOSITY = 1.42e-5 # sea level at 10 degrees C
N_CRIT = 7 # lower numbers -> more turbulent conditions
CL_MULTIPLIER = 1 / 1.5
CD_MULTIPLIER = 1.5
CM_MULTIPLIER = 1

def get_reynolds(chord, v):
    return int(chord * v / (KINEMATIC_VISCOSITY))

def airfoil_full_name(baseName, flap_x, flap_y, flap_ang):
    """
    flap_x = position of hinge along chord (0-1)
    flap_y = height of hinge between bottom and top surface of wing (0-1)
    """
    if flap_ang == 0:
        return baseName
    return baseName + "_f" + str(int(100 * flap_x)) + "_" + str(int(100 * flap_y)) + "_" + str(int(flap_ang))

def get_airfoil_cache_file_path(fullname, reynolds):
    return "xfoilcache/" + fullname + "_Re=" + str(int(reynolds)) + "_nCrit=" + str(int(N_CRIT)) + ".csv"

def xFoil_generate_csv_naca(nacaCode, reynolds, d=0.1):
    xf = xfoil.XFoil()
    xf.naca(nacaCode)
    return xFoil_generate_csv(xf, "NACA" + str(nacaCode), reynolds, d=d)

def xFoil_generate_csv_selig(name, reynolds, d=0.1):
    x = []
    y = []
    file = open("airfoils/" + name + ".dat", "r")
    for line in file.read().split("\n")[1:]:
        if line != "":
            line = line.split()
            x.append(line[0])
            y.append(line[1])
    file.close()
    
    xf = xfoil.XFoil()
    xf.airfoil = xfoil.Airfoil(np.array(x), np.array(y))
    return xFoil_generate_csv(xf, name, reynolds, d=d)

def xFoil_generate_csv(xf, name, reynolds, d=0.1):
    print(" [ running xFoil on airfoil:", name, " with Re =", reynolds, end=" ] ... ")
    xf.repanel(n_nodes=300)
    xf.Re = reynolds
    xf.max_iter = 500
    xf.mach = 0
    xf.n_crit = N_CRIT

    results = []

    xf.reset_bls()  
    AoA = 0 - d
    numFails = 0
    while numFails < XFOIL_FAIL_TOLERANCE:
        cl, cd, cm, cp = xf.a(AoA)
        failed = False
        if math.isnan(cl):
            cl = 0
            failed = True
        if math.isnan(cd):
            cd = 10
            failed = True
        if math.isnan(cm):
            cm = 0
            failed = True

        if failed: numFails += 1
        else: numFails = 0

        if not failed:
            results.append([AoA, cl, cd, cm])
            #print(AoA, cl, cd, cm)
        AoA -= d

    results = results[::-1]


    xf.reset_bls() 
    AoA = 0
    numFails = 0
    while numFails < XFOIL_FAIL_TOLERANCE:
        cl, cd, cm, cp = xf.a(AoA)
        failed = False
        if math.isnan(cl):
            cl = 0
            failed = True
        if math.isnan(cd):
            cd = 10
            failed = True
        if math.isnan(cm):
            cm = 0
            failed = True

        if failed: numFails += 1
        else: numFails = 0

        if not failed:
            results.append([AoA, cl, cd, cm])
            #print(AoA, cl, cd, cm)
        AoA += d
    
    file = open(get_airfoil_cache_file_path(name, reynolds), "w")
    for line in results:
        AoA, cl, cd, cm = line
        AoA = str(AoA)
        cl = str(cl)
        cd = str(cd)
        cm = str(cm)
        file.write(",".join([AoA, cl, cd, cm]) + "\n")
    file.close()

    print("done")
    return results


xFoil_generate_csv_selig(airfoil_full_name("NACA4412", 0.8, 1, -30), 200000)
