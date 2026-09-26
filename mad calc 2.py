import xfoil
from matplotlib import pyplot as plt
import math

XFOIL_FAIL_TOLERANCE = 2
KINEMATIC_VISCOSITY = 1.42e-5 # sea level at 10 degrees C

def get_reynolds(chord, v):
    return int(chord * v / (KINEMATIC_VISCOSITY))


def xFoil_generate_csv_naca(nacaCode, reynolds, d=0.1):
    xf = xfoil.XFoil()
    xf.naca(nacaCode)
    return xFoil_generate_csv(xf, "NACA" + str(nacaCode), reynolds, d=d)

def get_airfoil_cache_file_path(name, reynolds):
    return "xfoilcache/" + name + "-Re=" + str(int(reynolds)) + ".csv"

def xFoil_generate_csv(xf, name, reynolds, d=0.1):
    print(" [ running xFoil on airfoil:", name, " with Re =", reynolds, end=" ] ... ")
    xf.Re = reynolds
    xf.max_iter = 500
    xf.mach = 0

    results = []

    
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


xFoil_generate_csv_naca(4412, 200000)
