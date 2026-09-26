import xfoil
from matplotlib import pyplot as plt
import math

XFOIL_FAIL_TOLERANCE = 2

def xFoil_generate_csv_naca(nacaCode, path, reynolds, d=0.1):
    xf = xfoil.XFoil()
    xf.Re = reynolds
    xf.max_iter = 500
    xf.mach = 0
    xf.naca(nacaCode)

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
            print(AoA, cl, cd, cm)
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
            print(AoA, cl, cd, cm)
        AoA += d
    
    file = open(path, "w")
    for line in results:
        AoA, cl, cd, cm = line
        AoA = str(AoA)
        cl = str(cl)
        cd = str(cd)
        cm = str(cm)
        file.write(",".join([AoA, cl, cd, cm]) + "\n")
    file.close()


print("start")
xFoil_generate_csv_naca(4412, "xfoilcache/4412.csv", 200000)
print("end")
