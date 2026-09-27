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

REYNOLDS_MAX_INTERPOLATION_DISTANCE = 10000

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

def get_selig_data(name):
    x = []
    y = []
    file = open("airfoils/" + name + ".dat", "r")
    for line in file.read().split("\n")[1:]:
        if line != "":
            line = line.split()
            x.append(float(line[0]))
            y.append(float(line[1]))
    file.close()
    return [x, y]
    
def xFoil_generate_csv_selig(name, reynolds, d=0.1):
    x, y = get_selig_data(name)
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
        if math.isnan(cl) or math.isnan(cd) or math.isnan(cm):
            numFails += 1
        else:
            numFails = 0
            results.append([AoA, cl, cd, cm])
        AoA -= d

    results = results[::-1]


    xf.reset_bls() 
    AoA = 0
    numFails = 0
    while numFails < XFOIL_FAIL_TOLERANCE:
        cl, cd, cm, cp = xf.a(AoA)
        
        if math.isnan(cl) or math.isnan(cd) or math.isnan(cm):
            numFails += 1
        else:
            numFails = 0
            results.append([AoA, cl, cd, cm])
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


def get_available_flap_angles(name, flapx, flapy): # returns list of angles for which selig .dat coordinate file was provided
    angles = []
    for i in range(-90, 90):
        try:
            file = open("airfoils/" + airfoil_full_name(name, flapx, flapy, i) + ".dat", "r")
            file.close()
            angles.append(i)
        except:
            pass
    return angles


def get_xfoil_data_no_interpolation(name, flap_x, flap_y, flap_ang, reynolds):
    name = airfoil_full_name(name, flap_x, flap_y, flap_ang)

    try:
        data = []
        file = open(get_airfoil_cache_file_path(name, reynolds), "r")
        for line in file.read().split("\n"):
            if line != "":
                line = line.split(",")
                for i in range(0, 4):
                    line[i] = float(line[i])
                data.append(line)
        file.close()
        return data
    
    except:
        return xFoil_generate_csv_selig(name, reynolds)

print(get_xfoil_data_no_interpolation("NACA4412", 0.8, 0.5, -20, 200000))
print(get_xfoil_data_no_interpolation("NACA4412", 0.8, 0.5, 0, 200000))
print(get_available_flap_angles("NACA4412", 0.8, 0.5))
print(get_available_flap_angles("NACA4412", 0.8, 0.4))
print(get_available_flap_angles("NACA4412", 0.9, 0.5))
print(get_available_flap_angles("NACA4419", 0.8, 0.5))

model_diagram_fig = plt.figure()
model_diagram = model_diagram_fig.add_subplot(projection = "3d")

def clear_model_diagram():
    model_diagram.cla()
    model_diagram.set_axis_off()
    model_diagram.set_aspect("equal")

clear_model_diagram()

def model_diagram_add_airfoil(name, x, y, z, chord):
    y_data, z_data = get_selig_data(name)
    x_data = []
    for i in range(0, len(z_data)):
        y_data[i] = y + chord * y_data[i]
        z_data[i] = z + chord * z_data[i]
        x_data.append(x)
    
    model_diagram.plot(x_data, y_data, z_data, c="blue")


# testing plots:
"""
for i in range(0, 50):
    model_diagram.scatter([i], [i], [i])
    plt.show(block=False)
    plt.pause(1)
    if (i % 10 == 0):
        clear_model_diagram()
"""

model_diagram_add_airfoil(airfoil_full_name("NACA4412", 0.8, 0.5, -20), 0, 0, 0, 1)
model_diagram_add_airfoil(airfoil_full_name("NACA4412", 0.8, 1, 0), 1, 0, 0, 1)
model_diagram_add_airfoil(airfoil_full_name("NACA4412", 0.8, 1, 0), 2, -0.5, 0, 1.5)
model_diagram_add_airfoil(airfoil_full_name("NACA4412", 0.8, 0.5, -20), 0, 0, 1, 1)
plt.show(block=False)
