import xfoil
from matplotlib import pyplot as plt
import math
import numpy as np

# axis dirextions:
# positive x: towards tail
# positive y: starboard
# positive z: up


XFOIL_FAIL_TOLERANCE = 2
KINEMATIC_VISCOSITY = 1.42e-5 # sea level at 10 degrees C
N_CRIT = 7 # lower numbers -> more turbulent conditions
CL_MULTIPLIER = 1 / 1.5
CD_MULTIPLIER = 1.5
CM_MULTIPLIER = 1

REYNOLDS_MAX_INTERPOLATION_DISTANCE = 10000

def get_reynolds(chord, v):
    return int(chord * v / (KINEMATIC_VISCOSITY))

def airfoil_full_name(baseName, flapx, flapy, flapang):
    """
    flapx = position of hinge along chord (0-1)
    flapy = height of hinge between bottom and top surface of wing (0-1)
    """
    if flapang == 0:
        return baseName
    return baseName + "_f" + str(int(100 * flapx)) + "_" + str(int(100 * flapy)) + "_" + str(int(flapang))

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


def get_available_flapangles(name, flapx, flapy): # returns list of angles for which selig .dat coordinate file was provided
    angles = []
    for i in range(-90, 90):
        try:
            file = open("airfoils/" + airfoil_full_name(name, flapx, flapy, i) + ".dat", "r")
            file.close()
            angles.append(i)
        except:
            pass
    return angles


def get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, reynolds):
    name = airfoil_full_name(name, flapx, flapy, flapang)

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


def rotate_point(x, y, angle_degrees):
    theta = math.radians(angle_degrees)
    c = math.cos(theta)
    s = math.sin(theta)
    
    x_new = x * c + y * s
    y_new = -x * s + y * c
    return x_new, y_new

class AeroElement:
    def __init__(self, name, chord, incidence, width, flapx, flapy, flapang, x, y, z):
        self.incidence = incidence
        self.name = name
        self.chord = chord
        self.A = chord * width
        self.flapx = flapx
        self.flapy = flapy
        self.flapang = flapang
        self.position = [x, y, z]
        self.availableFlapAngles = get_available_flapangles(name, flapx, flapy)
        self.interpolFlapAngles = self.choose_interpol_flap_angles()
        self.xfoilData = None

    def choose_interpol_flap_angles(self):
        if len(self.availableFlapAngles) == 0:
            print("ERROR: No coordinate file available for '", self.name, "'")
            return [None, None]

        if self.flapang in self.availableFlapAngles:
            return [self.flapang, self.flapang]

        if self.flapang < self.availableFlapAngles[0]:
            print("Warning: ", self.flapang, "outside of flap data range for airfoil '", self.name, "'")
            return [self.availableFlapAngles[0], self.availableFlapAngles[0]]

        if self.flapang > self.availableFlapAngles[-1]:
            print("Warning: ", self.flapang, "outside of flap data range for airfoil '", self.name, "'")
            return [self.availableFlapAngles[-1], self.availableFlapAngles[-1]]

        for i in range(0, len(self.availableFlapAngles)):
            if self.flapang < self.availableFlapAngles[i]:
                return [self.availableFlapAngles[i-1], self.availableFlapAngles[i]]

    def nearest_flap_angle(self):
        if self.interpolFlapAngles == [None, None]:
            self.interpolFlapAngles = self.chooseInterpolFlapAngles()
        if abs(self.interpolFlapAngles[1] - self.flapang) < abs(self.interpolFlapAngles[0] - self.flapang):
            return self.interpolFlapAngles[1]
        return self.interpolFlapAngles[0]

    def display_outline(self, parentModel):
        fullName = airfoil_full_name(self.name, self.flapx, self.flapy, self.nearest_flap_angle())
        x, y, z = self.position
        parentModel.diagram_draw_airfoil(fullName, x, y, z, self.chord, self.incidence)
        
        


class StaticModel:
    def __init__(self):
        self.model_diagram_fig = plt.figure(figsize=(20,9), layout="constrained")
        self.model_diagram = self.model_diagram_fig.add_subplot(projection = "3d")
        
        self.clear_model_diagram()
        self.aeroElements = []
        self.massElements = []
        

    def clear_model_diagram(self):
        self.model_diagram.cla()
        self.model_diagram.set_axis_off()
        self.model_diagram.set_aspect("equal")
        #if self.model_diagram_fig.canvas is not None:
        #    self.model_diagram_fig.canvas.draw_idle()
        #    self.model_diagram_fig.canvas.flush_events()
        


    def diagram_draw_airfoil(self, name, x, y, z, chord, incidence=0):
        x_data, z_data = get_selig_data(name)
        y_data = []
        for i in range(0, len(z_data)):
            x_data[i] = x + chord * x_data[i]
            z_data[i] = z + chord * z_data[i]

            x_data[i], z_data[i] = rotate_point(x_data[i], z_data[i], incidence)
            y_data.append(y)

        self.model_diagram.plot(x_data, y_data, z_data, c="blue", clip_on=False)

    def display_airfoil_outlines(self):
        for a in self.aeroElements:
            a.display_outline(self)

    


# testing plots:
test = StaticModel()

width = 0.2

#tried cool animation, matplotlib doesn't seem to like it (resetting zoom between each frame)
for a in range(0, 10):
    for dihedral in range(20, 30):
        dihedral = dihedral / 1000
        test.aeroElements = []
        for i in range(0, 21):
            test.aeroElements.append(AeroElement("NACA4412", 2, i/3, width, 0.8, 0.5, -i, 0, i*width, i*dihedral))
            test.aeroElements.append(AeroElement("NACA4412", 2, i/3, width, 0.8, 0.5, -i, 0, -i*width, i*dihedral))
        
        test.clear_model_diagram()
        test.display_airfoil_outlines()
        #plt.show()
        plt.show(block=False)
        #test.model_diagram_fig.canvas.draw()
        plt.pause(0.1)

test.display_airfoil_outlines()
plt.show(block=False)
