import xfoil
from matplotlib import pyplot as plt
import math
import numpy as np
from pathlib import Path

# axis dirextions:
# positive x: towards tail
# positive y: starboard
# positive z: up


XFOIL_FAIL_TOLERANCE = 2
KINEMATIC_VISCOSITY = 1.42e-5 # sea level at 10 degrees C
AIR_DENSITY = 1.225
N_CRIT = 7 # lower numbers -> more turbulent conditions
CL_MULTIPLIER = 1 / 1.5
CD_MULTIPLIER = 1.5
CM_MULTIPLIER = 1

REYNOLDS_MAX_INTERPOLATION_DISTANCE = 100000

ARROW_SIZE_MULTIPLIER = 0.02

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
    print(" [ running xFoil on '", name, "' with Re =", reynolds, end=" ] ... ")
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


def choose_interpol_flap_angles(flapang, availableAngles):
    if len(availableAngles) == 0:
        print("ERROR: No coordinate file available for '", self.name, "'")
        return [None, None]

    if flapang in availableAngles:
        return [flapang, flapang]

    if flapang < availableAngles[0]:
        print("Warning: ", flapang, "outside of flap data range for airfoil")
        return [availableAngles[0], availableAngles[0]]

    if flapang > availableAngles[-1]:
        print("Warning: ", flapang, "outside of flap data range for airfoil")
        return [availableAngles[-1], availableAngles[-1]]

    for i in range(0, len(availableAngles)):
        if flapang < availableAngles[i]:
            return [availableAngles[i-1], availableAngles[i]]


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



def interpolate_xfoil(data1, data2, ratio):
    xfoilData = []
    i = 0
    j = 0
    while i < len(data1) and j < len(data2):
        point1 = data1[i]
        point2 = data2[j]
        if point1[0] == point2[0]:
            newPoint = [data1[i][0]]
            for k in range(1, 4):
                newPoint.append(point1[k] + ratio * (point2[k] - point1[k]))
            xfoilData.append(newPoint)
            i += 1
            j += 1

        elif point1[0] > point2[0]:
            j += 1
        else:
            i += 1
    
    return xfoilData

def get_xfoil_data_no_reynolds_interpolation(name, flapx, flapy, flapang, reynolds, availableAngles=None): # allow interpolation between flap angles
    if availableAngles == None:
        availableAngles = get_available_flapangles(name, flapx, flapy)
    
    a1, a2 = choose_interpol_flap_angles(flapang, availableAngles)
    if a1 == a2:
        return get_xfoil_data_no_interpolation(name, flapx, flapy, a1, reynolds)

    data1 = get_xfoil_data_no_interpolation(name, flapx, flapy, a1, reynolds)
    data2 = get_xfoil_data_no_interpolation(name, flapx, flapy, a2, reynolds)
    ratio = (flapang - a1) / (a2 - a1)
    return interpolate_xfoil(data1, data2, ratio)
    


def get_already_calculated_reynolds(name, flapx, flapy, flapang):
    folder = Path("./xfoilcache")
    files = [str(f) for f in folder.iterdir() if f.is_file()]
    available = []
    for f in files:
        if f.endswith("_nCrit=" + str(N_CRIT) + ".csv") and f.startswith("xfoilcache/" + airfoil_full_name(name, flapx, flapy, flapang) + "_Re="):
            r = f[f.find("Re=") + 3 : f.find("_nCrit=")]
            available.append(int(r))
    return sorted(available)


def get_xfoil_data_no_flap_interpolation(name, flapx, flapy, flapang, reynolds):  # interpolates between reynolds numbers if within tolerance
    available = get_already_calculated_reynolds(name, flapx, flapy, flapang)
    print(name, flapx, flapy, flapang, reynolds)
    
    if reynolds < available[0]:
        if available[0] - reynolds > REYNOLDS_MAX_INTERPOLATION_DISTANCE:
            return get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, reynolds)
        
        if available[0] - REYNOLDS_MAX_INTERPOLATION_DISTANCE < 0:
            get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, reynolds * 0.8)
        else:
            get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, available[0] - REYNOLDS_MAX_INTERPOLATION_DISTANCE)
        available = get_already_calculated_reynolds(name, flapx, flapy, flapang)
    
    if reynolds > available[-1]:
        if reynolds - available[-1] > REYNOLDS_MAX_INTERPOLATION_DISTANCE:
            return get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, reynolds)
        get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, available[-1] + REYNOLDS_MAX_INTERPOLATION_DISTANCE)
        available = get_already_calculated_reynolds(name, flapx, flapy, flapang)


    
    for i in range(0, len(available)):
        if reynolds == available[i]:
            return get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, reynolds)
        
        if reynolds < available[i]:
            r1 = available[i-1]
            r2 = available[i]

            if (r2 - r1) > REYNOLDS_MAX_INTERPOLATION_DISTANCE:
                
                if r2 - reynolds < REYNOLDS_MAX_INTERPOLATION_DISTANCE:
                    if r2 - REYNOLDS_MAX_INTERPOLATION_DISTANCE > 0:
                        get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, r2 - REYNOLDS_MAX_INTERPOLATION_DISTANCE)
                    else:
                        get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, reynolds * 0.8)
                    return get_xfoil_data_no_flap_interpolation(name, flapx, flapy, flapang, reynolds)

                if reynolds - r1 < REYNOLDS_MAX_INTERPOLATION_DISTANCE:
                    get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, r1 + REYNOLDS_MAX_INTERPOLATION_DISTANCE)
                    return get_xfoil_data_no_flap_interpolation(name, flapx, flapy, flapang, reynolds)
                
                return get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, reynolds)
                
            ratio = (reynolds - r1) / (r2 - r1)
            data1 = get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, r1)
            data2 = get_xfoil_data_no_interpolation(name, flapx, flapy, flapang, r2)
            return interpolate_xfoil(data1, data2, ratio)


def get_xfoil_data(name, flapx, flapy, flapang, reynolds, availableAngles=None):
    if availableAngles == None:
        availableAngles = get_available_flapangles(name, flapx, flapy)

    a1, a2 = choose_interpol_flap_angles(flapang, availableAngles)
    
    if a1 == a2:
        return get_xfoil_data_no_flap_interpolation(name, flapx, flapy, flapang, reynolds)

    data1 = get_xfoil_data_no_flap_interpolation(name, flapx, flapy, a1, reynolds)
    data2 = get_xfoil_data_no_flap_interpolation(name, flapx, flapy, a2, reynolds)
    ratio = (flapang - a1) / (a2 - a1)
    return interpolate_xfoil(data1, data2, ratio)


def rotate_point(x, y, angle_degrees):
    theta = math.radians(angle_degrees)
    c = math.cos(theta)
    s = math.sin(theta)
    
    x_new = x * c - y * s
    y_new = x * s + y * c
    return [x_new, y_new]

def scale_vector(vec, scale):
    return [i * scale for i in vec]

class AeroElement:
    def __init__(self, name, chord, incidence, width, flapx, flapy, flapang, x, y, z):
        self.incidence = incidence
        self.name = name
        self.chord = chord
        self.S = chord * width
        self.flapx = flapx
        self.flapy = flapy
        self.flapang = flapang
        self.position = [x, y, z]
        self.availableFlapAngles = get_available_flapangles(name, flapx, flapy)
        self.interpolFlapAngles = choose_interpol_flap_angles(flapang, self.availableFlapAngles)
        self.xfoilData = None
        self.lift = 0
        self.drag = 0
        self.moment = 0
        self.inverted = False
        self.stalled = True
        self.lastCalcModelIncidence = 0 # incidence of model when forces were last updated

    def nearest_flap_angle(self):
        if self.interpolFlapAngles == [None, None]:
            self.interpolFlapAngles = self.chooseInterpolFlapAngles()
        if abs(self.interpolFlapAngles[1] - self.flapang) < abs(self.interpolFlapAngles[0] - self.flapang):
            return self.interpolFlapAngles[1]
        return self.interpolFlapAngles[0]

    def display_outline(self, parentModel):
        fullName = airfoil_full_name(self.name, self.flapx, self.flapy, self.nearest_flap_angle())
        x, y, z = self.position
        if self.stalled:
            parentModel.diagram_draw_airfoil(fullName, x, y, z, self.chord, self.incidence, color="orange")
        else:
            parentModel.diagram_draw_airfoil(fullName, x, y, z, self.chord, self.incidence, color="blue")

    def display_forces(self, parentModel):
        x, y, z = self.position
        x += self.chord / 4
        lift, drag = self.get_force_vectors()
        u, v, w = scale_vector(lift, ARROW_SIZE_MULTIPLIER)
        parentModel.model_diagram.quiver(x, y, z, u, v, w, color="green", arrow_length_ratio=0.05, clip_on=False)
        u, v, w = scale_vector(drag, ARROW_SIZE_MULTIPLIER)
        parentModel.model_diagram.quiver(x, y, z, u, v, w, color="red", arrow_length_ratio=0.05, clip_on=False)

    def display(self, parentModel):
        self.display_outline(parentModel)
        self.display_forces(parentModel)

    def get_max_dimension(self):
        x, y, z = self.position
        return max([ abs(x), abs(y), abs(z), abs(x + self.chord) ])

    def choose_interpol_AoA_indices(self, AoA):
        if AoA < self.xfoilData[0][0]:
            return [None, None]
        if AoA > self.xfoilData[-1][0]:
            return [None, None]
        for i in range(0, len(self.xfoilData)):
            a = self.xfoilData[i][0]
            if a == AoA:
                return [i, i]
            if a > AoA:
                return [i-1, i]
            

    def update_forces(self, modelIncidence, airspeed):
        self.lastCalcModelIncidence = modelIncidence
        self.stalled = False
        AoA = modelIncidence + self.incidence
        if self.inverted: AoA = -AoA

        i1, i2 = self.choose_interpol_AoA_indices(AoA)
        
        if i1 == None: # stall
            self.lift = 0
            self.moment = 0
            self.drag = 0 # use largest cd in whole polar
            for dataPoint in self.xfoilData:
                d = 0.5 * AIR_DENSITY * airspeed**2 * self.S * dataPoint[2]
                if d > self.drag:
                    self.drag = d
            
            self.stalled = True
            return
            
        dataPoint = [AoA, None, None, None]
        
        if i1 == i2:
            dataPoint = self.xfoilData[i1]
        else:
            data1 = self.xfoilData[i1]
            data2 = self.xfoilData[i2]
            ratio = (AoA - data1[0]) / (data2[0] - data1[0])
            
            for i in range(1, 4):
                dataPoint[i] = data1[i] + ratio * (data2[i] - data1[i])

        self.lift = 0.5 * AIR_DENSITY * airspeed**2 * self.S * dataPoint[1]
        self.drag = 0.5 * AIR_DENSITY * airspeed**2 * self.S * dataPoint[2]
        self.moment = 0.5 * AIR_DENSITY * airspeed**2 * self.S * self.chord * dataPoint[3]

        if self.inverted:
            self.lift = -self.lift

    def get_force_vectors(self):
        lift = rotate_point(0, self.lift, self.lastCalcModelIncidence)
        drag = rotate_point(self.drag, 0, self.lastCalcModelIncidence)

        lift = [lift[0], 0, lift[1]]
        drag = [drag[0], 0, drag[1]]
        return [lift, drag]

    def update_xfoil_data_no_interpolation(self, airspeed):
        reynolds = get_reynolds(self.chord, airspeed)
        self.xfoilData = get_xfoil_data_no_interpolation(self.name, self.flapx, self.flapy, self.flapang, reynolds)

    def update_xfoil_data_no_reynolds_interpolation(self, airspeed):
        reynolds = get_reynolds(self.chord, airspeed)
        self.xfoilData = get_xfoil_data_no_reynolds_interpolation(self.name, self.flapx, self.flapy, self.flapang, reynolds, self.availableFlapAngles)

    def update_xfoil_data_no_flap_interpolation(self, airspeed):
        reynolds = get_reynolds(self.chord, airspeed)
        self.xfoilData = get_xfoil_data_no_flap_interpolation(self.name, self.flapx, self.flapy, self.flapang, reynolds)

    def update_xfoil_data(self, airspeed):
        reynolds = get_reynolds(self.chord, airspeed)
        self.xfoilData = get_xfoil_data(self.name, self.flapx, self.flapy, self.flapang, reynolds, self.availableFlapAngles)
        
        


class StaticModel:
    def __init__(self):
        self.model_diagram_fig = plt.figure(figsize=(20,9), layout="constrained")
        self.model_diagram = self.model_diagram_fig.add_subplot(projection = "3d")
        self.aeroElements = []
        self.massElements = []
        self.clear_model_diagram()
        

    def clear_model_diagram(self):
        self.model_diagram.cla()
        self.model_diagram.set_axis_off()
        #self.model_diagram.set_aspect("equal")
        r = self.get_max_dimension() * 0.8
        self.model_diagram.set_xlim3d(-r, r)
        self.model_diagram.set_ylim3d(-r, r)
        self.model_diagram.set_zlim3d(-r, r)
        


    def diagram_draw_airfoil(self, name, x, y, z, chord, incidence=0, color="blue"):
        x_data, z_data = get_selig_data(name)
        y_data = []
        for i in range(0, len(z_data)):
            x_data[i] = x + chord * x_data[i]
            z_data[i] = z + chord * z_data[i]

            x_data[i], z_data[i] = rotate_point(x_data[i], z_data[i], -incidence)
            y_data.append(y)

        self.model_diagram.plot(x_data, y_data, z_data, c=color, clip_on=False)

    def display_airfoil_outlines(self):
        for a in self.aeroElements:
            a.display_outline(self)
        plt.show(block=False)
    
    def display_airfoils(self):
        for a in self.aeroElements:
            a.display(self)
        plt.show(block=False)

    def get_max_dimension(self):
        maxDim = 0
        for element in self.aeroElements:
            r = element.get_max_dimension()
            if r > maxDim:
                maxDim = r

        if maxDim == 0:
            return 1
        return maxDim

    def update_xfoil_data_no_interpolation(self, airspeed):
        for element in self.aeroElements:
            element.update_xfoil_data_no_interpolation(airspeed)

    def update_xfoil_data_no_reynolds_interpolation(self, airspeed): # interpolates for flap angle only
        for element in self.aeroElements:
            element.update_xfoil_data_no_reynolds_interpolation(airspeed)

    def update_xfoil_data_no_flap_interpolation(self, airspeed): # interpolates reynolds number only
        for element in self.aeroElements:
            element.update_xfoil_data_no_flap_interpolation(airspeed)

    def update_xfoil_data(self, airspeed): # interpolates for flaps and reynolds number
        for element in self.aeroElements:
            element.update_xfoil_data(airspeed)

    def update_forces(self, incidence, airspeed):
        for element in self.aeroElements:
            element.update_forces(incidence, airspeed)

    







# testing plots:
test = StaticModel()

width = 0.2

dihedral = 2 / 100
test.aeroElements = []
for i in range(0, 21):
    test.aeroElements.append(AeroElement("NACA4412", 2 * (1-0.02*i), i/2, width, 0.8, 0.5, -i, 0, i*width, i*dihedral))
    test.aeroElements.append(AeroElement("NACA4412", 2 * (1-0.02*i), i/2, width, 0.8, 0.5, -20, 0, -i*width, i*dihedral))

test.clear_model_diagram()
test.display_airfoils()
plt.pause(1)


airspeed = 20
test.update_xfoil_data(airspeed)
test.update_forces(5, airspeed)
test.clear_model_diagram()
test.display_airfoils()
plt.show(block=False)
