import prettytable, math, statistics, itertools
from statistics import mean
from math import sin, cos, tan, pi, sqrt, asin, acos, atan
from fractions import gcd
from operator import mul, add, sub, truediv
from functools import partial


# ----------------- Global Variables -------------------

bar_dia_list = {3: 0.375, 4: 0.5, 5: 0.625, 6: 0.75, 7: 0.875, 8: 1,
           9: 1.128, 10: 1.27, 11: 1.41, 14: 1.693, 18: 2.257}

bar_area_list = {3: 0.11, 4: 0.2, 5: 0.31, 6: 0.44, 7: 0.6, 8: 0.79,
           9: 1, 10: 1.27, 11: 1.56, 14: 2.25, 18: 4}

bar_weight_list = {3: 0.376, 4: 0.668, 5: 1.044, 6: 1.503, 7: 2.046, 8: 2.673,
           9: 3.4, 10: 4.311, 11: 5.313, 14: 7.660, 18: 13.614}


# --------------- Class Definitions ------------------

class Bar:
    """
    Bar(size)

    Class representing a piece of rebar.

    ---Attributes---
    Bar.size - the rebar size
    Bar.dia - the rebar diameter in inches
    Bar.area - the rebar area in square inches
    Bar.weight - the weight of rebar in lb per linear foot
    Bar.bend_rad - basic bend radius for a certain bar size
    Bar.stirrup_bend_rad - stirrup bend radius for a certain bar size
    Bar.bend_ext_180 - extension following bend for a 180 deg hook
    Bar.bend_ext_90 - extension following bend for a 90 deg hook
    Bar.stirrup_bend_ext_135 - extension following bend for a 135 deg stirrup hook
    Bar.stirrup_bend_ext_90 - extension following bend for a 90 deg stirrup hook
    Bar.dev(fc=5) - basic tension development length with default fc value = 5
    Bar.splice(class_="b",fc=5) - class A or B splice length with default class = "b" and fc = 5
    Bar.info(fc=5) - returns info for a bar size with default fc = 5
    """
    def __init__(self, size):
        """
        Initializes a piece of rebar with its size and assigns it the
        following attributes:
        self.dia - bar diameter in inches
        self.area - bar area in square inches
        self.weight - bar weight in lb per linear foot
        """
        valid_sizes = [3,4,5,6,7,8,9,10,11,14,18]

        if size not in valid_sizes:
            raise ValueError("Not a valid bar size!")

        self.size = size

        bars = {'dia': bar_dia_list,
                'area': bar_area_list,
                'weight': bar_weight_list}

        self.dia = bars['dia'][self.size]
        self.area = bars['area'][self.size]
        self.weight = bars['weight'][self.size]

        self.dia = bars['dia'][self.size]
        self.area = bars['area'][self.size]
        self.weight = bars['weight'][self.size]

        # bend radii
        if size <= 8:
            self.bend_rad = self.dia * 6 / 2
        elif size == 9 or size == 10 or size == 11:
            self.bend_rad = self.dia * 8 / 2
        else:
            self.bend_rad = self.dia * 10 / 2

        # stirrup bend radii
        self.stirrup_bend_rad = self.dia * 4 / 2

        # bend extension
        self.bend_ext_180 = max(2.5/12, self.dia * 4)
        self.bend_ext_90 = self.dia * 12

        # stirrup bend extension
        if size < 9:
            self.stirrup_bend_ext_135 = self.dia * 4
        else:
            self.stirrup_bend_ext_135 = "Invalid"
        if size <= 5:
            self.stirrup_bend_ext_90 = self.dia * 6
        elif size == 6 or size == 7 or size == 8:
            self.stirrup_bend_ext_90 = self.dia * 12
        else:
            self.stirrup_bend_ext_90 = "Invalid"
    
    # implement string representation
    def __str__(self):
        return "#{0!s} bar".format(self.size)


    # implement equality behaviors
    def __eq__(self, other):
        return self.size == other.size
    

    def __ne__(self, other):
        return self.size != other.size
    

    def __gt__(self, other):
        return self.size > other.size
    

    def __lt__(self, other):
        return self.size < other.size
    

    def __ge__(self, other):
        return self.size >= other.size
    

    def __le__(self, other):
        return self.size <= other.size


    #development & splices
    def dev_old(self, fc=5):
        """
        dev_old(self,fc=5)

        Returns the basic tension development length for a bar size, with the default fc = 5

        AASHTO 5.11.2.1.1:

        For bars #11 and under, the development length is 1.25 * A_bar * fy / sqrt(f`c) but not less than 0.4 * Dia_bar * fy

        For #14 bars, the development length is 2.7 * fy / sqrt(f`c)

        For #18 bars, the development length is 3.5 * fy / sqrt(f`c)
        """
        if self.size <= 11:
            return max(1.25*self.area*60/sqrt(fc), 0.4*self.dia*60)
        elif self.size == 14:
            return 2.7*60/sqrt(fc)
        elif self.size == 18:
            return 3.5*60/sqrt(fc)
        else:
            return "Invalid"
    

    def dev(self, fc=5):
        """
        dev(self,fc)

        Returns the basic development length for a piece of rebar.

        AASHTO Equation 5.11.2.1.1-2
        ldb = 2.4*db*fy/sqrt(f'c)
        """
        return 2.4*self.dia*60/sqrt(fc)
    

    def splice_old(self, class_="c", fc=5):
        """
        splice_old(self, class_="c", fc=5)

        Returns the required splice length per AASHTO given the splice class and the concrete strength f`c.

        AASHTO 5.11.5.3.1:

        Class A splice - 1.0 * ld
        Class B splice - 1.3 * ld
        Class C splice - 1.7 * ld
        """
        if class_ == "a":
            return self.dev(fc=fc)
        elif class_ == "b":
            return self.dev(fc=fc) * 1.3
        elif class_ == "c":
            return self.dev(fc=fc) * 1.7
        else:
            return "Invalid class."


    def splice(self, class_="b", fc=5):
        """
        splice(self, class_="b", fc=5)

        Returns the required splice length per AASHTO given the splice class and the concrete strength f`c.

        AASHTO 5.11.5.3.1:

        Class A splice - 1.0 * ld
        Class B splice - 1.3 * ld
        """
        if class_ == "a":
            return self.dev(fc=fc)
        elif class_ == "b":
            return self.dev(fc=fc) * 1.3
        else:
            return "Invalid class."


    def info(self, fc=5):
        print("#{0!s} bar data:".format(self.size))
        print("Diameter:          {0} in".format(self.dia))
        print("                   {0} ft".format(round(self.dia/12,4)))
        print("Radius:            {0} in".format(self.dia/2))
        print("                   {0} ft".format(round(self.dia/24,4)))
        print("Area:              {0} sq in".format(self.area))
        print("Weight:            {0} lb/ft".format(self.weight))
        print("Bend radius:       {0} in".format(self.bend_rad))
        print("                   {0} ft".format(round(self.bend_rad/12,4)))
        if type(self.stirrup_bend_rad) is not str:
            print("(stirrup):         {0} in".format(self.stirrup_bend_rad))
            print("                   {0} ft".format(round(self.stirrup_bend_rad/12,4)))
        print("Extension (180):   {0} in".format(self.bend_ext_180))
        print("                   {0} ft".format(round(self.bend_ext_180/12,4)))
        print("Extension (90):    {0} in".format(self.bend_ext_90))
        print("                   {0} ft".format(round(self.bend_ext_90/12,4)))
        if type(self.stirrup_bend_ext_135) is not str:
            print("Stirrup ext (135): {0} in".format(self.stirrup_bend_ext_135))
            print("                   {0} ft".format(round(self.stirrup_bend_ext_135/12,4)))
        if type(self.stirrup_bend_ext_90) is not str:
            print("Stirrup ext (90):  {0} in".format(self.stirrup_bend_ext_90))
            print("                   {0} ft".format(round(self.stirrup_bend_ext_90/12,4)))
        if type(self.dev(fc=fc)) is not str:
            print("Basic dev. length: {0} in".format(round(self.dev(fc=fc),3)))
            print("                   {0} ft".format(round(self.dev(fc=fc)/12,3)))
        print("Class A splice:    {0} in".format(round(self.splice(class_="a",fc=fc),3)))
        print("                   {0} ft".format(round(self.splice(class_="a",fc=fc)/12,3)))
        print("Class B splice:    {0} in".format(round(self.splice(class_="b",fc=fc),3)))
        print("                   {0} ft".format(round(self.splice(class_="b",fc=fc)/12,3)))

           
class Prestress:
    def __init__(self,strand=0.6, fpu=270):
        if strand == 0.6:
            self.area = 0.217
        elif strand == 0.5:
            self.area = 0.153
        else:
            raise ValueError("Strand diameter must be 0.5 or 0.6!")
        self.strand = strand
        self.fpu = fpu

           
class Concrete:
    """
    Concrete(fc=4000)

    Class representing a concrete material. Default f'c = 4000 psi,
    can be specified.

    ---Attributes---
    Concrete.fc - f'c (compressive strength of the concrete) default 
        in psi
    Concrete.fc_units - f'c units (default psi)
    Concrete.unit_weight - 150 pcf
    Concrete.e - modulus of elasticity - 33*unit_wt^1.5*sqrt(fc)/1000

    ---Functions---
    Concrete.set_fc(new_fc) - changes the f'c of the concrete (in psi)
    """
    def __init__(self, fc=4000):
        self.fc = fc
        self.fc_units = 'psi'
        self.unit_weight = 150
        self.e = 33*self.unit_weight**1.5*math.sqrt(self.fc)/1000

    def set_fc(self, new_fc):
        """
        set_fy(new_fc)

        Changes the f'c of the concrete and recalculates the modulus
        of elasticity.
        """
        if type(new_fc) is int or type(new_fc) is float:
            self.fc = new_fc
        else:
            raise ValueError("The argument of the 'set_fc' function must be a valid integer or float.")
        self.e = 33*self.unit_weight**1.5*math.sqrt(new_fc)/1000


class Steel:
    """
    Steel(fy=60)

    Class representing a steel material. Default fy = 60 ksi,
    can be specified.

    ---Attributes---
    Steel.fy - fy (yield strength) of the steel - default in ksi
    Steel.fy_units - default is ksi
    Steel.e - modulus of elasticity - 29000 ksi
    Steel.unit_weight - 490 pcf

    ---Functions---
    Steel.set_fy(new_fy) - changes the fy of the steel (in ksi)
    """
    def __init__(self, fy=60):
        self.fy = fy
        self.fy_units = 'ksi'
        self.e = 29000
        self.unit_weight = 490

    def set_fy(self, new_fy):
        """
        set_fy(new_fy)

        Changes the fy of the steel.
        """
        if type(new_fy) is int or type(new_fy) is float:
            self.fy = new_fy
        else:
            raise ValueError("The argument of the 'set_fy' function must be a valid integer or float.")


class Circle:
    """
    Circle(r)

    Circle object. Initialized with r (radius).

    ---Attributes---
    Circle.r - radius
    Circle.area - area of the circle - pi * r^2
    Circle.i - Moment of inertia of the circle - pi/4 * r^4
    Circle.s - section modulus of the circle - I/r
    Circle.circumf - circumference of the circle - 2*pi*r
    Circle.dia - diameter of the circle - 2 * r
    Circle.rg - radius of gyration - sqrt(I/A)
    """
    def __init__(self, r):
        self.r = r
        self.area = pi * r**2
        self.i = pi/4 * r**4
        self.s = self.i / self.r
        self.circumf = 2*pi*self.r
        self.dia = self.r * 2
        self.rg = sqrt(self.i/self.area)


class Tube:
    """
    Tube(r1, r2)

    Tube object. Initialized with radius 1 and radius 2 (doesn't matter
    which one is bigger).

    ---Attributes---
    Tube.r1 - max of r1, r2 arguments
    Tube.r2 - min of r1, r2 arguments
    Tube.area - area - pi * (r1^2 - r2^2)
    Tube.i - moment of inertia - pi/4 * (r1^4 - r2^4)
    Tube.s - section modulus - I/r1
    Tube.circumf - circumference - 2*pi*r1
    Tube.dia - diameter - 2*r1
    Tube.rg - radius of gyration - sqrt(I/A)
    """
    def __init__(self, r1, r2):
        self.r1 = max(r1,r2)
        self.r2 = min(r1,r2)
        self.area = pi * (r1**2 - r2**2)
        self.i = pi/4 * (r1**4 - r2**4)
        self.s = self.i / self.r1
        self.circumf = 2*pi*self.r1
        self.dia = self.r1 * 2
        self.rg = sqrt(self.i/self.area)


class Rectangle:
    """
    Rectangle(b, h)

    Rectangle object. Initialized with base and height.

    ---Attributes---
    Rectangle.b - base
    Rectangle.h - height
    Rectangle.area - area - b * h
    Rectangle.perim - perimeter - 2*b + 2*h
    Rectangle.i - moment of inertia - b*h^3/12
    Rectangle.s - section modulus - I/(h/2)
    Rectangle.rg - radius of gyration - sqrt(I/A)
    """
    def __init__(self, b, h):
        self.b = b
        self.h = h
        self.area = self.b * self.h
        self.perim = 2 * self.b + 2 * self.h
        self.i = self.b * self.h ** 3 / 12
        self.s = self.i / (self.h/2)
        self.rg = sqrt(self.i/self.area)

# ------------------- Rebar Functions --------------------

def rebar():
    """
    rebar()

    Print a rebar table.
    """
    x = prettytable.PrettyTable()

    x.field_names = ["Bar Size", "Diameter, in", "Area, sq in", "Weight, plf"]

    for bar in [3,4,5,6,7,8,9,10,11,14,18]:
        b = Bar(bar)
        x.add_row([b.size, b.dia, b.area, b.weight])

    print(x)


def rebarxl():
    """
    rebarxl()

    Outputs a simple rebar table delimited with spaces. Used for
    pasting a rebar table into excel for use with VLOOKUP etc.
    """
    print("size dia area wt")
    for bar in [3,4,5,6,7,8,9,10,11,14,18]:
        b = Bar(bar)
        print("{0} {1} {2} {3}".format(b.size, b.dia, b.area, b.weight))


def eff_d(h, cover, bars):
    """
    eff_d(h, cover, bars)

    Return the effective depth of a concrete beam given h, the height of the
    section, cover, the concrete cover to the edge of the outermost rebar, and
    bars, a list of bar sizes from the outside in, with the flexural bars, the
    ones you are calculating d for, are the last ones in the list. Most lists
    will be one or two elements.
    """
    bar_dia = bar_dia_list
    bar_dist = sum([bar_dia[bar] for bar in bars]) - 0.5 * bar_dia[bars[-1]]
    return h - cover - bar_dist


def rebar_cent(cover, bars):
    """
    rebar_cent(cover, bars)

    Return the distance from the outside of a concrete beam to the centroid of a
    piece of rebar given cover, the concrete cover to the edge of the outermost
    rebar, and bars, a list of bar sizes from the outside in, with the flexural
    bars, the ones you are calculating d for, are the last ones in the list.
    Most lists will be one or two elements.
    """
    bar_dia = bar_dia_list
    bar_dist = sum([bar_dia[bar] for bar in bars]) - 0.5 * bar_dia[bars[-1]]
    return cover + bar_dist


# ------------------- Unit Conversions ----------------------

def feet(i):
    """
    feet(i)

    Return i (in inches) converted to feet.
    """
    return i / 12


def inches(f):
    """
    inches(f)

    Return f (in feet) converted to inches.
    """
    return f * 12


def dms_to_dec(deg, min_, sec):
    """
    dms_to_dec(deg, min_, sec)

    Return an angle in degrees, decimal form, given degrees, minutes, and seconds.
    """
    min_ = min_ + sec/60
    return deg + min_/60


def rad_deg(rad):
    """
    rad_deg(rad)

    Return an angle in degrees given the angle in radians.
    """
    return rad*180/pi


def deg_rad(deg):
    """
    deg_rad(deg)

    Return an angle in radians given the angle in degrees.
    """
    return deg*pi/180


def feetdisp(num, precision=8, roundup=False):
    """
    feetdisp(num, precision)

    Print a value in "feet - inches" with fractions of an inch.
    """
    feet = math.floor(num)
    inches_tot = (num - feet)*12
    whole_inches = math.floor(inches_tot)
    if roundup:
        inches_frac = int(math.ceil((inches_tot - whole_inches)*precision))
    else:
        inches_frac = int(round((inches_tot - whole_inches)*precision, 0))
    gcd_ = gcd(precision, inches_frac)
    if gcd_ > 1:
        inches_frac = int(inches_frac / gcd_)
        precision = int(precision / gcd_)
    if inches_frac == precision:
        whole_inches += 1
        inches_frac = 0
    if whole_inches == 12:
        feet += 1
        whole_inches = 0
    if inches_frac != 0:
        print("{0}'-{1} {2}/{3}\"".format(feet, whole_inches, inches_frac, precision))
    else:
        print("{0}'-{1}\"".format(feet, whole_inches))


def sigfigs(num, n):
    """
    sig_figs(num, n)
    
    Converts num to a float with n significant figures.
    """
    return round(num,n-len(str(int(num))))

# --------------------- Section Properties ----------------------

def cent_grav(areas, ybars):
    """
    cent_grav(areas, ybars)

    Return the center of gravity of an object given the areas and centers
    of gravity of its individual parts.
    """
    if len(areas) != len(ybars):
        raise ValueError("The two input lists must be the same length!")

    return sum(map(mul, areas, ybars)) / sum(areas)


def area(shape):
    """
    area(shape)

    Calculates the area of a shape. shape can be 'circle', 'tube', or
    'rect'. This function is deprecated in favor of the shape classes.
    """
    shape_table = {'circle': {'dims': ['r: '],
                              'eq': lambda d: pi * d[0]**2},
                   'tube': {'dims': ['r1: ', 'r2: '],
                            'eq': lambda d: pi * abs(d[0]**2 - d[1]**2)},
                   'rect': {'dims': ['b: ', 'h: '],
                            'eq': lambda d: d[0] * d[1]}
    }

    dims = list(map(lambda dim: eval(input(dim)), shape_table[shape]['dims']))
    return shape_table[shape]['eq'](dims)


def i(shape):
    """
    i(shape)

    Calculates the moment of inertia of a shape. shape can be 'circle',
    'tube', or 'rect'. This function is deprecated in favor of the
    shape classes.
    """
    shape_table = {'circle': {'dims': ['r: '],
                              'eq': lambda d: pi * d[0]**4 / 4},
                   'tube': {'dims': ['r1: ', 'r2: '],
                            'eq': lambda d: (pi/4) * abs(d[0]**4 - d[1]**4)},
                   'rect': {'dims': ['b: ', 'h: '],
                            'eq': lambda d: d[0] * d[1]**3 / 12}
    }

    dims = list(map(lambda dim: eval(input(dim)), shape_table[shape]['dims']))
    return shape_table[shape]['eq'](dims)


def rg(shape):
    """
    rg(shape)

    Caclulates the radius of gyration of a shape. shape can be 'circle',
    'tube', or 'rect'. This function is deprecated in favor of the
    shape classes.
    """
    shape_table = {'circle': {'dims': ['r: '],
                              'i': lambda d: pi * d[0]**4 / 4,
                              'area': lambda d: pi * d[0]**2},
                   'tube': {'dims': ['r1: ', 'r2: '],
                            'i': lambda d: (pi/4) * abs(d[0]**4 - d[1]**4),
                            'area': lambda d: pi * abs(d[0]**2 - d[1]**2)},
                   'rect': {'dims': ['b: ', 'h: '],
                            'i': lambda d: d[0] * d[1]**3 / 12,
                            'area': lambda d: d[0] * d[1]}
    }
    dims = list(map(lambda dim: eval(input(dim)), shape_table[shape]['dims']))
    return sqrt(shape_table[shape]['i'](dims) / shape_table[shape]['area'](dims))

# ------------------ Trig (Custom sin, cos, tan with degrees)-----------------

def cosd(x):
    """
    cosd(x)

    Return the cosine of x (measured in degrees).
    """
    return cos(x * pi / 180)


def sind(x):
    """
    sind(x)

    Return the sine of x (measured in degrees).
    """
    return sin(x * pi / 180)


def tand(x):
    """
    tand(x)

    Return the tangent of x (measured in degrees).
    """
    return tan(x * pi / 180)


def acosd(x):
    """
    acosd(x)

    Return the arc cosine in degrees.
    """
    return acos(x) * 180 / pi


def asind(x):
    """
    asind(x)

    Return the arc sine in degrees.
    """
    return asin(x) * 180 / pi


def atand(x):
    """
    atand(x)

    Return the arc tangent in degrees.
    """
    return atan(x) * 180 / pi


# -------------------- Engineering Calcs ---------------------

def norm_stress(p, a):
    """
    norm_stress(p, a)

    Returns the normal stress given the force, p, and the area, a.
    """
    return p / a


def bend_stress(mom, y, i):
    """
    bend_stress(mom, y, i)

    Returns the bending stress given the applied moment, distance to the
    centroid, and moment of inertia.
    """
    return mom * y / i


def e_conc(fc, wc):
    """
    e_conc(fc, wc)

    Returns the modulus of elasticity of concrete in psi, given the strength, fc
    in psi, and the unit weight, wc in pcf.
    """
    return 33 * wc**1.5 * (sqrt(fc))


def simple_moment(w, l):
    """
    simple_moment(w, l)

    Returns the simple span moment given w and l.
    """
    return w * l**2 / 8

# ---------------------- Statistics -----------------------

def pct_diff(x, y):
    """
    pct_diff(x, y)

    Returns the percent difference between two numbers, base number x, and
    "new" number y.
    """
    pct = round((abs(y - x) / x) * 100, 2)
    print(str(pct) + '%')
    return pct/100


def spaces(x, spa):
    """
    spaces(x, spa)

    Returns an array of the spacings from the original point given the distance
    and the number of spaces.
    """
    return list(map(lambda y: y * (x/spa), range(1, spa + 1)))


def sumlist(x):
    """
    sumlist(x)

    Returns the running sum of a list of numbers, x. For example:
    sumlist([1,2,3]) would return [1,3,6].
    """
    try:
        return list(itertools.accumulate(x))
    except TypeError:
        raise TypeError("unsupported argument type for 'sumlist'")


def shift_list(l, x):
    """
    shift_list(l, x)

    Shifts all values in a list l, by x.
    """
    return [i+x for i in l]


def scale_list(l, x):
    """
    scale_list(l, x)

    Scales all values in a list l, by x.
    """
    return [i*x for i in l]


def lin_interp(x,x0,x1,y0,y1):
    """
    lin_interp(x,x0,x1,y0,y1)
    
    Linearly interpolates for y between y0 and y1, given x0, x1, and x.
    """
    
    return y0+(y1-y0)*(x-x0)/(x1-x0)

def sum_range(x):
    """
    sum_range(x)
    
    Sums a range from 1 to x.
    """
    return sum(range(1, x+1))

def avg(*args):
    """
    avg(*args)

    Basically a wrapper around statistics.mean that allows the user to not put the
    values in a list.
    """
    return mean(args)

# -------------------- Geometry ----------------------------


def spiral(h, pitch, dia):
    """
    sprial(h, pitch, dia)

    Returns the length of a spiral (helix), given the total height, pitch, and
    diameter of the spiral (does not include full turns at bottom and top, if
    present).
    """
    circumf = sqrt((pi*dia)**2 + pitch**2)
    num_spirals = h/pitch
    tot_length = circumf * num_spirals
    return tot_length
