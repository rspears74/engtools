# Structural Engineering Tools
## Description
A structural/civil engineering tools python module. It contains structural/civil engineering-specific tools along with miscellaneous utilities for working with engineering data sets.

## Installation
```
pip install engtools
```

## Feature Highlights
### Bar class
This class represents a piece of rebar, and has attributes that store its properties like diameter, area, weight, etc. It also has attributes and modules that contain bend radii, bend extensions, and development and splice lengths from AASHTO code equations. There is also a `rebar` function that outputs a formatted rebar table.
#### Example
```python
>>> b5 = Bar(5)
>>> b5.dia
0.625
>>> b5.area
0.31
>>> b5.dev(fc=4)
45.0
>>> rebar()
+----------+--------------+-------------+-------------+
| Bar Size | Diameter, in | Area, sq in | Weight, plf |
+----------+--------------+-------------+-------------+
|    3     |    0.375     |     0.11    |    0.376    |
|    4     |     0.5      |     0.2     |    0.668    |
|    5     |    0.625     |     0.31    |    1.044    |
|    6     |     0.75     |     0.44    |    1.503    |
|    7     |    0.875     |     0.6     |    2.046    |
|    8     |      1       |     0.79    |    2.673    |
|    9     |    1.128     |      1      |     3.4     |
|    10    |     1.27     |     1.27    |    4.311    |
|    11    |     1.41     |     1.56    |    5.313    |
|    14    |    1.693     |     2.25    |     7.66    |
|    18    |    2.257     |      4      |    13.614   |
+----------+--------------+-------------+-------------+
```
### Shape classes
Rectangle, Cirle, and Tube shapes that have attributes and methods defining their properties like area, moment of inertia (i), section modulus (s), and radius of gyration (rg).
#### Example
```python
>>> r = Rectangle(10,12)
>>> r.area
120
>>> r.i
1440.0
>>> r.s
240.0
>>> r.rg
3.4641016151377544
```
### feetdisp function
Takes a decimal value in feet and prints a feet-inches representation with user specified precision (default 1/8 inch).
#### Example
```python
>>> feetdisp(num)
4'-4 5/8"
>>> feetdisp(num, precision=2)
4'-4 1/2"
```
### Miscellaneous utilities highlights
- `sigfigs(num, n)` Converts num to a float with n significant figures
- `cent_grav(areas, ybars)` Calculates the center of gravity of something given areas (or weights) and centers of individual pieces
- `sind`, `cosd`, `tand`, `asind`, etc. - Trig functions for degrees
- `spaces(x, spa)` Returns an array of x divided into spa spaces
- `sumlist(x)` Returns running sum of a list. Basically an alias for `list(itertools.accumulate(x))`
- `shift_list(l, x)` and `scale_list(l, x)` Shift/scale list `l` by `x`
