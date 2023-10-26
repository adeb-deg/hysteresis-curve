# Hysteresis Curve Analysis Tool

This Python module provides tools for analyzing hysteresis curves using various mathematical and visualization techniques.

## Features
- Generate and analyze hysteresis curves with functionalities to:
    - Resample data points.
    - Estimate initial slopes.
    - Identify unloading and reloading branches.
    - Calculate Effective Dissipation Ratio (EDR).
    - Visualize hysteresis polygons and data plots.
    - Find intersection points and extend LineString objects.
    - Generate hysteresis polygons.

## Dependencies
- `numpy`
- `matplotlib`
- `shapely`
- `scipy`

# HysteresisCurve Class
Provides methods to handle and analyze hysteresis curves.

## Methods
- `__init__(self, d, f, labels=('D', 'F'))`: Constructor method to initialize the curve with d and f arrays.
- `resample(self, n)`: Resamples the data points of the curve.
- `get_init_slope_estimate(self)`: Estimates the initial slope.
- `unloading_branches(self)`: Identifies the unloading branches.
- `reloading_branches(self)`: Identifies the reloading branches.
- `hysteresis_polygons(self)`: Generates polygons representing the hysteresis.
- `edr(self)`: Calculates the Effective Dissipation Ratio (EDR).
- `plot(self, plot_cycle_num=None, **kwargs)`: Visualizes the hysteresis curve, supporting optional arguments like color and custom figure axes.

## Quick Start
To use the functionalities provided in the module, you can start by creating an instance of the `HysteresisCurve` class with data arrays for D and F values:

```python
from hysteresis_module import HysteresisCurve  # Assuming this module is named 'hysteresis_module.py'

# Example data
d_values = [ ... ]  # Array of D values
f_values = [ ... ]  # Array of F values

curve = HysteresisCurve(d_values, f_values)
curve.plot()

![Hysteresis Curve Analysis](https://github.com/adeb-deg/hysteresis-curve/raw/main/plot_example.png)