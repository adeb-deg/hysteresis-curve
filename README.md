# Hysteresis Curve Analysis Tool

This Python module provides tools for analyzing hysteresis curves using various mathematical and visualization techniques.

## Features
- Generate and analyze hysteresis curves with functionalities to:
    - Resample data points.
    - Estimate initial slopes.
    - Identify unloading and reloading branches.
    - Calculate Energy Dissipation Ratio (EDR), κ<sub>e</sub>
    - Visualize hysteresis polygons and data plots.
    - Find intersection points and extend LineString objects.
    - Generate hysteresis polygons.

## Dependencies
- `numpy`
- `matplotlib`
- `shapely`
- `scipy`
- `piecewise-regression` https://piecewise-regression.readthedocs.io/en/latest/

# HysteresisCurve Class
Provides methods to handle and analyze hysteresis curves.

## Methods
- `__init__(self, d, f, n_seg_env, labels=('D', 'F'))`: Constructor method to initialize the curve with d and f arrays.
- `resample(self, n)`: Resamples the data points of the curve to have specified number of points in each unloading/reloading branch (Total number of points in 1 cycle = 2 * n).
- `get_init_slope_estimate(self)`: Estimates the initial slope.
- `unloading_branches(self)`: Identifies the unloading branches.
- `reloading_branches(self)`: Identifies the reloading branches.
- `hysteresis_polygons(self)`: Generates polygons representing the hysteresis.
- `edr(self)`: Calculates the Energy Dissipation Ratio (EDR).
- `plot(self, plot_cycle_num=None, **kwargs)`: Visualizes the hysteresis curve, supporting optional arguments like color and whether or not to plot envelope.

## Quick Start
To use the functionalities provided in the module, you can start by creating an instance of the `HysteresisCurve` class with data arrays for D and F values.
You will also need to provide `n_seg_env` which is the number of piecewise segments in the piecewise linear envelope.

```python
from hysteresis import HysteresisCurve  # Assuming this module is named 'hysteresis.py'

# Example data
d_values = [ ... ]  # Array of D values
f_values = [ ... ]  # Array of F values

curve = HysteresisCurve(d_values, f_values, n_seg_env=4)
```

You will need to provide a `prominence` value for the `HysteresisCurve` to avoid noise or small data fluctuations from being identified as peaks. For this you can use `curve.see_peaks_and_valleys(prominence=p)` for different values of p to determine prominence.

```python
curve.see_peaks_and_valleys(prominence=0.55)
```
![Hysteresis Curve Analysis](https://github.com/adeb-deg/hysteresis-curve/blob/main/README_figures/PeaksAndValleys.PNG)

Then set the prominence attribute of the `HysteresisCurve` object to the desired value

```python
curve.prominence = 0.55
```

Finally, you can plot the curve with the following arguments:

- `plot_cycle_num`: `None, int, list of int, default=None`  
  Specifies which cycle numbers (0-indexed) to plot. If `None` only the hysteresis curve is plotted. If an integer, plots the specified cycle. If a list of integers, plots the specified cycles.

- `axs`: `None, list of matplotlib.axes, default=None`  
  Allows you to pass a list of matplotlib axes to plot on. If `None`, a new set of (2 or 3) axes will be created and returned based on whether `plot_cycle_num` was `None` or not.  

- `color`: `str, default='b'`  
  Defines the color of the plot. The default color is blue (`'b'`). You can specify any matplotlib color code or name.

- `plot_env`: `bool, default=False`  
  If `True`, includes the plot of the envelope along with the cycle data. The default is `False`, meaning only the cycle data is plotted.


```python
curve.plot()
```
![Hysteresis Curve Analysis](https://github.com/adeb-deg/hysteresis-curve/blob/main/README_figures/fig1.png)
![Hysteresis Curve Analysis](https://github.com/adeb-deg/hysteresis-curve/blob/main/README_figures/fig2.png)

```python
curve.plot([2, 3], plot_env=True)
```
![Hysteresis Curve Analysis](https://github.com/adeb-deg/hysteresis-curve/blob/main/README_figures/fig3.png)
![Hysteresis Curve Analysis](https://github.com/adeb-deg/hysteresis-curve/blob/main/README_figures/fig4.png)
![Hysteresis Curve Analysis](https://github.com/adeb-deg/hysteresis-curve/blob/main/README_figures/fig5.png)
