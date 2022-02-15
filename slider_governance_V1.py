# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 07:56:03 2018

@author: amador
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.integrate import odeint


# define function
def resilience(X, t, parameters):
    x, y = X
    A,B,C,D = parameters
    derivatives = [A*x*y-B*x, C*x*y-D*y] 
    return derivatives

def PE(parameters):
    A,B,C,D = parameters
    Geq = D/C
    Beq = B/A
    eq = Geq, Beq
    return eq

def UV(Xn, Yn, parameters):
    A,B,C,D = parameters
    u = A*Xn*Yn-B*Xn
    v = C*Xn*Yn-D*Yn
    return u,v

# define initial conditions and parameters
Amin = 0
Amax = 1    
A0 = 1
Bmin = 0
Bmax = 1
B0 = 0.5
Cmin = 0
Cmax = 1    
C0 = 1
Dmin = 0
Dmax = 1
D0 = 0.5
xmin = 0
xmax = 1
ymin = 0
ymax = 1
x0 = 0.5
y0 = 0.5

# ODE optiones
parameters0 = [A0,B0,C0,D0]
X0 =[x0,y0]
tmin = 0
tmax = 8
tstepsize = 0.0001
t = np.arange(tmin, tmax, tstepsize)
xs, ys = odeint(resilience, X0, t, args=(parameters0,)).T   # solve equations

# Create figures to plot in
fig = plt.figure(figsize=(10, 5))   # define the size of plotting area
ax1 = fig.add_axes([0.1, 0.4, 0.45, 0.5]) #  [left, bottom, width, height]
x_plot = ax1.plot(t, xs, label = "Governance")[0]
y_plot = ax1.plot(t, ys, label = "Nature's benefits")[0]
ax1.set_xlabel("Time")
ax1.legend(bbox_to_anchor = (0.1, 1.02, 0.8, .102), loc = 4,ncol = 2, mode = "expand",borderaxespad = 1.)
ax1.grid(color = "b", alpha = 0.5, linestyle = "dashed", linewidth = 0.5)
ax1.set_ylim([np.min(ys)-0.1, np.max(ys)+0.1])
ax1.set_xlim([tmin-1, tmax+1])
ax2 = fig.add_axes([0.65, 0.4, 0.3, 0.5]) #  [left, bottom, width, height]
xy_plot = ax2.plot(xs, ys,'r')[0]
xy_plot0 = ax2.plot(xs[0], ys[0],'bo')[0]
xy_ploteq = ax2.plot(PE(parameters0)[0],PE(parameters0)[1],'ro')[0]
ax2.set_xlabel("Governance")
ax2.set_ylabel("Nature's benefit")
#ax2.grid(color = "b", alpha = 0.5, linestyle = "dashed", linewidth = 0.5)
ax2.set_xlim([0,1])
ax2.set_ylim([0,1])


# generate initial vector field 
X = np.linspace(0, 1, 30)
Y = np.linspace(0, 1, 30)
Xn, Yn = np.meshgrid(X, Y)
u, v = UV(Xn, Yn, parameters0)
vector_field = ax2.quiver(X, Y, u, v, color = '#4d4d4d', scale = 4, headwidth=4,  headlength=5, linewidths=5, scale_units='x')


# create a space in the figure to place the n sliders
axcolor = 'lightgoldenrodyellow'
axis_A = fig.add_axes([0.10, 0.25, 0.3, 0.03], facecolor=axcolor)  #[left, bottom, width, height
axis_B = fig.add_axes([0.10, 0.2, 0.3, 0.03], facecolor=axcolor)
axis_C = fig.add_axes([0.10, 0.15, 0.3, 0.03], facecolor=axcolor)
axis_D = fig.add_axes([0.10, 0.1, 0.3, 0.03], facecolor=axcolor)
axis_x = fig.add_axes([0.50, 0.2, 0.3, 0.03], facecolor=axcolor)
axis_y = fig.add_axes([0.50, 0.15, 0.3, 0.03], facecolor=axcolor)

# create each slider 
slider_A = Slider(axis_A, 'A', Amin, Amax, valinit=A0)
slider_B = Slider(axis_B, 'B', Bmin, Bmax, valinit=B0)
slider_C = Slider(axis_C, 'C', Cmin, Cmax, valinit=C0)
slider_D = Slider(axis_D, 'D', Dmin, Dmax, valinit=D0)
slider_x = Slider(axis_x, 'x', xmin, xmax, valinit=x0)
slider_y = Slider(axis_y, 'y', ymin, ymax, valinit=y0)

def update(val):
    parameters = [slider_A.val, slider_B.val, slider_C.val, slider_D.val]
    x = [slider_x.val, slider_y.val]
    # recalculate the function values
    xs, ys = odeint(resilience, x, t, args=(parameters,)).T
    u, v = UV(Xn, Yn, parameters)
    # update the value on the graph
    x_plot.set_ydata(xs)
    y_plot.set_ydata(ys)
    xy_plot.set_data(xs,ys)
    xy_plot0.set_data(xs[0],ys[0])
    xy_ploteq.set_data(PE(parameters)[0],PE(parameters)[1])
    vector_field.set_UVC( u, v)
    # redraw the graph
    fig.canvas.draw_idle()
    ax1.set_ylim([-0.1, 1.1])
    ax1.set_xlim([tmin-1, tmax+1])
    ax2.set_xlim([np.min(xs)-0.5, np.max(xs)+0.5])
    ax2.set_ylim([np.min(ys)-0.5, np.max(ys)+0.5])
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    
# set all sliders to call update when their value is changed
slider_A.on_changed(update)
slider_B.on_changed(update)
slider_C.on_changed(update)
slider_D.on_changed(update)
slider_x.on_changed(update)
slider_y.on_changed(update)

# create the reset button axis (where its drawn)
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
# and the button itself
button = Button(resetax, 'Reset', color='yellow', hovercolor='0.575')

def reset(event):
    slider_A.reset()
    slider_B.reset()
    slider_C.reset()
    slider_D.reset()
    slider_x.reset()
    slider_y.reset()

button.on_clicked(reset)