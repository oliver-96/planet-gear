import matplotlib
matplotlib.use("Qt5Agg")  # works for vscode

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.animation import FuncAnimation


# --- Initial values ---
r = 1
a_initial = 1
k_initial = 5
R_initial = k_initial * r
mode_initial = "epicycloid"

# --- Figure setup ---
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25) 
plt.subplots_adjust(left=0.2)
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
axis_factor = 1.2

# --- Cycloid Equations ---
def cycloid_curve(t, R, r, a, mode):
    """Return x, y for epicycloid or hypocycloid."""

    if mode == "epicycloid":
        x = (R + r) * np.cos(t) - a * np.cos(((R + r) / r) * t)
        y = (R + r) * np.sin(t) - a * np.sin(((R + r) / r) * t)
    elif mode == "hypocycloid":
        x = (R - r) * np.cos(t) + a * np.cos(((R - r) / r) * t)
        y = (R - r) * np.sin(t) - a * np.sin(((R - r) / r) * t)
    else:
        raise ValueError("Mode must be 'epicycloid' or 'hypocycloid'")
    return x, y

# --- Cycloid Points ---
t = np.linspace(0, 2 * np.pi * 10, 2000)
x_inital, y_initial = cycloid_curve(t, R_initial, r, a_initial, mode_initial)

# --- Gears ---
planet_x = (R_initial - r) if mode_initial == "hypocycloid" else (R_initial + r)

sun_circle = plt.Circle((0, 0), R_initial, color='black', fill=False, linestyle='-')
planet_circle = plt.Circle((planet_x, 0), r, color='black', fill=False, linestyle='-')

centre_point = plt.Circle((planet_x, 0), 0.1, color='black', fill=True)
cir_point = plt.Circle((planet_x - r, 0), 0.1, color='black', fill=True)
cir_point_outer = plt.Circle((planet_x + r, 0), 0.1, color='black', fill=True)
trace_point = plt.Circle((planet_x - r,0), 0.1, color='blue', fill=True)
rotating_line = plt.Line2D([], [], color='black', lw=1)

# --- Slider setup ---
k_ax_slider = plt.axes([0.25, 0.1, 0.5, 0.03])  # [left, bottom, width, height]
a_ax_slider = plt.axes([0.25, 0.05, 0.5, 0.03])  # [left, bottom, width, height]
speed_ax_slider = plt.axes([0.25, 0.01, 0.5, 0.03])  # [left, bottom, width, height]

k_slider = Slider(k_ax_slider, 'k = R/r', 0, 10, valinit=k_initial, valstep=0.1)
a_slider = Slider(a_ax_slider, 'a', 0, 10, valinit=a_initial, valstep=0.1)
speed_slider = Slider(speed_ax_slider, 'Speed', 0, 1, valinit=0.5, valstep=0.05)

# --- Radio button setup ---
mode_ax = plt.axes([0.01, 0.4, 0.2, 0.15], facecolor='lightgray')
mode_radio = RadioButtons(mode_ax, ('epicycloid', 'hypocycloid'), active=0)

# --- Create Plot ---
ax.add_patch(sun_circle)
ax.add_patch(planet_circle)
[line] = ax.plot(x_inital, y_initial, 'r-', lw=1)

ax.add_patch(centre_point)
ax.add_patch(cir_point)
ax.add_patch(cir_point_outer)
ax.add_patch(trace_point)
ax.add_line(rotating_line)

# --- Update for sliders ---
def update(val):
    k = k_slider.val
    R = k * r
    a = a_slider.val
    mode = mode_radio.value_selected

    # Update curve
    x, y = cycloid_curve(t, R, r, a, mode)
    line.set_data(x, y)

    # Update circles
    sun_circle.set_radius(R)

    # Calc axis limits
    if mode == "epicycloid":
        ax.set_xlim( (-R - max(2*r, a + r)) * axis_factor, (R + max(2*r, a + r)) * axis_factor )
        ax.set_ylim( (-R - max(2*r, a + r)) * axis_factor, (R + max(2*r, a + r)) * axis_factor )
    else:
        ax.set_xlim( -max(R, R-r+a) * axis_factor, max(R, R-r+a) * axis_factor )
        ax.set_ylim( -max(R, R-r+a) * axis_factor, max(R, R-r+a) * axis_factor )

    fig.canvas.draw_idle()


# --- Connect widgets ---
k_slider.on_changed(update)
a_slider.on_changed(update)
mode_radio.on_clicked(update)

# --- Animation update ---
def animate_update(frame):
    k = k_slider.val
    R = k * r
    a = a_slider.val
    mode = mode_radio.value_selected
    speed_input = speed_slider.val

    theta_step = (2 * np.pi * speed_input) * (1/80)
    theta = theta_step * frame

    static_centre = (R - r) if mode == "hypocycloid" else (R + r)

    # rotation around sun circle
    cx = static_centre * np.cos(theta)
    cy = static_centre * np.sin(theta)
    planet_circle.center = (cx, cy)

    # cycloid point at given theta
    x_point, y_point = cycloid_curve(theta, R, r, a, mode)

    centre_point.center = (cx, cy)

    # rotation of planet circle
    phi = ((R - r) / r) * theta * -1 if mode == "hypocycloid" else ((R + r) / r) * theta
    
    cir_point.center = (cx - r * np.cos(phi), cy - r * np.sin(phi))
    cir_point_outer.center = (cx + r * np.cos(phi), cy + r * np.sin(phi))

    trace_point.center = (x_point, y_point)
    rotating_line.set_data([cx, x_point], [cy, y_point])

    return planet_circle, centre_point, cir_point, cir_point_outer, trace_point, rotating_line

def frame_gen():
    n = 0
    while True:
        yield n
        n += 1

ani = FuncAnimation(fig, animate_update, frames=frame_gen, interval=30, blit=False)

plt.show()
