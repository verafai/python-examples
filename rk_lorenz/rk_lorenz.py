import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

font_family = 'Myriad Pro'
title_font = fm.FontProperties(family=font_family, style='normal', size=20, weight='normal', stretch='normal')
# system parameters
rho = 28.0
sigma = 10.0
beta = 8.0 / 3.0
# initial system state ( x, y, z positions in space)
initial_state = [1.0, 1.0, 1.0] 

# Lorenz system 
def F(current_state, t):
  #positions of x, y, z in space at the current time point
  x, y, z = current_state 
  #3 ordinary differential equations that describe the system
  rhs = sigma * (y - x), x * (rho - z) - y, x * y - beta * z 
  return np.array(rhs)

#function F computes the derivatives:
# dx/dt = sigma * (y - x)
# dy/dt = x * (rho - z) - y
# dz/dt = x * y - beta * z 

#solve using 4th order Runge-Kutta scheme 

# Runge-Kutta step
def step_runge_kutta(F, vn, tn, dt):
    k1 = dt * F(vn, tn)
    k2 = dt * F(vn + dt/2.0, tn + k1/2.0)
    k3 = dt * F(vn + dt/2.0, tn + k2/2.0)
    k4 = dt * F(vn + dt, tn + k3)
    return vn + ( k1 + 2.0*k2 + 2.0*k3 + k4)/6.0

def integrate_runge_kutta(f, initial_state, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    states = [initial_state]
    for tn in t:
        vn = np.array(states[-1])
        step = step_runge_kutta(F, vn, tn, dt)
        states.append(step)
    return np.array(states)


t0 = 0.0
tf = 40.0
dt1 = 0.01 #step size


states_dt1 = integrate_runge_kutta(F, initial_state, t0, tf, dt1)

fig = plt.figure(figsize=(10, 7))
ax = fig.gca(projection="3d")
ax.xaxis.set_pane_color((1,1,1,1))
ax.yaxis.set_pane_color((1,1,1,1))
ax.zaxis.set_pane_color((1,1,1,1))
ax.plot(states_dt1[:, 0], states_dt1[:, 1], states_dt1[:, 2],color='c', alpha=0.7, linewidth=0.7, label='dt=0.01')
ax.legend()
ax.set_xlabel("X axis") 
ax.set_ylabel("Y axis")
ax.set_zlabel("Z axis")
ax.set_title("Lorenz system attractor",fontproperties=title_font)
plt.savefig('rk_lorenz.png',dpi=180, bbox_inches='tight')


# plot two-dimensional cuts of the three-dimensional phase space
fig, ax = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(17, 6))

# plot x values vs y values
ax[0].plot(states_dt1[:, 0], states_dt1[:, 1], color='y', alpha=0.7, linewidth=0.3)
ax[0].set_title('x-y phase plane', fontproperties=title_font)

# plot x values vs z values
ax[1].plot(states_dt1[:, 0], states_dt1[:, 2], color='m', alpha=0.7, linewidth=0.3)
ax[1].set_title('x-z phase plane', fontproperties=title_font)

# plot y values vs  z values
ax[2].plot(states_dt1[:, 1], states_dt1[:, 2], color='b', alpha=0.7, linewidth=0.3)
ax[2].set_title('y-z phase plane', fontproperties=title_font)

fig.savefig('lorenz-attractor-phase-plane.png', dpi=180, bbox_inches='tight')

plt.show()

