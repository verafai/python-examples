"""For the flame propagation model we use the implicit Euler scheme.
In each step we will need to solve the non linear equation x_n+1=x_n+h*(x_n+1**2-x_n+1**3)
using the Newton-Raphson method to get an aproximation for x_n+1.
For example: in the first step n=0 => x_1=x_0+h*(x_1**2-x_1**3)
we are given x_0 and we will find an aproximation for x_1 using NR
and so on."""

import numpy as np
import matplotlib.pyplot as plt

def NewtonRaphson(f,f_der,x0):
    x_new = x0  
    for i in range (200):
        x1 = x_new-(f(x_new)/f_der(x_new)) 
        if abs(x1-x_new)< 0.0001: break
        x_new = x1
    return([x_new,i])

#test function
f = lambda x : 2*x**3-9.5*x+7.5
f_der = lambda x : 6*x**2-9.5
x0 = 5
NR = NewtonRaphson(f,f_der,x0)
print ('with initial value',x0,',in',NR[1],'iterations, we get the root r= {:.8f}'.format(NR[0]))


def ImplicitEuler(h,x0,t0,T):
    xn = [x0]
    t = [t0]
    while t[-1]<T:
          def g(x):
              return (x - h*(x**2 - x**3) - xn[-1])
          def g_der(x):
              return (1 - h*(2*x - 3*x**2))
          x_NR = NewtonRaphson(g,g_der,xn[-1])
          tN = t[-1] + h
          xn.append(x_NR[0])
          t.append(tN)
    return([t,xn])
    #print (xn)
    #print (t)

h = 2
x0 = 0.01
t0 = 0
T = 200
E_NR = ImplicitEuler(h,x0,t0,T)
#print (E_NR[0])
#print (E_NR[1])


plt.plot(E_NR[0],E_NR[1], c='c')
plt.title('Flame propagation')
plt.savefig('ImplicitEuler.png')
plt.show()

