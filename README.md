# Python examples
Here are some computational exercises I solved using Python for a intro course in  numerical analysis.
<br/>

## [Lorenz Attractor](https://github.com/verafai/python-examples/blob/main/rk_lorenz/rk_lorenz.py)
In 1963, Edward Lorenz developed a simple mathematical model of the way air moves around in the atmosphere:

![\begin{equation*}
    \begin{cases}
     \frac{dx}{dt} = \sigma (y-x) \\
     \frac{dy}{dt} =  x (\rho-z) - y\\
      \frac{dz}{dt} =  xy - \beta z
      \end{cases}  
\end{equation*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0A++++%5Cbegin%7Bcases%7D%0A+++++%5Cfrac%7Bdx%7D%7Bdt%7D+%3D+%5Csigma+%28y-x%29+%5C%5C%0A+++++%5Cfrac%7Bdy%7D%7Bdt%7D+%3D++x+%28%5Crho-z%29+-+y%5C%5C%0A++++++%5Cfrac%7Bdz%7D%7Bdt%7D+%3D++xy+-+%5Cbeta+z%0A++++++%5Cend%7Bcases%7D++%0A%5Cend%7Bequation%2A%7D%0A)

This system of three ordinary differential equations is a nonlinear dynamical system that has a very sensitive dependence on initial conditions and certain parameter values.
Now known as the Lorenz system, it describes the movement of a point in a three-dimensional space over time. Used the ![\begin{align*}4^{th}\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D4%5E%7Bth%7D%5Cend%7Balign%2A%7D%0A) order Runge Kutta method to solve the equations system and matplotlib to draw the 
attractor in two and three dimensions.


## [A model of flame propagation](https://github.com/verafai/python-examples/blob/main/Euler_NR/Euler_NR.py)
If you light a match, the ball of flame grows rapidly until it reaches a critical size. Then it remains at that size because the amount of oxygen being consumed by the combustion in the interior of the ball balances the amount available through the surface. 

![\begin{equation*}
    \begin{cases}
     \dot{x} = x^{2} - x^{3} \\
      x(0) = \delta \\
       0 \leq t \leq \tfrac{2}{\delta}
      \end{cases}  
\end{equation*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0A++++%5Cbegin%7Bcases%7D%0A+++++%5Cdot%7Bx%7D+%3D+x%5E%7B2%7D+-+x%5E%7B3%7D+%5C%5C%0A++++++x%280%29+%3D+%5Cdelta+%5C%5C%0A+++++++0+%5Cleq+t+%5Cleq+%5Ctfrac%7B2%7D%7B%5Cdelta%7D%0A++++++%5Cend%7Bcases%7D++%0A%5Cend%7Bequation%2A%7D%0A)

The scalar variable ![\begin{align*}x(t)\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7Dx%28t%29%5Cend%7Balign%2A%7D%0A) represents the radius of the ball. The ![\begin{align*}x^2\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7Dx%5E2%5Cend%7Balign%2A%7D%0A) and ![\begin{align*}x^3\end{align*}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7Dx%5E3%5Cend%7Balign%2A%7D%0A) terms come from the surface area and the volume. The critical parameter is the small initial radius ![\begin{align*}\delta\end{align*}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%5Cdelta%5Cend%7Balign%2A%7D). We seek the solution over a length of time that is inversely proportional to ![\begin{align*}\delta\end{align*}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%5Cdelta%5Cend%7Balign%2A%7D). The plot starts at ![\begin{align*}x = 0.01\end{align*}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7Dx+%3D+0.01%5Cend%7Balign%2A%7D), grows at a modestly increasing rate until t approaches 100, which is ![\begin{align*}1 / \delta\end{align*}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D1+%2F+%5Cdelta%5Cend%7Balign%2A%7D), then grows rapidly until it reaches a value close to 1, where it remains.

## [Iterative methods comparison](https://github.com/verafai/python-examples/blob/main/Jacobi-Gauss%20Seidel/JvsGS.py)
This program applies the iterative methods of Jacobi and Gauss Seidel to randomly generated matrices and vectors and counts the number of steps that must be taken to reach the specified tolerance level.



