# Simulation of the Heat Equation
This repository will contain a simulation of the 2D heat equation in a rectangular room. 
The temperatures of each wall will be chosen in beforehand, and a simulation using the Finite Difference Method will show how the temperature develops in the room over time.

## The math behind the finite-difference-method
The heat equation is a partial differential equation that describes how the temperature $u$ changes with respect to time and position. The equation in 2D can be states as:

$$
u_t = c^2 (u_{xx}+u_{yy})
$$

In order to simulate a system that satisifes the heat equation, we can make use of a numerical method called the finite-difference-method (FDM for short). We recall that the definition of the derivative comes from:

$$ u'(x) = \lim_{h\rightarrow 0} \frac{u(x+h)-u(x)}{h} $$

For small $h$, we can approximate this derivative as 

$$ u'(x) \approx \frac{u(x+h)-u(x)}{h}. $$

For the second derivative we utilize the above approximation again to get:

$$ u''(x) \approx \frac{u'(x+h)-u'(x)}{h} = \frac{u(x+h)-2u(x)+u(x-h)}{h^2} $$

In the 2D case, we can discretize our spatial coordinates with increments of $\Delta x$ and $\Delta y$ respectively. This just means replacing the $h$'s in the equation above.

Let $u^m_{i,j}$ be the temperature $u$ at position $i,j$ in our discretization at iteration ("time") $m$. Substituting the ideas developed above into the heat equation, we get:

$$\frac{u^{m+1}_{i,j}-u^{m}_{i,j}}{\Delta t} = c^2 \left(\frac{u^{m}_{i+1,j}-2u^{m}_{i,j}+u^{m}_{i-1,j}}{\Delta x^2} + \frac{u^{m}_{i,j+1}-2u^{m}_{i,j}+u^{m}_{i,j-1}}{\Delta y^2}\right) $$

For which we the ncan solve the temperature at position $i,j$ for time $m+1$ as:

$$u^{m+1}_{i,j} = u^{m}_{i,j} + \Delta t c^2 \left(\frac{u^{m}_{i+1,j}-2u^{m}_{i,j}+u^{m}_{i-1,j}}{\Delta x^2} + \frac{u^{m}_{i,j+1}-2u^{m}_{i,j}+u^{m}_{i,j-1}}{\Delta y^2}\right) $$


Which is precisely what is utilized in the simulation.  Notice that in the code, $\Delta x$ and $\Delta y$ are assumed to be the same for simplicity. 

# Simulation of the Wave Equation

This part is not done yet!

## User instructions

In order to use the code, use the following image as a template.

<p align="center">
    <img src="https://user-images.githubusercontent.com/121384892/215316805-f38d2595-0808-485a-b122-0011875258bb.png" width="500" height="400">
    
    The BC[i] corresponds to the i:th boundary condition in correspondance to the image above. rect[i] stands for the $i:th$ coordinate where for i=0,1 is used for x-coordinates and for i=2,3 is used for the y-coordinates. 
</p>

Below is an image from the code itself:

<p align="center">
    <img src="https://user-images.githubusercontent.com/121384892/215294698-77c83a8f-ed9a-4985-8e77-0414350c6bfc.png" width="500" height="400">
    
    This is an image from the simulation where the initial conditions for the heat equation is given by u(x,y,0) = 3-x^2-y^2 in a 
    10x10 room with boundary conditions given by 0 degrees across all walls.
</p>
