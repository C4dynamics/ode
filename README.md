![](https://github.com/C4dynamics/ode/blob/main/rocketsim%203.gif)

The starting point to any dynamical system, namely a system that evolves in time, is the mathematical representation of the evolution in time. This representation is called the model of the system.

The mathematical representation of most systems is a set of differential equations. To solve this set of equations, matlab, python, and other programs, provide tools which allow to write the equations and to solve them within the desired times.

The program here simulates the dynamics of a rocket. The model is known as 'short-period', because it represents the behavior of a rocket or other flying object in the short time after any change in the flight conditions.

The program is written in matlab and the function that holds the equations is called dx_nonlinear(). The concept of object oriented programming that used here means that we first generate an entity, let it be a rocket. 
In the matlab command window write:

rckt = rocket_dynamics();

Now the variable rckt is a rocket object that has properties and functionality of rocket:

rckt.v 

Returns the velocity of the rocket.
The command: 

rckt.x0(3) = 15 * pi / 180;

Means to start the simulation with the fins of the rocket deflected by 15 degrees. 

rckt.run_sim(2);

Runs the rocket for 2 seconds.

rckt.drawstate();

Plots the properties of the system in the time of simulation.

rckt.ranimate();

Animates the motion of the system in the time of simulation. 
gif edited with: https://ezgif.com/cut/ezgif-4-1e68375c6e.gif

# ode45
# dynamical systems
# object oriented 
# short period
# ordinary differential equations 
# nonlinear dynamics 



