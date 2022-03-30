# robotic_scripts_learning

`coin_rolling.py` - super naive script that move a coin to any C-Space point 
from any C-space by rolling it.

`grubler_calculator/` is a script that takes on a file that describe a mechanism using
its links and joints. We specify a list of joints with its 2 links and the program
compute the number of degrees of freedom.

`rigid_body_motion` is an inverse kinematic library in python still in progress.
Once I have it in python, the goal is to port it to c++ so it can run on an
stm32.
