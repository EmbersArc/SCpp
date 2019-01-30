
# SCpp
Implementation of "Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time" in C++

**Features:**

 - Basic Successive Convexification algorithm
 - JIT code generation with CppAD
 - Intuitive interface to implement custom models
 - Rocket landing models in 2D and 3D
 
**Dependencies:**

 - Eigen
 - Boost (odeint and ptree)
 - CppAD and CppADCodeGen
 - fmt (included as submodule)
 - ECOS (included as submodule)

Rocket trajectory model with free-final-time:

![](https://thumbs.gfycat.com/MenacingThornyGrackle-small.gif)

![](https://i.imgur.com/W6E0rgL.png)

[Video example of generated fixed-time trajectories](https://gfycat.com/InbornCoarseArcticseal)

- 2D rocket landing problem
[Feed-forward input tested in a box2d physics simulation](https://gfycat.com/DaringPortlyBlacklab)
