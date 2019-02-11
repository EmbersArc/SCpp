
# SCpp
A C++ implementation of the Successive Convexification algorithm.

**Features:**

 - Efficient Successive Convexification algorithm
 - JIT derivative code generation with CppAD
 - Intuitive interface to implement custom models
 - Rocket landing models in 2D and 3D
 
**Dependencies:**

 - C++17
 - Eigen
 - Boost (odeint and ptree)
 - fmt (included as submodule)
 - ECOS (included as submodule)
 
 **Papers:**
 - [Successive Convexification: A Superlinearly Convergent Algorithm for Non-convex Optimal Control Problems
](https://arxiv.org/abs/1804.06539)
 - [Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time
](https://arxiv.org/abs/1802.03827)
 
**Examples:**

Rocket trajectory model with free-final-time:

![](https://thumbs.gfycat.com/DeliriousCandidAldabratortoise-size_restricted.gif)

![](https://i.imgur.com/W6E0rgL.png)

[Video example of generated fixed-time trajectories](https://gfycat.com/InbornCoarseArcticseal)

- 2D rocket landing problem
[Feed-forward input tested in a box2d physics simulation](https://gfycat.com/DaringPortlyBlacklab)
