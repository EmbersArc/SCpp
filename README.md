A C++ library implementing various optimal control algorithms.

## Features

 * Efficient Successive Convexification algorithm
 * Generic SOCP MPC algorithm
 * LQR algorithm
 * JIT derivative code generation with CppAD
 * Intuitive interface to implement custom models
 

## Current Models

 * Generic Rocket Model
 * SpaceX Starship Landing Model
 

## Dependencies

 * C++17
 * Eigen
 * Boost (odeint and ptree)
 * fmt (included as submodule)
 * ECOS (included as submodule)

## Instructions

### Install

``` 
git clone --recurse-submodules https://github.com/EmbersArc/SCpp.git
cd SCpp
mkdir build
cd build
cmake ..
make
```

### Run

Available executables are:

* **MPC_sim** to simulate a trajectory with the classic MPC controller

* **SC_oneshot** to calculate one trajectory with Successive Convexification

* **SC_sim** to simulate a trajectory with Successive Convexification

Calculated trajectories are written to the `output/<modelname>` directory.

### Create a Custom Model

See existing models in the `socp_mpc/models` folder for some examples.

## Papers

* [Successive Convexification: A Superlinearly Convergent Algorithm for Non-convex Optimal Control Problems](https://arxiv.org/abs/1804.06539)

* [Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time](https://arxiv.org/abs/1802.03827)

## Examples

### Rocket trajectory model with free-final-time

![](https://thumbs.gfycat.com/DeliriousCandidAldabratortoise-size_restricted.gif)

### SpaceX Starship landing trajectory

![StarshipLanding](https://user-images.githubusercontent.com/1352472/66057427-f736be00-e538-11e9-8078-727282910f54.png)

### 2D rocket landing problem

[Feed-forward input tested in a box2d physics simulation](https://gfycat.com/DaringPortlyBlacklab)

### Cartpole

![](https://thumbs.gfycat.com/KnobbyFlatCanvasback-small.gif)

## Contributing

I'm looking forward to contributions, both problem formulations and improvements to the core library.
