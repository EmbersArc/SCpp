import numpy as np
import os
import random
import bpy

# set current path
abspath = os.path.abspath(__file__)
dname = os.path.dirname(os.path.dirname(abspath))
os.chdir(dname)

scn = bpy.context.scene

FPS = scn.render.fps

# set output folder and get highest index 
data_folder = '../../output/Cartpole/'
data_folder += sorted(os.listdir(data_folder))[-1]

# load data
X = np.loadtxt(open(f"{data_folder}/X.txt", "rb"))
U = np.loadtxt(open(f"{data_folder}/U.txt", "rb"))
t = np.loadtxt(open(f"{data_folder}/t.txt", "rb"))

# get objects
cart = bpy.data.objects["Cart"]
pole = bpy.data.objects["Pole"]

# get timesteps, set total frames and timestep
K = X.shape[1]
trajectory_frames = FPS * t
step_size = trajectory_frames / K

# clear animation data
cart.animation_data_clear()
pole.animation_data_clear()

current_frame = 1
scn.frame_current = 1
# for each timestep in trajectory
for i in range(K):
    current_frame += step_size
    scn.frame_current = int(current_frame)

    x = X[:, i]

    # location and rotation
    cart.location[0] = x[0]
    pole.rotation_euler[1] =  x[2]
    cart.keyframe_insert(data_path='location')
    pole.keyframe_insert(data_path='rotation_euler')

scn.frame_current += FPS

# set frame range
scn.frame_start = 1
scn.frame_end = scn.frame_current
# go back to start
scn.frame_current = 1

# select all objects
bpy.context.area.type = 'VIEW_3D'
bpy.ops.object.select_all(action='SELECT')

# set all to linear interpolation
bpy.context.area.type = 'GRAPH_EDITOR'
bpy.ops.graph.select_all(action='SELECT')
bpy.ops.graph.interpolation_type(type='LINEAR')

# deselect all objects
bpy.context.area.type = 'VIEW_3D'
bpy.ops.object.select_all(action='DESELECT')

bpy.context.area.type = 'TEXT_EDITOR'
