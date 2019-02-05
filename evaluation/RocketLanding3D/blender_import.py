import numpy as np
import os
import random
import bpy

# set current path
abspath = os.path.abspath(__file__)
dname = os.path.dirname(os.path.dirname(abspath))
os.chdir(dname)

scn = bpy.context.scene

# maximum thrust for scaling
thrust_magnitude = 800000

FPS = scn.render.fps

# set output folder and get highest index
data_folder = '../../output/RocketLanding3D'

# load data
last_output = str(max(int(t[:3]) for t in os.listdir(data_folder))).zfill(3)
X = np.loadtxt(open(f"{data_folder}/{last_output}_X.txt", "rb"))
U = np.loadtxt(open(f"{data_folder}/{last_output}_U.txt", "rb"))
t = np.loadtxt(open(f"{data_folder}/{last_output}_t.txt", "rb"))

# get objects
rck = bpy.data.objects["rck"]
eng = bpy.data.objects["eng"]
fir = bpy.data.objects["fir"]
eng_light = bpy.data.lights.get("Point")

# get timesteps, set total frames and timestep
K = X.shape[1]
trajectory_frames = FPS * t
step_size = int(trajectory_frames / K)

scn.frame_current = 1
# for each timestep in trajectory
for i in range(K):
    x = X[:, i]
    u = U[:, i]

    # location and rotation
    rck.location = np.array((x[1], x[2], x[3]))
    rck.rotation_quaternion = (x[7], x[8], x[9], x[10])
    rck.keyframe_insert(data_path='location')
    rck.keyframe_insert(data_path='rotation_quaternion')

    # thrust magnitude and angles
    rx = np.arctan(-u[1] / u[2])
    ry = np.arctan(u[0] / u[2])
    eng.rotation_euler = (rx, ry, 0)
    fir.scale[2] = np.linalg.norm(u) / thrust_magnitude
    eng_light.energy = fir.scale[2] * 10

    eng_light.keyframe_insert(data_path='energy')
    eng.keyframe_insert(data_path='rotation_euler')
    fir.keyframe_insert(data_path='scale')

    legs = ['legA', 'legB', 'legC', 'legD']
    # legs default position
    if i == int(K / 3 * 2):
        for l in legs:
            leg = bpy.data.objects[l]
            leg.rotation_euler[0] = 0.
            leg.keyframe_insert(data_path='rotation_euler')
    # legs deployed position
    if i == int(K / 5 * 4):
        for l in legs:
            leg = bpy.data.objects[l]
            leg.rotation_euler[0] = 120 / 180 * np.pi
            leg.keyframe_insert(data_path='rotation_euler')

    scn.frame_current += step_size

scn.frame_current += FPS

# shutdown engine
fir.scale[2] = 0
eng_light.energy = 0
fir.keyframe_insert(data_path='scale')
eng_light.keyframe_insert(data_path='energy')

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
