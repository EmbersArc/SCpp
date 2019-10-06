import numpy as np
import bpy
from os import listdir

scene = bpy.context.scene

PATH = bpy.path.abspath("//../../output/Starship/")
FPS = scene.render.fps


scene.frame_current = 0


path = PATH + sorted(listdir(PATH))[-1]
print(path)
try:
    X = np.genfromtxt(path+"/X.txt", delimiter=',')
    U = np.genfromtxt(path+"/U.txt", delimiter=',')
    t = np.genfromtxt(path+"/t.txt", delimiter=',')
except OSError:
    print("Data not found.")
    #continue

K = X.shape[0]
dt = int(FPS * t / K)
scene.frame_end = int(t * FPS)

body_ob = bpy.data.objects.get("Starship")
eng1_ob = bpy.data.objects.get("Engine1")
eng2_ob = bpy.data.objects.get("Engine2")
eng3_ob = bpy.data.objects.get("Engine3")
fir1_ob = bpy.data.objects.get("Fire1")
fir2_ob = bpy.data.objects.get("Fire2")
fir3_ob = bpy.data.objects.get("Fire3")
body_ob.animation_data_clear()
eng1_ob.animation_data_clear()
eng2_ob.animation_data_clear()
eng3_ob.animation_data_clear()
fir1_ob.animation_data_clear()
fir2_ob.animation_data_clear()
fir3_ob.animation_data_clear()

T_max = np.max(np.linalg.norm(U, axis=1))

for k in range(K):
    scene.frame_current = dt * k
    x = X[k]
    u = U[k]
    body_ob.location = x[1:4] / 100
    body_ob.rotation_quaternion = x[7:11]
    body_ob.keyframe_insert(data_path='location')
    body_ob.keyframe_insert(data_path='rotation_quaternion')

    rx = np.arctan(-u[1] / u[2])
    ry = np.arctan(u[0] / u[2])
    min_l = 0.6
    l = min_l + (1.-min_l) * (np.linalg.norm(u) / T_max)
    for eng in [eng1_ob, eng2_ob, eng3_ob]:
        eng.rotation_euler = (rx, ry, 0)
        eng.keyframe_insert(data_path='rotation_euler')
    for fir in [fir1_ob, fir2_ob, fir3_ob]:
        fir.scale[2] = l
        fir.keyframe_insert(data_path='scale')
        
scene.frame_current += 5
for fir in [fir1_ob, fir2_ob, fir3_ob]:
    fir.scale[2] = 0
    fir.keyframe_insert(data_path='scale')

scene.frame_current = 0
