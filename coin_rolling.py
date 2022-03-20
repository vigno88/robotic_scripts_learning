import math

r = 1
d_rot = 0
d_orientation = 0

pos = [0,0,0,0] # x y angle_rot angle_orientation
req_pos = [6,7,130,50]
print('Starting position: ', pos)
print('Required position: ', pos)
# get to the x coord
# first rotate the coin
d_orientation = 0 - pos[3]
pos[3] += d_orientation
# Move to the required x
change_x = req_pos[0]-pos[0]
d_rot = change_x/r
pos[0] += d_rot*r

# get to the y coord
# first rotate the coin
d_orientation = 90 - pos[3]
pos[3] += d_orientation
# Move to the required x
change_y = req_pos[1]-pos[1]
d_rot = change_y/r
pos[1] += d_rot*r

# get to the required rotation 
req_rot = req_pos[2] - pos[2]
# move forward
pos[2] += req_rot/2
# rotate
pos[3] = (pos[3] + 180)%360
# move "backward"
pos[2] += req_rot/2

# get the require orientation
change_orientation = req_pos[3] - pos[3]
pos[3] += change_orientation

print(pos)
