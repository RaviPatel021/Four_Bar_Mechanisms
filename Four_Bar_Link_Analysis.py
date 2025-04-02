import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

def get_second_point(point, angle, length):
    x = point[0]
    y = point[1]
    angle_rad = np.radians(angle)  # Convert angle to radians
    x2 = x + length * np.cos(angle_rad)
    y2 = y + length * np.sin(angle_rad)
    return [point, [x2, y2]]

def find_angle(p1, p2):
    x1 = p1[0] 
    y1 = p1[1] 
    x2 = p2[0] 
    y2 = p2[1]
    # Calculate the difference in coordinates
    delta_y = y2 - y1
    delta_x = x2 - x1
    
    # Calculate the angle in radians using atan2
    angle_radians = math.atan2(delta_y, delta_x)
    
    # Convert radians to degrees
    angle_degrees = math.degrees(angle_radians)
    
    return angle_degrees

def find_distance(p1, p2):
    x1 = p1[0] 
    y1 = p1[1] 
    x2 = p2[0] 
    y2 = p2[1]
    # Calculate the difference in coordinates
    delta_x = x2 - x1
    delta_y = y2 - y1
    
    # Calculate the Euclidean distance
    distance = math.sqrt(delta_x**2 + delta_y**2)
    
    return distance

def solve_fourbar(theta_2, l1, l2, l3, l4, delta):
    
    theta_2_rad = np.radians(theta_2)

    K1 = l1/l2
    K2 = l1/l4
    K3 = (l2**2 - l3**2 + l4**2 + l1**2)/(2*l2*l4)
    K4 = l1/l3
    K5 = (l4**2 - l1**2 - l2**2 - l3**2)/(2*l2*l3)

    A = -K1 - K2*np.cos(theta_2_rad) + K3 + np.cos(theta_2_rad)
    B = -2*np.sin(theta_2_rad)
    C = K1 - K2*np.cos(theta_2_rad) + K3 - np.cos(theta_2_rad)
    D = np.cos(theta_2_rad) - K1 + K4*np.cos(theta_2_rad) + K5
    E = -2*np.sin(theta_2_rad)
    F = K1 + (K4-1)*np.cos(theta_2_rad) + K5

    theta_3 = 2*np.arctan2(-E-delta*np.sqrt(np.power(E,2)-4*np.multiply(D,F)),2*D)
    theta_4 = 2*np.arctan2(-B-delta*np.sqrt(np.power(B,2)-4*np.multiply(A,C)),2*A)

    theta_3 = (np.degrees(theta_3)) % 360
    theta_4 = (np.degrees(theta_4)) % 360
    
    return (theta_2), theta_3, theta_4

# Define a function to plot lines between points
def plot_line(segment, label=None, linestyle = "solid"):
    x_values = [segment[0][0], segment[1][0]]
    y_values = [segment[0][1], segment[1][1]]
    plt.plot(x_values, y_values, marker='o', label=label, linestyle = linestyle)


def fourbarvel(th_1, th_2, th_3, th_4, l2, l3, l4, omega_2):
    th_1_vec = np.deg2rad(th_1)
    th_2_vec = np.deg2rad(th_2 - th_1)
    th_3_vec = np.deg2rad(th_3 - th_1)
    th_4_vec = np.deg2rad(th_4 - th_1)

    
    # solve for the angular velocities
    omega_3 = omega_2*l2/l3*(np.sin(th_4_vec-th_2_vec))/(np.sin(th_3_vec-th_4_vec))
    omega_4 = omega_2*l2/l4*(np.sin(th_2_vec-th_3_vec))/(np.sin(th_4_vec-th_3_vec))
    
    VBA = [l3*omega_3*-np.sin(th_3_vec + th_1_vec),l3*omega_3*np.cos(th_3_vec + th_1_vec)]
    VB = [l4*omega_4*-np.sin(th_4_vec + th_1_vec),l4*omega_4*np.cos(th_4_vec + th_1_vec)]



    return VBA, VB

loop1_theta2 = 0
loop2_theta2 = 0
loop3_theta2 = 0


fig, ax = plt.subplots(figsize=(16, 9))
ax.set_xlim(-20, 20)
ax.set_ylim(-5, 25)

def update(frame):
    # lowest loop
    loop1_origin = [0,0]
    loop1_l1 = 10
    loop1_l2 = 2
    loop1_l3 = 15
    loop1_l4 = 8
    loop1_theta1 = 0
    loop1_gear = 47
    loop1_increment = 3 # changes this to increase the speed. It controls how many degrees will the first gear move by every iteration

    # middle loop
    loop2_origin = [0,8]
    loop2_l2 = 4
    loop2_l3 = 15
    loop2_l4 = 10
    loop2_gear = 43
    loop2_increment = loop1_increment * loop1_gear/loop2_gear

    # upper loop
    loop3_origin = [0,20]
    loop3_l2 = 6
    loop3_l3 = 15
    loop3_l4 = 14
    loop3_gear = 40
    loop3_increment = loop2_increment * loop2_gear/loop3_gear
    loop3_omega_2 = 1
    global loop1_theta2
    global loop2_theta2
    global loop3_theta2


    ax.clear()

    loop1_theta2, loop1_theta3, loop1_theta4 = solve_fourbar(loop1_theta2-loop1_theta1, loop1_l1, loop1_l2, loop1_l3, loop1_l4, 1)

    loop1_theta2 = loop1_theta2 + loop1_theta1
    loop1_theta3 = loop1_theta3 + loop1_theta1
    loop1_theta4 = loop1_theta4 + loop1_theta1

    loop1_l1_cord = get_second_point(loop1_origin, loop1_theta1, loop1_l1)
    loop1_l2_cord = get_second_point(loop1_origin, loop1_theta2, loop1_l2)
    loop1_l3_cord = get_second_point(loop1_l2_cord[1], loop1_theta3, loop1_l3)
    loop1_l4_cord = get_second_point(loop1_l1_cord[1], loop1_theta4, loop1_l4)



    loop2_l1 = find_distance(loop2_origin, loop1_l4_cord[1])
    loop2_theta1 = find_angle(loop2_origin, loop1_l4_cord[1])

    loop2_theta2, loop2_theta3, loop2_theta4 = solve_fourbar(loop2_theta2 - loop2_theta1, loop2_l1, loop2_l2, loop2_l3, loop2_l4, 1)

    loop2_theta2 = loop2_theta2 + loop2_theta1
    loop2_theta3 = loop2_theta3 + loop2_theta1
    loop2_theta4 = loop2_theta4 + loop2_theta1

    loop2_l1_cord = get_second_point(loop2_origin, loop2_theta1, loop2_l1)
    loop2_l2_cord = get_second_point(loop2_origin, loop2_theta2, loop2_l2)
    loop2_l3_cord = get_second_point(loop2_l2_cord[1], loop2_theta3, loop2_l3)
    loop2_l4_cord = get_second_point(loop2_l1_cord[1], loop2_theta4, loop2_l4)



    loop3_l1 = find_distance(loop3_origin, loop2_l4_cord[1])
    loop3_theta1 = find_angle(loop3_origin, loop2_l4_cord[1])

    loop3_theta2, loop3_theta3, loop3_theta4 = solve_fourbar(loop3_theta2-loop3_theta1, loop3_l1, loop3_l2, loop3_l3, loop3_l4, 1)
    

    loop3_theta2 = loop3_theta2 + loop3_theta1
    loop3_theta3 = loop3_theta3 + loop3_theta1
    loop3_theta4 = loop3_theta4 + loop3_theta1

    VBA, VB = fourbarvel(loop3_theta1, loop3_theta2, loop3_theta3, loop3_theta4, loop3_l2, loop3_l3, loop3_l4, loop3_omega_2)


    loop3_l1_cord = get_second_point(loop3_origin, loop3_theta1, loop3_l1)
    loop3_l2_cord = get_second_point(loop3_origin, loop3_theta2, loop3_l2)
    loop3_l3_cord = get_second_point(loop3_l2_cord[1], loop3_theta3, loop3_l3)
    loop3_l4_cord = get_second_point(loop3_l1_cord[1], loop3_theta4, loop3_l4)

    

    # Plot each segment
    plot_line(loop1_l1_cord, label="Loop1 L1",  linestyle = "dashed")
    plot_line(loop1_l2_cord, label="Loop1 L2")
    plot_line(loop1_l3_cord, label="Loop1 L3")
    plot_line(loop1_l4_cord, label="Loop1 L4")
    plot_line(loop2_l1_cord, label="Loop2 L1", linestyle = "dashed")
    plot_line(loop2_l2_cord, label="Loop2 L2")
    plot_line(loop2_l3_cord, label="Loop2 L3")
    plot_line(loop2_l4_cord, label="Loop2 L4")
    plot_line(loop3_l1_cord, label="Loop3 L1", linestyle = "dashed")
    plot_line(loop3_l2_cord, label="Loop3 L2")
    plot_line(loop3_l3_cord, label="Loop3 L3")
    plot_line(loop3_l4_cord, label="Loop3 L4")

    # print(loop3_l4_cord[0])

    ax.arrow(loop3_l4_cord[1][0],loop3_l4_cord[1][1],VBA[0],VBA[1],head_width=0.1,width=0.001,length_includes_head=True,linestyle='dashed',color='blue',overhang=1.0)
    ax.arrow(loop3_l4_cord[1][0],loop3_l4_cord[1][1],VB[0],VB[1],head_width=0.1,width=0.001,length_includes_head=True,linestyle='dashed',color='green',overhang=1.0)


    # Display the plot
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.legend()
    ax.grid(True)
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_title("Segment Plot")
    ax.set_ylim(-5, 35)
    ax.set_xlim(-10, 25)
    ax.set_aspect(1, adjustable='datalim')

    loop1_theta2 += loop1_increment
    loop2_theta2 -= loop2_increment
    loop3_theta2 += loop3_increment




ani = FuncAnimation(fig, update, frames=100, interval=10)
ani.save("fourbar_mechanism.gif", writer=PillowWriter(fps=30))

# plt.show()




