import numpy as np
import matplotlib.pyplot as plt
np.seterr(divide='ignore', invalid='ignore')

#Setting the axes
plt.figure(1)
plt.axhline(y=0 , color="black")
plt.axvline(x=0, color = "black")

print(f'Enter the type of interaction - "a" for amensalism ; "cm" for commensalism; "cn" for competition; "m" for mutualism; "pm" for parasitism and "pn" for predation.')
interaction = input().lower()
interactions = {"a":[-1,0], "cm":[1,0], "cn":[-1,-1], "m":[1,1], "pm":[1,-1], "pn":[1,-1]}
while interaction not in interactions:
    interaction = input("Input the interaction correctly.   -   ").lower()

print()

r1 = float(input("Enter the intrinsic growth rate of population 1 = "))
a1 = float(input("Enter the interaction parameter of population 1 = "))
k1 = int(input("Enter the carrying capacity of population 1 = "))
r2 = float(input("Enter the intrinsic growth rate of population 2 = "))
a2 = float(input("Enter the interaction parameter of population 2 = "))
k2 = int(input("Enter the carrying capacity of population 2 = "))

def add_tuple(tuple1,tuple2):
    #The function takes two vectors in R^2 and returns their sum.
    return (tuple1[0]+tuple2[0],tuple1[1]+tuple2[1])

def add_n_tuple(list_tuple):
    #This function takes n  vectors in R^2 and returns their sum.
    temp1 = 0
    temp2 = 0
    for i in list_tuple:
        temp1 += i[0]
        temp2 += i[1]
    return (temp1,temp2)
    

def scalar_multiply(intuple,c):
    #This function takes a vector in R^2 and a scalar and returns the scalar multiple of that vector.
    return (intuple[0]*c,intuple[1]*c)

sign1 = interactions[interaction][0]
sign2 = interactions[interaction][1]

def g(intuple):
    # g is a function, g: R^2 --> R^2
    global r1,a1,k1,r2,a2,k2,sign1,sign2
    x = intuple[0]
    y = intuple [1]
    if sign1>0:
        sign1 = sign1*(1-x/k1)
    if sign2>0:
        sign2 = sign2*(1-y/k2)
    return (r1*x*(1-x/k1)+sign1*x*y*a1 , r2*y*(1 - y/k2) + sign2*x*y*a2)

def runge_kutta_solve(f,initial_condition,final_time):
    #f is a function defined from R^2 -> R^2.
    #initial = ((x(t_0),y(t_0)),t_0)
    #If X(t) is the solution of the system satisfying the initial condition then the function outputs X(final_t).
    out = initial_condition[0]
    t=initial_condition[1]
    #finds a sufficiently small time interval so that approximation holds true..
    delta_t= (final_time - t)/10**2
    while t+delta_t<final_time:
        #We use the runge-kutta method in the following steps.
        k1 = scalar_multiply(f(out),delta_t)
        k2 = scalar_multiply(f(add_tuple(out , scalar_multiply(k1,0.5))),delta_t)
        k3 = scalar_multiply(f(add_tuple(out , scalar_multiply(k2,0.5))),delta_t)
        k4 = scalar_multiply(f(add_tuple(out , k3)),delta_t)
        out = add_tuple(out , scalar_multiply(add_n_tuple([k1 , scalar_multiply(k2,2) , scalar_multiply(k3,2) , k4]),1/6))
        t+= delta_t
    delta_t = final_time - t
    k1 = scalar_multiply(f(out),delta_t)
    k2 = scalar_multiply(f(add_tuple(out , scalar_multiply(k1,0.5))),delta_t)
    k3 = scalar_multiply(f(add_tuple(out , scalar_multiply(k2,0.5))),delta_t)
    k4 = scalar_multiply(f(add_tuple(out , k3)),delta_t)
    out = add_tuple(out , scalar_multiply(add_n_tuple([k1 , scalar_multiply(k2,2) , scalar_multiply(k3,2) , k4]),1/6))
    return out

def trajectory(f,initial_condition,final_time,lower_limit, upper_limit):
    # The function graphs the trajectory in the phase plane.
    T = [initial_condition[1]]
    X = [initial_condition[0][0]]
    Y = [initial_condition[0][1]]
    t = T[0] 
    while t<=final_time:
        n = len(T) -1
        temp = runge_kutta_solve(f,((X[n],Y[n]),T[n]),t+0.01)
        X.append(temp[0])
        Y.append(temp[1])
        if temp[0]>upper_limit or temp[1]>upper_limit or temp[0]<lower_limit or temp[1]<lower_limit : #controls blow-up
            break
        T.append(T[len(T)-1]+0.01)
        t = T[len(T)-1]
    return [X,Y]


N1_0 = int(input("Enter the initial value of population 1 = "))
N2_0 = int(input("Enter the initial value of population 2 = "))
final_time = float(input("Enter the final time upto which the population growth has to be simulated = "))

def graph_N1_t(f,initial_condition,final_time):
    # The function graphs the x-t curve.
    t= initial_condition[1]
    X = []
    T = []
    while t<=final_time:
        temp = runge_kutta_solve(f,initial_condition,t)
        T.append(t)
        X.append(temp[0])
        t+=0.01
    plt.figure(1)
    plt.plot(T,X)
    plt.xlabel("t")
    plt.ylabel("N_1")
    plt.title("N_1 - t graph")
    plt.show()

graph_N1_t(g,((N1_0,N2_0),0),final_time)

def graph_N2_t(f,initial_condition,final_time):
    # The function graphs the x-t curve.
    t= initial_condition[1]
    Y = []
    T = []
    while t<=final_time:
        temp = runge_kutta_solve(f,initial_condition,t)
        T.append(t)
        Y.append(temp[1])
        t+=0.01
    plt.figure(2)
    plt.plot(T,Y)
    plt.xlabel("t")
    plt.ylabel("N_2")
    plt.title("N_2 - t graph")
    plt.show()

graph_N2_t(g,((N1_0,N2_0),0),final_time)

sign1 = interactions[interaction][0]
sign2 = interactions[interaction][1]

plt.figure(3)
plt.axhline(y=0 , color="black")
plt.axvline(x=0, color = "black")
def graph_phase_portrait(f):
    global k1,k2,sign1,sign2
    max1 = int(1.5*max([k1,k2]))
    Xn , Yn = np.meshgrid(np.linspace(0,max1,25), np.linspace(0,max1,25))
    if sign1>0:
        sign11 = sign1*(1-Xn/k1)
    else:
        sign11 = sign1
    if sign2>0:
        sign21 = sign2*(1-Yn/k2)
    else:
        sign21 = sign2
    U = r1*Xn*(1-Xn/k1)+sign11*Xn*Yn*a1
    V = r2*Yn*(1 - Yn/k2) + sign21*Xn*Yn*a2

    N = np.sqrt(U**2 + V**2)
    U = U/N
    V = V/N

    plt.quiver(Xn,Yn,U,V)
    initials = [] #List of initial conditions of the form (x(0),y(0)).
    #If you want, you can manually add the required initial conditions to the above list, if you are interested in trajectories with that initial condition.
    for i in range(1,max1+1,int(max1/6.5)):
        for j in range(1,max1+1, int(max1/6.5)):
            initials.append((i,j))
    for i in initials:
        out = trajectory(f,(i,0),5,0,max1)
        plt.plot(out[0],out[1])
    plt.title("Phase portrait")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

graph_phase_portrait(g)