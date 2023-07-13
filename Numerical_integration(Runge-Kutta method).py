import matplotlib.pyplot as plt
import math

def runge_kutta_integrate(initial_condition, final_time, f):
    #Let x(t) be the solution to differential equation dx/dt = f(x). The function outputs x(final_time).
    #initial_condition is a tuple (x0, t0) that is, x0 = x(t0)
    out = initial_condition[0]
    t=initial_condition[1]
    #finds a sufficiently small time interval so that  all higher order terms of delta_t goes negligible.
    delta_t= (final_time - t)/10**2
    while t+delta_t<final_time:
        k1 = f(out)*delta_t
        k2 = f(out + 0.5*k1)*delta_t
        k3 = f(out + 0.5*k2)*delta_t
        k4 = f(out + k3)*delta_t
        out+= (k1 + 2*k2 + 2*k3 + k4)/6
        t+= delta_t
    delta_t = final_time - t
    k1 = f(out)*delta_t
    k2 = f(out + 0.5*k1)*delta_t
    k3 = f(out + 0.5*k2)*delta_t
    k4 = f(out + k3)*delta_t
    out+= (k1 + 2*k2 + 2*k3 + k4)/6
    return out

def g(x):
    return x + math.exp( -1*x)

print(runge_kutta_integrate((0,0),1,g))

def graph_soln(infn,initial_condition):
    # given the differential equation dx/dt = f(x) and initial condition x(t_0) = x_0 - 
    # infn = f ; initial_condition = (x_0, t_0)
    #The function graphs the trajectory at (x_0,t_0).
    X = []
    start = initial_condition[1]
    while start<=initial_condition[1]+10:
        # adds different time values at a difference of 0.1 to the list X.
        X.append(start)
        start += 0.1
    Y=[]
    for i in X:
        # adds x(t) corresponding to each t in X, obtained using runge-kutta method,  to the list Y
        Y.append(runge_kutta_integrate(initial_condition,i,infn))
    plt.plot(X,Y)
    plt.show()

graph_soln(g,(0,0))