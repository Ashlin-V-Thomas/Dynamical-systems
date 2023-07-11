import matplotlib.pyplot as plt
import math

def integrate(initial_condition, final_time, f):
    #Let x(t) be the solution to differential equation dx/dt = f(x). The function outputs X(final_time).
    #f=dx/dt
    #initial_condition is a tuple (x0, t0) that is, x0 = x(t0)
    out = initial_condition[0]
    t=initial_condition[1]
    #finds a sufficiently small time interval so that  all higher order terms of delta_t goes negligible.
    delta_t= (final_time - t)/10**2
    while t+delta_t<final_time:
        # Xn+1 = Xn + f(Xn)*delta_t
        out+= f(out)*delta_t
        t+= delta_t
    out+= f(out)*(final_time - t)
    return out

def g(x):
    return (1-x**2)**0.5

def graph_soln(infn,initial_condition):
    X = []
    start = initial_condition[1]
    while start<=initial_condition[1]+10:
        X.append(start)
        start += 0.1
    Y=[]
    for i in X:
        Y.append(integrate(initial_condition,i,infn))
    plt.plot(X,Y)
    plt.show()

graph_soln(g,(0,0))