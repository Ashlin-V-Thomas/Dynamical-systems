import math

def integrate(initial_condition, final_time, f):
    #Let x(t) be the solution to differential equation dx/dt = f(x). The function outputs x(final_time).
    #initial_condition is a tuple (x0, t0) that is, x0 = x(t0)
    out = initial_condition[0]
    t=initial_condition[1]
    #finds a sufficiently small time interval so that our approximation holds .
    delta_t= (final_time - t)/10**2
    while t+delta_t<final_time:
        # Xn+1 = Xn + f(Xn)*delta_t
        out+= f(out)*delta_t
        t+= delta_t
    out+= f(out)*(final_time - t)
    return out

def g(x):
    return x + math.exp(-x)

print(integrate((0,0), 1, g))