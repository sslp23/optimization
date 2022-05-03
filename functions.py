import numpy as np
PI = np.pi

def ackley(x, a = 20, b = 2, c = 2*PI, is_graphic = False): 
    # x -> vetor
    

    d = len(x)
    first_term = -a*np.exp(-b*np.sqrt(sum(x*x)/d))
    second_term = -np.exp(sum(np.cos(c*x)/d))

        
    ackley = first_term + second_term + a + np.exp(1)
    return ackley

def ackley_plot(x, a = 20, b = 2, c = 2*PI): 
    # x -> vetor
    
    X = x[0]
    Y = x[1]
    first_term = -a*np.exp(-b*np.sqrt(X*X + Y*Y)/2)
    second_term = -np.exp((np.cos(c*X) + np.cos(c*Y))/2)

    ackley = first_term + second_term + a + np.exp(1)
    return ackley

def rastrigin(x): 
    # x -> vetor
    d = len(x)
    second_term = sum((x*x) - 10*np.cos(2*PI*x))
    rastrigin = 10*d + second_term
    return rastrigin

def rastrigin_plot(x): 
    # x -> vetor
    X = x[0]
    Y = x[1]
    second_term = ((X*X) - 10*np.cos(2*PI*X)) + ((Y*Y) - 10*np.cos(2*PI*Y))
    rastrigin = 20 + second_term
    return rastrigin

def schwefel(x):
    # x -> vetor
    d = len(x)
    second_term = sum(x*np.sin(np.sqrt(abs(x))))
    schwefel = 418.9829*d - second_term
    return schwefel

def schwefel_plot(x):
    # x -> vetor
    X = x[0]
    Y = x[1]
    
    second_term = X*np.sin(np.sqrt(abs(X))) + Y*np.sin(np.sqrt(abs(Y)))
    schwefel = 418.9829*2 - second_term
    return schwefel

def rosenbrock(x):
    # x -> vetor
    d = len(x)
    x1 = np.roll(x, -1)[:d-1]
    xt = x[:d-1]
    first_term = 100*np.power((x1 - xt*xt), 2)
    second_term = np.power((xt-1), 2)
    rosenbrock = sum(first_term + second_term)
    return rosenbrock

def rosenbrock_plot(x):
    # x -> vetor
    X = x[0]
    Y = x[1]
    
    d = len(X)
    x1 = np.roll(X, -1)[:d-1]
    xt = X[:d-1]
    d = len(Y)
    y1 = np.roll(Y, -1)[:d-1]
    yt = Y[:d-1]
    
    first_term = 100*np.power((x1 - xt*xt), 2) + 100*np.power((y1 - yt*yt), 2)
    second_term = np.power((xt-1), 2) + np.power((yt-1), 2)
    rosenbrock = first_term + second_term
    return rosenbrock

def levy_plot(x):
    # x -> vetor
    X = x[0]
    Y = x[1]
    
    d = len(X)
    w = 1 + (X-1)/4
    wi = w[:d-1]
    wd = w[d-1]
    w1 = w[0]
    
    w = 1 + (Y-1)/4
    wiy = w[:d-1]
    wdy = w[d-1]
    w1y = w[0]
    
    
    first_term = np.power(np.sin(PI*w1), 2) + np.power(np.sin(PI*w1y), 2)
    sum_aux_1 = np.power((wi-1), 2) + np.power((wiy-1), 2)
    sum_aux_2 = (1 + 10*np.power(np.sin(PI*wi + 1), 2)) + (1 + 10*np.power(np.sin(PI*wiy + 1), 2))
    second_term = sum_aux_1*sum_aux_2
    third_term = 1 + np.power(np.sin(2*PI*wd), 2) + 1 + np.power(np.sin(2*PI*wdy), 2)
    
    levy = first_term + second_term + third_term
    return levy

def levy(x):
    # x -> vetor
    d = len(x)
    w = 1 + (x-1)/4
    wi = w[:d-1]
    wd = w[d-1]
    w1 = w[0]
    
    first_term = np.power(np.sin(PI*w1), 2)
    sum_aux_1 = np.power((wi-1), 2)
    sum_aux_2 = (1 + 10*np.power(np.sin(PI*wi + 1), 2))
    second_term = sum(sum_aux_1*sum_aux_2)
    third_term = 1 + np.power(np.sin(2*PI*wd), 2)
    
    levy = first_term + second_term + third_term
    return levy
    
