from scipy.optimize import minimize, rosen, rosen_der
from scipy.optimize import OptimizeResult,approx_fprime
from matplotlib import pyplot as plt
import numpy as np
from copy import deepcopy
import random


# for this from of minimize, y is not used at all
def adam(fun, x0, args=(),beta1=0.9,beta2=0.999,eps=1E-8,eta=0.001,
         maxiter=1000000,tol=1E-5, n_data=1,batch=1, **options):

     m0 = 0 
     v0 = 0

     x = x0
     mt = m0
     vt = v0
     y = fun(x0)
     eps_grad = 1E-8
     funcalls = 1
     niter = 0
     while niter < maxiter:
        niter += 1
        random.seed(niter) 
        random_ind = random.sample(range(0,n_data),batch)
        g = approx_gradient(x, fun, eps_grad,random_ind = random_ind)
        mt = beta1*mt + (1-beta1)*g
        vt = beta2*vt + (1-beta2)*(g**2)
        mhat = mt/(1-beta1**niter)
        vhat = vt/(1-beta2**niter)
        update = -eta*mhat/(np.sqrt(vhat)+eps)
        if np.max(abs(update)) < tol: 
           y = fun(x,*args)
           return OptimizeResult(fun=y, x=x, nit=niter,
                           nfev=funcalls, success=True)
        x += update
     y = fun(x,*args)
     
     return OptimizeResult(fun=y, x=x, nit=niter,
                           nfev=funcalls, success=False)
def gd(fun, x0, args=(),  eta=0.1,
         maxiter=100,tol=1E-5,  **options):
     x = x0
     y = fun(x0)
     funcalls = 1
     eps_grad = 1E-8
     niter = 0
     while niter < maxiter:
        niter += 1
        g = approx_gradient(x, fun, eps_grad)
        if np.max(abs(g)) < tol: 
           y = fun(x,*args)
           return OptimizeResult(fun=y, x=x, nit=niter,
                           nfev=funcalls, success=True)
        x -= eta*g 

     y = fun(x,*args)
     return OptimizeResult(fun=y, x=x, nit=niter,
                    nfev=funcalls, success=False)

def gd_minibatches(fun, x0, args=(),  eta=0.1,
         maxiter=100,tol=1E-5,n_data=1,batch=1,  **options):
     x = x0
     funcalls = 1
     eps_grad = 1E-8
     niter = 0
     plot_xhi = [fun(x)] 
     plot_iter = [0]
     
     # random seed has to be changed in every iteration
     # when evaluating gradient(g) and objective function(y), same random_ind has to be used
     while niter < maxiter:
        niter += 1
        list_random_ind = create_minibatch_ind(n_data,batch,seed=niter)
        for i in list_random_ind:
            g = approx_gradient(x, fun, eps_grad,random_ind = i)
            if np.max(abs(g)) < tol: 
               y = fun(x,*args)
               return OptimizeResult(fun=y, x=x, nit=niter,
                           nfev=funcalls, success=True)
            x -= eta*g 

     y = fun(x)
     return OptimizeResult(fun=y, x=x, nit=niter,
                    nfev=funcalls, success=False)

def RMSProp(fun,x0,args=(), eta=0.1,eps=1E-6,beta=0.9,
             maxiter=1E4,tol=1E-5,n_data=1,batch=1,**options):

     x = x0
     funcalls = 1
     eps_grad = 1E-8
     A = np.zeros(len(x))
     eps = np.array([eps]*len(x))
     niter = 0
     while niter < maxiter:
        niter += 1
        random.seed(niter)
        random_ind = random.sample(range(0,n_data),batch)
        g = approx_gradient(x, fun, eps_grad,random_ind = random_ind)
        A = beta*A + (1-beta)*g**2.0
        update = -eta*np.dot(np.diag(1.0/(A**(0.5)+eps)),g)
        if np.max(abs(update)) < tol:
           y = fun(x,*args)
           return OptimizeResult(fun=y, x=x, nit=niter,
                           nfev=funcalls, success=True)
        x += update
     y = fun(x,*args)
     return OptimizeResult(fun=y, x=x, nit=niter,
                    nfev=funcalls, success=False)
       
def AdaGrad(fun,x0,args=(),eta=0.1,eps=1E-6,maxiter=1E4,tol=1E-5,n_data=1,batch=1,**options): 
     
     x = x0
     funcalls = 1
     eps_grad = 1E-8
     A = np.zeros(len(x))
     eps = np.array([eps]*len(x))
     niter = 0
     while niter < maxiter:
        niter += 1
        random.seed(niter)
        random_ind = random.sample(range(0,n_data),batch)
        g = approx_gradient(x, fun, eps_grad,random_ind = random_ind)
        A = A + g**2.0
        update = -eta*np.dot(np.diag(1.0/(A**(0.5) + eps)),g)
        if np.max(abs(update)) < tol:
           y = fun(x,*args)
           return OptimizeResult(fun=y, x=x, nit=niter,
                           nfev=funcalls, success=True)
        x += update 
     y = fun(x,*args)
     return OptimizeResult(fun=y, x=x, nit=niter,
                    nfev=funcalls, success=False)
        


def sgd(fun, x0, args=(),  eta=0.1,
         maxiter=100,tol=1E-5,n_data=1,batch=1,  **options):
     x = x0
     funcalls = 1
     eps_grad = 1E-8
     niter = 0
     
     # random seed has to be changed in every iteration
     # when evaluating gradient(g) and objective function(y), same random_ind has to be used
     while niter < maxiter:
        niter += 1
        random.seed(niter) 
        random_ind = random.sample(range(0,n_data),batch)
        g = approx_gradient(x, fun, eps_grad,random_ind = random_ind)
        if np.max(abs(eta*g)) < tol: 
           y = fun(x,*args)
           return OptimizeResult(fun=y, x=x, nit=niter,
                          nfev=funcalls, success=True)
        x -= eta*g 

     y = fun(x)
     return OptimizeResult(fun=y, x=x, nit=niter,
                    nfev=funcalls, success=False)
# Stochastic Gradient Descent + Nesterov
def sgd_mom(fun, x0, args=(),  eta=0.1,b=0.9,
         maxiter=100,tol=1E-5,n_data=1,batch=1,  **options):
     x = x0
     funcalls = 1
     eps_grad = 1E-8
     m = np.zeros(len(x))
     niter = 0
     
     # random seed has to be changed in every iteration
     while niter < maxiter:
        niter += 1
        random.seed(niter) 
        random_ind = random.sample(range(0,n_data),batch)
        g = approx_gradient(x+b*m, fun, eps_grad,random_ind = random_ind)
        m = b*m - eta*g
        if np.max(abs(m)) < tol: 
           y = fun(x,*args)
           return OptimizeResult(fun=y, x=x, nit=niter,
                          nfev=funcalls, success=True)
        x += m 

     y = fun(x)
     return OptimizeResult(fun=y, x=x, nit=niter,
                    nfev=funcalls, success=False)
# this function creates random_ind for minibatches
# it return n random_ind list where n is the number of batches
def create_minibatch_ind(ndata,batchsize,seed=0):
    np.random.seed(seed)
    ind = np.arange(0,ndata)
    np.random.shuffle(ind)
    list_random_ind = []
    n_minibatches = ndata // batchsize
    for i in range(n_minibatches):
         list_random_ind.append(list(ind[i*batchsize:(i+1)*batchsize]))
    if ndata % batchsize != 0:
         list_random_ind.append(list(ind[batchsize*n_minibatches:])) 
    return list_random_ind
         
     

# finite difference gradient that ca give argument to anonymous function
# this is something approx_fprime can't do
# if random_ind is not given, this function behaves identical to approx_fprime
def approx_gradient(x,func,eps_grad=1E-8,random_ind=[]):
    #initialize empty list
    fprime = np.zeros(len(x))
    for count_i,i in enumerate(x):
        tmp_x = deepcopy(x)
        tmp_x[count_i] = i+eps_grad
        if len(random_ind) == 0:
            fprime[count_i] = (func(tmp_x)-func(x))/eps_grad 
        else:
            fprime[count_i] = (func(tmp_x,random_ind = random_ind)-func(x,random_ind=random_ind))/eps_grad 
    return fprime


