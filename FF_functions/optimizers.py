"""
"""

# Authors:

import numpy as np

class Optimizer(object):
    """ Basic stochastic gradient descent optimizer

    Parameters:


    """

    def __init__(self, Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch):
        self.Pairs = Pairs
        self.All_Pairs = All_Pairs
        self.VDW_dict = dict(VDW_dict)
        self.VDW_init = dict(VDW_init)
        self.params_dict = params_dict 
        self.epoch = epoch
        self.lambd1 = L2e
        self.lambd2 = L2s
        self.mixing_rule = mixing_rule

        self.bs = len(y_batch) # batch_size, can be different for the last batch
        self.X_batch = X_batch # Matrix consisting of (r**(-12), -r**6) pairs: 2*len(all_Pairs) X len(y_batch)
        self.y_batch = y_batch # DFT calculated energy
        # Parameters conversion from epsilon and sigma to A = a**2 B=b**2
        
        if mixing_rule == 'none': 
            Self_Terms = Pairs
        if mixing_rule in ['wh', 'lb']:
            Self_Terms = [ i for i in Pairs if i[0] == i[1] ]
        self.Self_Terms = Self_Terms # used for get_grads

        params_fit = np.zeros(len(Self_Terms)*2)
        for count_p, p in enumerate(Self_Terms):
            eps = self.VDW_dict[p][0]
            sigma = self.VDW_dict[p][1]
            params_fit[count_p*2 + 0] = 4 * eps * sigma**(12)
            params_fit[count_p*2 + 1] = 4 * eps * sigma**(6)
        self.params_fit = params_fit
        
        if mixing_rule == 'wh':
            params = np.zeros(len(All_Pairs)*2)
            for count_p, p in enumerate(All_Pairs):
                if p in Self_Terms:
                    eps = self.VDW_dict[p][0]
                    sigma = self.VDW_dict[p][1]
                else:
                    sigma = ( (self.VDW_dict[(p[0], p[0])][1]**(6.0) + self.VDW_dict[(p[1], p[1])][1]**(6.0)) / 2.0 )**(1.0/6.0)
                    eps = ( self.VDW_dict[(p[0],p[0])][0]*self.VDW_dict[(p[0],p[0])][1]**(6.0) * self.VDW_dict[(p[1],p[1])][0]*self.VDW_dict[(p[1],p[1])][1]**(6.0) )**0.5 / sigma**(6.0)
                params[count_p*2 + 0] = 4 * eps * sigma**(12.0)
                params[count_p*2 + 1] = 4 * eps * sigma**(6.0)
        
        elif mixing_rule == 'lb':
            params = np.zeros(len(All_Pairs)*2)
            for count_p, p in enumerate(All_Pairs):
                if p in Self_Terms:
                    eps = self.VDW_dict[p][0]
                    sigma = self.VDW_dict[p][1]
                else:
                    eps = ( self.VDW_dict[(p[0],p[0])][0] * self.VDW_dict[(p[1],p[1])][0] )**(0.5)
                    sigma = ( self.VDW_dict[(p[0],p[0])][1] + self.VDW_dict[(p[1],p[1])][1] ) / 2.0
                params[count_p*2 + 0] = 4 * eps * sigma**(12.0)
                params[count_p*2 + 1] = 4 * eps * sigma**(6.0)
        
        elif mixing_rule == 'none': params = np.copy(params_fit)

        # LJ energy converison using all possible pairs
        self.E_batch = self.X_batch * params


    def get_grads(self):
        
        # Calculate gradients for LJ potential
        deltaE = self.y_batch - np.sum(self.E_batch, axis=1)

        if self.mixing_rule == 'wh':
            Fit_Atoms = [i[0] for i in self.Self_Terms]
            mask_a = np.zeros([len(self.All_Pairs), len(self.Self_Terms)])
            mask_b = np.zeros([len(self.All_Pairs)*2, len(self.Self_Terms)])
            for count_a, fit_atom in enumerate(Fit_Atoms):
                ai = self.params_fit[count_a*2 + 0]**0.5
                bi = self.params_fit[count_a*2 + 1]**0.5
                for count_p, p in enumerate(self.All_Pairs):
                    if (fit_atom == p[0] and p[0] != p[1]) or (fit_atom == p[1] and p[0] != p[1]):
                        if fit_atom == p[0]: other = p[1]
                        elif fit_atom == p[1]: other = p[0]
                        
                        if (other, other) in self.Pairs:
                            aj = self.params_fit[2*Fit_Atoms.index(other) + 0]**0.5
                            bj = self.params_fit[2*Fit_Atoms.index(other) + 1]**0.5
                        elif (other, other) in self.params_dict.keys():
                            aj = self.params_dict[(other,other)][0]**0.5
                            bj = self.params_dict[(other,other)][1]**0.5
                        else:
                            print('Parameter not found')
                            quit()

                        mask_a[count_p, count_a] = (-ai / bi) * bj
                        mask_b[2*count_p + 0, count_a] = (ai**2 / bi**2 * bj - aj**2 / bj) / 2.0
                        mask_b[2*count_p + 1, count_a] = -bj
                    
                    elif fit_atom == p[0] and fit_atom == p[1]:
                        mask_a[count_p, count_a] = -2 * ai 
                        mask_b[2*count_p+1, count_a] = -2 * bi 
            
            grads_a = np.mean(2 * deltaE[:, np.newaxis] * np.matmul(self.X_batch[:, 0::2], mask_a), axis=0) # np.dot(sum(1/r12), mask_a).shape = (bs_size, len(Self_Terms) = 128, 4)
            grads_b = np.mean(2 * deltaE[:, np.newaxis] * np.matmul(self.X_batch, mask_b), axis=0) # np.dot(self.X_batch, mask_b).shape = (bs_size, len(Self_Terms) = 128, 4)
            grads = np.zeros(len(self.Self_Terms)*2)
            
            for count_a in range(len(self.Self_Terms)):
                grads[count_a*2 + 0] = grads_a[count_a]
                grads[count_a*2 + 1] = grads_b[count_a]
        
        elif self.mixing_rule == 'none':
            grads = np.mean(-4 * (self.X_batch * self.params_fit[np.newaxis,:]**0.5) * deltaE[:,np.newaxis], axis=0)
        
        # L2 regularization terms with respect to the UFF values
        grads_l2 = np.zeros(len(self.Self_Terms)*2)
        for count_p, p in enumerate(self.Self_Terms):
            
            # a and b values used for regularization 
            old_a = self.params_fit[count_p*2 + 0]**0.5
            old_b = self.params_fit[count_p*2 + 1]**0.5
            
            # L2 regularization terms with respect to the UFF values
            grads_l2[count_p*2 + 0] +=  self.lambd1 * (self.VDW_dict[p][0] - self.VDW_init[p][0]) * (-1/2 * old_b**4 * old_a**(-3)) #/ self.bs
            grads_l2[count_p*2 + 0] +=  self.lambd2 * (self.VDW_dict[p][1] - self.VDW_init[p][1]) * (1/3 * old_a**(-2/3) * old_b**(-1/3)) #/ self.bs
            grads_l2[count_p*2 + 1] +=  self.lambd1 * (self.VDW_dict[p][0] - self.VDW_init[p][0]) * (old_b**3 * old_a**(-2)) #/ self.bs
            grads_l2[count_p*2 + 1] +=  self.lambd2 * (self.VDW_dict[p][1] - self.VDW_init[p][1]) * (-4/3 * old_a**(1/3) * old_b**(-4/3)) #/ self.bs
        
        grads += grads_l2
        
        return grads


    def convert_params(self, update, lr): # Accepts lr for learning rate decay
        
        # Initialize a new VDW dictionary
        VDW_new = {}

        # Parameters conversion and update
        for p in self.params_dict.keys():
            VDW_new[p] = self.VDW_dict[p]
        
        for count_p, p in enumerate(self.Self_Terms):
            # Parameter update
            new_a = self.params_fit[count_p*2 + 0]**0.5 + update[count_p*2 + 0]
            new_b = self.params_fit[count_p*2 + 1]**0.5 + update[count_p*2 + 1]
            
            # Parameter conversion
            eps = (new_b**4.0)/(4.0*new_a**2.0) # epsilon is equal to B^2/4*A
            sigma = (new_a**2/new_b**2)**(1.0/6.0) # sigma is equal to (A/B)^(1/6)

            VDW_new[p] = (eps, sigma)
            #VDW_new[(p[1],p[0])] = (eps, sigma)
        
        for count_p, p in enumerate(self.All_Pairs):
            if (p in self.Self_Terms) == False:
                
                sigma = ( (VDW_new[(p[0], p[0])][1]**(6.0) + VDW_new[(p[1], p[1])][1]**(6.0)) / 2.0 )**(1.0/6.0)
                eps = ( VDW_new[(p[0],p[0])][0]*VDW_new[(p[0],p[0])][1]**(6.0) * VDW_new[(p[1],p[1])][0]*VDW_new[(p[1],p[1])][1]**(6.0) )**0.5 / sigma**(6.0)
                VDW_new[p] = (eps, sigma)
                VDW_new[(p[1], p[0])] = (eps, sigma)
        
        return VDW_new


class SGD(Optimizer):
    """ Stochastic gradient descent optimizer.
        SGD with momentum + learning rate decay 
        
        # Inputs
            learning_rate: Initial learning rate
            beta: Decay constant for momeuntum, recovers vanilla SGD if beta = 0
            decay_rate: Decay constant for learning rate
        # Variables
            grads: Gradient
            m: Momentum, has size equal to the size of the gradients (stored in cache1)
    """

    def __init__(self, Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch, moments, learning_rate, beta, decay_rate=0):
        super().__init__(Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch)
        self.learning_rate = learning_rate
        self.beta = beta
        self.m = moments
        self.decay = decay_rate

    def get_update(self):
        
        lr = self.learning_rate
        # Learning rate decay
        lr = lr * (1. / (1. + self.decay * self.epoch))
        # Get gradients and update momentum
        grads = self.get_grads()
        m = self.beta * self.m - lr * grads 
        update = m
        
        VDW_new = self.convert_params(update, lr)
        
        return VDW_new, m, [], []

        
class Adagrad(Optimizer):
    """ Adaptive subgradient optimizer.

        # Inputs
            learning_rate: Initial learning rate
            eps: Avoid dividing by zero (1e-4 ~ 1e-8)
        # Variables
            grads: Gradient
            cache: Has size equal to the size of the gradient and keeps track of per-parameter sum of square gradients
    
        Reference: Duchi et al., 2011
    """

    def __init__(self, Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch, cache, learning_rate, eps=1E-8):
        super().__init__(Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch)
        self.learning_rate = learning_rate
        self.cache = cache
        self.eps = eps #np.array([eps] * (len(Pairs)*2))
        
    def get_update(self):
        
        lr = self.learning_rate

        # Get gradients
        grads = self.get_grads()
        
        # Update the second momentum term
        cache = self.cache 
        cache += grads**(2.0) # second moment
        update = -lr * np.dot(np.diag(1.0/(cache**(0.5) + self.eps)), grads)
        #update = -lr * grads / (np.sqrt(cache) + self.eps)

        VDW_new = self.convert_params(update, lr)

        return VDW_new, cache, [], []


class RMSProp(Optimizer):
    """ Root Mean Square Prop

        # Inputs
            learning_rate: Initial learning rate
            beta: Decay constant for the norm
            eps: Avoid dividing by zero
        # Variables
            grads: Gradient
            cache: Has size equal to the size of the gradient and keeps track of per-parameter sum of square gradients
    
        Reference: Tieleman and Hinton, 2012
    """

    def __init__(self, Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch, cache, learning_rate, beta=0.9, eps=1E-8):
        super().__init__(Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch)
        self.learning_rate = learning_rate
        self.cache = cache

        self.beta = beta
        self.eps = eps #np.array([eps] * (len(Pairs)*2))
        
    def get_update(self):
        
        lr = self.learning_rate

        # Get gradients
        grads = self.get_grads()
        
        # Update the second momentum term
        cache = self.beta * self.cache + (1.0 - self.beta) * grads**(2.0)
        update = -lr * np.dot(np.diag(1.0/(cache**(0.5) + self.eps)), grads)
        #update = -lr * grads / (np.sqrt(cache) + self.eps)

        VDW_new = self.convert_params(update, lr)

        return VDW_new, cache, [], []


class Adam(Optimizer):
    """ Adaptive moment estimation optimizer.
        RMSProp + Momentum

        # Inputs
            learning_rate: Initial learning rate (Default: 0.001)
            beta1 : Exponential decay rate for moving average of first moment of gradient (Default: 0.9)
            beta2 : Exponential decay rate for moving average of second moment of gradient (Default: 0.999)
            eps: Avoid dividing by zero (Default: 1E-8)
        # Variables
            grads: Gradient
            m: Moving average of first moment estimate, stored in cache1
            mt: Bias-corrected first moment
            v: Moving average of second moment estimate, stored in cache2
            vt: Bias-corrected second moment
   
        Reference: Kingma and Ba, 2015
    """

    def __init__(self, Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch, cache1, cache2, learning_rate, beta1 = 0.9, beta2 = 0.999, eps = 1E-8):
        super().__init__(Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch)
        # Hyperparameters
        self.learning_rate = learning_rate
        self.beta1 = beta1 
        self.beta2 = beta2 

        self.m = cache1 
        self.v = cache2 
        self.eps = eps 

    def get_update(self):
        
        lr = self.learning_rate

        # Get gradients
        grads = self.get_grads()
        
        # Main parts
        m = self.beta1 * self.m + (1-self.beta1) * grads # update first momentum
        mt = m / (1-self.beta1 ** self.epoch) # bias corrected m
        v = self.beta2 * self.v + (1-self.beta2) * (grads**2) # update second moemntum
        vt = v / (1-self.beta2 ** self.epoch) # bias corrected v
        update = -lr * mt / (np.sqrt(vt) + self.eps)
         
        VDW_new = self.convert_params(update, lr)
        #print(VDW_new)

        return VDW_new, m, v, []


class Nadam(Optimizer):
    """ Nesterov-accelerated adaptive moment estimation optimizer.

        # Inputs
            learning_rate: Initial learning rate
            beta1 : Decay rate for moving average of first moment of gradient (Default: 0.99)
            beta2 : Decay rate for moving average of second moment of gradient (Default: 0.999)
            eps: Avoid dividing by zero (referred to Adam paper: 1E-8)
        # Variables
            m_schedule: [Sutskever et al., "On the importance of initialization and momentum in deep learning", 2015]
            mtb: 
   
        Reference: Dozat, 2015
                   Keras documentation for optimizers.py
    """

    def __init__(self, Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch, cache1, cache2, cache3, learning_rate, beta1 = 0.99, beta2 = 0.999, eps = 1E-8):
        super().__init__(Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch)
        # Hyperparameters
        self.learning_rate = learning_rate
        self.beta1 = beta1 
        self.beta2 = beta2 

        self.m = cache1 
        self.v = cache2
        self.m_schedule = cache3
        if self.epoch is 1:
            self.m_schedule = np.ones(len(self.Self_Terms)*2)
        self.eps = eps 

    def get_update(self):
        
        lr = self.learning_rate

        # Get gradients 
        grads = self.get_grads()

        # Warming momentum schedule. m_schedule suggested in the original paper, originally in Sutskever et al, 2013
        mu_t = self.beta1 * (1.- 0.5*0.96**(self.epoch / 250.0)) 
        mu_tp1 = self.beta1 * (1 - 0.5*0.96**((self.epoch +1) / 250.0)) 
        m_schedule = self.m_schedule * mu_t
        m_schedule_next = self.m_schedule * mu_t * mu_tp1
        
        # Main parts
        gradst = grads / (1. - m_schedule)
        m = self.beta1 * self.m + (1. - self.beta1) * grads
        mt = m / (1. - m_schedule_next) 
        v = self.beta2 * self.v + (1. - self.beta2) * (grads**2) 
        vt = v / (1. - self.beta2**self.epoch) 
        mtb = (1. - mu_t) *gradst + mu_tp1 * mt        
        update = -lr * mtb / (np.sqrt(vt) + self.eps)

        VDW_new = self.convert_params(update, lr)

        return VDW_new, m, v, m_schedule


class AMSGrad(Optimizer):
    """ Exponential moving average variant with guaranteed convergence

        # Inputs
            learning_rate: Initial learning rate (Default: 0.001)
            beta1 : Exponential decay rate for moving average of first moment of gradient (Default: 0.9)
            beta2 : Exponential decay rate for moving average of second moment of gradient (Default: 0.99 or 0.999)
            eps: Avoid dividing by zero (Default: 1E-8)
        # Variables
            grads: Gradient
            m: Moving average of first moment estimate, stored in cache1
            v: Moving average of second moment estimate, stored in cache2
   
        Reference: Reddi et al., 2018
    """

    def __init__(self, Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch, cache1, cache2, cache3, learning_rate, beta1 = 0.9, beta2 = 0.99, eps = 1E-8):
        super().__init__(Pairs, All_Pairs, VDW_dict, VDW_init, params_dict, epoch, L2e, L2s, mixing_rule, X_batch, y_batch)
        # Hyperparameters
        self.learning_rate = learning_rate
        self.beta1 = beta1 
        self.beta2 = beta2 

        self.m = cache1 
        self.v = cache2
        self.v_hat = cache3
        self.eps = eps 

    def get_update(self):
        
        lr = self.learning_rate

        # Get gradients
        grads = self.get_grads()
        
        # Main parts
        m = self.beta1 * self.m + (1-self.beta1) * grads
        v = self.beta2 * self.v + (1-self.beta2) * (grads**2)
        v_hat = np.maximum(self.v_hat, v)

        update = -lr * m / (np.sqrt(v_hat) + self.eps)

        VDW_new = self.convert_params(update, lr)

        return VDW_new, m, v, v_hat
