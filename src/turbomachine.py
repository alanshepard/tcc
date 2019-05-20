import numpy as np
from scipy.optimize import root
from collections import OrderedDict
import random
import itertools
import matplotlib.pyplot as plt

def gridmap(t, x_grid, pr_grid, x_param='MFP1', plot=False):
    """ 
    t: Turbomachine
    """
    assert np.shape(x_grid) == np.shape(pr_grid)


    # Initialize storage
    success = np.full_like(x_grid, False, dtype=bool)
    nfev = np.zeros_like(x_grid, dtype=int)
    tries = np.zeros_like(x_grid, dtype=int)
    params = dict([(k, np.full_like(x_grid, np.nan)) for k in t.initial_guess.keys()])


    # Find a random point which converges
    ijgrid = list(itertools.product(*map(range, x_grid.shape)))
    random.shuffle(ijgrid)
    for i,j in ijgrid:
        sol = t.general_explicit_map({x_param: x_grid[i,j], 'P0_ratio': pr_grid[i,j]})
        if sol.success: 
            print(i,j)
            break

    ### Breadth-first search ###
    assert sol.success, "Make sure initial solution is successful"
    
    # List to store search frontier 
    queue = OrderedDict() #substitute this by a pool to paralelize code (use concurrent.futures)
    queue[(i,j)] = sol
    
    def visit(ii,jj,sol):
        if 0<=ii<x_grid.shape[0] and 0<=jj<x_grid.shape[1]: # Only visit valid addresses
            if not success[ii,jj]:
                
                if (ii,jj) in queue:
                    if queue[(ii,jj)].success:
                        return
                    else:
                        del queue[(ii,jj)]
                        
                assert (ii,jj) not in queue
                sol=t.general_explicit_map(params={x_param: x_grid[ii,jj], 'P0_ratio': pr_grid[ii,jj]},
                                           initial_guesses=sol.params)
                queue[(ii, jj)] = sol
                tries[ii,jj]+=1
                
    while queue:
        (i,j), sol = queue.popitem(last=False)
        
        success[i,j] = sol.success
        nfev[i,j] = sol.nfev
        
        
        for k in sol.params.keys():
            params[k][i,j] = sol.params[k] if sol.success else np.nan
        
        if sol.success:
            visit(i+1,j  ,sol)
            visit(i-1,j  ,sol)
            visit(i  ,j+1,sol)
            visit(i  ,j-1,sol)
            visit(i+1,j+1,sol)
            visit(i+1,j-1,sol)
            visit(i-1,j-1,sol)
            visit(i-1,j+1,sol)

    params['success'] = success
    params['nfev'] = nfev
    params['tries'] = tries

    if plot:
        map_imshow(x_grid, pr_grid, params)

    return params

def map_imshow(x_grid, pr_grid, params):
    fig, axes = plt.subplots(3,len(params)//3+1,squeeze=False)
    for (name, data), ax in zip(params.items(), axes.flatten()):
        im=ax.imshow(data, origin='lower')
        ax.set_title(name)
        fig.colorbar(im, ax=ax)

class Turbomachine():

    def __init__(self):
        self.initial_guess = self.__class__.DEFAULT_PARAMS.copy()

    def general_explicit_map(self, params, initial_guesses=None):
        
        if initial_guesses is None:
            initial_guesses = self.initial_guess.copy()
        
        initial_guesses = {k: v for k,v in initial_guesses.items() if k not in params.keys()}
        
        assert len(params) == self.__class__.N_FREE_PARAMS, "exactly {} parameters are required".format(self.__class__.N_FREE_PARAMS)
        assert len(initial_guesses) == len(self.__class__.DEFAULT_PARAMS)-self.__class__.N_FREE_PARAMS, "exactly {} initial guesses are required".format(len(self.__class__.DEFAULT_PARAMS)-self.__class__.N_FREE_PARAMS)
        assert params.keys() | initial_guesses.keys() == self.__class__.DEFAULT_PARAMS.keys(), "missing value for {}".format(default_params.keys() - (params.keys() | initial_guesses.keys()))
        assert params.keys().isdisjoint(initial_guesses.keys()), "A parameter cannot be supplied as both a initial guess and a fixed parameter"
        
        def fsv(x):
            kwargs = params
            kwargs.update(dict(zip(initial_guesses.keys(), x)))
            res = self.implicit_map(**kwargs)
            
            return res
          
        sol = root(fsv, list(initial_guesses.values()),method='hybr', options={'eps':1e-10,})  
        params.update(dict(zip(initial_guesses.keys(), sol.x)))
        sol.params=params
        
        return sol
