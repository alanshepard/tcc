import abc
from scipy.optimize import root

class Turbomachine(abc.ABC):

    @abc.abstractmethod
    def implicit_map(MFP1, MFP2, Mb, T0_ratio, P0_ratio):
        pass

    def general_explicit_map(self, params, initial_guesses=None):
        
        default_params = {'MFP1':0.3, 'MFP2':0.3, 'Mb':0.5, 'T0_ratio':1, 'P0_ratio':1}
        
        if initial_guesses is None:
            initial_guesses = {k: v for k,v in default_params.items() if k not in params.keys()}
        
        assert len(params) == 2, "exactly two parameters are required"
        assert len(initial_guesses) == 3, "exactly three initial guesses are required"
        assert params.keys() | initial_guesses.keys() == default_params.keys(), "missing value for {}".format(default_params.keys() - (params.keys() | initial_guesses.keys()))
        assert params.keys().isdisjoint(initial_guesses.keys()), "A parameter cannot be supplied as both a initial guess and a fixed parameter"
        
        def fsv(x):
            kwargs = params
            kwargs.update(dict(zip(initial_guesses.keys(), x)))
            res = self.implicit_map(**kwargs)
            
            return res
          
        sol = root(fsv, list(initial_guesses.values()))  
        params.update(dict(zip(initial_guesses.keys(), sol.x)))
        sol.params=params
        
        return sol
