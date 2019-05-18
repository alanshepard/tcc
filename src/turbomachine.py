from scipy.optimize import root

class Turbomachine():

    def general_explicit_map(self, params, initial_guesses=None):
        
        if initial_guesses is None:
            initial_guesses = {k: v for k,v in self.__class__.DEFAULT_PARAMS.items() if k not in params.keys()}
        
        assert len(params) == self.__class__.N_FREE_PARAMS, "exactly {} parameters are required".format(self.__class__.N_FREE_PARAMS)
        assert len(initial_guesses) == len(self.__class__.DEFAULT_PARAMS)-self.__class__.N_FREE_PARAMS, "exactly {} initial guesses are required".format(len(self.__class__.DEFAULT_PARAMS)-self.__class__.N_FREE_PARAMS)
        assert params.keys() | initial_guesses.keys() == self.__class__.DEFAULT_PARAMS.keys(), "missing value for {}".format(default_params.keys() - (params.keys() | initial_guesses.keys()))
        assert params.keys().isdisjoint(initial_guesses.keys()), "A parameter cannot be supplied as both a initial guess and a fixed parameter"
        
        def fsv(x):
            kwargs = params
            kwargs.update(dict(zip(initial_guesses.keys(), x)))
            res = self.implicit_map(**kwargs)
            
            return res
          
        sol = root(fsv, list(initial_guesses.values()),method='lm', options={'eps':1e-10,})  
        params.update(dict(zip(initial_guesses.keys(), sol.x)))
        sol.params=params
        
        return sol
