import json
import numpy as np
from scipy.special import gamma, gammaincc
from cobaya.likelihood import Likelihood

def inc_gamma(a, x):
    if a <= 0:
        raise ValueError("a must be positive")
    
    first = gamma(a)

    if np.isinf(first): return 0
    
    second = gammaincc(a, x)
    prod = first * second
    return prod

def weibull_pdf(x, k, r):
    first = k/x
    second = (x*r*gamma(1 + k**(-1)))**k

    if np.isinf(second):  return 0
    
    third = np.exp(-second)
    prod = first*second*third
    return prod

def weibull_ccdf(x, k, r):
    return np.exp(-(x*r*gamma(1 + k**(-1)))**k)

class WeibullLikelihood(Likelihood):
    data_file: str

    def initialize(self):
        self.obs_start_list = []
        self.obs_end_list = []
        self.burst_mjds_list = []
        
        with open(self.data_file, 'r') as f:
            data = json.load(f)

        for epoch in data:
            self.obs_start_list.append(epoch['start'])
            self.obs_end_list.append(epoch['end'])
            self.burst_mjds_list.append(np.array(epoch['bursts']))

        # Data Validation
        for i in range(len(self.obs_start_list)):
            start = self.obs_start_list[i]
            end = self.obs_end_list[i]
            if i < len(self.burst_mjds_list):
                 bursts = self.burst_mjds_list[i]
                 if len(bursts) > 0:
                      if np.any(bursts < start) or np.any(bursts > end):
                           bad_bursts = bursts[(bursts < start) | (bursts > end)]
                           raise ValueError(f"Bursts found outside observation interval in observation {i+1}. "
                                            f"Observation range: [{start:.5f}, {end:.5f}]. "
                                            f"Out of bounds bursts: {bad_bursts}")

    def logp(self, k, logr, **kwargs):
        r = 10**logr
        lnposterior = 0.0

        for i in range(len(self.obs_start_list)):
            obs_start = self.obs_start_list[i]
            obs_end = self.obs_end_list[i]
            burst_mjds = self.burst_mjds_list[i]

            obs_len = (obs_end - obs_start)*24 # in hours
            mjds_sorted = np.sort(burst_mjds)
            mjds_sorted = (mjds_sorted - obs_start)*24 # in hours
            
            n = len(mjds_sorted)
            
            logL = 0.0
            
            if n == 0:
                try:
                    numerator = inc_gamma(1/k, (obs_len*r*gamma(1 + k**(-1)))**k)
                    if numerator == 0:
                        logL = -1e32
                    else:
                        denominator = k*gamma(1 + k**(-1))
                        prod = numerator/denominator
                        if prod <= 0 or np.isnan(prod):
                            logL = -1e32
                        else:
                            logL = np.log(prod)
                except Exception:
                     logL = -1e32
            
            elif n == 1:
                t1 = mjds_sorted[0]
                ccdf1 = weibull_ccdf(t1, k, r)
                ccdf2 = weibull_ccdf(obs_len - t1, k, r)
                
                if ccdf1 <= 0 or ccdf2 <= 0:
                    logL = -1e32
                else:
                    logL = np.log(r) + np.log(ccdf1) + np.log(ccdf2)
            
            elif n > 1:
                t1 = mjds_sorted[0]
                tn = mjds_sorted[-1]
                
                ccdf1 = weibull_ccdf(t1, k, r)
                ccdf2 = weibull_ccdf(obs_len - tn, k, r)
                
                if ccdf1 <= 0 or ccdf2 <= 0:
                     logL = -1e32
                else:
                    logL += np.log(r)
                    logL += np.log(ccdf1)
                    logL += np.log(ccdf2)
                    
                    diffs = np.diff(mjds_sorted)
                    for diff in diffs:
                        pdf = weibull_pdf(diff, k, r)
                        if pdf <= 0 or np.isnan(pdf):
                            logL = -1e32
                            break
                        logL += np.log(pdf)
            
            lnposterior += logL
            
        return lnposterior
