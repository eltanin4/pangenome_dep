import timeit
import numpy as np
from collections import Counter
from tqdm import tqdm

def get_sample(arr, n_iter=None, sample_size=10, 
               fast=True):
    """Get random sample from arr.
    
    Parameters
    ----------
    arr: np.array
        array to sample from.
    n_iter: int
        current iteration number.
    sample_size: int
        sample size
    fast: bool
        use sampling optimized for fast consecutive samples 
        from the same array.
    
    Returns
    -------
    sample: np.array
        sample from arr of length n_iter.
    """
    if fast:
        # find the index we last sampled from
        start_idx = (n_iter * sample_size) % n
        if start_idx + sample_size >= n:
            # shuffle array if we have reached the end and repeat again
            np.random.shuffle(arr)
            
        return arr[start_idx:start_idx+sample_size] 
    else:
        return np.random.choice(arr, sample_size, replace=False)
    
def collect_samples(arr,
                    sample_size,
                    n_samples,
                    fast=False):
    """
    Collect several samples from arr.
    
    Parameters
    ----------
    arr: np.array
        array to sample from.
    sample_size: int
        sample size.
    n_samples: int
        number of samples to take.
    fast: bool
        use sampling optimized for fast consecutive samples 
        from the same array.
    
    Returns
    -------
    samples: np.ndarray
        sample matrix of shape (n_samples, sample_size)
    """
    samples = np.zeros((n_samples + 1, sample_size), np.int32)
    
    for sample_n in range(0, n_samples):
        sample = get_sample(arr, 
                            n_iter=sample_n,
                            sample_size=sample_size,
                            fast=fast)
        samples[sample_n] = sample
        
    return samples


