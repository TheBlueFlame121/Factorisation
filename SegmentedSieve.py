# ----------------------------------------------------------------------------
#  primes.py 
#
#  This is a modified version of the script found at:
#  https://www.kaggle.com/saswatpadhi/python-sieve-48m-primes-in-4-sec
#
#  Instead of printing primes to a file, this version just returns a
#  numpy array containing all the primes. Not recommended to run for
#  n greater than 10**9.
#
# ---------------------------------------------------------------------------


import numpy as np
from math import isqrt
from multiprocessing.pool import ThreadPool


def sieve_primes(n):

    # Find all primes n > prime > 2 using the Sieve of Eratosthenes 
    # For efficiency, track only odd numbers (evens are nonprime)

    sieve = np.ones(n//2, dtype=np.bool) 
    limit = isqrt(n) + 1 

    for i in range(3, limit, 2): 
        if sieve[i//2]:
            sieve[i*i//2 :: i] = False

    prime_indexes = np.nonzero(sieve)[0][1::]
    primes  = 2 * prime_indexes.astype(int) + 1 
    return primes


def is_odd( i): 
    return i % 2

def make_odd(i, delta):
    assert delta in (-1, 1)
    return i if is_odd(i) else i + delta  


def seg_sieve_primes(seg_range):
    # Segmented Sieve of Sieve of Eratosthenes finds primes in a range.

    # As in sieve_primes(), only odd numbers are tracked (evens are nonprime)
    # So first adjust the start/end of the segement so they're odd numbers, then
    # map into the sieve as follows: [seg_start, seg_start+1*2...seg_start+n*2]
    seg_start, seg_end = seg_range
    seg_start = make_odd(seg_start, +1)
    seg_end   = make_odd(seg_end,   -1)
    seg_len   = seg_end - seg_start + 1 
    sieve_len = (seg_len + 1) // 2      # only track odds; evens nonprime
    sieve = np.ones(sieve_len, dtype=np.bool) 

    # Find a short list of primes used to strike out non-primes in the segment
    root_limit = isqrt(seg_end) + 1 
    root_primes = sieve_primes(root_limit)
    assert seg_len > root_limit

    for root_prime in root_primes:

        # find the first odd multiple of root_prime within the segment 
        prime_multiple = seg_start - seg_start % root_prime
        while not( is_odd(prime_multiple) and (prime_multiple >= seg_start) ):
            prime_multiple += root_prime

        # strike all multiples of the prime in the range...
        sieve_start = (prime_multiple - seg_start) // 2
        sieve[sieve_start : sieve_len : root_prime] = False

        # ...except for the prime itself
        if seg_start <= root_prime <= seg_end:
            ix = (root_prime - seg_start) // 2
            sieve[ix] = True

    prime_indexes = np.nonzero(sieve)[0]  
    primes  = 2 * prime_indexes.astype(int) + seg_start 
    return primes

def multiprocess_calc_primes(n, n_processes):
    # First, create seperate non-overlapping ranges so multiple 
    # processes can work on sieve segments independently 
    seg_size = n // n_processes
    seg_ranges = [(s, min(n, s+seg_size)) for s in range(2, n, seg_size)]
    if len(seg_ranges) > n_processes:
        # merge the last 2 ranges (if there's a bit left over)
        range1_start, range1_end = seg_ranges.pop()
        range2_start, range2_end = seg_ranges.pop()
        range_merged = (range2_start, range1_end)
        seg_ranges.append(range_merged)

    # Launch the processes to work on each sieve segment.
    # Each returns string of primes to write to the .CSV file
    processes = ThreadPool(n_processes)
    primes_pieces = processes.map(seg_sieve_primes, seg_ranges)

    primes = np.concatenate(primes_pieces)
    print(primes)
    return primes
