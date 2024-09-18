"""
Function to calculate EOS parameters using critical parameter values
"""
import numpy as np

def calculate_a_CS(p_c, T_c, R):
    return 0.4963 * (R*T_c)**2 / p_c

def calculate_b_CS(p_c, T_c, R):
    return 0.18727 * (R*T_c) / p_c

def calculate_a_PR(p_c, T_c, R):
    return 0.45724 * (R*T_c)**2 / p_c

def calculate_b_PR(p_c, T_c, R):
    return 0.0778 * (R*T_c) / p_c

def calculate_a_RK(p_c, T_c, R):
    return 0.42748 * R**2 * T_c**2.5 / p_c

def calculate_b_RK(p_c, T_c, R):
    return 0.08664 * (R*T_c) / p_c

def calculate_a_RKS(p_c, T_c, R):
    return 0.42748 * (R*T_c)**2 / p_c

def calculate_b_RKS(p_c, T_c, R):
    return 0.08664 * (R*T_c) / p_c

def calculate_a_b_vdW(p_c, T_c, R):
    b = (R * T_c) / (8 * p_c)
    a = (27/8) * R * b * T_c
    return a, b

def generate_intervals(f, start=0, end=5, num_intervals=100):
    """
    Generate three distinct intervals where the function f(x) changes sign.
    
    Parameters:
    - f: the function for which to find intervals with sign changes
    - start: start of the range to search for intervals
    - end: end of the range to search for intervals
    - num_intervals: number of points to sample between start and end
    
    Returns:
    - intervals: a list of three intervals where the function changes sign
    """
    # Generate evenly spaced points
    x_vals = np.linspace(start, end, num_intervals)
    
    # Initialize list to store intervals
    intervals = []
    
    # Loop through points to find intervals where f(x) changes sign
    for i in range(len(x_vals) - 1):
        x1, x2 = x_vals[i], x_vals[i + 1]
        f1, f2 = f(x1), f(x2)
        # Check if the function changes sign between x1 and x2
        if f1 * f2 < 0:
            intervals.append((x1, x2))
        
        # If we already have 3 intervals, stop
        if len(intervals) == 3:
            break
    
    # # If fewer than 3 intervals found, raise an error
    # if len(intervals) < 3:
    #     raise ValueError("Could not find 3 intervals with sign changes.")
    
    return intervals
