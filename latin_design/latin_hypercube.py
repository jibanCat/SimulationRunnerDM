"""
The original file is from git:sbird/lyaemu

This file contains functions which pick a set of samples from a parameter space
which will allow a Gaussian process to best interpolate the samples to new positions in parameter space.
Several schemes for this are possible.

We use rejection-sampled latin hypercubes.
"""

import numpy as np

def remove_single_parameter(center, prior_points):
    """Remove all values within cells covered by prior samples for a particular parameter.
    Arguments:
    center contains the central values of each (evenly spaced) bin.
    prior_points contains the values of each already computed point."""
    #Find which bins the previously computed points are in
    already_taken = np.array([np.argmin(np.abs(center - pp)) for pp in prior_points])
    #Find the indices of points not in already_taken
    not_taken = np.setdiff1d(range(np.size(center)), already_taken)
    new_center = center[not_taken]
    assert np.size(new_center) == np.size(center) - np.size(prior_points)
    return new_center,not_taken

def lhscentered(n, samples, prior_points = None):
    """
    Generate a latin hypercube design where all samples are
    centered on their respective cells. Can specify an already
    existing set of points using prior_points; these must also
    be a latin hypercube on a smaller sample, but need not be centered.
    """
    #Set up empty prior points if needed.
    if prior_points is None:
        prior_points = np.empty([0,n])

    npriors = np.shape(prior_points)[0]
    # Generate the intervals
    cut = np.linspace(0, 1, samples + 1)

    # Fill points uniformly in each interval
    # Number of stratified layers used is samples desired + prior_points.
    a = cut[:samples]
    b = cut[1:samples + 1]
    #Get list of central values
    _center = (a + b)/2
    # Choose a permutation so each sample is in one bin for each factor.
    H = np.zeros((samples, n))
    for j in range(n):
        #Remove all values within cells covered by prior samples for this parameter.
        #The prior samples must also be a latin hypercube!
        if npriors > 0:
            H[:,j] = _center
            new_center, not_taken = remove_single_parameter(_center, prior_points[:,j])
            H[not_taken, j] = np.random.permutation(new_center)
        else:
            H[:, j] = np.random.permutation(_center)
    assert np.shape(H) == (samples, n)
    return H
