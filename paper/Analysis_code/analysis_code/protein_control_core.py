"""Label-profile randomization restricted to retained training proteins."""
import numpy as np


def permute_profiles(Y, rng):
    y = np.asarray(Y, dtype=float)
    if y.ndim != 2 or not np.all(np.isfinite(y)) or not np.all((y == 0) | (y == 1)):
        raise ValueError('Expected binary finite protein-label incidence matrix')
    training = np.flatnonzero(y.sum(1)); donor = rng.permutation(training)
    altered = np.zeros_like(y); altered[training] = y[donor]
    if not np.array_equal(y.sum(0), altered.sum(0)): raise ValueError('Candidate counts changed')
    if not np.array_equal(np.sort(y.sum(1)), np.sort(altered.sum(1))): raise ValueError('Profile-size multiset changed')
    if sorted(row.tobytes() for row in y[training]) != sorted(row.tobytes() for row in altered[training]):
        raise ValueError('Profile multiset changed')
    return altered, {'recipient_indices': training.tolist(), 'donor_indices': donor.tolist(),
                     'fixed_points': int(np.sum(training == donor)), 'training_proteins': len(training)}
