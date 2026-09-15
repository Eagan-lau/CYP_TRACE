"""Feature construction and data-driven two-residual conditional coefficients."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from conditional_core import normalized_features, group_weights, select_lambda


def local_descriptions(onehot, present, global_features):
    onehot = np.asarray(onehot, dtype=float); present = np.asarray(present, dtype=bool)
    if onehot.ndim != 3 or onehot.shape[-1] != 21 or onehot.shape[:2] != present.shape:
        raise ValueError('Invalid ordered feature shape')
    if not np.allclose(onehot.sum(-1), present.astype(float)) or not np.all((onehot == 0) | (onehot == 1)):
        raise ValueError('Invalid residue encoding or missingness')
    absent = (~present).astype(float)
    ordered = np.concatenate([onehot, absent[..., None]], axis=2).reshape(len(onehot), -1)
    missing = np.stack([onehot[..., 20], absent], axis=2).reshape(len(onehot), -1)
    return {'global': normalized_features(global_features), 'ordered': normalized_features(ordered),
            'missing': normalized_features(missing), 'composition': normalized_features(onehot.sum(1))}


def joint_objective(logprior, deltas, targets, records):
    lp = np.asarray(logprior, dtype=float); d = np.asarray(deltas, dtype=float)
    if d.shape != (len(lp), lp.shape[1], 2): raise ValueError('Expected two conditional residuals')
    weights = group_weights(records)
    true_prior = np.array([lp[i, t].mean() for i, t in enumerate(targets)])
    true_delta = np.array([d[i, t].mean(0) for i, t in enumerate(targets)])
    def value_gradient(coefficients):
        contribution = np.einsum('nck,k->nc', d, coefficients)
        score = lp + contribution; z = logsumexp(score, axis=1)
        probability = np.exp(score - z[:, None])
        value = weights @ (z - true_prior - true_delta @ coefficients)
        gradient = weights @ (np.einsum('nc,nck->nk', probability, d) - true_delta)
        return float(value), gradient
    return value_gradient


def select_joint(logprior, global_delta, local_delta, targets, records):
    if len({r['group'] for r in records}) < 2:
        return {'lambdas': [0., 0.], 'reason': 'fewer_than_two_inner_groups', 'joint_optimizer_accepted': False}
    d = np.stack([global_delta, local_delta], axis=2); objective = joint_objective(logprior, d, targets, records)
    g = select_lambda(logprior, global_delta, targets, records)['lambda']
    l = select_lambda(logprior, local_delta, targets, records)['lambda']
    starts = [np.array([0., 0.]), np.array([g, 0.]), np.array([0., l])]
    start = min(starts, key=lambda x: (objective(x)[0], float(x.sum()), float(x[0])))
    best_loss = objective(start)[0]; base = objective(np.zeros(2))[0]
    opt = minimize(objective, start, jac=True, method='L-BFGS-B', bounds=[(0., 1024.), (0., 1024.)],
                   options={'gtol': 1e-9, 'ftol': 1e-13, 'maxiter': 500, 'maxls': 100})
    value, gradient = objective(opt.x)
    projected = np.where(opt.x <= 1e-8, np.minimum(gradient, 0), gradient)
    accepted = bool(opt.success and np.all(np.isfinite(opt.x)) and np.max(opt.x) < 1023.999 and
                    np.max(np.abs(projected)) <= 1e-6 and value <= best_loss + 1e-10)
    chosen = opt.x if accepted else start
    actual, final_gradient = objective(chosen)
    design = d.reshape(-1, 2); gram = design.T @ design
    eigen = np.linalg.eigvalsh(gram); ratio = float(max(eigen[0], 0) / max(eigen[-1], 1e-30))
    return {'lambdas': list(map(float, chosen)), 'reason': 'convex_joint_optimum' if accepted else 'flagged_best_single_or_zero_fallback',
            'joint_optimizer_accepted': accepted, 'optimizer_message': str(opt.message), 'optimizer_success': bool(opt.success),
            'optimizer_coefficients': list(map(float, opt.x)), 'optimizer_projected_gradient_max': float(np.max(np.abs(projected))),
            'selected_loss': actual, 'zero_loss': base, 'best_single_or_zero_loss': best_loss,
            'selected_gradient': list(map(float, final_gradient)), 'residual_design_eigenvalue_ratio': ratio,
            'separate_coefficient_identifiability_warning': ratio < 1e-8, 'inner_groups': len({r['group'] for r in records}),
            'inner_rows': len(records)}
