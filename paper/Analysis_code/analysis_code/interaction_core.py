"""Exact separable-ridge ordered-protein x candidate-centre interaction."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.sparse import csr_matrix


def normalized_rows(values):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 2 or not np.all(np.isfinite(values)):
        raise ValueError("features must be a finite matrix")
    return values / np.maximum(np.linalg.norm(values, axis=1, keepdims=True), 1e-12)


def fit_interaction(protein_features, labels, chemical_kernel, reaction_transform):
    """Fit one exact bilinear ridge model with train-only GCV.

    The reaction SVD is a fixed catalogue transform.  No outer labels select a
    rank or penalty, and a zero-interaction solution is an explicit candidate.
    """
    p = np.asarray(protein_features, dtype=np.float64)
    y = np.asarray(labels, dtype=np.float64)
    kernel = np.asarray(chemical_kernel, dtype=np.float64)
    training = np.flatnonzero(y.sum(axis=1))
    if not len(training):
        raise ValueError("empty interaction training pool")
    balanced = y[training] / y[training].sum(axis=1, keepdims=True)
    target = csr_matrix(balanced) @ kernel
    target = np.asarray(target, dtype=np.float64)
    prior = target.mean(axis=0)
    probability = np.maximum(prior / prior.sum(), 1e-12)
    probability /= probability.sum()
    residual = target - prior
    residual_rms = float(np.sqrt(np.mean(residual * residual)))

    x = p[training]
    x_center = x.mean(axis=0)
    xc = x - x_center
    x_scale = float(np.trace(xc @ xc.T) / len(training))
    z_u = np.asarray(reaction_transform["u"], dtype=np.float64)
    z_s = np.asarray(reaction_transform["singular_values"], dtype=np.float64)
    kernel_u = np.asarray(reaction_transform["chemical_kernel_times_u"], dtype=np.float64)
    base = {
        "training": training,
        "prior_log": np.log(probability),
        "residual_rms": residual_rms,
        "x_center": x_center,
        "x_scale": x_scale,
        "x_singular_values": np.zeros(0),
        "x_vt": np.zeros((0, p.shape[1])),
        "reaction_rank": int(len(z_s)),
        "core_coefficients": np.zeros((0, len(z_s))),
        "alpha": None,
        "zero_interaction": True,
        "zero_reason": None,
        "x_rank_tolerance": None,
        "gcv_zero": None,
        "gcv_selected": None,
        "gcv_curve_alpha": np.zeros(0),
        "gcv_curve_value": np.zeros(0),
        "alpha_bound_adjacent": False,
    }
    if len(training) < 3 or residual_rms < 1e-12 or x_scale < 1e-12 or not len(z_s):
        base["zero_reason"] = "insufficient_training_or_feature_variation"
        return base

    xs = xc / np.sqrt(x_scale)
    ux, sx, vtx = np.linalg.svd(xs, full_matrices=False)
    tolerance = float(np.finfo(np.float64).eps * max(xs.shape) * sx[0]) if len(sx) else 0.0
    keep = sx > tolerance
    ux, sx, vtx = ux[:, keep], sx[keep], vtx[keep]
    base.update({"x_singular_values": sx, "x_vt": vtx, "x_rank_tolerance": tolerance})
    if not len(sx):
        base["zero_reason"] = "zero_numerical_protein_rank"
        return base

    # R @ U_z = balanced(Y) @ (K @ U_z) - its across-protein mean.
    target_u = np.asarray(csr_matrix(balanced) @ kernel_u, dtype=np.float64)
    residual_u = target_u - target_u.mean(axis=0)
    projected = ux.T @ residual_u
    energy = projected * projected
    spectrum = (sx * sx)[:, None] * (z_s * z_s)[None, :]
    total_residual_ss = float(np.sum(residual * residual))
    observations = int(residual.size)
    intercept_fraction = 1.0 / len(training)

    def gcv_alpha(alpha):
        if alpha == 0.0:
            shrink = np.ones_like(spectrum)
        else:
            shrink = spectrum / (spectrum + alpha)
        rss = total_residual_ss - 2.0 * float(np.sum(energy * shrink)) + float(np.sum(energy * shrink * shrink))
        denominator = max(1.0 - intercept_fraction - float(np.sum(shrink)) / observations, 1e-12)
        return max(rss, 0.0) / observations / (denominator * denominator)

    gcv_zero = total_residual_ss / observations / max(1.0 - intercept_fraction, 1e-12) ** 2
    log_grid = np.linspace(-14.0, 14.0, 81)
    alpha_grid = np.concatenate(([0.0], np.exp(log_grid)))
    value_grid = np.asarray([gcv_alpha(float(a)) for a in alpha_grid])
    optimization = minimize_scalar(lambda value: gcv_alpha(float(np.exp(value))), bounds=(-14.0, 14.0), method="bounded", options={"xatol": 1e-6})
    candidates = [(float(value_grid[0]), 0.0)]
    if optimization.success and np.isfinite(optimization.fun):
        candidates.append((float(optimization.fun), float(np.exp(optimization.x))))
    selected_gcv, alpha = min(candidates, key=lambda item: (item[0], item[1]))
    base.update({
        "gcv_zero": float(gcv_zero),
        "gcv_selected": float(selected_gcv),
        "gcv_curve_alpha": alpha_grid,
        "gcv_curve_value": value_grid,
        "alpha_bound_adjacent": bool(alpha > 0 and abs(np.log(alpha)) > 13.99),
    })
    if not np.isfinite(selected_gcv) or selected_gcv >= gcv_zero * (1.0 - 1e-10):
        base["zero_reason"] = "gcv_favors_zero_interaction"
        return base
    denominator = spectrum + alpha
    core = projected * (sx[:, None] * z_s[None, :]) / denominator
    base.update({"alpha": float(alpha), "zero_interaction": False, "core_coefficients": core})
    return base


def predict_interaction(model, protein_query, reaction_transform):
    protein_query = np.asarray(protein_query, dtype=np.float64)
    candidates = len(reaction_transform["reaction_keys"])
    if model["zero_interaction"]:
        return np.zeros((len(protein_query), candidates), dtype=np.float64)
    x_latent = ((protein_query - model["x_center"]) / np.sqrt(model["x_scale"])) @ model["x_vt"].T
    z_u = np.asarray(reaction_transform["u"], dtype=np.float64)
    z_s = np.asarray(reaction_transform["singular_values"], dtype=np.float64)
    prediction = ((x_latent @ model["core_coefficients"]) * z_s[None, :]) @ z_u.T
    return prediction / model["residual_rms"]


def interaction_receipt(model):
    return {
        "alpha": model["alpha"],
        "alpha_bound_adjacent": model["alpha_bound_adjacent"],
        "zero_interaction": model["zero_interaction"],
        "zero_reason": model["zero_reason"],
        "residual_rms": model["residual_rms"],
        "x_scale": model["x_scale"],
        "x_rank": int(len(model["x_singular_values"])),
        "reaction_rank": int(model["reaction_rank"]),
        "x_rank_tolerance": model["x_rank_tolerance"],
        "gcv_zero": model["gcv_zero"],
        "gcv_selected": model["gcv_selected"],
        "gcv_curve_alpha": model["gcv_curve_alpha"].tolist(),
        "gcv_curve_value": model["gcv_curve_value"].tolist(),
    }
