"""Small missing-aware joint probabilistic model; no biological data I/O.

For each *connected co-observation block*, fit a finite latent-class mixture:
  p(x_observed) = sum_k pi_k product_{j observed} p_j(x_j | k).
Continuous factors are diagonal Gaussians and categorical factors are tables.
EM uses only observed entries. A symmetric Dirichlet(alpha>1) prior regularizes
mixture weights/tables, and variances have a floor in training-standardized
units. The optimized objective is observed-data log likelihood + log prior;
both are retained, rather than claiming a smoothed step maximizes plain ML.

Disconnected blocks are fitted independently. Their component labels cannot
be aligned from these data, so evidence in one never conditions another.
Connected overlap is necessary, not a guarantee of identifiable conditionals.
Missingness itself is not modeled; informative missingness can bias results.
Intervals describe the fitted distribution, not parameter uncertainty or
validated coverage. Associations/predictive distributions are not causes.
"""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
import math
import numbers

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp, ndtr, ndtri


MAX_COMPONENTS = 16
MAX_ABS_VALUE = 1e100


def _missing(value):
    return value is None or isinstance(value, (float, np.floating)) and bool(np.isnan(value))


def _continuous(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError("Continuous observations must be finite real numbers, not booleans.")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError("Continuous value is outside the finite numerical bound 1e100.") from error
    if not math.isfinite(result) or abs(result) > MAX_ABS_VALUE:
        raise ValueError("Continuous value is outside the finite numerical bound 1e100.")
    return result


def _category(value):
    # Type tags keep True, 1, and '1' distinct; output labels remain JSON scalars.
    if isinstance(value, (bool, np.bool_)):
        return ("bool", bool(value))
    if isinstance(value, str):
        return ("str", value)
    if isinstance(value, numbers.Integral):
        result = int(value)
        if abs(result) > MAX_ABS_VALUE:
            raise ValueError("Categorical numeric label exceeds the numerical bound.")
        return ("int", result)
    if isinstance(value, numbers.Real):
        return ("float", _continuous(value))
    raise ValueError("Categorical labels must be finite JSON scalar strings/numbers/booleans.")


def _int_option(name, value, minimum, maximum):
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return int(value)


def _real_option(name, value, minimum, maximum, *, open_min=False):
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a finite number.")
    result = float(value)
    if not math.isfinite(result) or result > maximum or (result <= minimum if open_min else result < minimum):
        raise ValueError(f"{name} is outside its allowed bounds.")
    return result


@dataclass
class _Parameters:
    weights: np.ndarray
    means: dict
    variances: dict
    tables: dict


@dataclass
class _Block:
    name: str
    fields: tuple
    values: dict
    params: _Parameters
    n_rows: int
    n_components: int
    log_likelihood: float
    objective: float
    trace: list
    converged: bool
    iterations: int


class JointModel:
    """Bounded mixed-data mixture fitted only on supplied training mappings.

    ``max_components`` requests fixed K for each supported block, capped by
    that block's informative row count; it does NOT perform BIC/model search.
    Choose K outside the test set. Rows may omit fields or use None/float NaN.
    Unknown fields are rejected to make accidental hidden-label use explicit.
    """

    def __init__(self, field_types, *, max_components=3, seed=0, n_init=5,
                 max_iter=200, tol=1e-6, variance_floor=1e-3,
                 dirichlet_alpha=1.1):
        if not isinstance(field_types, Mapping) or not field_types:
            raise ValueError("field_types must be a nonempty mapping.")
        if any(not isinstance(f, str) or not f for f in field_types):
            raise ValueError("Field names must be nonempty strings.")
        if any(t not in ("continuous", "categorical") for t in field_types.values()):
            raise ValueError("Field types must be 'continuous' or 'categorical'.")
        self.field_types = dict(sorted(field_types.items()))
        self.max_components = _int_option("max_components", max_components, 1, MAX_COMPONENTS)
        self.seed = _int_option("seed", seed, 0, 2**32 - 1)
        self.n_init = _int_option("n_init", n_init, 1, 50)
        self.max_iter = _int_option("max_iter", max_iter, 1, 10000)
        self.tol = _real_option("tol", tol, 0., 1.)
        self.variance_floor = _real_option("variance_floor", variance_floor, 1e-12, 1e6)
        self.dirichlet_alpha = _real_option("dirichlet_alpha", dirichlet_alpha, 1., 1e6, open_min=True)
        self._fitted = False

    def _validate_rows(self, rows):
        rows = list(rows)
        if not rows:
            raise ValueError("At least one training row is required.")
        result = []
        for row in rows:
            if not isinstance(row, Mapping):
                raise ValueError("Each row must be a mapping.")
            unknown = set(row) - set(self.field_types)
            if unknown:
                raise ValueError(f"Undeclared fields: {sorted(unknown)}")
            converted = {}
            for f, t in self.field_types.items():
                value = row.get(f)
                if not _missing(value):
                    converted[f] = _continuous(value) if t == "continuous" else _category(value)
            result.append(converted)
        return result

    def fit(self, training_rows):
        """Fit observed-data EM, returning self. Never accepts query labels."""
        self._fitted = False
        rows = self._validate_rows(training_rows)
        informative = [r for r in rows if r]
        if not informative:
            raise ValueError("All training rows are entirely missing; nothing can be fitted.")
        fields = tuple(self.field_types)
        counts = {f: sum(f in r for r in informative) for f in fields}
        supported = [f for f in fields if counts[f]]
        coobs = {f: {g: sum(f in r and g in r for r in informative) for g in fields} for f in fields}
        remaining, components = set(supported), []
        while remaining:
            stack, component = [min(remaining)], []
            remaining.remove(stack[0])
            while stack:
                f = stack.pop()
                component.append(f)
                linked = sorted(g for g in remaining if coobs[f][g])
                remaining.difference_update(linked)
                stack.extend(linked)
            components.append(tuple(sorted(component)))
        components.sort()
        self.field_metadata_ = {}
        self._category_indices = {}
        for f, t in self.field_types.items():
            meta = {"type": t, "observed_count": counts[f], "missing_count": len(rows) - counts[f],
                    "supported": bool(counts[f]), "block": None,
                    "coobserved_counts": coobs[f].copy()}
            if not counts[f]:
                meta["status"] = "unsupported_all_missing"
            elif t == "continuous":
                vals = np.asarray([r[f] for r in informative if f in r])
                # Center before summation: even many identical large values
                # can otherwise acquire a fictitious nonzero sample variance.
                centered = vals - vals[0]
                mean = float(vals[0] + np.mean(centered))
                std = float(np.std(centered))
                # Constant/near-underflow columns retain finite unit scaling.
                scale = std if std >= 1e-100 else 1.
                meta.update({"training_mean": mean, "training_scale": scale,
                             "training_std": std, "constant": bool(np.all(vals == vals[0])),
                             "standardized_variance_floor": self.variance_floor})
            else:
                categories = sorted({r[f] for r in informative if f in r}, key=lambda x: (x[0], repr(x[1])))
                self._category_indices[f] = {cat: i for i, cat in enumerate(categories)}
                meta["categories"] = [cat[1] for cat in categories]
                meta["category_types"] = [cat[0] for cat in categories]
            self.field_metadata_[f] = meta
        self.blocks_ = []
        for index, fs in enumerate(components):
            name = f"block_{index}"
            for f in fs:
                self.field_metadata_[f]["block"] = name
            block_rows = [r for r in informative if any(f in r for f in fs)]
            values = {}
            for f in fs:
                if self.field_types[f] == "continuous":
                    meta = self.field_metadata_[f]
                    values[f] = np.asarray([(r[f] - meta["training_mean"]) / meta["training_scale"]
                                            if f in r else np.nan for r in block_rows])
                else:
                    values[f] = np.asarray([self._category_indices[f][r[f]] if f in r else -1
                                            for r in block_rows], dtype=int)
            k = min(self.max_components, len(block_rows))
            block = self._fit_block(name, fs, values, k, index)
            self.blocks_.append(block)
        self._blocks_by_name = {b.name: b for b in self.blocks_}
        self.training_metadata_ = {"rows_supplied": len(rows), "informative_rows": len(informative),
                                   "ignored_all_missing_rows": len(rows) - len(informative),
                                   "requested_components": self.max_components,
                                   "component_selection": "fixed_requested_K_capped_by_block_informative_rows",
                                   "seed": self.seed, "n_init": self.n_init,
                                   "variance_floor": self.variance_floor,
                                   "dirichlet_alpha": self.dirichlet_alpha,
                                   "missingness_model": "not_modeled",
                                   "connected_overlap_guarantees_identifiability": False}
        self._fitted = True
        return self

    def _initial(self, fs, values, k, rng):
        n = len(next(iter(values.values())))
        seeds = rng.choice(n, size=k, replace=False)
        means, variances, tables = {}, {}, {}
        for f in fs:
            vals = values[f]
            if self.field_types[f] == "continuous":
                observed = vals[np.isfinite(vals)]
                # Missing seed entries initialize a parameter from another
                # observed value; no row entries are imputed for EM.
                means[f] = np.array([vals[s] if np.isfinite(vals[s]) else rng.choice(observed) for s in seeds])
                variances[f] = np.full(k, max(float(np.var(observed)), self.variance_floor))
            else:
                c = len(self.field_metadata_[f]["categories"])
                marginal = np.bincount(vals[vals >= 0], minlength=c).astype(float) + self.dirichlet_alpha - 1
                marginal /= marginal.sum()
                table = np.tile(marginal, (k, 1))
                for j, s in enumerate(seeds):
                    if vals[s] >= 0:
                        table[j, vals[s]] += 1.
                tables[f] = table / table.sum(axis=1, keepdims=True)
        return _Parameters(np.full(k, 1. / k), means, variances, tables)

    def _e_step(self, fs, values, params):
        n = len(next(iter(values.values())))
        joint = np.broadcast_to(np.log(params.weights), (n, len(params.weights))).copy()
        for f in fs:
            vals = values[f]
            if self.field_types[f] == "continuous":
                mask = np.isfinite(vals)
                delta = vals[mask, None] - params.means[f]
                joint[mask] += -.5 * (math.log(2 * math.pi) + np.log(params.variances[f]) + delta**2 / params.variances[f])
            else:
                mask = vals >= 0
                joint[mask] += np.log(params.tables[f][:, vals[mask]].T)
        maxima = joint.max(axis=1, keepdims=True)
        shifted = joint - maxima
        shifted_normalizers = logsumexp(shifted, axis=1)
        normalizers = maxima[:, 0] + shifted_normalizers
        responsibilities = np.exp(shifted - shifted_normalizers[:, None])
        likelihood = float(normalizers.sum())
        if not np.isfinite(likelihood) or not np.isfinite(responsibilities).all():
            raise FloatingPointError("Non-finite observed-data likelihood.")
        return responsibilities, likelihood

    def _m_step(self, fs, values, responsibilities, previous):
        smoothing = self.dirichlet_alpha - 1.
        weights = responsibilities.sum(axis=0) + smoothing
        weights /= weights.sum()
        means, variances, tables = {}, {}, {}
        for f in fs:
            vals = values[f]
            if self.field_types[f] == "continuous":
                mask = np.isfinite(vals)
                r, x = responsibilities[mask], vals[mask]
                mass = r.sum(axis=0)
                alive = mass > np.finfo(float).tiny
                means[f] = previous.means[f].copy()
                variances[f] = previous.variances[f].copy()
                means[f][alive] = (r[:, alive] * x[:, None]).sum(axis=0) / mass[alive]
                differences = x[:, None] - means[f][alive]
                variances[f][alive] = np.maximum((r[:, alive] * differences**2).sum(axis=0) / mass[alive], self.variance_floor)
            else:
                mask = vals >= 0
                c = len(self.field_metadata_[f]["categories"])
                table = np.full((len(weights), c), smoothing)
                for category in range(c):
                    table[:, category] += responsibilities[mask & (vals == category)].sum(axis=0)
                tables[f] = table / table.sum(axis=1, keepdims=True)
        return _Parameters(weights, means, variances, tables)

    def _objective(self, likelihood, params):
        prior = np.log(params.weights).sum()
        prior += sum(np.log(table).sum() for table in params.tables.values())
        return float(likelihood + (self.dirichlet_alpha - 1.) * prior)

    def _fit_block(self, name, fs, values, k, index):
        best = None
        for init in range(self.n_init):
            rng = np.random.default_rng(np.random.SeedSequence([self.seed, index, init]))
            params = self._initial(fs, values, k, rng)
            responsibilities, likelihood = self._e_step(fs, values, params)
            objective = self._objective(likelihood, params)
            trace = [{"log_likelihood": likelihood, "objective": objective}]
            converged = False
            for iteration in range(1, self.max_iter + 1):
                updated = self._m_step(fs, values, responsibilities, params)
                updated_r, updated_likelihood = self._e_step(fs, values, updated)
                updated_objective = self._objective(updated_likelihood, updated)
                improvement = updated_objective - objective
                if improvement < -1e-8 * (1 + abs(objective)):
                    raise FloatingPointError("Regularized EM objective decreased beyond numerical tolerance.")
                if improvement < 0:
                    converged = True
                    break  # retain the previous numerically better parameters
                params, responsibilities = updated, updated_r
                likelihood, objective = updated_likelihood, updated_objective
                trace.append({"log_likelihood": likelihood, "objective": objective})
                if improvement <= self.tol * (1 + abs(objective)):
                    converged = True
                    break
            result = _Block(name, fs, values, params, len(next(iter(values.values()))), k,
                            likelihood, objective, trace, converged, iteration)
            if best is None or result.objective > best.objective:
                best = result
        return best

    def _require_fitted(self):
        if not self._fitted:
            raise RuntimeError("Fit the model on training rows before inference.")

    def _prepare_observations(self, observations):
        if not isinstance(observations, Mapping):
            raise ValueError("Observations must be a mapping.")
        unknown = set(observations) - set(self.field_types)
        if unknown:
            raise ValueError(f"Undeclared observed fields: {sorted(unknown)}")
        observed, usable, warnings = {}, {}, []
        missing_mask = {}
        for f, t in self.field_types.items():
            raw = observations.get(f)
            missing_mask[f] = _missing(raw)
            if missing_mask[f]:
                continue
            value = _continuous(raw) if t == "continuous" else _category(raw)
            output_value = value if t == "continuous" else value[1]
            status = "observed"
            if not self.field_metadata_[f]["supported"]:
                status = "unsupported_training_field"
            elif t == "categorical" and value not in self._category_indices[f]:
                status = "unknown_category"
            else:
                usable[f] = value
            observed[f] = {"value": output_value, "status": status}
            if status != "observed":
                warnings.append(f"{f}: {status}; this observed fact is retained but not used for conditioning.")
        return observed, usable, missing_mask, warnings

    def _posterior(self, block, usable):
        joint = np.log(block.params.weights).copy()
        fields = []
        for f in block.fields:
            if f not in usable:
                continue
            fields.append(f)
            if self.field_types[f] == "continuous":
                meta = self.field_metadata_[f]
                x = (usable[f] - meta["training_mean"]) / meta["training_scale"]
                if not np.isfinite(x) or abs(x) > MAX_ABS_VALUE:
                    raise ValueError("Standardized observation exceeds safe numerical bounds.")
                joint += -.5 * (math.log(2 * math.pi) + np.log(block.params.variances[f]) +
                                (x - block.params.means[f])**2 / block.params.variances[f])
            else:
                joint += np.log(block.params.tables[f][:, self._category_indices[f][usable[f]]])
        shifted = joint - joint.max()
        weights = np.exp(shifted - logsumexp(shifted))
        if not np.isfinite(weights).all():
            raise FloatingPointError("Non-finite conditional component probabilities.")
        return weights, fields

    @staticmethod
    def _quantile(probability, weights, means, variances):
        if not 0 < probability < 1:
            raise ValueError("Mixture quantile probability must be representable strictly inside (0, 1).")
        scales = np.sqrt(variances)
        individual = means + float(ndtri(probability)) * scales
        lower, upper = float(individual.min()), float(individual.max())
        if lower == upper:
            return lower
        # Each component CDF is <=p at the smallest component p-quantile,
        # and >=p at the largest. Solve the exact mixture CDF in that bracket.
        def cdf_minus_p(x):
            return float(np.dot(weights, ndtr((x - means) / scales)) - probability)
        fl, fu = cdf_minus_p(lower), cdf_minus_p(upper)
        if fl >= 0:
            return lower
        if fu <= 0:
            return upper
        return float(brentq(cdf_minus_p, lower, upper, xtol=1e-10, rtol=1e-12))

    def infer(self, observations, *, interval_level=.90):
        """Return only missing-field predictions; supplied facts remain separate."""
        self._require_fitted()
        level = _real_option("interval_level", interval_level, 0., 1., open_min=True)
        if level >= 1:
            raise ValueError("interval_level must be strictly between zero and one.")
        lower_probability, upper_probability = (1 - level) / 2, (1 + level) / 2
        if not 0 < lower_probability < upper_probability < 1:
            raise ValueError("interval_level has unrepresentable floating-point tail probabilities.")
        observed, usable, mask, warnings = self._prepare_observations(observations)
        posteriors, block_records = {}, {}
        for b in self.blocks_:
            weights, conditioned = self._posterior(b, usable)
            posteriors[b.name] = weights
            block_records[b.name] = {"fields": list(b.fields), "n_components": b.n_components,
                                     "posterior_weights": weights.tolist(),
                                     "prior_weights": b.params.weights.tolist(),
                                     "conditioning_fields": conditioned,
                                     "informative_training_rows": b.n_rows,
                                     "connected_overlap_guarantees_identifiability": False}
        predictions = {}
        for f, t in self.field_types.items():
            if not mask[f]:
                continue
            meta = self.field_metadata_[f]
            if not meta["supported"]:
                predictions[f] = {"type": t, "status": "unsupported_all_missing", "distribution": None,
                                  "flags": ["No training observations; no estimate is invented."]}
                continue
            block = self._blocks_by_name[meta["block"]]
            weights = posteriors[block.name]
            conditioned = block_records[block.name]["conditioning_fields"]
            outside = sorted(set(usable) - set(block.fields))
            result = {"type": t, "status": "conditional" if conditioned else "marginal",
                      "block": block.name, "conditioning_fields": conditioned,
                      "ignored_cross_block_fields": outside, "flags": []}
            if outside:
                result["flags"].append("Cross-block dependence is not learned: using this block alone.")
            if t == "continuous":
                # Quantiles in standardized units avoid overflow/cancellation;
                # then transform exactly once into original training units.
                means, variances = block.params.means[f], block.params.variances[f]
                mean = float(np.dot(weights, means))
                variance = float(np.dot(weights, variances + (means - mean)**2))
                scale, offset = meta["training_scale"], meta["training_mean"]
                qlo = self._quantile(lower_probability, weights, means, variances)
                qhi = self._quantile(upper_probability, weights, means, variances)
                result.update({"mean": offset + scale * mean, "variance": scale**2 * variance,
                               "interval": {"level": level, "lower": offset + scale * qlo,
                                            "upper": offset + scale * qhi, "calibrated": False},
                               "components": [{"weight": float(w), "mean": offset + scale * float(mu),
                                               "variance": scale**2 * float(var)}
                                              for w, mu, var in zip(weights, means, variances)],
                               "uncertainty_scope": "fitted_distribution_only_not_parameter_uncertainty"})
                if result["interval"]["lower"] == result["interval"]["upper"]:
                    result["flags"].append("Interval endpoints coincide at floating-point resolution; inspect the positive model variance.")
            else:
                probabilities = np.dot(weights, block.params.tables[f])
                result["probabilities"] = [{"value": label, "probability": float(p)}
                                            for label, p in zip(meta["categories"], probabilities)]
            predictions[f] = result
        if len(self.blocks_) > 1:
            warnings.append("Disconnected co-observation blocks are independent in this model; cross-block conditionals are unidentifiable here.")
        warnings.extend(["Connected co-observation does not guarantee identifiability or calibrated predictions.",
                         "Missingness is not modeled; predictions assume the observed-data model is appropriate."])
        return {"observed": observed, "missing_mask": mask, "predictions": predictions,
                "blocks": block_records, "field_metadata": deepcopy(self.field_metadata_),
                "training_metadata": deepcopy(self.training_metadata_), "warnings": warnings}

    def conditional_log_prob(self, field, value, observations):
        """Held-out scalar log density/mass; target must not appear in evidence.

        Unknown categorical targets have zero modeled support (log p = -inf).
        Unsupported all-missing training targets return None, not a fake score.
        The helper does not refit or use the target to compute responsibilities.
        """
        self._require_fitted()
        if field not in self.field_types:
            raise ValueError("Target field is undeclared.")
        if not isinstance(observations, Mapping):
            raise ValueError("Observations must be a mapping.")
        if field in observations and not _missing(observations[field]):
            raise ValueError("Target value must be hidden from conditioning observations.")
        if _missing(value):
            raise ValueError("A held-out target value is required.")
        _, usable, _, _ = self._prepare_observations(observations)
        t = self.field_types[field]
        target = _continuous(value) if t == "continuous" else _category(value)
        meta = self.field_metadata_[field]
        if not meta["supported"]:
            return None
        block = self._blocks_by_name[meta["block"]]
        weights, _ = self._posterior(block, usable)
        if t == "categorical":
            if target not in self._category_indices[field]:
                return -math.inf
            probability = float(np.dot(weights, block.params.tables[field][:, self._category_indices[field][target]]))
            return math.log(probability)
        x = (target - meta["training_mean"]) / meta["training_scale"]
        if not np.isfinite(x) or abs(x) > MAX_ABS_VALUE:
            raise ValueError("Standardized target exceeds safe numerical bounds.")
        # log weights from a normalized log posterior avoid log(0) warnings;
        # zero underflowed component mass legitimately contributes -inf.
        with np.errstate(divide="ignore"):
            log_weights = np.log(weights)
        log_density = -.5 * (math.log(2 * math.pi) + np.log(block.params.variances[field]) +
                            (x - block.params.means[field])**2 / block.params.variances[field])
        return float(logsumexp(log_weights + log_density) - math.log(meta["training_scale"]))

    def diagnostics(self):
        """Training-only fit/overlap metadata; no predictions or hidden truth."""
        self._require_fitted()
        records = []
        for b in self.blocks_:
            jacobian = sum(int(np.isfinite(b.values[f]).sum()) * math.log(self.field_metadata_[f]["training_scale"])
                           for f in b.fields if self.field_types[f] == "continuous")
            records.append({"name": b.name, "fields": list(b.fields), "n_components": b.n_components,
                            "informative_rows": b.n_rows, "converged": b.converged,
                            "iterations": b.iterations, "likelihood_space": "original_observed_units",
                            "log_likelihood": b.log_likelihood - jacobian,
                            "regularized_objective": b.objective - jacobian,
                            "standardized_log_likelihood": b.log_likelihood,
                            "standardized_regularized_objective": b.objective,
                            "continuous_log_jacobian": jacobian,
                            "trace": [{"log_likelihood": r["log_likelihood"] - jacobian,
                                       "objective": r["objective"] - jacobian} for r in b.trace]})
        return {"training": deepcopy(self.training_metadata_), "fields": deepcopy(self.field_metadata_),
                "blocks": records}

    def sample_missing(self, observations, *, n_samples=1, seed=0):
        """Draw coherent missing fields with one shared latent per block/draw.

        Supplied facts (including flagged unknown categories) are preserved.
        Unsupported missing fields remain None. Disconnected blocks receive
        independent latent draws; this does not assert biological independence.
        """
        self._require_fitted()
        n = _int_option("n_samples", n_samples, 1, 100000)
        seed = _int_option("seed", seed, 0, 2**32 - 1)
        observed, usable, mask, warnings = self._prepare_observations(observations)
        rng = np.random.default_rng(seed)
        samples = [{f: observed[f]["value"] if f in observed else None for f in self.field_types}
                   for _ in range(n)]
        latents = {}
        unsupported = [f for f in self.field_types if mask[f] and not self.field_metadata_[f]["supported"]]
        for block in self.blocks_:
            weights, _ = self._posterior(block, usable)
            choices = rng.choice(block.n_components, size=n, p=weights)
            latents[block.name] = choices.tolist()
            for f in block.fields:
                if not mask[f]:
                    continue
                meta = self.field_metadata_[f]
                if self.field_types[f] == "continuous":
                    values = rng.normal(block.params.means[f][choices], np.sqrt(block.params.variances[f][choices]))
                    values = meta["training_mean"] + meta["training_scale"] * values
                    if not np.isfinite(values).all():
                        raise FloatingPointError("Generated values exceed finite numerical range.")
                    for row, value in zip(samples, values):
                        row[f] = float(value)
                else:
                    labels = meta["categories"]
                    for row, component in zip(samples, choices):
                        row[f] = labels[int(rng.choice(len(labels), p=block.params.tables[f][component]))]
        return {"samples": samples, "latent_components": latents, "observed": observed,
                "unsupported_fields": unsupported, "warnings": warnings,
                "uncertainty_scope": "fitted_distribution_only_not_parameter_uncertainty"}
