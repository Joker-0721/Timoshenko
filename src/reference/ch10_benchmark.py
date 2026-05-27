import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / 'data'
FIG_DIR = ROOT / 'figures'
DATA_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

GAUSS5_PTS, GAUSS5_WTS = leggauss(5)
GAUSS16_PTS, GAUSS16_WTS = leggauss(16)

VIBRATION_SCAN_THRESHOLD = 5.0e-2
VIBRATION_ACCEPT_TOL = 1.0e-2
ROOT_CLUSTER_TOL = 5.0e-2


def element_matrices(E, nu, L, h, kappa, n_elem, rho=1.0):
    G = E / (2.0 * (1.0 + nu))
    A = h
    I = h**3 / 12.0
    ndof = 2 * (n_elem + 1)
    K = np.zeros((ndof, ndof), dtype=float)
    M = np.zeros((ndof, ndof), dtype=float)
    KG = np.zeros((ndof, ndof), dtype=float)
    le = L / n_elem

    EI = E * I
    kGA = kappa * G * A

    Kb = np.zeros((4, 4), dtype=float)
    Kb[1, 1] = 1.0
    Kb[1, 3] = -1.0
    Kb[3, 1] = -1.0
    Kb[3, 3] = 1.0
    Kb *= EI / le

    Bs = np.array([-1.0 / le, -0.5, 1.0 / le, -0.5], dtype=float)
    Ks = (kGA * le) * np.outer(Bs, Bs)
    Kloc = Kb + Ks

    Mloc = np.zeros((4, 4), dtype=float)
    m_t = rho * A * le
    m_r = rho * I * le
    Mloc[0, 0] = m_t / 3.0
    Mloc[2, 2] = m_t / 3.0
    Mloc[0, 2] = m_t / 6.0
    Mloc[2, 0] = m_t / 6.0
    Mloc[1, 1] = m_r / 3.0
    Mloc[3, 3] = m_r / 3.0
    Mloc[1, 3] = m_r / 6.0
    Mloc[3, 1] = m_r / 6.0

    KGloc = np.zeros((4, 4), dtype=float)
    KGloc[0, 0] = 1.0 / le
    KGloc[0, 2] = -1.0 / le
    KGloc[2, 0] = -1.0 / le
    KGloc[2, 2] = 1.0 / le

    for e in range(n_elem):
        i = e + 1
        dofs = np.array([2 * i - 2, 2 * i - 1, 2 * i, 2 * i + 1], dtype=int)
        K[np.ix_(dofs, dofs)] += Kloc
        M[np.ix_(dofs, dofs)] += Mloc
        KG[np.ix_(dofs, dofs)] += KGloc

    return K, M, KG, A, I


def free_dofs(n_elem, bc):
    ndof = 2 * (n_elem + 1)
    all_dofs = np.arange(ndof)
    if bc == 'cantilever':
        fixed = np.array([0, 1], dtype=int)
    elif bc == 'fixed-fixed':
        fixed = np.array([0, 1, ndof - 2, ndof - 1], dtype=int)
    elif bc == 'simply-supported':
        fixed = np.array([0, ndof - 2], dtype=int)
    elif bc == 'pinned-pinned':
        fixed = np.array([0, ndof - 2], dtype=int)
    else:
        raise ValueError(f'Unknown bc: {bc}')
    mask = np.ones(ndof, dtype=bool)
    mask[fixed] = False
    return all_dofs[mask]


def output_scale(E, I, A, L, rho):
    return (L**2) * math.sqrt(rho * A / (E * I))


def solve_vibration_eigenpairs(E, nu, L, h, kappa, n_elem, bc, n_modes=15, rho=1.0, output='beta2'):
    K, M, _, A, I = element_matrices(E, nu, L, h, kappa, n_elem, rho=rho)
    dofs = free_dofs(n_elem, bc)
    vals, vecs = eig(K[np.ix_(dofs, dofs)], M[np.ix_(dofs, dofs)])
    vals = np.real(vals)
    vecs = np.real(vecs)
    keep = np.isfinite(vals) & (vals > 1e-12)
    vals = vals[keep]
    vecs = vecs[:, keep]
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    omegas = np.sqrt(vals)
    scale = output_scale(E, I, A, L, rho)
    outputs = omegas * scale
    if output == 'beta':
        outputs = np.sqrt(outputs)
    elif output != 'beta2':
        raise ValueError(f'Unknown vibration output: {output}')
    full_modes = []
    for i in range(min(n_modes, vecs.shape[1])):
        full = np.zeros(2 * (n_elem + 1), dtype=float)
        full[dofs] = vecs[:, i]
        full_modes.append(full)
    return np.asarray(outputs[:n_modes], dtype=float), full_modes


def solve_vibration(E, nu, L, h, kappa, n_elem, bc, n_modes=15, rho=1.0, output='beta2'):
    vals, _ = solve_vibration_eigenpairs(E, nu, L, h, kappa, n_elem, bc, n_modes=n_modes, rho=rho, output=output)
    return vals


def solve_vibration_modes(E, nu, L, h, kappa, n_elem, bc, n_modes=4, rho=1.0):
    omega_bar, full_modes = solve_vibration_eigenpairs(E, nu, L, h, kappa, n_elem, bc, n_modes=n_modes, rho=rho, output='beta')
    x = np.linspace(0.0, L, n_elem + 1)
    mode_shapes = []
    for full in full_modes:
        w = full[0::2].copy()
        m = np.max(np.abs(w))
        if m > 0:
            w = w / m
        if w[len(w) // 2] < 0:
            w = -w
        mode_shapes.append((x, w))
    return omega_bar, mode_shapes


def solve_buckling_eigenpairs(E, nu, L, h, kappa, n_elem, bc, n_modes=4):
    K, _, KG, _, _ = element_matrices(E, nu, L, h, kappa, n_elem, rho=1.0)
    dofs = free_dofs(n_elem, bc)
    vals, vecs = eig(K[np.ix_(dofs, dofs)], KG[np.ix_(dofs, dofs)])
    vals = np.real(vals)
    vecs = np.real(vecs)
    keep = np.isfinite(vals) & (vals > 1e-9)
    vals = vals[keep]
    vecs = vecs[:, keep]
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    full_modes = []
    for i in range(min(n_modes, vecs.shape[1])):
        full = np.zeros(2 * (n_elem + 1), dtype=float)
        full[dofs] = vecs[:, i]
        full_modes.append(full)
    return vals[:n_modes], full_modes


def solve_buckling(E, nu, L, h, kappa, n_elem, bc):
    vals, _ = solve_buckling_eigenpairs(E, nu, L, h, kappa, n_elem, bc, n_modes=1)
    return float(vals[0])


def solve_buckling_modes(E, nu, L, h, kappa, n_elem, bc, n_modes=4):
    vals, full_modes = solve_buckling_eigenpairs(E, nu, L, h, kappa, n_elem, bc, n_modes=n_modes)
    x = np.linspace(0.0, L, n_elem + 1)
    mode_shapes = []
    for full in full_modes:
        w = full[0::2].copy()
        m = np.max(np.abs(w))
        if m > 0:
            w = w / m
        if np.mean(w[1:-1]) < 0:
            w = -w
        mode_shapes.append((x, w))
    return vals, mode_shapes


def solve_static(E, nu, L, h, kappa, n_elem, bc, q=1.0):
    nodes, U = solve_static_state(E, nu, L, h, kappa, n_elem, bc, q=q)
    w = U[0::2]
    if bc == 'simply-supported':
        return float(np.max(np.abs(w)))
    return float(abs(w[-1]))


def solve_static_state(E, nu, L, h, kappa, n_elem, bc, q=1.0):
    G = E / (2.0 * (1.0 + nu))
    A = h
    I = h**3 / 12.0
    EI = E * I
    kGA = kappa * G * A
    ndof = 2 * (n_elem + 1)
    K = np.zeros((ndof, ndof), dtype=float)
    F = np.zeros(ndof, dtype=float)
    le = L / n_elem
    nodes = np.linspace(0.0, L, n_elem + 1)

    Kb = np.zeros((4, 4), dtype=float)
    Kb[1, 1] = 1.0
    Kb[1, 3] = -1.0
    Kb[3, 1] = -1.0
    Kb[3, 3] = 1.0
    Kb *= EI / le
    Bs = np.array([-1.0 / le, -0.5, 1.0 / le, -0.5], dtype=float)
    Ks = (kGA * le) * np.outer(Bs, Bs)
    Kloc = Kb + Ks
    Floc = np.array([q * le / 2.0, 0.0, q * le / 2.0, 0.0], dtype=float)

    for e in range(n_elem):
        dofs = np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3], dtype=int)
        K[np.ix_(dofs, dofs)] += Kloc
        F[dofs] += Floc

    all_dofs = np.arange(ndof)
    if bc == 'simply-supported':
        fixed = np.array([0, ndof - 2], dtype=int)
    elif bc == 'cantilever':
        fixed = np.array([0, 1], dtype=int)
    else:
        raise ValueError(f'Unknown bc for static solve: {bc}')
    free = np.setdiff1d(all_dofs, fixed)
    U = np.zeros(ndof, dtype=float)
    U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    return nodes, U


def shape_l2(xi):
    return 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)


def mode_element_fields(xi, Ue):
    n1, n2 = shape_l2(xi)
    w = n1 * Ue[0] + n2 * Ue[2]
    phi = n1 * Ue[1] + n2 * Ue[3]
    return float(w), float(phi)


def w_fem_element(xi, Ue):
    return mode_element_fields(xi, Ue)[0]


def w_fem_at_x(nodes, U, x_query):
    if x_query <= nodes[0]:
        return float(U[0])
    if x_query >= nodes[-1]:
        return float(U[-2])
    e = np.searchsorted(nodes, x_query, side='right') - 1
    e = max(0, min(e, len(nodes) - 2))
    x1 = nodes[e]
    x2 = nodes[e + 1]
    le = x2 - x1
    xi = 2.0 * (x_query - 0.5 * (x1 + x2)) / le
    Ue = U[np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3], dtype=int)]
    return float(w_fem_element(xi, Ue))


def w_exact_ss(x, E, I, kappa, G, A, L, q):
    wb = q * x * (L**3 - 2.0 * L * x**2 + x**3) / (24.0 * E * I)
    ws = q * x * (L - x) / (2.0 * kappa * G * A)
    return wb + ws


def w_exact_cf(x, E, I, kappa, G, A, L, q):
    wb = q * x**2 * (6.0 * L**2 - 4.0 * L * x + x**2) / (24.0 * E * I)
    ws = q * x * (2.0 * L - x) / (2.0 * kappa * G * A)
    return wb + ws


def field_l2_static(nodes, U, exact_w_fn, exact_args):
    err_sq = 0.0
    ref_sq = 0.0
    for e in range(len(nodes) - 1):
        x1 = nodes[e]
        x2 = nodes[e + 1]
        le = x2 - x1
        det_j = le / 2.0
        Ue = U[np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3], dtype=int)]
        for xi, wt in zip(GAUSS5_PTS, GAUSS5_WTS):
            x = 0.5 * (x1 + x2) + det_j * xi
            w_h = w_fem_element(xi, Ue)
            w_exact = exact_w_fn(x, *exact_args)
            err_sq += (w_h - w_exact) ** 2 * wt * det_j
            ref_sq += (w_exact ** 2) * wt * det_j
    return math.sqrt(err_sq / ref_sq) * 100.0


def vibration_spatial_roots(E, G, A, I, kappa, rho, omega):
    kGA = kappa * G * A
    EI = E * I
    a = kGA * EI
    b = omega**2 * (kGA * rho * I + rho * A * EI)
    c = rho * A * rho * I * omega**4 - rho * A * kGA * omega**2
    y_roots = np.roots([a, b, c])
    roots = []
    for y in y_roots:
        s = np.sqrt(complex(y))
        roots.extend([s, -s])
    return roots


def vibration_phi_ratio(s, E, G, A, I, kappa, rho, omega):
    kGA = kappa * G * A
    if abs(s) < 1.0e-12:
        den = E * I * s * s - kGA + rho * I * omega**2
        if abs(den) < 1.0e-12:
            return 0.0
        return -(kGA * s) / den
    return (kGA * s * s + rho * A * omega**2) / (kGA * s)


def output_to_omega(output_value, output_kind, E, I, A, L, rho):
    scale = output_scale(E, I, A, L, rho)
    if output_kind == 'beta2':
        omega_bar = output_value
    elif output_kind == 'beta':
        omega_bar = output_value**2
    else:
        raise ValueError(f'Unknown output kind: {output_kind}')
    return omega_bar / scale


def omega_to_output(omega, output_kind, E, I, A, L, rho):
    omega_bar = omega * output_scale(E, I, A, L, rho)
    if output_kind == 'beta2':
        return omega_bar
    if output_kind == 'beta':
        return math.sqrt(omega_bar)
    raise ValueError(f'Unknown output kind: {output_kind}')


def vibration_boundary_matrix(output_value, E, nu, L, h, kappa, rho, bc, output_kind):
    G = E / (2.0 * (1.0 + nu))
    A = h
    I = h**3 / 12.0
    omega = output_to_omega(output_value, output_kind, E, I, A, L, rho)
    roots = vibration_spatial_roots(E, G, A, I, kappa, rho, omega)
    B = np.zeros((4, 4), dtype=complex)
    for j, s in enumerate(roots):
        r = vibration_phi_ratio(s, E, G, A, I, kappa, rho, omega)
        eL = np.exp(s * L)
        if bc == 'cantilever':
            B[0, j] = 1.0
            B[1, j] = r
            B[2, j] = s * r * eL
            B[3, j] = (s - r) * eL
        elif bc == 'fixed-fixed':
            B[0, j] = 1.0
            B[1, j] = r
            B[2, j] = eL
            B[3, j] = r * eL
        else:
            raise ValueError(f'Unsupported vibration reference bc: {bc}')
    return B, roots, omega, G, A, I


def vibration_sigma_min(output_value, E, nu, L, h, kappa, rho, bc, output_kind):
    B, _, _, _, _, _ = vibration_boundary_matrix(output_value, E, nu, L, h, kappa, rho, bc, output_kind)
    return float(np.linalg.svd(B, compute_uv=False)[-1])


def cluster_sorted_roots(candidates, tol=ROOT_CLUSTER_TOL):
    roots = []
    for x, sig in sorted(candidates, key=lambda item: item[0]):
        if not roots or abs(x - roots[-1][0]) > tol:
            roots.append((x, sig))
        elif sig < roots[-1][1]:
            roots[-1] = (x, sig)
    return roots


_VIBRATION_ROOT_CACHE = {}


def find_vibration_reference_roots(E, nu, L, h, kappa, rho, bc, output_kind, n_modes, upper_hint):
    cache_key = (E, nu, L, h, kappa, rho, bc, output_kind, n_modes, upper_hint)
    if cache_key in _VIBRATION_ROOT_CACHE:
        return _VIBRATION_ROOT_CACHE[cache_key]

    lower = 0.1 if output_kind == 'beta2' else 0.5
    upper = max(upper_hint, lower * 10.0)
    grid_count = 30000 if bc == 'cantilever' else 60000
    xs = np.linspace(lower, upper, grid_count)
    vals = np.array([
        vibration_sigma_min(x, E, nu, L, h, kappa, rho, bc, output_kind)
        for x in xs
    ])

    candidates = []
    for i in range(1, len(xs) - 1):
        if vals[i] < vals[i - 1] and vals[i] < vals[i + 1] and vals[i] < VIBRATION_SCAN_THRESHOLD:
            lo = xs[i - 1]
            hi = xs[i + 1]
            res = minimize_scalar(
                lambda x: vibration_sigma_min(x, E, nu, L, h, kappa, rho, bc, output_kind),
                bounds=(lo, hi),
                method='bounded',
                options={'xatol': 1.0e-12},
            )
            candidates.append((float(res.x), float(res.fun)))

    roots = [
        (x, sig) for x, sig in cluster_sorted_roots(candidates)
        if sig <= VIBRATION_ACCEPT_TOL
    ]
    if len(roots) < n_modes:
        raise RuntimeError(
            f'Could not identify {n_modes} vibration reference roots for {bc}, h={h}. '
            f'Found {len(roots)} roots up to {upper}.'
        )
    _VIBRATION_ROOT_CACHE[cache_key] = roots
    return roots


def make_exponential_vibration_reference(output_value, E, nu, L, h, kappa, rho, bc, output_kind):
    B, roots, omega, G, A, I = vibration_boundary_matrix(output_value, E, nu, L, h, kappa, rho, bc, output_kind)
    _, _, vh = np.linalg.svd(B)
    coeffs = vh.conj().T[:, -1]
    phase_x = 0.37 * L
    sample = 0.0j
    for s, c in zip(roots, coeffs):
        sample += c * np.exp(s * phase_x)
    if abs(sample) > 1.0e-12:
        coeffs *= np.exp(-1j * np.angle(sample))

    ratios = [vibration_phi_ratio(s, E, G, A, I, kappa, rho, omega) for s in roots]

    def evaluator(x):
        exps = [np.exp(s * x) for s in roots]
        w_val = sum(c * ex for c, ex in zip(coeffs, exps))
        phi_val = sum(r * c * ex for r, c, ex in zip(ratios, coeffs, exps))
        w_real = np.real_if_close(w_val, tol=1000.0)
        phi_real = np.real_if_close(phi_val, tol=1000.0)
        return float(np.real(w_real)), float(np.real(phi_real))

    return evaluator


def simply_supported_frequency_scalar(mode, E, nu, L, h, kappa, rho=1.0, output_kind='beta'):
    G = E / (2.0 * (1.0 + nu))
    A = h
    I = h**3 / 12.0
    lam = mode * math.pi / L
    kGA = kappa * G * A
    EI = E * I
    a = rho * A * rho * I
    b = -(rho * A * EI * lam**2 + rho * A * kGA + kGA * rho * I * lam**2)
    c = kGA * EI * lam**4
    disc = max(b * b - 4.0 * a * c, 0.0)
    omega_sq = (-b - math.sqrt(disc)) / (2.0 * a)
    omega = math.sqrt(omega_sq)
    if output_kind == 'beta2':
        return omega * output_scale(E, I, A, L, rho)
    if output_kind == 'beta':
        return math.sqrt(omega * output_scale(E, I, A, L, rho))
    raise ValueError(f'Unknown output kind: {output_kind}')


def make_simply_supported_vibration_reference(mode, E, nu, L, h, kappa, rho=1.0):
    G = E / (2.0 * (1.0 + nu))
    A = h
    I = h**3 / 12.0
    lam = mode * math.pi / L
    kGA = kappa * G * A
    EI = E * I
    a = rho * A * rho * I
    b = -(rho * A * EI * lam**2 + rho * A * kGA + kGA * rho * I * lam**2)
    c = kGA * EI * lam**4
    disc = max(b * b - 4.0 * a * c, 0.0)
    omega_sq = (-b - math.sqrt(disc)) / (2.0 * a)
    r_n = (kGA * lam**2 - rho * A * omega_sq) / (kGA * lam)

    def evaluator(x):
        return math.sin(lam * x), r_n * math.cos(lam * x)

    return evaluator


def make_vibration_reference(mode, E, nu, L, h, kappa, rho, bc, output_kind, upper_hint, n_modes):
    if bc == 'simply-supported':
        scalar = simply_supported_frequency_scalar(mode, E, nu, L, h, kappa, rho=rho, output_kind=output_kind)
        return scalar, make_simply_supported_vibration_reference(mode, E, nu, L, h, kappa, rho=rho)

    roots = find_vibration_reference_roots(
        E, nu, L, h, kappa, rho, bc, output_kind, n_modes=n_modes, upper_hint=upper_hint
    )
    scalar = roots[mode - 1][0]
    evaluator = make_exponential_vibration_reference(scalar, E, nu, L, h, kappa, rho, bc, output_kind)
    return scalar, evaluator


def make_buckling_reference(case_name, E, nu, L, h, kappa):
    G = E / (2.0 * (1.0 + nu))
    A = h
    I = h**3 / 12.0
    if case_name == 'pinned-pinned':
        lam = math.pi / L

        def w_fn(x):
            return math.sin(lam * x)
    elif case_name == 'fixed-fixed':
        lam = 2.0 * math.pi / L

        def w_fn(x):
            return 1.0 - math.cos(lam * x)
    else:
        raise ValueError(f'Unknown buckling case: {case_name}')

    Pcr = (E * I * lam**2) / (1.0 + (E * I * lam**2) / (kappa * G * A))
    rot_factor = 1.0 - Pcr / (kappa * G * A)

    def evaluator(x):
        w_val = w_fn(x)
        if case_name == 'pinned-pinned':
            dw = lam * math.cos(lam * x)
        else:
            dw = lam * math.sin(lam * x)
        return w_val, rot_factor * dw

    return Pcr, evaluator


def modal_field_l2(nodes, mode_full, reference_eval):
    norm_w_h_sq = 0.0
    norm_phi_h_sq = 0.0
    norm_w_ref_sq = 0.0
    norm_phi_ref_sq = 0.0

    for e in range(len(nodes) - 1):
        x1 = nodes[e]
        x2 = nodes[e + 1]
        le = x2 - x1
        det_j = le / 2.0
        Ue = mode_full[np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3], dtype=int)]
        for xi, wt in zip(GAUSS16_PTS, GAUSS16_WTS):
            x = 0.5 * (x1 + x2) + det_j * xi
            w_h, phi_h = mode_element_fields(xi, Ue)
            w_ref, phi_ref = reference_eval(x)
            norm_w_h_sq += (w_h**2) * wt * det_j
            norm_phi_h_sq += (phi_h**2) * wt * det_j
            norm_w_ref_sq += (w_ref**2) * wt * det_j
            norm_phi_ref_sq += (phi_ref**2) * wt * det_j

    norm_w_h = math.sqrt(max(norm_w_h_sq, 1.0e-30))
    norm_phi_h = math.sqrt(max(norm_phi_h_sq, 1.0e-30))
    norm_w_ref = math.sqrt(max(norm_w_ref_sq, 1.0e-30))
    norm_phi_ref = math.sqrt(max(norm_phi_ref_sq, 1.0e-30))

    inner_w = 0.0
    inner_phi = 0.0
    for e in range(len(nodes) - 1):
        x1 = nodes[e]
        x2 = nodes[e + 1]
        le = x2 - x1
        det_j = le / 2.0
        Ue = mode_full[np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3], dtype=int)]
        for xi, wt in zip(GAUSS16_PTS, GAUSS16_WTS):
            x = 0.5 * (x1 + x2) + det_j * xi
            w_h, phi_h = mode_element_fields(xi, Ue)
            w_ref, phi_ref = reference_eval(x)
            inner_w += (w_h / norm_w_h) * (w_ref / norm_w_ref) * wt * det_j
            inner_phi += (phi_h / norm_phi_h) * (phi_ref / norm_phi_ref) * wt * det_j

    sign = 1.0 if abs(inner_w) >= 1.0e-10 and inner_w >= 0.0 else -1.0
    if abs(inner_w) < 1.0e-10:
        sign = 1.0 if inner_phi >= 0.0 else -1.0

    err_w_sq = 0.0
    err_phi_sq = 0.0
    for e in range(len(nodes) - 1):
        x1 = nodes[e]
        x2 = nodes[e + 1]
        le = x2 - x1
        det_j = le / 2.0
        Ue = mode_full[np.array([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3], dtype=int)]
        for xi, wt in zip(GAUSS16_PTS, GAUSS16_WTS):
            x = 0.5 * (x1 + x2) + det_j * xi
            w_h, phi_h = mode_element_fields(xi, Ue)
            w_ref, phi_ref = reference_eval(x)
            w_h_hat = w_h / norm_w_h
            phi_h_hat = phi_h / norm_phi_h
            w_ref_hat = sign * w_ref / norm_w_ref
            phi_ref_hat = sign * phi_ref / norm_phi_ref
            err_w_sq += (w_h_hat - w_ref_hat) ** 2 * wt * det_j
            err_phi_sq += (phi_h_hat - phi_ref_hat) ** 2 * wt * det_j

    return math.sqrt(err_w_sq) * 100.0, math.sqrt(err_phi_sq) * 100.0


def write_csv(path, headers, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)


def mean_col(rows, idx, *, where=None):
    vals = [float(r[idx]) for r in rows if where is None or where(r)]
    return float(np.mean(vals))


def plot_table_10_1(rows_10_1):
    h_vals = [0.001, 0.01, 0.1]
    ss_rows = [r for r in rows_10_1 if r[0] == 'SS']
    cf_rows = [r for r in rows_10_1 if r[0] == 'CF']
    ss_l2_pct = [float(r[4]) for r in ss_rows]
    cf_l2_pct = [float(r[4]) for r in cf_rows]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(h_vals, ss_l2_pct, 'o-', label='SS')
    ax.plot(h_vals, cf_l2_pct, 's-', label='CF')
    ax.set_xscale('log')
    ax.set_xlabel('h/L')
    ax.set_ylabel('Field L2 (%)')
    ax.set_title('Table 10.1 Field L2 (%)')
    ax.grid(True, alpha=0.35)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'ch10_l2_10_1.png', dpi=180)
    plt.close(fig)


def plot_table_10_2_compare(rows_10_2):
    x = np.array([1, 2, 5, 10, 50], dtype=float)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for idx, mode in enumerate([1, 2]):
        subset = [r for r in rows_10_2 if int(r[1]) == mode]
        formula_vals = [r[2] for r in subset]
        fem_vals = [r[3] for r in subset]
        ax[idx].plot(x, formula_vals, 's--', label='Formula')
        ax[idx].plot(x, fem_vals, 'o-', label='FEM')
        ax[idx].set_title(f'Table 10.2 Mode {mode}')
        ax[idx].set_xlabel('n_elements')
        ax[idx].set_ylabel('omega_bar')
        ax[idx].grid(True, alpha=0.3)
        ax[idx].legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'ch10_compare_10_2.png', dpi=180)
    plt.close(fig)


def plot_table_10_2_field_l2(rows_10_2):
    x = np.array([1, 2, 5, 10, 50], dtype=float)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
    for mode in [1, 2]:
        subset = [r for r in rows_10_2 if int(r[1]) == mode]
        ax[0].plot(x, [r[4] for r in subset], 'o-', label=f'Mode {mode}')
        ax[1].plot(x, [r[5] for r in subset], 'o-', label=f'Mode {mode}')
    ax[0].set_title('Field L2_w')
    ax[1].set_title('Field L2_phi')
    for axis in ax:
        axis.set_xlabel('n_elements')
        axis.grid(True, alpha=0.3)
        axis.legend()
    ax[0].set_ylabel('Field L2 (%)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'ch10_field_l2_10_2.png', dpi=180)
    plt.close(fig)


def plot_mode_comparison(rows, outname, title, y_label):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharex=True)
    for idx, hL in enumerate([0.002, 0.01, 0.1]):
        subset = [r for r in rows if abs(float(r[1]) - hL) < 1.0e-12]
        modes = [r[0] for r in subset]
        formula_vals = [r[2] for r in subset]
        fem_vals = [r[3] for r in subset]
        axes[idx].plot(modes, formula_vals, 's--', ms=3, label='Formula')
        axes[idx].plot(modes, fem_vals, 'o-', ms=3, label='FEM')
        axes[idx].set_title(f'h/L={hL}')
        axes[idx].set_xlabel('mode')
        axes[idx].grid(True, alpha=0.3)
        if idx == 0:
            axes[idx].set_ylabel(y_label)
            axes[idx].legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(FIG_DIR / outname, dpi=180)
    plt.close(fig)


def plot_field_l2_by_h(rows, outname, title):
    labels = ['0.002', '0.01', '0.1']
    h_vals = [0.002, 0.01, 0.1]
    w_means = [mean_col(rows, 4, where=lambda r, h=h: abs(float(r[1]) - h) < 1.0e-12) for h in h_vals]
    phi_means = [mean_col(rows, 5, where=lambda r, h=h: abs(float(r[1]) - h) < 1.0e-12) for h in h_vals]
    pos = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(6.4, 4))
    ax.bar(pos - width / 2, w_means, width=width, label='Field L2_w')
    ax.bar(pos + width / 2, phi_means, width=width, label='Field L2_phi')
    ax.set_xticks(pos, labels)
    ax.set_xlabel('h/L')
    ax.set_ylabel('Field L2 (%)')
    ax.set_title(title)
    ax.grid(True, axis='y', alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / outname, dpi=180)
    plt.close(fig)


def plot_mode_shapes():
    E = 1.0
    nu = 0.3
    L = 1.0
    kappa = 5.0 / 6.0
    n_elem = 40
    for bc, prefix, title_prefix in [
        ('fixed-fixed', '10_3', 'Table 10.3 fixed-fixed'),
        ('simply-supported', '10_4', 'Table 10.4 simply-supported'),
    ]:
        for hL in [0.002, 0.01, 0.1]:
            h = hL * L
            _, modes = solve_vibration_modes(E, nu, L, h, kappa, n_elem, bc, n_modes=4, rho=1.0)
            fig, axes = plt.subplots(4, 1, figsize=(8, 5), sharex=True)
            for i, (xv, wv) in enumerate(modes):
                axes[i].set_facecolor('#f2f2f2')
                axes[i].plot(xv / L, wv * 5.0, 'o-', ms=2, lw=1.0, color='#1f4aa8')
                axes[i].grid(True, alpha=0.55, linestyle=':')
                axes[i].set_xlim(0.0, 1.0)
                axes[i].set_ylim(-5.2, 5.2)
                axes[i].set_ylabel(f'M{i + 1}', fontsize=8)
            axes[-1].set_xlabel('x/L')
            fig.suptitle(f'{title_prefix} mode shapes (h/L={hL})', fontsize=10)
            fig.tight_layout()
            fig.savefig(FIG_DIR / f'ch10_mode_shapes_{prefix}_h{str(hL).replace(".", "p")}.png', dpi=180)
            plt.close(fig)


def plot_buckling_mode_shapes():
    E = 1.0e7
    nu = 1.0 / 3.0
    L = 1.0
    kappa = 5.0 / 6.0
    n_elem = 40
    Lh = 100
    h = L / Lh
    for bc, suffix, title_prefix in [
        ('pinned-pinned', 'pp', 'Table 10.5A pinned-pinned'),
        ('fixed-fixed', 'ff', 'Table 10.5B fixed-fixed'),
    ]:
        _, modes = solve_buckling_modes(E, nu, L, h, kappa, n_elem, bc, n_modes=4)
        fig, axes = plt.subplots(4, 1, figsize=(8, 5), sharex=True)
        for i, (xv, wv) in enumerate(modes):
            axes[i].set_facecolor('#f2f2f2')
            axes[i].plot(xv / L, wv * 0.2, 'o-', ms=2, lw=1.0, color='#1f4aa8')
            axes[i].grid(True, alpha=0.55, linestyle=':')
            axes[i].set_xlim(0.0, 1.0)
            axes[i].set_ylim(-0.22, 0.22)
            axes[i].set_ylabel(f'M{i + 1}', fontsize=8)
        axes[-1].set_xlabel('x/L')
        fig.suptitle(f'{title_prefix} buckling modes (L/h={Lh})', fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f'ch10_mode_shapes_10_5_{suffix}_Lh{Lh}.png', dpi=180)
        plt.close(fig)


def plot_table_10_5_compare(rows_10_5):
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.array([10, 100, 1000], dtype=float)
    pp_rows = [r for r in rows_10_5 if r[1] == 'pinned-pinned']
    ff_rows = [r for r in rows_10_5 if r[1] == 'fixed-fixed']
    ax.plot(x, [r[2] for r in pp_rows], 'o--', label='Formula pinned-pinned')
    ax.plot(x, [r[3] for r in pp_rows], 'o-', label='FEM pinned-pinned')
    ax.plot(x, [r[2] for r in ff_rows], 's--', label='Formula fixed-fixed')
    ax.plot(x, [r[3] for r in ff_rows], 's-', label='FEM fixed-fixed')
    ax.set_xscale('log')
    ax.set_xlabel('L/h')
    ax.set_ylabel('P_cr')
    ax.set_title('Table 10.5 Buckling comparison')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'ch10_compare_10_5.png', dpi=180)
    plt.close(fig)


def plot_table_10_5_field_l2(rows_10_5):
    x = np.array([10, 100, 1000], dtype=float)
    pp_rows = [r for r in rows_10_5 if r[1] == 'pinned-pinned']
    ff_rows = [r for r in rows_10_5 if r[1] == 'fixed-fixed']
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
    ax[0].plot(x, [r[4] for r in pp_rows], 'o-', label='PP')
    ax[0].plot(x, [r[4] for r in ff_rows], 's-', label='FF')
    ax[1].plot(x, [r[5] for r in pp_rows], 'o-', label='PP')
    ax[1].plot(x, [r[5] for r in ff_rows], 's-', label='FF')
    ax[0].set_title('Field L2_w')
    ax[1].set_title('Field L2_phi')
    for axis in ax:
        axis.set_xscale('log')
        axis.set_xlabel('L/h')
        axis.grid(True, alpha=0.3)
        axis.legend()
    ax[0].set_ylabel('Field L2 (%)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'ch10_field_l2_10_5.png', dpi=180)
    plt.close(fig)


def build_static_table():
    rows_10_1 = []
    E = 1.0e8
    nu = 0.3
    L = 1.0
    kappa = 5.0 / 6.0
    G = E / (2.0 * (1.0 + nu))
    for hL in [0.001, 0.01, 0.1]:
        h = hL * L
        A = h
        I = h**3 / 12.0
        q = 1.0
        nodes_ss, U_ss = solve_static_state(E, nu, L, h, kappa, 100, 'simply-supported', q=q)
        nodes_cf, U_cf = solve_static_state(E, nu, L, h, kappa, 100, 'cantilever', q=q)
        w_ss_formula = w_exact_ss(L / 2.0, E, I, kappa, G, A, L, q)
        w_cf_formula = w_exact_cf(L, E, I, kappa, G, A, L, q)
        w_ss_fem = w_fem_at_x(nodes_ss, U_ss, L / 2.0)
        w_cf_fem = w_fem_at_x(nodes_cf, U_cf, L)
        w_ss_field_l2 = field_l2_static(nodes_ss, U_ss, w_exact_ss, (E, I, kappa, G, A, L, q))
        w_cf_field_l2 = field_l2_static(nodes_cf, U_cf, w_exact_cf, (E, I, kappa, G, A, L, q))
        rows_10_1.append(['SS', hL, w_ss_formula, w_ss_fem, w_ss_field_l2])
        rows_10_1.append(['CF', hL, w_cf_formula, w_cf_fem, w_cf_field_l2])
    write_csv(
        DATA_DIR / 'ch10_table10_1_compare.csv',
        ['case', 'h_over_L', 'formula_same_point', 'fem_same_point', 'field_l2_vs_formula_pct'],
        rows_10_1,
    )
    return rows_10_1


def build_vibration_tables():
    rows_10_2 = []
    rows_10_3 = []
    rows_10_4 = []

    E = 1.0e8
    nu = 0.3
    L = 1.0
    h = 0.001
    kappa = 5.0 / 6.0
    rho = 1.0
    formula_m1, ref_eval_m1 = make_vibration_reference(
        1, E, nu, L, h, kappa, rho, 'cantilever', 'beta2', upper_hint=35.0, n_modes=2
    )
    formula_m2, ref_eval_m2 = make_vibration_reference(
        2, E, nu, L, h, kappa, rho, 'cantilever', 'beta2', upper_hint=35.0, n_modes=2
    )
    for n in [1, 2, 5, 10, 50]:
        nodes = np.linspace(0.0, L, n + 1)
        fem_vals, fem_modes = solve_vibration_eigenpairs(
            E, nu, L, h, kappa, n, 'cantilever', n_modes=2, rho=rho, output='beta2'
        )
        for mode, formula_scalar, ref_eval in [
            (1, formula_m1, ref_eval_m1),
            (2, formula_m2, ref_eval_m2),
        ]:
            fem_scalar = float(fem_vals[mode - 1])
            l2_w, l2_phi = modal_field_l2(nodes, fem_modes[mode - 1], ref_eval)
            eig_rel = abs(fem_scalar - formula_scalar) / abs(formula_scalar) * 100.0
            rows_10_2.append([n, mode, formula_scalar, fem_scalar, l2_w, l2_phi, eig_rel])

    write_csv(
        DATA_DIR / 'ch10_table10_2_compare.csv',
        ['n_elements', 'mode', 'formula', 'fem', 'field_l2_w_pct', 'field_l2_phi_pct', 'eig_rel_error_pct'],
        rows_10_2,
    )

    E = 1.0
    nu = 0.3
    L = 1.0
    kappa = 5.0 / 6.0
    rho = 1.0
    n_elem = 40
    selected_modes = [1, 2, 5, 10, 15]
    for hL in [0.002, 0.01, 0.1]:
        h = hL * L
        nodes = np.linspace(0.0, L, n_elem + 1)
        fem_fixed_vals, fem_fixed_modes = solve_vibration_eigenpairs(
            E, nu, L, h, kappa, n_elem, 'fixed-fixed', n_modes=15, rho=rho, output='beta'
        )
        fem_simply_vals, fem_simply_modes = solve_vibration_eigenpairs(
            E, nu, L, h, kappa, n_elem, 'simply-supported', n_modes=15, rho=rho, output='beta'
        )
        upper_fixed = max(90.0, 1.6 * float(fem_fixed_vals[-1]))
        for mode in selected_modes:
            formula_fixed, ref_fixed_eval = make_vibration_reference(
                mode, E, nu, L, h, kappa, rho, 'fixed-fixed', 'beta', upper_hint=upper_fixed, n_modes=15
            )
            formula_simply, ref_simply_eval = make_vibration_reference(
                mode, E, nu, L, h, kappa, rho, 'simply-supported', 'beta', upper_hint=0.0, n_modes=15
            )
            fem_fixed = float(fem_fixed_vals[mode - 1])
            fem_simply = float(fem_simply_vals[mode - 1])
            l2_fixed_w, l2_fixed_phi = modal_field_l2(nodes, fem_fixed_modes[mode - 1], ref_fixed_eval)
            l2_simply_w, l2_simply_phi = modal_field_l2(nodes, fem_simply_modes[mode - 1], ref_simply_eval)
            rows_10_3.append([
                mode, hL, formula_fixed, fem_fixed, l2_fixed_w, l2_fixed_phi,
                abs(fem_fixed - formula_fixed) / abs(formula_fixed) * 100.0,
            ])
            rows_10_4.append([
                mode, hL, formula_simply, fem_simply, l2_simply_w, l2_simply_phi,
                abs(fem_simply - formula_simply) / abs(formula_simply) * 100.0,
            ])

    write_csv(
        DATA_DIR / 'ch10_table10_3_compare.csv',
        ['mode', 'h_over_L', 'formula', 'fem', 'field_l2_w_pct', 'field_l2_phi_pct', 'eig_rel_error_pct'],
        rows_10_3,
    )
    write_csv(
        DATA_DIR / 'ch10_table10_4_compare.csv',
        ['mode', 'h_over_L', 'formula', 'fem', 'field_l2_w_pct', 'field_l2_phi_pct', 'eig_rel_error_pct'],
        rows_10_4,
    )
    return rows_10_2, rows_10_3, rows_10_4


def build_buckling_table():
    rows_10_5 = []
    E = 1.0e7
    nu = 1.0 / 3.0
    L = 1.0
    kappa = 5.0 / 6.0
    n_elem = 40
    for Lh in [10, 100, 1000]:
        h = L / Lh
        nodes = np.linspace(0.0, L, n_elem + 1)
        for case_name in ['pinned-pinned', 'fixed-fixed']:
            fem_vals, fem_modes = solve_buckling_eigenpairs(E, nu, L, h, kappa, n_elem, case_name, n_modes=1)
            formula_scalar, ref_eval = make_buckling_reference(case_name, E, nu, L, h, kappa)
            fem_scalar = float(fem_vals[0])
            l2_w, l2_phi = modal_field_l2(nodes, fem_modes[0], ref_eval)
            rows_10_5.append([
                Lh, case_name, formula_scalar, fem_scalar, l2_w, l2_phi,
                abs(fem_scalar - formula_scalar) / abs(formula_scalar) * 100.0,
            ])

    write_csv(
        DATA_DIR / 'ch10_table10_5_compare.csv',
        ['L_over_h', 'case', 'formula', 'fem', 'field_l2_w_pct', 'field_l2_phi_pct', 'load_rel_error_pct'],
        rows_10_5,
    )
    return rows_10_5


def write_summary_files(rows_10_1, rows_10_2, rows_10_3, rows_10_4, rows_10_5):
    l2_rows = [
        ['Table 10.1 (field L2, SS mean)', mean_col(rows_10_1, 4, where=lambda r: r[0] == 'SS')],
        ['Table 10.1 (field L2, CF mean)', mean_col(rows_10_1, 4, where=lambda r: r[0] == 'CF')],
        ['Table 10.2 (Field L2_w mean)', mean_col(rows_10_2, 4)],
        ['Table 10.2 (Field L2_phi mean)', mean_col(rows_10_2, 5)],
        ['Table 10.3 (Field L2_w mean)', mean_col(rows_10_3, 4)],
        ['Table 10.3 (Field L2_phi mean)', mean_col(rows_10_3, 5)],
        ['Table 10.4 (Field L2_w mean)', mean_col(rows_10_4, 4)],
        ['Table 10.4 (Field L2_phi mean)', mean_col(rows_10_4, 5)],
        ['Table 10.5 (Field L2_w mean)', mean_col(rows_10_5, 4)],
        ['Table 10.5 (Field L2_phi mean)', mean_col(rows_10_5, 5)],
    ]
    write_csv(DATA_DIR / 'ch10_l2_summary.csv', ['case', 'metric_value'], l2_rows)

    scalar_rows = [
        ['Table 10.2 (eig_rel_error mean)', mean_col(rows_10_2, 6)],
        ['Table 10.3 (eig_rel_error mean)', mean_col(rows_10_3, 6)],
        ['Table 10.4 (eig_rel_error mean)', mean_col(rows_10_4, 6)],
        ['Table 10.5 (load_rel_error mean)', mean_col(rows_10_5, 6)],
    ]
    write_csv(DATA_DIR / 'ch10_scalar_error_summary.csv', ['case', 'metric_value'], scalar_rows)


def build_plots(rows_10_1, rows_10_2, rows_10_3, rows_10_4, rows_10_5):
    plot_table_10_1(rows_10_1)
    plot_table_10_2_compare(rows_10_2)
    plot_table_10_2_field_l2(rows_10_2)
    plot_mode_comparison(rows_10_3, 'ch10_compare_10_3.png', 'Table 10.3 fixed-fixed', 'beta')
    plot_mode_comparison(rows_10_4, 'ch10_compare_10_4.png', 'Table 10.4 simply-supported', 'beta')
    plot_field_l2_by_h(rows_10_3, 'ch10_field_l2_10_3_by_h.png', 'Table 10.3 mean Field L2 by h/L')
    plot_field_l2_by_h(rows_10_4, 'ch10_field_l2_10_4_by_h.png', 'Table 10.4 mean Field L2 by h/L')
    plot_mode_shapes()
    plot_buckling_mode_shapes()
    plot_table_10_5_compare(rows_10_5)
    plot_table_10_5_field_l2(rows_10_5)


def main():
    rows_10_1 = build_static_table()
    rows_10_2, rows_10_3, rows_10_4 = build_vibration_tables()
    rows_10_5 = build_buckling_table()
    write_summary_files(rows_10_1, rows_10_2, rows_10_3, rows_10_4, rows_10_5)
    build_plots(rows_10_1, rows_10_2, rows_10_3, rows_10_4, rows_10_5)

    print('Generated CSV files in', DATA_DIR)
    print('Generated figures in', FIG_DIR)
    print('Field L2 10.1 SS mean =', mean_col(rows_10_1, 4, where=lambda r: r[0] == 'SS'))
    print('Field L2 10.1 CF mean =', mean_col(rows_10_1, 4, where=lambda r: r[0] == 'CF'))
    print('Field L2 10.2 w mean =', mean_col(rows_10_2, 4))
    print('Field L2 10.2 phi mean =', mean_col(rows_10_2, 5))
    print('Field L2 10.3 w mean =', mean_col(rows_10_3, 4))
    print('Field L2 10.3 phi mean =', mean_col(rows_10_3, 5))
    print('Field L2 10.4 w mean =', mean_col(rows_10_4, 4))
    print('Field L2 10.4 phi mean =', mean_col(rows_10_4, 5))
    print('Field L2 10.5 w mean =', mean_col(rows_10_5, 4))
    print('Field L2 10.5 phi mean =', mean_col(rows_10_5, 5))


if __name__ == '__main__':
    main()
