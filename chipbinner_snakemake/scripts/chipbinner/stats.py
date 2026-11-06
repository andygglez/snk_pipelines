"""Numerical helpers ported from R: hexbin::hexbin, MASS::kde2d, base::pretty."""
from __future__ import annotations

import math
import sys

import numpy as np
import pandas as pd

SQRT3 = math.sqrt(3.0)


# --------------------------------------------------------------------------- hexbin
def hexbin(x: np.ndarray, y: np.ndarray, xbins: int = 30, shape: float = 1.0,
           xbnds=None, ybnds=None):
    """Port of hexbin::hexbin(x, y, xbins, shape, IDs = TRUE) + hcell2xy().

    Returns (cells, cell_id_per_point) where ``cells`` is a DataFrame with the
    non-empty cells (cell, x, y, count) in increasing cell order and the second
    array gives the cell of every input point (hex@cID).
    Pointy-top hexagons; geometry for drawing is given by :func:`hex_geometry`.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xbnds = (x.min(), x.max()) if xbnds is None else xbnds
    ybnds = (y.min(), y.max()) if ybnds is None else ybnds
    jmax = math.floor(xbins + 1.5001)
    c1 = xbins / (xbnds[1] - xbnds[0])
    c2 = xbins * shape / ((ybnds[1] - ybnds[0]) * SQRT3)
    iinc = 2 * jmax
    lat = jmax + 1

    sx = c1 * (x - xbnds[0])
    sy = c2 * (y - ybnds[0])
    j1 = np.floor(sx + 0.5).astype(np.int64)  # Fortran int() truncates; sx >= 0
    i1 = np.floor(sy + 0.5).astype(np.int64)
    dist1 = (sx - j1) ** 2 + 3.0 * (sy - i1) ** 2
    j2 = np.floor(sx).astype(np.int64)
    i2 = np.floor(sy).astype(np.int64)
    dist2 = (sx - j2 - 0.5) ** 2 + 3.0 * (sy - i2 - 0.5) ** 2
    l_even = i1 * iinc + j1 + 1
    l_odd = i2 * iinc + j2 + lat
    cid = np.where(dist1 < 0.25, l_even,
                   np.where(dist1 > 1.0 / 3.0, l_odd,
                            np.where(dist1 <= dist2, l_even, l_odd)))

    cells, counts = np.unique(cid, return_counts=True)
    c3 = (xbnds[1] - xbnds[0]) / xbins
    c4 = (ybnds[1] - ybnds[0]) * SQRT3 / (2 * shape * xbins)
    c = cells - 1
    i = c // jmax
    j = c % jmax
    cy = c4 * i + ybnds[0]
    cx = c3 * np.where(i % 2 == 0, j, j + 0.5) + xbnds[0]
    df = pd.DataFrame({"cell": cells, "x": cx, "y": cy, "count": counts})
    df.attrs.update(dx=c3, dy=c4)
    return df, cid


def hex_geometry(dx: float, dy: float):
    """Vertex offsets of a pointy-top hexagon of width dx and row spacing dy."""
    ry = dy / 1.5  # circumradius in y units (row spacing = 1.5 R)
    hx = dx / 2.0
    return np.array([[0, ry], [hx, ry / 2], [hx, -ry / 2], [0, -ry], [-hx, -ry / 2], [-hx, ry / 2]])


# --------------------------------------------------------------------------- kde2d
def bandwidth_nrd(x: np.ndarray) -> float:
    """MASS::bandwidth.nrd."""
    q1, q3 = np.quantile(x, [0.25, 0.75])
    h = (q3 - q1) / 1.34
    return 4 * 1.06 * min(np.std(x, ddof=1), h) * len(x) ** (-1 / 5)


def kde2d(x: np.ndarray, y: np.ndarray, n: int = 100, chunk: int = 20000):
    """MASS::kde2d(x, y, n) with the default bandwidths; returns (gx, gy, z[nx, ny])."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    gx = np.linspace(x.min(), x.max(), n)
    gy = np.linspace(y.min(), y.max(), n)
    hx = bandwidth_nrd(x) / 4
    hy = bandwidth_nrd(y) / 4
    z = np.zeros((n, n))
    norm = 1.0 / math.sqrt(2 * math.pi)
    for k in range(0, len(x), chunk):
        ax = (gx[:, None] - x[None, k:k + chunk]) / hx
        ay = (gy[:, None] - y[None, k:k + chunk]) / hy
        z += (norm * np.exp(-0.5 * ax ** 2)) @ (norm * np.exp(-0.5 * ay ** 2)).T
    z /= len(x) * hx * hy
    return gx, gy, z


# --------------------------------------------------------------------------- pretty
def pretty(values, n: int = 5) -> np.ndarray:
    """base::pretty(x, n) (R_pretty with the default high.u.bias = 1.5, u5.bias)."""
    lo, up = float(np.nanmin(values)), float(np.nanmax(values))
    h, h5 = 1.5, 0.5 + 1.5 * 1.5
    min_n = n // 3
    shrink = 0.75
    eps = 1e-10
    dx = up - lo
    if dx == 0 and up == 0:
        cell, i_small = 1.0, True
    else:
        cell = max(abs(lo), abs(up))
        u = 1 + (1 / (1 + h) if h5 >= 1.5 * h + 0.5 else 1.5 / (1 + h5))
        u *= max(1, n) * sys.float_info.epsilon
        i_small = dx < cell * u * 3
    if i_small:
        if cell > 10:
            cell = 9 + cell / 10
        cell *= shrink
        if min_n > 1:
            cell /= min_n
    else:
        cell = dx
        if n > 1:
            cell /= n
    base = 10.0 ** math.floor(math.log10(cell))
    unit = base
    if 2 * base - cell < h * (cell - unit):
        unit = 2 * base
        if 5 * base - cell < h5 * (cell - unit):
            unit = 5 * base
            if 10 * base - cell < h * (cell - unit):
                unit = 10 * base
    ns = math.floor(lo / unit + eps)
    nu = math.ceil(up / unit - eps)
    while ns * unit > lo + eps * unit:
        ns -= 1
    while nu * unit < up - eps * unit:
        nu += 1
    k = int(0.5 + nu - ns)
    if k < min_n:
        k = min_n - k
        if ns >= 0:
            nu += k // 2
            ns -= k // 2 + k % 2
        else:
            ns -= k // 2
            nu += k // 2 + k % 2
        ndiv = min_n
    else:
        ndiv = k
    return np.linspace(ns * unit, nu * unit, ndiv + 1)


def quantile_filter(df: pd.DataFrame, lo: float, hi: float, cols=("x", "y")) -> pd.DataFrame:
    """dplyr::filter(x > quantile(x, lo), x < quantile(x, hi), y > ..., y < ...)."""
    keep = np.ones(len(df), bool)
    for c in cols:
        v = df[c].to_numpy(float)
        ql, qh = np.quantile(v, [lo, hi])  # R type-7 quantiles == numpy "linear"
        keep &= (v > ql) & (v < qh)
    return df[keep].reset_index(drop=True)


def expand_range(lo: float, hi: float, mul: float = 0.05):
    """scales::expand_range(c(lo, hi), mul)."""
    d = (hi - lo) * mul
    return lo - d, hi + d
