import math
import itertools as it
from collections import defaultdict

import numpy as np


# ============================================================
# Distribution Utilities
# ============================================================

def cbd_pmf(eta):
    pmf = defaultdict(float)
    for a in range(eta + 1):
        pa = math.comb(eta, a) / (2 ** eta)
        for b in range(eta + 1):
            pb = math.comb(eta, b) / (2 ** eta)
            pmf[a - b] += pa * pb
    return dict(sorted(pmf.items()))


def law_convolution(A, B):
    C = {}
    for a in A:
        for b in B:
            c = a + b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def clean_dist(A):
    """Preserve original pruning exactly."""
    B = {}
    for x, y in A.items():
        if y > 2 ** (-300):
            B[x] = y
    return B


def iter_law_convolution(A, i):
    """Exact original double-and-add structure."""
    D = {0: 1.0}
    i_bin = bin(i)[2:]
    for ch in i_bin:
        D = law_convolution(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A)
            D = clean_dist(D)
    return D


# ============================================================
# Modular Arithmetic
# ============================================================

def mod_switch(x, q, rq):
    return int(round(1.0 * rq * x / q) % rq)


def mod_centered(x, q):
    a = x % q
    if a < q / 2:
        return a
    return a - q


def mod_sym(x, q):
    r = x % q
    if r > q // 2:
        r -= q
    return r


# ============================================================
# Error Distributions
# ============================================================

def build_mod_switching_error_law(q, rq):
    D = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch(y, rq, q)
        d = mod_centered(x - z, q)
        D[d] = D.get(d, 0) + 1.0 / q
    return D


# ============================================================
# Core Distribution Computation: Compute the Distribution D_i= Pr{<P_i,g-g'>} in our paper and [11]
# *Pr(\Delta n_e)=P_1^{*kn/\ell}*P_2^{*kn/\ell}*P_3
# *Pr(<\Delta n_e, g-g'>)=D_1^{*kn/\ell}*D_2^{*kn/\ell}*D_3
# ============================================================

#Computing D_1 and D_2
def compute_distribution_S(B1, B2, d):
    N = len(d)
    result = defaultdict(float)

    for b in it.product(B2.keys(), repeat=N):
        p_b = math.prod(B2[x] for x in b)

        alpha = [0] * N
        for j in range(N):
            total = 0
            for k in range(N):
                idx = (j + k) % N
                sign = -1 if (j + k) >= N else 1
                total += b[k] * d[idx] * sign
            alpha[j] = total

        dist_sum = {0: 1.0}
        for j in range(N):
            if alpha[j] == 0:
                continue

            aj_scaled = defaultdict(float)
            for x, p_x in B1.items():
                val = x * alpha[j]
                aj_scaled[val] += p_x

            temp = defaultdict(float)
            for v1, p1 in dist_sum.items():
                for v2, p2 in aj_scaled.items():
                    temp[v1 + v2] += p1 * p2

            dist_sum = temp

        for val, p in dist_sum.items():
            result[val] += p * p_b

    return dict(result)

def scale_law(D, w):
    if w == 1:
        return dict(D)
    if w == -1:
        return {-x: p for x, p in D.items()}
    return {w * x: p for x, p in D.items()}

#Computing D_3
def weighted_sum_distribution(D, d):
    result = {0: 1.0}

    weight_counts = defaultdict(int)
    for w in d:
        if w != 0:
            weight_counts[w] += 1

    for w, m in weight_counts.items():
        scaled = scale_law(D, w)
        conv = iter_law_convolution(scaled, m)
        result = law_convolution(result, conv)
        result = clean_dist(result)

    return result


# ============================================================
# MAIN 
# ============================================================


def main():
    print('Computing DFR of 4D-GTD Code')
    
if __name__ == "__main__":
    main()


    Q = 3329
    N = 4   # GTC code dimension (or block length)
    eta = 2
    k = 4
    p = 6

    dv = 5
    du = 11

    print("KYBER1024:", "dv =", dv, ",", "du =", du)

    knN = int(k * 256 / N)

    # 4D-GTD Codebook (16 codewords selected from D4 lattice)
    D4_16_points = {
        0: (0, 0, 0, 0),
        1: (4, 2, 4, 0),
        2: (3, 3, 3, 3),
        3: (2, 0, 4, 2),
        4: (2, 4, 2, 0),
        5: (3, 1, 1, 1),
        6: (1, 1, 3, 5),
        7: (1, 5, 1, 3),
        8: (0, 2, 2, 2),
        9: (3, 5, 5, 5),
        10: (1, 3, 5, 1),
        11: (0, 4, 4, 4),
        12: (4, 0, 2, 4),
        13: (4, 4, 0, 2),
        14: (5, 1, 5, 3),
        15: (5, 5, 3, 1),
    }

    size_codebook = len(D4_16_points)

    # Distributions
    B1 = cbd_pmf(eta)
    B2 = B1

    cu = build_mod_switching_error_law(Q, 2 ** du)
    C = law_convolution(B1, cu)

    cv = build_mod_switching_error_law(Q, 2 ** dv)
    e2 = cbd_pmf(eta)
    e2cv = law_convolution(cv, e2)

    sum_of_values = 0  # IMPORTANT: preserved global accumulation

    for j in range(size_codebook):
        code_copy = list(range(size_codebook))
        del code_copy[j]

        index = 1

        for i in code_copy:
            diff = [
                mod_sym(vi - ci, p)
                for vi, ci in zip(D4_16_points[j], D4_16_points[i])
            ]

            dist_squared = sum(x ** 2 for x in diff)

            D1 = compute_distribution_S(B1, B2, diff)
            Tempt_D1 = iter_law_convolution(D1, knN)

            D2 = compute_distribution_S(C, B1, diff)
            Tempt_D2 = iter_law_convolution(D2, knN)

            Tempt_D3 = weighted_sum_distribution(e2cv, diff)

           
            T = law_convolution(Tempt_D1, Tempt_D2)
            T = law_convolution(T, Tempt_D3)

            for key, value in T.items():
                if key > Q / (p * 2) * dist_squared:
                    sum_of_values += value

            if sum_of_values > 0:
                print("Code", j, "Itr", index,
                      "Sum DFR (in log2):", math.log(sum_of_values) / math.log(2))

            index += 1

    print("\nFinal Total DFR (in log2):", (math.log(sum_of_values)+math.log(256/N/(2**N))) / math.log(2))


# ============================================================
# ENTRY POINT
# ============================================================




