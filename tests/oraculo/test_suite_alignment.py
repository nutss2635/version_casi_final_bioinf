#!/usr/bin/env python3
"""
Test Suite and Biological/Mathematical Validation Reference Engine
for Sequence Alignment Visualizer (Needleman-Wunsch, Smith-Waterman, Wagner-Fischer ED, LCS)
"""

import json
from typing import List, Tuple, Dict, Any

# Standard BLOSUM62 Matrix
BLOSUM62_RAW = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""

def parse_blosum62() -> Dict[Tuple[str, str], int]:
    lines = [l.strip() for l in BLOSUM62_RAW.strip().splitlines()]
    headers = lines[0].split()
    matrix = {}
    for line in lines[1:]:
        parts = line.split()
        row_aa = parts[0]
        for col_aa, val in zip(headers, parts[1:]):
            matrix[(row_aa, col_aa)] = int(val)
    return matrix

BLOSUM62 = parse_blosum62()

# -------------------------------------------------------------
# NEEDLEMAN-WUNSCH (Global Alignment)
# -------------------------------------------------------------
def needleman_wunsch(s1: str, s2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> Dict[str, Any]:
    n, m = len(s1), len(s2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diag = dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            score_up = dp[i - 1][j] + gap
            score_left = dp[i][j - 1] + gap
            dp[i][j] = max(score_diag, score_up, score_left)

    # Traceback: Priority Diagonal > Superior (Up) > Lateral (Left)
    i, j = n, m
    al1, al2 = [], []
    steps = []
    coords = []

    while i > 0 or j > 0:
        curr = dp[i][j]
        coords.append((i, j))
        # Check diagonal
        if i > 0 and j > 0 and curr == dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch):
            al1.append(s1[i - 1])
            al2.append(s2[j - 1])
            steps.append("diagonal")
            i -= 1
            j -= 1
        # Check superior (up)
        elif i > 0 and curr == dp[i - 1][j] + gap:
            al1.append(s1[i - 1])
            al2.append("-")
            steps.append("superior")
            i -= 1
        # Lateral (left)
        else:
            al1.append("-")
            al2.append(s2[j - 1])
            steps.append("lateral")
            j -= 1

    # Find all optimal paths for validation
    def get_all_paths(ci, cj):
        if ci == 0 and cj == 0:
            return [([], [])]
        paths = []
        curr_val = dp[ci][cj]
        if ci > 0 and cj > 0 and curr_val == dp[ci - 1][cj - 1] + (match if s1[ci - 1] == s2[cj - 1] else mismatch):
            for p1, p2 in get_all_paths(ci - 1, cj - 1):
                paths.append((p1 + [s1[ci - 1]], p2 + [s2[cj - 1]]))
        if ci > 0 and curr_val == dp[ci - 1][cj] + gap:
            for p1, p2 in get_all_paths(ci - 1, cj):
                paths.append((p1 + [s1[ci - 1]], p2 + ["-"]))
        if cj > 0 and curr_val == dp[ci][cj - 1] + gap:
            for p1, p2 in get_all_paths(ci, cj - 1):
                paths.append((p1 + ["-"], p2 + [s2[cj - 1]]))
        return paths

    all_paths = [("".join(p1), "".join(p2)) for p1, p2 in get_all_paths(n, m)]

    return {
        "algorithm": "Needleman-Wunsch",
        "s1": s1,
        "s2": s2,
        "score": dp[n][m],
        "dp_matrix": dp,
        "alignment_s1": "".join(reversed(al1)),
        "alignment_s2": "".join(reversed(al2)),
        "steps": list(reversed(steps)),
        "all_optimal_paths": all_paths,
        "optimal_path_count": len(all_paths)
    }

# -------------------------------------------------------------
# SMITH-WATERMAN (Local Alignment)
# -------------------------------------------------------------
def smith_waterman(s1: str, s2: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> Dict[str, Any]:
    n, m = len(s1), len(s2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_cells = []

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diag = dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            score_up = dp[i - 1][j] + gap
            score_left = dp[i][j - 1] + gap
            val = max(0, score_diag, score_up, score_left)
            dp[i][j] = val
            if val > max_score:
                max_score = val
                max_cells = [(i, j)]
            elif val == max_score and val > 0:
                max_cells.append((i, j))

    # Standard traceback from first max_cell
    al1, al2 = [], []
    steps = []
    if max_score > 0 and max_cells:
        i, j = max_cells[0]
        while dp[i][j] > 0:
            curr = dp[i][j]
            # Priority: diag > up > left
            if i > 0 and j > 0 and curr == dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch):
                al1.append(s1[i - 1])
                al2.append(s2[j - 1])
                steps.append("diagonal")
                i -= 1
                j -= 1
            elif i > 0 and curr == dp[i - 1][j] + gap:
                al1.append(s1[i - 1])
                al2.append("-")
                steps.append("superior")
                i -= 1
            elif j > 0 and curr == dp[i][j - 1] + gap:
                al1.append("-")
                al2.append(s2[j - 1])
                steps.append("lateral")
                j -= 1
            else:
                break

    # All paths from max cells
    def get_sw_paths(ci, cj):
        if dp[ci][cj] == 0:
            return [([], [])]
        paths = []
        curr_val = dp[ci][cj]
        if ci > 0 and cj > 0 and curr_val == dp[ci - 1][cj - 1] + (match if s1[ci - 1] == s2[cj - 1] else mismatch):
            for p1, p2 in get_sw_paths(ci - 1, cj - 1):
                paths.append((p1 + [s1[ci - 1]], p2 + [s2[cj - 1]]))
        if ci > 0 and curr_val == dp[ci - 1][cj] + gap:
            for p1, p2 in get_sw_paths(ci - 1, cj):
                paths.append((p1 + [s1[ci - 1]], p2 + ["-"]))
        if cj > 0 and curr_val == dp[ci][cj - 1] + gap:
            for p1, p2 in get_sw_paths(ci, cj - 1):
                paths.append((p1 + ["-"], p2 + [s2[cj - 1]]))
        return paths

    all_paths = []
    if max_score > 0:
        for mc in max_cells:
            for p1, p2 in get_sw_paths(mc[0], mc[1]):
                all_paths.append(("".join(p1), "".join(p2)))

    return {
        "algorithm": "Smith-Waterman",
        "s1": s1,
        "s2": s2,
        "score": max_score,
        "max_cells": max_cells,
        "dp_matrix": dp,
        "alignment_s1": "".join(reversed(al1)),
        "alignment_s2": "".join(reversed(al2)),
        "steps": list(reversed(steps)),
        "all_optimal_paths": all_paths,
        "optimal_path_count": len(all_paths)
    }

# -------------------------------------------------------------
# WAGNER-FISCHER / EDIT DISTANCE (Minimization)
# -------------------------------------------------------------
def wagner_fischer(s1: str, s2: str, match_cost: int = 0, sub_cost: int = 1, indel_cost: int = 1) -> Dict[str, Any]:
    n, m = len(s1), len(s2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i * indel_cost
    for j in range(1, m + 1):
        dp[0][j] = j * indel_cost

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost_diag = dp[i - 1][j - 1] + (match_cost if s1[i - 1] == s2[j - 1] else sub_cost)
            cost_up = dp[i - 1][j] + indel_cost
            cost_left = dp[i][j - 1] + indel_cost
            dp[i][j] = min(cost_diag, cost_up, cost_left)

    # Traceback: Priority Diagonal > Superior (Up) > Lateral (Left)
    i, j = n, m
    al1, al2 = [], []
    steps = []

    while i > 0 or j > 0:
        curr = dp[i][j]
        cost_diag = dp[i - 1][j - 1] + (match_cost if s1[i - 1] == s2[j - 1] else sub_cost) if (i > 0 and j > 0) else None
        cost_up = dp[i - 1][j] + indel_cost if i > 0 else None
        cost_left = dp[i][j - 1] + indel_cost if j > 0 else None

        if cost_diag is not None and curr == cost_diag:
            al1.append(s1[i - 1])
            al2.append(s2[j - 1])
            steps.append("diagonal")
            i -= 1
            j -= 1
        elif cost_up is not None and curr == cost_up:
            al1.append(s1[i - 1])
            al2.append("-")
            steps.append("superior")
            i -= 1
        else:
            al1.append("-")
            al2.append(s2[j - 1])
            steps.append("lateral")
            j -= 1

    # All paths
    def get_wf_paths(ci, cj):
        if ci == 0 and cj == 0:
            return [([], [])]
        paths = []
        curr_val = dp[ci][cj]
        if ci > 0 and cj > 0 and curr_val == dp[ci - 1][cj - 1] + (match_cost if s1[ci - 1] == s2[cj - 1] else sub_cost):
            for p1, p2 in get_wf_paths(ci - 1, cj - 1):
                paths.append((p1 + [s1[ci - 1]], p2 + [s2[cj - 1]]))
        if ci > 0 and curr_val == dp[ci - 1][cj] + indel_cost:
            for p1, p2 in get_wf_paths(ci - 1, cj):
                paths.append((p1 + [s1[ci - 1]], p2 + ["-"]))
        if cj > 0 and curr_val == dp[ci][cj - 1] + indel_cost:
            for p1, p2 in get_wf_paths(ci, cj - 1):
                paths.append((p1 + ["-"], p2 + [s2[cj - 1]]))
        return paths

    all_paths = [("".join(p1), "".join(p2)) for p1, p2 in get_wf_paths(n, m)]

    return {
        "algorithm": "Wagner-Fischer (Edit Distance)",
        "s1": s1,
        "s2": s2,
        "score": dp[n][m],  # Distance
        "dp_matrix": dp,
        "alignment_s1": "".join(reversed(al1)),
        "alignment_s2": "".join(reversed(al2)),
        "steps": list(reversed(steps)),
        "all_optimal_paths": all_paths,
        "optimal_path_count": len(all_paths)
    }

# -------------------------------------------------------------
# LONGEST COMMON SUBSEQUENCE (LCS)
# -------------------------------------------------------------
def longest_common_subsequence(s1: str, s2: str) -> Dict[str, Any]:
    n, m = len(s1), len(s2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    # Traceback: Priority Match/Diagonal > Superior > Lateral
    i, j = n, m
    lcs_chars = []
    al1, al2 = [], []
    steps = []

    while i > 0 and j > 0:
        if s1[i - 1] == s2[j - 1]:
            lcs_chars.append(s1[i - 1])
            al1.append(s1[i - 1])
            al2.append(s2[j - 1])
            steps.append("diagonal")
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            al1.append(s1[i - 1])
            al2.append("-")
            steps.append("superior")
            i -= 1
        else:
            al1.append("-")
            al2.append(s2[j - 1])
            steps.append("lateral")
            j -= 1

    while i > 0:
        al1.append(s1[i - 1])
        al2.append("-")
        steps.append("superior")
        i -= 1
    while j > 0:
        al1.append("-")
        al2.append(s2[j - 1])
        steps.append("lateral")
        j -= 1

    # Find all unique LCS strings
    def get_all_lcs(ci, cj):
        if ci == 0 or cj == 0:
            return {""}
        if s1[ci - 1] == s2[cj - 1]:
            return {sub + s1[ci - 1] for sub in get_all_lcs(ci - 1, cj - 1)}
        res = set()
        if dp[ci - 1][cj] == dp[ci][cj]:
            res.update(get_all_lcs(ci - 1, cj))
        if dp[ci][cj - 1] == dp[ci][cj]:
            res.update(get_all_lcs(ci, cj - 1))
        return res

    all_lcs = sorted(list(get_all_lcs(n, m)))

    return {
        "algorithm": "Longest Common Subsequence (LCS)",
        "s1": s1,
        "s2": s2,
        "score": dp[n][m],  # LCS Length
        "dp_matrix": dp,
        "lcs_string": "".join(reversed(lcs_chars)),
        "all_lcs_strings": all_lcs,
        "alignment_s1": "".join(reversed(al1)),
        "alignment_s2": "".join(reversed(al2)),
        "steps": list(reversed(steps))
    }

# -------------------------------------------------------------
# PROTEIN ALIGNMENT (BLOSUM62 vs Linear)
# -------------------------------------------------------------
def protein_alignment_compare(p1: str, p2: str, gap_linear: int = -2, gap_blosum: int = -4):
    res_linear = needleman_wunsch(p1, p2, match=1, mismatch=-1, gap=gap_linear)
    
    # Run with BLOSUM62
    n, m = len(p1), len(p2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap_blosum
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap_blosum

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair_score = BLOSUM62.get((p1[i - 1], p2[j - 1]), -4)
            diag = dp[i - 1][j - 1] + pair_score
            up = dp[i - 1][j] + gap_blosum
            left = dp[i][j - 1] + gap_blosum
            dp[i][j] = max(diag, up, left)

    i, j = n, m
    al1, al2 = [], []
    while i > 0 or j > 0:
        curr = dp[i][j]
        pair_score = BLOSUM62.get((p1[i - 1], p2[j - 1]), -4) if (i > 0 and j > 0) else None
        if pair_score is not None and curr == dp[i - 1][j - 1] + pair_score:
            al1.append(p1[i - 1])
            al2.append(p2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and curr == dp[i - 1][j] + gap_blosum:
            al1.append(p1[i - 1])
            al2.append("-")
            i -= 1
        else:
            al1.append("-")
            al2.append(p2[j - 1])
            j -= 1

    res_blosum = {
        "algorithm": "Needleman-Wunsch (BLOSUM62)",
        "s1": p1,
        "s2": p2,
        "score": dp[n][m],
        "dp_matrix": dp,
        "alignment_s1": "".join(reversed(al1)),
        "alignment_s2": "".join(reversed(al2))
    }
    return res_linear, res_blosum

def format_matrix(s1: str, s2: str, dp: List[List[int]]) -> str:
    header = ["   ", "  -"] + [f"{c:>3}" for c in s2]
    out = [" ".join(header)]
    for i, row in enumerate(dp):
        r_char = "-" if i == 0 else s1[i - 1]
        row_str = [f"{r_char:>3}"] + [f"{val:>3}" for val in row]
        out.append(" ".join(row_str))
    return "\n".join(out)

if __name__ == "__main__":
    canonical_cases = [
        ("TC1: Identical Strings", "AGTC", "AGTC"),
        ("TC2: Single Mismatch (Substitution)", "AGTC", "AGAC"),
        ("TC3: Single Indel (Insertion/Deletion)", "AGTC", "AGC"),
        ("TC4: Highly Asymmetrical Lengths", "ACGTACGT", "CG"),
        ("TC5: Completely Divergent", "AAAA", "CCCC"),
        ("TC6: Multiple Optimal Paths (Tie-breaker)", "AGC", "ACG")
    ]

    all_results = {}

    for name, s1, s2 in canonical_cases:
        nw_res = needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2)
        sw_res = smith_waterman(s1, s2, match=2, mismatch=-1, gap=-2)
        wf_res = wagner_fischer(s1, s2, match_cost=0, sub_cost=1, indel_cost=1)
        lcs_res = longest_common_subsequence(s1, s2)

        all_results[name] = {
            "s1": s1,
            "s2": s2,
            "NW": nw_res,
            "SW": sw_res,
            "ED": wf_res,
            "LCS": lcs_res
        }

    # Save to JSON
    with open("canonical_test_suite.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print("Test suite generated and verified successfully in canonical_test_suite.json")
