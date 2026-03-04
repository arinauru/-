import itertools, numpy as np, pandas as pd

G = np.array([
    [1,1,0,1,0],
    [0,1,1,0,1],
    [1,0,1,1,0]
], dtype=int)

k, n = G.shape
print("Given G (k x n):")
print(G)
print(f"k={k}, n={n}, R=k/n={k}/{n}={k/n:.3f}")
print()

codewords = []
for m in itertools.product([0,1], repeat=k):
    m_vec = np.array(m, dtype=int)
    c = (m_vec @ G) % 2
    codewords.append((tuple(m_vec.tolist()), tuple(c.tolist())))

def weight(v):
    return int(sum(v))

print("All 2^k codewords (m -> c) and weights:")
rows = []
for m, c in codewords:
    rows.append({"m": "".join(map(str,m)), "c": "".join(map(str,c)), "w(c)": weight(c)})
df = pd.DataFrame(rows)
print(df.to_string(index=False))
print()

nonzero_weights = [weight(c) for m,c in codewords if any(c)]
dmin = min(nonzero_weights) if nonzero_weights else 0
print(f"d_min = min_{'{c!=0}'} w(c) = {dmin}")
t = (dmin - 1)//2
print(f"t = floor((d_min-1)/2) = {t}  (guaranteed correctable errors)")
print()

C = [np.array(list(c), dtype=int) for _,c in codewords]
D = np.zeros((len(C), len(C)), dtype=int)
for i in range(len(C)):
    for j in range(len(C)):
        D[i,j] = weight((C[i] ^ C[j]).tolist())
print("Distance matrix D (rows/cols in the same order as above table):")
print(D)
print()

def gf2_rank(A):
    A = A.copy() % 2
    r = 0
    m, nA = A.shape
    for col in range(nA):
        pivot = None
        for i in range(r, m):
            if A[i,col] == 1:
                pivot = i
                break
        if pivot is None:
            continue
        A[[r,pivot]] = A[[pivot,r]]
        for i in range(m):
            if i != r and A[i,col] == 1:
                A[i] ^= A[r]
        r += 1
        if r == m:
            break
    return r

solutions = []
for h in itertools.product([0,1], repeat=n):
    h_vec = np.array(h, dtype=int)
    if np.all((G @ h_vec) % 2 == 0):
        solutions.append(h_vec)
solutions = np.array(solutions, dtype=int)

print(f"Number of vectors h with G h^T=0: {len(solutions)} (expected 2^(n-k)={2**(n-k)})")
print("All such h vectors:")
print(["".join(map(str,h.tolist())) for h in solutions])
print()

basis = []
for h in solutions:
    if not np.any(h):
        continue
    candidate = basis + [h]
    M = np.array(candidate, dtype=int)
    if gf2_rank(M) == len(candidate):
        basis.append(h)
    if len(basis) == n-k:
        break

H = np.array(basis, dtype=int)
print("One possible parity-check matrix H ((n-k) x n):")
print(H)
print()

print("Check G * H^T mod 2 (should be all zeros):")
print((G @ H.T) % 2)
print()

def syndrome(v):
    v = np.array(v, dtype=int)
    return tuple(((v @ H.T) % 2).tolist())

errs = []
for e in itertools.product([0,1], repeat=n):
    e_vec = np.array(e, dtype=int)
    s = syndrome(e_vec)
    errs.append((e_vec, s, weight(e_vec)))

leader = {}
for e_vec, s, w in sorted(errs, key=lambda x: (x[2], "".join(map(str,x[0].tolist())))):
    if s not in leader:
        leader[s] = e_vec

print("Syndrome -> coset leader e (minimal weight, then lexicographically smallest):")
for s in sorted(leader.keys()):
    print(f"S={''.join(map(str,s))}  ->  e={''.join(map(str,leader[s].tolist()))}  (w={weight(leader[s])})")
print()

print("Syndromes for single-bit errors (weight 1):")
for i in range(n):
    e = np.zeros(n, dtype=int)
    e[i] = 1
    print(f"e at pos {i+1}: {''.join(map(str,e.tolist()))} -> S={''.join(map(str,syndrome(e)))}")
print()

def decode(y_bits):
    y = np.array(list(map(int, y_bits)), dtype=int)
    s = syndrome(y)
    e = leader[s]
    c_hat = y ^ e
    return s, e, c_hat

example_y = "11111"
s, e, c_hat = decode(example_y)
print("Example decoding:")
print("y =", example_y)
print("S =", "".join(map(str,s)))
print("chosen e =", "".join(map(str,e.tolist())))
print("decoded c =", "".join(map(str,c_hat.tolist())))