from itertools import product

'''Строим проверочную матрицу H'''
def build_hamming_parity_check(r: int):
    n = 2**r - 1
    columns = []

    for x in range(1, n + 1):
        bits = [(x >> i) & 1 for i in range(r)]
        columns.append(bits)

    H = [[columns[j][i] for j in range(n)] for i in range(r)]
    return H

'''Умножение вектора на матрицу'''
def vec_mat_mul(v, M):
    k = len(v)
    n = len(M[0])
    res = [0] * n

    for j in range(n):
        s = 0
        for i in range(k):
            s ^= (v[i] & M[i][j])
        res[j] = s

    return res


def hamming_weight(v):
    return sum(v)


def hamming_distance(a, b):
    return sum(x ^ y for x, y in zip(a, b))


def all_binary_vectors(length):
    return [list(v) for v in product([0, 1], repeat=length)]


def print_matrix(M, name):
    print(f"{name} =")
    for row in M:
        print(" ", row)
    print()


def verify_simplex_dual_hamming(r: int):
    print("=" * 60)
    print(f"Для r = {r}")
    print("=" * 60)

    n = 2**r - 1
    k_hamming = n - r

    print(f"Код Хэмминга имеет параметры: (n, k) = ({n}, {k_hamming})")
    print(f"Тогда дуальный код имеет параметры: ({n}, {r})")
    print()

    H = build_hamming_parity_check(r)
    print_matrix(H, "Проверочная матрица H кода Хэмминга")
    print(f"Размер H: {len(H)} x {len(H[0])}")
    print()

    messages = all_binary_vectors(r)
    codewords = []

    for m in messages:
        c = vec_mat_mul(m, H)
        codewords.append((m, c))

    print("Все кодовые слова дуального кода:")
    for m, c in codewords:
        print(f"m = {m} -> c = {c}, вес = {hamming_weight(c)}")
    print()

    print(f"Число кодовых слов: {len(codewords)}")
    print(f"Ожидаемое число: 2^{r} = {2**r}")
    print()

    nonzero_weights = []
    for m, c in codewords:
        if any(c):
            nonzero_weights.append(hamming_weight(c))

    print("Веса всех ненулевых кодовых слов:")
    print(nonzero_weights)
    print()

    expected_weight = 2**(r - 1)
    print(f"Ожидаемый вес для симплекс-кода: 2^({r}-1) = {expected_weight}")
    print()

    if all(w == expected_weight for w in nonzero_weights):
        print("Проверка весов пройдена:")
        print(f"все ненулевые кодовые слова имеют вес {expected_weight}.")
    else:
        print("Проверка весов не пройдена.")
    print()

    distances = []
    for i in range(len(codewords)):
        for j in range(i + 1, len(codewords)):
            c1 = codewords[i][1]
            c2 = codewords[j][1]
            d = hamming_distance(c1, c2)
            distances.append(d)

    print("Все расстояния между различными парами кодовых слов:")
    print(distances)
    print()

    if all(d == expected_weight for d in distances):
        print("Проверка расстойний пройдена:")
        print(f"все расстояния между различными кодовыми словами равны {expected_weight}.")
    else:
        print("Проверка расстойний не пройдена.")
    print()

    if all(w == expected_weight for w in nonzero_weights) and all(d == expected_weight for d in distances):
        print("ИТОГ:")
        print("Дуальный коду Хэмминга код действительно является симплекс-кодом.")
    else:
        print("ИТОГ:")
        print("Проверка не подтвердила свойство симплекс-кода.")
    print()


verify_simplex_dual_hamming(3)
verify_simplex_dual_hamming(4)