from math import comb


def eta_inverse_coefficients(power: int, order: int) -> list[int]:
    coeffs = [1] + [0] * order
    for n in range(1, order + 1):
        next_coeffs = [0] * (order + 1)
        for i, a_i in enumerate(coeffs):
            if a_i == 0:
                continue
            j = 0
            while i + n * j <= order:
                next_coeffs[i + n * j] += a_i * comb(power + j - 1, j)
                j += 1
        coeffs = next_coeffs
    return coeffs


def convolve(a: list[int], b: list[int], order: int) -> list[int]:
    return [sum(a[j] * b[i - j] for j in range(i + 1)) for i in range(order + 1)]


def test_half_k3_unit_specialisation_gives_eta_inverse_four():
    assert eta_inverse_coefficients(4, 4) == [1, 4, 14, 40, 105]


def test_e4_eta_inverse_twelve_coefficients():
    e4 = [1, 240, 2160, 6720, 17520]
    assert convolve(e4, eta_inverse_coefficients(12, 4), 4) == [
        1,
        252,
        5130,
        54760,
        419895,
    ]
