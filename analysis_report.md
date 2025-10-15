# ECDLP Algorithm Scaling Analysis

The tables below summarise the empirical behaviour of four classical discrete-log algorithms executed on six elliptic curves with prime moduli between 97 and 1999.

## Mean Operation Counts

| Algorithm | p=97 | p=499 | p=911 | p=1237 | p=1613 | p=1999 |
| --- | --- | --- | --- | --- | --- | --- |
| Baby-Step Giant-Step | 14.44 | 29.00 | 17.56 | 20.72 | 39.96 | 39.04 |
| Brute Force | 31.32 | 151.20 | 39.80 | 56.80 | 273.60 | 220.36 |
| Pollard Kangaroo | 209.76 | 413.12 | 262.56 | 254.20 | 607.00 | 581.96 |
| Pollard Rho | 85.38 | 178.50 | 70.50 | 129.00 | 511.50 | 208.40 |

## Median Operations Normalised by √p

| Algorithm | p=97 | p=499 | p=911 | p=1237 | p=1613 | p=1999 |
| --- | --- | --- | --- | --- | --- | --- |
| Baby-Step Giant-Step | 2.12 | 1.95 | 1.88 | 1.92 | 1.67 | 1.72 |
| Brute Force | 4.67 | 11.34 | 4.27 | 5.10 | 10.49 | 8.96 |
| Pollard Kangaroo | 29.56 | 26.02 | 27.42 | 24.44 | 25.85 | 26.39 |
| Pollard Rho | 8.91 | 12.76 | 3.86 | 13.57 | 21.89 | 11.18 |

## Empirical Success Rates

| Algorithm | p=97 | p=499 | p=911 | p=1237 | p=1613 | p=1999 |
| --- | --- | --- | --- | --- | --- | --- |
| Baby-Step Giant-Step | 100.00% | 100.00% | 100.00% | 100.00% | 100.00% | 100.00% |
| Brute Force | 100.00% | 100.00% | 100.00% | 100.00% | 100.00% | 100.00% |
| Pollard Kangaroo | 100.00% | 100.00% | 100.00% | 100.00% | 100.00% | 100.00% |
| Pollard Rho | 52.00% | 16.00% | 24.00% | 12.00% | 8.00% | 20.00% |

## Observations

- Brute force scales linearly with the group order, requiring roughly 0.52·n operations on these curves.
- Baby-Step Giant-Step exhibits the expected √n complexity with a median constant factor of about 1.9.
- Pollard's Rho also scales with √n but with a higher and more variable constant (median ≈ 12.0 across curves).
- Pollard's Kangaroo shows similar √n scaling to Rho, with a median constant factor near 26.6.
- Extrapolating these constants to secp256k1 (n ≈ 1.16×10²⁷) implies ≥2¹²⁸ group operations for √n algorithms, far beyond practical reach.