
# GF2ⁿ Polynomial Arithmetic & Table Generator

A Python toolkit for working with polynomials over GF(2) and constructing finite‐field tables for GF(2ⁿ). Defines two core classes and produces 16×16 addition and multiplication tables for GF(2⁴) (modulus x⁴ + x³ + 1).

---

## 🚀 Features

- **`Polynomial2`**  
  - Represents arbitrary-degree polynomials over GF(2).  
  - Bitwise-XOR addition/subtraction.  
  - Carry-less multiplication with optional modular reduction.  
  - Long division yielding quotient & remainder.

- **`GF2N`**  
  - Wraps an integer as a field element in GF(2ⁿ).  
  - Automatic reduction modulo a given irreducible polynomial.  
  - Field operations (`add`, `sub`, `mul`) delegate to `Polynomial2`.  
  - (Optional) multiplicative inverse via Extended Euclidean Algorithm.  
  - AES-style affine map for n = 8 (unused in GF(2⁴) examples).

- **Table Generator**  
  - Enumerates all 16 elements of GF(2⁴).  
  - Outputs 16×16 **Addition** and **Multiplication** tables (mod x⁴ + x³ + 1) to `table1.txt`.

---

## 📋 Installation

1. **Clone this repo**  
   ```bash
   git clone https://github.com/YourUser/GF2n-Polynomial-Arithmetic.git
   cd GF2n-Polynomial-Arithmetic
