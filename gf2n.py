# 50.042 FCS Lab 5 Modular Arithmetic
# Year 2025

class Polynomial2:
    def __init__(self, coeffs):
        # coeffs[i] is the 0/1 bit for x^i
        self.coeffs = coeffs.copy()
        self._trim()

    def _trim(self):
        # remove leading zeros
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()

    def add(self, p2):
        max_len = max(len(self.coeffs), len(p2.coeffs))
        result = []
        for i in range(max_len):
            a = self.coeffs[i] if i < len(self.coeffs) else 0
            b = p2.coeffs[i]    if i < len(p2.coeffs)    else 0
            result.append(a ^ b)
        return Polynomial2(result)

    def sub(self, p2):
        # same as add in GF(2)
        return self.add(p2)

    def mul(self, p2, modp=None):
        # carry-less convolution
        res_len = len(self.coeffs) + len(p2.coeffs) - 1
        result = [0] * res_len
        for i in range(len(self.coeffs)):
            if self.coeffs[i] == 1:
                for j in range(len(p2.coeffs)):
                    if p2.coeffs[j] == 1:
                        result[i + j] ^= 1
        # optional reduction
        if modp is not None:
            m = len(modp.coeffs) - 1
            while len(result) - 1 >= m:
                shift = len(result) - 1 - m
                for k, bit in enumerate(modp.coeffs):
                    if bit:
                        result[k + shift] ^= 1
                while len(result) > 1 and result[-1] == 0:
                    result.pop()
        return Polynomial2(result)

    def div(self, p2):
        # long division: returns quotient, remainder
        r = self.coeffs.copy()
        q_len = len(self.coeffs) - len(p2.coeffs) + 1
        q = [0] * max(1, q_len)
        deg_b = len(p2.coeffs) - 1
        while len(r) - 1 >= deg_b and any(r):
            deg_r = len(r) - 1
            lead_r = r[-1]
            shift = deg_r - deg_b
            q[shift] = lead_r
            for i, bit in enumerate(p2.coeffs):
                if bit and lead_r:
                    r[i + shift] ^= 1
            while len(r) > 1 and r[-1] == 0:
                r.pop()
        while len(q) > 1 and q[-1] == 0:
            q.pop()
        return Polynomial2(q), Polynomial2(r)

    def __str__(self):
        if not any(self.coeffs):
            return '0'
        terms = []
        for i in range(len(self.coeffs)-1, -1, -1):
            if self.coeffs[i]:
                if i == 0:
                    terms.append('1')
                elif i == 1:
                    terms.append('x')
                else:
                    terms.append(f'x^{i}')
        return '+'.join(terms)

    def getInt(self):
        value = 0
        for i, bit in enumerate(self.coeffs):
            if bit:
                value |= (1 << i)
        return value


class GF2N:
    # AES affine transform matrix
    affinemat = [
        [1,0,0,0,1,1,1,1],
        [1,1,0,0,0,1,1,1],
        [1,1,1,0,0,0,1,1],
        [1,1,1,1,0,0,0,1],
        [1,1,1,1,1,0,0,0],
        [0,1,1,1,1,1,0,0],
        [0,0,1,1,1,1,1,0],
        [0,0,0,1,1,1,1,1]
    ]
    # AES constant vector (0x63)
    const = [1,1,0,0,0,1,1,0]

    def __init__(self, x, n=8, ip=None):
        self.n = n
        self.ip = ip
        coeffs = []
        for i in range(x.bit_length()):
            coeffs.append((x >> i) & 1)
        poly = Polynomial2(coeffs)
        if self.ip is not None:
            _, rem = poly.div(self.ip)
            self.p = rem
        else:
            self.p = poly
        self.x = self.p.getInt()

    def add(self, g2):
        s = self.p.add(g2.p)
        return GF2N(s.getInt(), self.n, self.ip)

    def sub(self, g2):
        d = self.p.sub(g2.p)
        return GF2N(d.getInt(), self.n, self.ip)

    def mul(self, g2):
        m = self.p.mul(g2.p, self.ip)
        return GF2N(m.getInt(), self.n, self.ip)

    def div(self, g2):
        q, r = self.p.div(g2.p)
        return GF2N(q.getInt(), self.n, self.ip), GF2N(r.getInt(), self.n, self.ip)

    def getPolynomial2(self):
        return self.p

    def __str__(self):
        return str(self.p)

    def getInt(self):
        return self.p.getInt()

    def mulInv(self):
        if self.ip is None:
            raise ZeroDivisionError("No modulus; cannot invert in polynomial ring.")
        a, b = self.ip, self.p
        s0, s1 = Polynomial2([1]), Polynomial2([0])
        t0, t1 = Polynomial2([0]), Polynomial2([1])
        while any(b.coeffs):
            q, r = a.div(b)
            a, b = b, r
            s0, s1 = s1, s0.sub(q.mul(s1))
            t0, t1 = t1, t0.sub(q.mul(t1))
        if not (len(a.coeffs) == 1 and a.coeffs[0] == 1):
            raise ZeroDivisionError("Element has no inverse.")
        inv = t0
        return GF2N(inv.getInt(), self.n, self.ip)

    def affineMap(self):
        if self.n != 8:
            raise ValueError("Affine map only defined for GF(2^8)")
        in_bits = self.p.coeffs + [0]*(8 - len(self.p.coeffs))
        out = [0]*8
        for i in range(8):
            s = 0
            for j in range(8):
                s ^= (self.affinemat[i][j] & in_bits[j])
            s ^= self.const[i]
            out[i] = s
        return GF2N(Polynomial2(out).getInt(), 8, self.ip)


# (Tests follow exactly as given in the prompt)

print('\nTest 1')
print('======')
print('p1=x^5+x^2+x')
print('p2=x^3+x^2+1')
p1 = Polynomial2([0,1,1,0,0,1])
p2 = Polynomial2([1,0,1,1])
p3 = p1.add(p2)
print('p3= p1+p2 =', p3)

print('\nTest 2')
print('======')
print('p4=x^7+x^4+x^3+x^2+x')
print('modp=x^8+x^7+x^5+x^4+1')
p4 = Polynomial2([0,1,1,1,1,0,0,1])
modp = Polynomial2([1,0,0,0,1,1,0,1,1])
p5 = p1.mul(p4, modp)
print('p5=p1*p4 mod (modp)=', p5)

print('\nTest 3')
print('======')
print('p6=x^12+x^7+x^2')
print('p7=x^8+x^4+x^3+x+1')
p6 = Polynomial2([0,0,1,0,0,0,0,1,0,0,0,0,1])
p7 = Polynomial2([1,1,0,1,1,0,0,0,1])
p8q, p8r = p6.div(p7)
print('q for p6/p7=', p8q)
print('r for p6/p7=', p8r)

print('\nTest 4')
print('======')
g1 = GF2N(100)
g2 = GF2N(5)
print('g1 =', g1.getPolynomial2())
print('g2 =', g2.getPolynomial2())
g3 = g1.add(g2)
print('g1+g2 =', g3)

print('\nTest 5')
print('======')
ip = Polynomial2([1,1,0,0,1])
print('irreducible polynomial', ip)
g4 = GF2N(0b1101,4,ip)
g5 = GF2N(0b110,4,ip)
print('g4 =', g4.getPolynomial2())
print('g5 =', g5.getPolynomial2())
g6 = g4.mul(g5)
print('g4 x g5 =', g6.p)

print('\nTest 6')
print('======')
g7 = GF2N(0b1000010000100,13,None)
g8 = GF2N(0b100011011,13,None)
print('g7 =', g7.getPolynomial2())
print('g8 =', g8.getPolynomial2())
q, r = g7.div(g8)
print('g7/g8 =')
print('q =', q.getPolynomial2())
print('r =', r.getPolynomial2())

print('\nTest 7')
print('======')
ip = Polynomial2([1,1,0,0,1])
print('irreducible polynomial', ip)
g9 = GF2N(0b101,4,ip)
print('g9 =', g9.getPolynomial2())
print('inverse of g9 =', g9.mulInv().getPolynomial2())

print('\nTest 8')
print('======')
ip = Polynomial2([1,1,0,1,1,0,0,0,1])
print('irreducible polynomial', ip)
g10 = GF2N(0xc2,8,ip)
print('g10 = 0xc2')
g11 = g10.mulInv()
print('inverse of g10 = g11 =', hex(g11.getInt()))
g12 = g11.affineMap()
print('affine map of g11 =', hex(g12.getInt()))

# --- GF(2^4) table generation appended here ---

# initialize GF(2^4) with modulus x^4 + x^3 + 1
table_mod = Polynomial2([1,0,0,1,1])
elems = [GF2N(i, 4, table_mod) for i in range(16)]

with open('table1.txt', 'w') as out:
    # write addition table header
    out.write("Addition Table (16×16)\n")
    # each row: element a added to all b
    for a in elems:
        line = " ".join(str(a.add(b).getInt()) for b in elems)
        out.write(line + "\n")

    # write multiplication table header with a leading blank line
    out.write("\nMultiplication Table (16×16)\n")
    for a in elems:
        line = " ".join(str(a.mul(b).getInt()) for b in elems)
        out.write(line + "\n")

print("\nGenerated table1.txt with GF(2^4) addition and multiplication tables.")
