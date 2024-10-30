import p_adic as pa
import algebraic_extension as ae
from math import lcm, gcd

class Alpha:
    poly = ae.Polynomial[int](0,0,0,0,1,1)
    def __init__(self, p: pa.pAdic):
        self.p = p.p
        self.still_digits = list()
        self.repeating_digits = list()
        l = len(p.digits[0])
        k = len(p.digits[1])
        denom = ae.Polynomial[int](1) - ae.prod(*[Alpha.poly for i in range(l)])
        num_still = ae.Polynomial(0)
        for i in p.digits[1]:
            num_still *= Alpha.poly
            num_still += i
        num_still *= denom
        num_repeating = ae.Polynomial(0)
        for j in p.digits[0]:
            num_repeating *= Alpha.poly
            num_repeating += j
        num_repeating *= ae.prod(*[ae.Polynomial(0,1) for i in range(k)])
        num = num_still + num_repeating
        # print(str(num) + "/" + str(denom))
        seen: dict[ae.Polynomial[int], int] = {}
        written = list()
        a = num[0]
        # iter = 0
        written.append(a)
        while num not in seen:
            seen[num] = len(written)
            num -= a * denom
            num = num >> 1
            if num[0] != num[0]%self.p:
                num -= ae.Polynomial(num[0]//self.p*self.p) - ae.Polynomial(num[0]//self.p) * Alpha.poly
            a = num[0]
            written.append(a)
            # print(str(a) + ":" + str(num) + "/" + str(denom))
            # iter+= 1
            # if iter > 100:
            #     raise Exception("Endless")
        self.still_digits = written[:seen[num]]
        self.repeating_digits = written[seen[num]:]

    def to_polynomial(self):
        denom = ae.Polynomial[int](1) - ae.Polynomial(*([0] * (len(self.repeating_digits)) + [1]))
        num = ae.Polynomial[int](*self.still_digits) * denom \
              + ae.Polynomial(*([0] * (len(self.still_digits)) + [1])) * ae.Polynomial[int](*self.repeating_digits)
        return num, denom

    @staticmethod
    def construct(still, infinite, p=3):
        a = Alpha(pa.pAdic.zero(p))
        a.still_digits = still
        a.repeating_digits = infinite
        return a

    def __str__(self):
        return f"{"".join(map(str, self.still_digits))}[{"".join(map(str, self.repeating_digits))}]"

    def _spew(self, i):
        loops, remainder = divmod(i, len(self.repeating_digits))
        self.still_digits.extend([x for x in self.repeating_digits for j in range(loops)] + self.repeating_digits[:remainder])
        self.repeating_digits = self.repeating_digits[remainder:] + self.repeating_digits[:remainder]

    def __rshift__(self, other):
        still_digits = [0]*other + self.still_digits
        return Alpha.construct(still_digits, self.repeating_digits)

    def __lshift__(self, other):
        i = self
        difference = other - len(self.still_digits)
        if difference > 0:
            i._spew(difference)
        return Alpha.construct(i.still_digits[other:], i.repeating_digits)

    def _multiply_loop(self, i):
        return Alpha.construct(self.still_digits, self.repeating_digits * i)


x = Alpha(pa.pAdic.zero(5))
y = Alpha(pa.pAdic.to_p_adic(5,5))<<1