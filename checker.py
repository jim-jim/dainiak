def write(x, m):
    t = type(x)
    if t is list:
        for i in x: write(i, m)
    elif t is dict:
        for k, v in x.items(): write(k, m); print(v)
    else:
        print(''.join(('1' if x & (1 << i) else '0') for i in range(m - 1, -1, -1)))
def bitc64(x):
    x -= ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F
    x = x + (x >> 8)
    x = x + (x >> 16)
    return (x + (x >> 32)) & 0x0000007F
def bitcL(x):
    r = 0
    while x > 0:
        r += bitc64(x & 0xFFFFFFFFFFFFFFFF)
        x >>= 64
    return r
def shift(x, m, sh):
    sh = (sh%m+m)%m
    mm = (1 << sh) - 1
    return ((x & mm) << (m-sh)) | (x >> sh)
def str_to_code(s):
    x = 0
    i = pos = 0
    for ch in s:
        ch = int(ch)
        if i % 2:
            x |= ((1 << ch) - 1) << pos
        pos += ch
        i += 1
    return (x, pos)
def get_ac(x, n):
    ac = []
    for i in range(1, n):
        y = shift(x, n, i)
        ac.append(int(bitcL(y^x) - n//2))
    return ac
"""
21143212114213111142
123442131121311112
2411411331122112
2113311411422112 //
1344112122311112
1223411421311112
1321124213411112
21133112134112
24411212311112
121233411112
134312211121
123421311112
221124411112
242211411112
211422411112
2431211112
1223411112
1342211112 //
13321121
23113112
23123111
12331121 //

"""
s = "211113431312412122111431"
x, n = str_to_code(s)
write(x, n)
print(get_ac(x, n))
print(n)