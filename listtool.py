import chemical as cm


def list_flat(lst: list):
    def single_flat(lst: list):
        return sum(lst, [])

    if len(lst) == 0: return []
    if type(lst[0]) == list:
        res = single_flat(lst)
        return list_flat(res)
    return lst


def isum(fs: list[cm.Formula]):
    res = fs[0]
    for i in range(1, len(fs)):
        res += fs[i]
    return res


class Multiplier:
    def __init__(self, obj: cm.Formula, num):
        self.obj = obj
        self.num = num

    def __repr__(self):
        return '(' + str(self.obj) + ')' + str(self.num)

    def eval(self):
        return self.obj * self.num

    def split(self, num: int):
        res = []
        n = len(self.obj)
        for i in range(1, n + 1):
            for j in range(0, i):
                res.append(split_formula_ij(self.obj, i, j, num))
        return res


class Region:
    def __init__(self, lst: list[Multiplier]):
        self.lst = lst

    def __repr__(self):
        return ''.join(list(map(str, self.lst)))

    def split(self, num: int):
        res = []
        for m in self.lst:
            a = m.split(num)
            if a is not None: res.append(a)
        return res

    def eval(self) -> cm.Formula:
        return isum(list(map(lambda x: x.eval(), self.lst)))


def formula_to_multiplier(f: cm.Formula, n: int):
    if len(f) == 0:
        return None
    else:
        return Multiplier(f, n)


# 使用时注意j<i
def split_formula_ij(f: cm.Formula, i: int, j: int, n: int):
    a = formula_to_multiplier(f[:j], 1)
    b = formula_to_multiplier(f[j:i], n)
    c = formula_to_multiplier(f[i:], 1)
    return Region(list(filter(lambda x: x is not None, [a, b, c])))


def formula_with_num_card(f: cm.Formula, ns: list[int]):
    rs = [Region([Multiplier(f, 1)])]
    for num in ns:
        s = []
        for r in rs:
            s.append(r.split(num))
        rs = list_flat(s)
    return rs


def mul_formula(f: cm.Formula, ns: list[int]) -> list[cm.Formula]:
    regions = formula_with_num_card(f, ns)
    res = list(map(lambda x: x.eval(), regions))
    for i in res:
        i.simplify()
    return res
