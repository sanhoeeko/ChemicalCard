import fractions as fr
import re
from collections.abc import Iterable

import numpy as np
import pandas as pd


class Atom:
    def __init__(self, name, num=1):
        self.name = name
        self.num = num

    def __repr__(self):
        return '(' + self.name + ',' + str(self.num) + ')'

    def __mul__(self, x: int):
        return Atom(self.name, self.num * x)

    def __add__(self, other):
        return Atom(self.name, self.num + other.num)

    def __eq__(self, o):
        return self.name == o.name and self.num == o.num

    def __hash__(self):
        return hash(self.name + str(self.num))

    def addable(self, o):
        return self.name == o.name


class Formula:
    def __init__(self, atoms: list[Atom], charge=0):
        self.atoms = atoms
        self.charge = charge

    def __repr__(self):
        s = str(self.atoms)
        return s + ' (' + str(self.charge) + ')' if self.charge != 0 else s

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, item):
        return Formula(self.atoms[item])

    def __mul__(self, x: int):
        return Formula(list(map(lambda y: y * x, self.atoms)))

    def __add__(self, o: 'Formula'):
        return Formula(self.atoms + o.atoms, self.charge + o.charge)

    def __eq__(self, o: 'Formula'):
        return set(self.atoms) == set(o.atoms)

    def simplify(self):
        flag = 0
        n = len(self.atoms)
        for i in range(n):
            for j in range(i):
                if self.atoms[i].addable(self.atoms[j]):
                    self.atoms[i] += self.atoms[j]
                    del self.atoms[j]
                    flag = 1
                    break
            if flag: break
        if flag: self.simplify()

    def elements(self):
        return set(list(map(lambda x: x.name, self.atoms)))

    def get_num_by_element(self, element: str):
        for a in self.atoms:
            if a.name == element:
                return a.num
        return 0


h2o = Formula([Atom('H', 2), Atom('O', 1)])

df = pd.read_csv('CMCtable.csv')


# 解析化学式
def ParseName(name: str):
    lst = []

    def subParseName(name):
        def matchElement(name):
            m = re.match('[A-Z][^A-Z().]*', name)
            return m.group() if m else None

        def parseAtomStr(atom_str):
            element = re.match('[a-zA-Z]*', atom_str).group()
            num_str = atom_str[len(element):]
            num = int(num_str) if num_str else 1
            return Atom(element, num)

        # 警告：这个正则表达式只能匹配一重括号！！
        def matchBracket(name):
            m = re.match(r'\(.*\)', name)
            return m.group() if m else None

        def matchDigit(name):
            m = re.match('[0-9]*', name)
            return m.group() if m else None

        atom_str = matchElement(name)
        if atom_str:
            # 如果是大写字母开头，则匹配元素和个数
            lst.append(parseAtomStr(atom_str))
            subname = name[len(atom_str):]
            if subname: subParseName(subname)
            return
        else:
            # 如果不是大写字母开头而是点号，则匹配结晶水
            if name[0] == '.':
                subname = name[1:]
                digit_str = matchDigit(subname)
                h2onum = int(digit_str) if digit_str else 1
                h2os = h2o * h2onum
                lst.extend(h2os.atoms)
                return  # 结晶水一定在末尾
            else:
                # 如果是括号
                bracket_str = matchBracket(name)
                sub_formula = ParseName(bracket_str[1:-1])
                # 匹配括号后面的数字
                subname = name[len(bracket_str):]
                digit_str = matchDigit(subname)
                bracket_num = int(digit_str)
                sub_formula *= bracket_num
                lst.extend(sub_formula.atoms)
                subname = subname[len(digit_str):]
                if subname: subParseName(subname)
                return

    # 匹配末尾电荷(+,-)
    def parseCharge(name: str):
        if '+' in name:
            idx = name.index('+')
            charge = int(name[idx + 1:])
            name = name[:idx]
            return name, charge
        elif '-' in name:
            idx = name.index('-')
            charge = -int(name[idx + 1:])
            name = name[:idx]
            return name, charge
        else:
            return name, 0

    try:
        name, charge = parseCharge(name)
        subParseName(name)
        res = Formula(lst)
        res.simplify()
        res.charge = charge
        return res
    except:
        return None


class Matter:
    def __init__(self, name: str, formula: Formula, state: str, gibbs: float):
        self.name = name
        self.formula = formula
        self.state = state
        self.gibbs = gibbs

    def __repr__(self):
        return self.name + self.state

    def msg(self):
        return '\n'.join([str(self.name), str(self.formula), str(self.gibbs), ''])

    def elements(self):
        return self.formula.elements()

    def get_num_by_element(self, element: str):
        return self.formula.get_num_by_element(element)


def LoadTable():
    res = []
    for i in range(len(df)):
        formula = ParseName(df['formula'][i])
        if formula and formula.charge == 0:  # 暂不考虑离子
            res.append(Matter(df['formula'][i], formula, df['state'][i], df['G'][i]))
    return res


def nullspace(A: np.ndarray, atol=1e-13, rtol=0) -> np.ndarray:
    A = np.atleast_2d(A)
    u, s, vh = np.linalg.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    ns[abs(ns) < 1e-6] = 0  # 防止后面求秩类操作出错
    return ns


# 配平反应方程式，若无法配平（电荷不守恒），返回None
# 参数element只是为了减少运算优化性能
def GetEquation(end: list[Matter], start: list[Matter], elements: set[str]) -> (np.ndarray, np.ndarray):
    n = len(start)
    elements = list(elements)
    matter_vec = start + end
    # 计算配平矩阵
    M = np.zeros((len(elements), len(matter_vec)), dtype=int)
    for i in range(len(matter_vec)):
        for j in range(len(elements)):
            M[j, i] = matter_vec[i].get_num_by_element(elements[j])
            if i < n: M[j, i] = -M[j, i]  # 这样所有配平系数必须是正值
    # 解方程MX=0，且X[0]>0（方便判断正负）
    vecs = nullspace(M)
    rank = vecs.shape[1]  # 独立解的个数
    for j in range(rank):
        if vecs[0, j] < 0:
            vecs[:, j] = -vecs[:, j]
    if rank == 0:
        return None, None
    elif rank == 1:
        vecs.reshape((-1,))
        # 判断解的正定性
        if not (vecs > 0).all():
            return None, None
        return vecs[:n], vecs[n:]
    else:
        # 挑出含0解作为冗余解
        verbose = []
        proper = []
        vs = list(vecs.T)
        for v in vs:
            if not (v > 0).all():
                verbose.append(v)
            else:
                proper.append(v)
        if len(proper) == 0:
            return None, None
        elif len(proper) == 1:
            vecs = proper[0]
            return vecs[:n], vecs[n:]
        else:
            raise ValueError


def CalGibbs(matters: list[Matter], coefs):
    gibbs_vec = np.array(list(map(lambda x: x.gibbs, matters)))
    return np.dot(gibbs_vec, coefs)


# 判断反应在元素守恒方面是否可行
def GetElements(ms: list[Matter]):
    elements = set()
    for m in ms:
        elements |= m.elements()
    return elements


def IsReactionByElement(end: list[Matter], start: set[str]):
    return GetElements(end) == start


# 判断反应在自由能方面是否可行
def IsReactionByGibbs(end: list[Matter], start: list[Matter], elements: set[str]):
    left, right = GetEquation(end, start, elements)
    if left is None: return False
    g_start = CalGibbs(start, left)
    g_end = CalGibbs(end, right)
    return g_end < g_start


# 判断生成物与反应物是否满足无交集
def IsReactionByName(end: list[Matter], start: list[Matter]):
    return not set(list(map(lambda x: x.name, end))) & set(list(map(lambda x: x.name, start)))


def IsReaction(end: list[Matter], start: list[Matter], elements: set[str]):
    return IsReactionByName(end, start) \
        and IsReactionByElement(end, elements) \
        and IsReactionByGibbs(end, start, elements)


matter_table = LoadTable()


def find_from_table(name: str):
    for m in matter_table:
        if m.name == name:
            return m
    return None


def result1(start: list[Matter]):
    res = []
    e = GetElements(start)
    matter_available = list(filter(lambda x: x.formula.elements().issubset(e), matter_table))
    for i in range(len(matter_available)):
        m1 = matter_available[i]
        if IsReaction([m1], start, e):
            res.append(Reaction(start, [m1]))
    return res


def result2(start: list[Matter]):
    res = []
    e = GetElements(start)
    matter_available = list(filter(lambda x: x.formula.elements().issubset(e), matter_table))
    for i in range(len(matter_available)):
        m1 = matter_available[i]
        for j in range(i):
            m2 = matter_available[j]
            if IsReaction([m1, m2], start, e):
                res.append(Reaction(start, [m1, m2]))
    return res


def result3(start: list[Matter]):
    res = []
    e = GetElements(start)
    matter_available = list(filter(lambda x: x.formula.elements().issubset(e), matter_table))
    for i in range(len(matter_available)):
        m1 = matter_available[i]
        for j in range(i):
            m2 = matter_available[j]
            for k in range(j):
                m3 = matter_available[k]
                if IsReaction([m1, m2, m3], start, e):
                    res.append(Reaction(start, [m1, m2, m3]))
    return res


class FrArray:  # 辅助类
    def __init__(self, arr: np.ndarray):
        def value(x):
            if isinstance(x, Iterable):
                return x[0]
            return x

        self.arr = list(
            map(lambda x: fr.Fraction(value(x)).limit_denominator(), list(arr))
        )

    def __getitem__(self, item):
        return self.arr[item]

    def __repr__(self):
        return str(self.arr)

    def float(self):
        return list(map(float, self.arr))

    def adequate(self):
        for f in self.arr:
            if f.denominator > 100:
                return False
        return True


class Reaction:
    def __init__(self, start: list[Matter], end: list[Matter]):
        self.start = start
        self.end = end
        self.left, self.right = GetEquation(end, start, GetElements(start))
        coef = 1 / self.left[0]
        self.left = FrArray(self.left * coef)  # 约束第一反应物的量为1mol
        self.right = FrArray(self.right * coef)
        self.dG = CalGibbs(self.end, self.right.float()) - CalGibbs(self.start, self.left.float())

    def __repr__(self):
        res = ''
        n = len(self.start)
        for i in range(n):
            res += str(self.left[i]) + ' ' + str(self.start[i])
            if i != n - 1: res += ' + '
        res += ' -> '
        n = len(self.end)
        for i in range(n):
            res += str(self.right[i]) + ' ' + str(self.end[i])
            if i != n - 1: res += ' + '
        res += ' dG='
        res += str('{:.1f}'.format(self.dG))
        return res

    def __lt__(self, o):
        return self.dG < o.dG

    def __gt__(self, o):
        return self.dG > o.dG

    def adequate(self):
        return self.right.adequate() and self.left.adequate()


# 判断化学式是否存在于数据库中
def isMatter(f: Formula):
    res = []
    for m in matter_table:
        if m.formula == f:
            res.append(m)
    return res


# 获取物质a与物质b所有可能的反应（最多3种产物），若没有，返回空列表
def find_reaction(a: Matter, b: Matter) -> list[Reaction]:
    start = [a, b]
    res = []
    res.extend(result1(start))
    res.extend(result2(start))
    res.extend(result3(start))
    # 删除系数不合理的反应式，并按自由能排序
    res = list(filter(lambda x: x.adequate(), res))
    res.sort()
    return res


def get_possible_reactions(a: Formula, bs: list[Formula]) -> list[Reaction]:
    res = []
    for b in bs:
        res.extend(find_reaction(a, b))
    return res
