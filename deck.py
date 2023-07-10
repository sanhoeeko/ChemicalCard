import copy as cp
import random as rd

import pandas as pd

import chemical as cm
import listtool as lt
import myio

console = myio.OutBuffer()


def init_prob(series: pd.Series):
    Z = series.sum()
    x = series / Z
    y = cp.deepcopy(x)
    s = 0
    for i in range(len(x)):
        s += x[i]
        y[i] = s
    return y


deckdf = pd.read_csv('cards.csv')
deckdf['r'] = init_prob(deckdf['weight'])


class Card:
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name


class ElementCard(Card):
    def __init__(self, name):
        super(ElementCard, self).__init__(name)
        self.atom = cm.Atom(name)


class NumCard(Card):
    def __init__(self, name):
        super(NumCard, self).__init__(name)
        self.num = int(name)


class SpecialCard(Card):
    def __init__(self, name):
        super(SpecialCard, self).__init__(name)

    def __le__(self, other):
        pass

    def __ge__(self, other):
        pass

    def act(self):
        pass


special_name = []


def getcard(name: str) -> Card:
    if name.isdigit():
        return NumCard(name)
    elif name in special_name:
        return SpecialCard(name)
    else:
        return ElementCard(name)


# 考虑所有数字牌的使用情况，返回所有可能的化学式，如果没有则返回空列表
# 注意数字牌能嵌套生效，但作用域不能交叉
def cards_to_matter(cards: list[Card]) -> list[cm.Matter]:
    def divide(cards: list[Card]):
        nums = list(filter(lambda x: type(x) == NumCard, cards))
        eles = list(filter(lambda x: type(x) == ElementCard, cards))
        return nums, eles

    nums, eles = divide(cards)
    nums = list(map(lambda x: x.num, nums))
    f = cm.Formula(list(map(lambda x: x.atom, eles)))
    fs = lt.mul_formula(f, nums)  # 仍有问题，如H2O2不能拼出来
    # 挑选出成立的化学式，并匹配物质
    res = []
    for f in fs:
        m = cm.isMatter(f)
        if m: res.extend(m)
    return res


class Hand:
    def __init__(self):
        self.cards = []

    def __repr__(self):
        return str(self.cards)

    def draw_card(self):
        r = rd.uniform(0, 1)
        for i in range(len(deckdf)):
            if deckdf['r'][i] > r:
                self.cards.append(getcard(deckdf['name'][i]))
                return
        self.cards.append(getcard(deckdf['name'][len(deckdf) - 1]))
        return

    # 开局按概率表随机发牌
    def start(self, num: int):
        for i in range(num): self.draw_card()

    # 输入牌的位置，获取出牌牌型
    def copy_push(self, cards: list[int]) -> list[Card]:
        return list(map(lambda idx: self.cards[idx], cards))

    # 输入牌的位置，丢弃这些牌
    def throw(self, cards: list[int]):
        for i in range(len(self.cards))[::-1]:
            if i in cards:
                del self.cards[i]


class Player:
    def __init__(self, name):
        self.name = name
        self.hand = Hand()
        self.draw_this_turn = False

    def start(self):
        self.hand.start(12)

    # 回合开始
    def turn(self):
        self.draw_this_turn = False
        self.show()

    # 展示手牌
    def show(self):
        console.print('现在轮到', self.name)
        console.print(*self.hand.cards)

    # 主动摸牌
    def draw_card(self):
        if not self.draw_this_turn:
            self.hand.draw_card()
            self.draw_this_turn = True
            return True
        else:
            console.print('这回合你已抽过牌了！')
            return False

    # 用数字牌和特殊牌换其他牌
    def exchange(self, idx: int):
        t = type(self.hand.cards[idx])
        if t == NumCard or t == SpecialCard:
            self.hand.draw_card()
            self.hand.throw([idx])
            return True
        else:
            return False

    # 放过
    def pass_out(self):
        pass

    # 出牌
    def push(self, cards: list[int]) -> list[cm.Matter]:
        def divide(cards: list[Card]):
            special = list(filter(lambda x: type(x) == SpecialCard, cards))
            others = list(filter(lambda x: type(x) != SpecialCard, cards))
            return special, others

        fs = self.hand.copy_push(cards)
        # 分离特殊牌，特殊牌先判定
        special, others = divide(fs)
        if special:
            special.sort()  # 按判定顺序排序
            for spcard in special:
                spcard.act()
        # 元素牌和数字牌生效，放出所有可能的物质
        ms = cards_to_matter(others)
        return ms


class Game:
    def __init__(self, players: list[Player]):
        self.players = players
        self.field = []  # 场上物质
        self.winside = -1
        self.side = 0  # 当前玩家
        self.state = myio.GameCache()

    def start(self):  # 游戏初始化
        # 发牌
        for p in self.players:
            p.start()
        self.players[self.side].turn()

    def next(self):
        self.side += 1
        if self.side == len(self.players): self.side = 0
        self.players[self.side].turn()
        console.print('场上的物质：', self.field)

    def show(self):
        self.players[self.side].show()
        console.print('场上的物质：', self.field)

    def push_success(self):
        self.players[self.side].hand.throw(self.state.keep['cards_cache'])

    def main_turn(self, cmd: str) -> bool:
        if self.winside == -1:
            p = self.players[self.side]
            # 读取玩家的操作
            args = cmd.split(' ')

            if self.state.state is None:
                if args[0] == 'draw':
                    p.draw_card()
                elif args[0] == 'pass':
                    p.pass_out()
                    self.field = []
                    return True  # 结束回合
                elif args[0] == 'exchange':
                    p.exchange(int(args[1]))
                elif args[0] == 'push':
                    nums = list(map(int, args[1:]))
                    fs = p.push(nums)
                    # 如果没有备选物质，则出牌失败
                    if not fs: return False
                    self.state.toState('PUSHING_CHOOSE_MATTER')  # 选择化学式的解释
                    self.state.vars['possible_formulas'] = fs
                    self.state.keep['cards_cache'] = nums
                    console.printlist(fs)
            # 选择化学式的解释
            elif self.state.state == 'PUSHING_CHOOSE_MATTER':
                formula = self.state['possible_formulas'][int(args[0])]
                # 如果场上没有物质，则出牌成功
                if not self.field:
                    self.field = [formula]
                    self.push_success()
                    self.state.toState(None)
                    return True
                else:
                    reactions = cm.get_possible_reactions(formula, self.field)
                    # 如果没有备选反应，则出牌失败，并复位
                    if not reactions:
                        self.state.toState(None)
                        return False
                    self.state.toState('PUSHING_CHOOSE_REACTION')
                    self.state['possible_reactions'] = reactions
                    console.printlist(reactions)
            # 选择反应方程式
            elif self.state.state == 'PUSHING_CHOOSE_REACTION':
                result = self.state['possible_reactions'][int(args[0])]
                self.field = result.end
                self.push_success()
                self.state.toState(None)
                return True  # 成功出牌，自动结束回合
            return False  # 未结束回合
        else:
            raise BaseException

    def win(self):
        pass


"""
if __name__ == '__main__':
    h = getcard('H')
    o = getcard('O')
    n2 = getcard('2')
    fs = cards_to_formula([h, n2, o])
    print(fs)
"""
