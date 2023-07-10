class OutBuffer:
    def __init__(self):
        self.buf = ''

    def print(self, *x):
        self.buf = list(map(str, x))
        print(' '.join(self.buf))

    def printlist(self, lst):
        self.buf = lst
        print('\n'.join(list(map(str, lst))))


class GameCache:
    def __init__(self):
        self.state = None
        self.vars = {}
        self.keep = {}

    def __getitem__(self, item):
        return self.vars[item]

    def __setitem__(self, key, value):
        self.vars[key] = value

    def get(self, state, item):
        if self.state == state:
            return self[item]
        else:
            raise ValueError

    def toState(self, state):
        self.state = state
        self.vars = {}
        if state is None:
            self.keep = {}
