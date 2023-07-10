from deck import *

# only main

p1 = Player('Alice')
p2 = Player('Bob')
g = Game([p1, p2])

g.start()
while g.winside == -1:
    cmd = input('>>> ')
    sign = g.main_turn(cmd)
    if sign:
        g.next()
    else:
        if g.state is None:
            g.show()
        else:
            console.print(g.state.state)

g.win()
