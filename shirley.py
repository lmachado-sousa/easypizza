import numpy as np

def shirley(y,initback, endback, reverse = False):
    background = np.zeros(y.size)
    background0 = np.zeros(y.size)

    if reverse == True:
        y = y[::-1]

    a = y[initback]
    b = y[endback]

    background0.fill(a)

    for nint in xrange(0,6):
        for k2 in xrange(endback, initback, -1):
            sum1 = 0
            sum2 = 0
            for k in xrange(endback, k2, -1):
                sum1 = sum1 + y[k] - background0[k]
            for k in xrange(endback, initback, -1):
                sum2 = sum2 + y[k] - background0[k]
            background[k2] = ((a-b)*sum1)/sum2+b
        background0 = background

    y2 = y-background0

    return y2
