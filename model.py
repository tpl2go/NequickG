import numpy as np
import matplotlib.pyplot as plt

from parameters import main


def electron_density(h, A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot, NmF2, H0):
    if h <= hmF2:  # bottomside
        # print "bot"
        return bottomside_electrondensity(h, A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot)
    else:
        return topside_electrondensity(h, NmF2, hmF2, H0)
        # return 0


def bottomside_electrondensity(h, A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot):
    """

    :param h:
    :param A1:
    :param A2:
    :param A3:
    :param hmF2:
    :param hmF1:
    :param hmE:
    :param B2bot:
    :param B1top:
    :param B1bot:
    :param BEtop:
    :param Bebot:
    :return:
    """
    #where are the 3 semi epstein layers?
    if h > hmE:
        BE = BEtop
    else:
        BE = BEbot

    if h > hmF1:
        BF1 = B1top
    else:
        BF1 = B1bot

    # Compute the exponential arguments for each layer
    if h < 100:
        h = 100
    alpha1 = (h - hmF2) / B2bot
    alpha2 = (h - hmF1) / BF1 * np.exp(10 / (1 + abs(h - hmF2)))
    alpha3 = (h - hmE) / BE * np.exp(10 / (1 + abs(h - hmF2)))
    alpha = np.array([alpha1, alpha2, alpha3])
    # print alpha
    A = np.array([A1, A2, A3])
    # print A

    S = np.empty(3)
    mask1 = np.abs(alpha) > 25
    mask2 = np.logical_not(mask1)
    S[mask1] = 0
    S[mask2] = A[mask2] * np.exp(alpha[mask2]) / np.power(1 + np.exp(alpha[mask2]), 2)
    # print S
    if h >= 100:
        N = (np.sum(S)) * 10 ** 11
    else: # compute corrective terms
        B = np.array([B2bot, BF1, BE])
        ds = np.empty(3)
        ds[mask1] = 0
        ds[mask2] = (1 - np.exp(alpha[mask2])) / (B[mask2] * (1 + np.exp(alpha[mask2])))

        # Chapman parameters
        BC = 1 - 10 * np.sum(S * ds) / np.sum(S)
        z = (h - 100) / 10.0

        # electron density
        N = np.sum(S) * np.exp(1 - BC * z - np.exp(-z)) * 10 ** 11

    return N


def topside_electrondensity(h, NmF2, hmF2, H0):
    """

    :param h:
    :param NmF2:
    :param hmF2:
    :param H0:
    :return:
    """
    g = 0.125
    r = 100

    deltah = h - hmF2
    z = deltah / (H0 * (1 + r * g * deltah / (r * H0 + g * deltah)))

    ea = np.exp(z)

    if ea > 10 ** 11:
        N = 4 * NmF2 / ea * 10 ** 11
    else:
        N = 4 * NmF2 * ea * 10 ** 11 / (1 + ea) ** 2

    return N


def run():
    bottomside_para, topside_para, Azr = main()
    [hmE, hmF1, hmF2, BEtop, BEbot, B1top, B1bot, B2bot, A1, A2, A3] = bottomside_para
    [NmF2, hmF2, H0] = topside_para
    h_array = np.arange(100, 1000)
    N = []
    for h in h_array:
        N.append(electron_density(h, A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot, NmF2, H0))
    plt.plot(h_array, N)
    plt.xlabel("height(km)")
    plt.ylabel("Electron Density (m^-3)")
    plt.title("Nequick-G:\nAzr=" + str(int(Azr)) + "Oct 12UT 30N 0E")
    plt.grid()
    plt.savefig("20Azr Oct 12UT 30N 0E")


def run2():
    [hmE, hmF1, hmF2, BEtop, BEbot, B1top, B1bot, B2bot, A1, A2, A3] = [120, 188.66777907194805, 257.33555814389609,
                                                                        34.333889535974023, 5.0, 20.600333721584413,
                                                                        34.333889535974023, 25.855305987785592,
                                                                        27.662095644483411, 2.288065199190279,
                                                                        4.0414666597141737]
    [NmF2, hmF2, H0] = [6.9155239111208529, 257.33555814389609, 72.319823726235867]
    print electron_density(253, A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot, NmF2, H0)
    print electron_density(254, A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot, NmF2, H0)


# run2()
run()
