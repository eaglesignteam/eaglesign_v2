from MSIS_security import MSIS_summarize_attacks, MSISParameterSet

class UniformEagleSignParameterSet(object):
    def __init__(self, n, l, eta_f, tau, td, tg, k, q, gamma1, gamma2):
        self.n = n
        self.k = k
        self.l = l
        self.eta_f = eta_f
        self.tau = tau
        self.td = td
        self.tg = tg
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.q = q

        self.beta = l*eta_f*tau
        self.alpha = max(2*tg*(self.gamma1 - self.beta), 2*self.gamma2)


UnifEagleSignMedium = UniformEagleSignParameterSet(
    1024, 1, 1, 16, 7, 7, 1, 2021377, 2**14, 2**17)
UnifEagleSignVeryHigh_I = UniformEagleSignParameterSet(
    2048, 1, 1, 30, 13, 14, 1 , 33292289, 2**16, 2**20)
UnifEagleSignVeryHigh_III = UniformEagleSignParameterSet(
    1024, 2, 1, 16, 3, 16, 3 , 7340033, 2**16, 2**20)


all_params_unif = [
    ("Uniform EagleSign Medium", UnifEagleSignMedium),
    ("Uniform EagleSign Very High I", UnifEagleSignVeryHigh_I),
    ("Uniform EagleSign Very High III", UnifEagleSignVeryHigh_III),
]


def EagleSign_to_MSIS(dps):
    return MSISParameterSet(dps.n, dps.k + dps.l + dps.l, dps.k, dps.alpha, dps.q, norm="linf")


if __name__ == "__main__":
    for (scheme, param) in all_params_unif:
        print("\n"+scheme)
        print(param.__dict__)

        print("")
        print("=== STRONG UF")

        MSIS_summarize_attacks(EagleSign_to_MSIS(param))