

def get_default_perturbation():
    """
    Get the default ice model variations imported from https://github.com/UTA-REST/multisim
    :returns: a tuple (parametrization, distribution) for use with icecube.snowstorm.Perturber
    """
    import numpy as np
    from icecube.dataclasses import I3Matrix
    from icecube.snowstorm import MultivariateNormal
    from .PlusModeParametrization import PlusModeParametrization

    amp_sigmas = np.asarray([0.00500100, 0.03900780, 0.04500900, 0.17903581, 0.07101420, 0.30306061, 0.14502901, 0.09501900, 0.16103221, 0.13302661, 0.15703141, 0.13302661])
    phase_sigmas = np.asarray([0.00000001, 0.01664937, 0.02708014, 0.43171273, 0.02351273, 2.33565571, 0.16767628, 0.05414841, 0.31355088, 0.04227052, 0.27955606, 4.02237848])
    modes_to_shift = np.arange(12)
    variance = np.concatenate((amp_sigmas,phase_sigmas))**2

    return PlusModeParametrization(modes_to_shift), MultivariateNormal(I3Matrix(np.diag(variance)), [0.]*variance.size)
