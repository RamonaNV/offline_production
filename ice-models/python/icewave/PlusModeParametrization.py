# -*- coding: utf-8 -*-
# Copyright (c) 2019
# Ben Jones <ben.jones@uta.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
# $Id$
#
# @file PlusModeParametrization.py
# @version $Revision$
# @date $Date$
# @author Ben Jones, Jakob van Santen

import numpy as np
import copy

from .FourierToolset import FourierSeries, PerturbPhases, PerturbAmplitudes
from icecube import clsim

class PlusModeParametrization:
    """
    An example parametriztion that varies the plus modes of the logarithmic
    abs+scat FFTs.
    """
    def __init__(self, modes_to_shift):
        """
        :param modes_to_shift: indices of the modes to perturb
        """
        self.modes_to_shift = modes_to_shift

    def transform(self, x, frame):
        """
        Transform the random variates `x` into ice layer perturbations
        """
        assert len(x) == 2*len(self.modes_to_shift)

        medium = frame['MediumProperties']
        scattering_coefficients = []
        absorption_coefficients = []
        for i in range(medium.GetLayersNum()):
            # CLSim changes the effective scattering coefficient to the
            # scattering coefficient by dividing by 1/(1-cos(dir))
            # where cos(dir) is 0.9
            scattering_coefficients.append(medium.GetScatteringLength(i).b400)
            # Note: This is just the absorption due to dust. There is another
            # term in PPC that is for the absorption of ice.
            # We'll just use the absorption of dust as an approximation.
            absorption_coefficients.append(medium.GetAbsorptionLength(i).aDust400)

        # The deepest layer represents the bedrock, and its absorption
        # coefficient is set to 999. This messes with the Fourier expansion, so
        # we treat it separately
        bedrock = [scattering_coefficients.pop(0), absorption_coefficients.pop(0)]
        sca = np.asarray(scattering_coefficients)
        abs = np.asarray(absorption_coefficients)

        # calculate central models via log prescription
        central_plus  = 0.5 * np.log10(abs*sca)
        central_minus = 0.5 * np.log10(abs/sca)
        # get central Fourier series
        z = medium.GetLayersZStart() + np.arange(1,medium.GetLayersNum())*medium.GetLayersHeight()
        central_fs_plus  = FourierSeries(z, central_plus)

        amp_shifts, phase_shifts = np.asarray(x).reshape((2,len(self.modes_to_shift)))
        fs_plus = PerturbPhases(PerturbAmplitudes(central_fs_plus, self.modes_to_shift, amp_shifts), self.modes_to_shift, phase_shifts)

        # convert from frequency space and add bedrock
        scattering_coefficients = np.concatenate(([bedrock[0]], 10**(fs_plus[1] - central_minus)))
        absorption_coefficients = np.concatenate(([bedrock[1]], 10**(fs_plus[1] + central_minus)))
        # force top layer to be identical
        scattering_coefficients[-1] = sca[-1]
        absorption_coefficients[-1] = abs[-1]

        # synthesize a new MediumProperties object
        medium = copy.deepcopy(medium)
        for i in range(0,medium.GetLayersNum()):
            oldScat = medium.GetScatteringLength(i)
            oldAbs  = medium.GetAbsorptionLength(i)

            newScat = clsim.I3CLSimFunctionScatLenIceCube(
                alpha = oldScat.alpha,
                b400  = float(scattering_coefficients[i])
            )

            newAbs = clsim.I3CLSimFunctionAbsLenIceCube(
                kappa    = oldAbs.kappa,
                A        = oldAbs.A,
                B        = oldAbs.B,
                D        = oldAbs.D,
                E        = oldAbs.E,
                aDust400    = float(absorption_coefficients[i]),
                deltaTau = oldAbs.deltaTau
            )

            medium.SetScatteringLength(i, newScat)
            medium.SetAbsorptionLength(i, newAbs)

        del frame['MediumProperties']
        frame['MediumProperties'] = medium
