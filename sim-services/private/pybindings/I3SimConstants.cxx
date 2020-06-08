//
//   Copyright (c) 2004, 2005, 2006, 2007   Troy D. Straszheim  
//   
//   This file is part of IceTray.
//
//   IceTray is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <boost/preprocessor.hpp>
#include <vector>

#include <sim-services/I3SimConstants.h>

using namespace boost::python;

#define I3_SIM_CONSTANT_NAMES			\
  (m_e)(m_mu)(m_tau)(m_p)(m_n)		\
  (m_W)(m_Z)(h)(hbar)(hbarc)      \
  (G_Fermi)(epsilon_0)(mu_0)      \
  (alpha)(k_B)(sinsq_theta_W)     \
  (cos_theta_W)(sin_theta_C)      \
  (cos_theta_C)                   \
  

#define I3SIMCONSTANT_DEF(r,data,t) \
  .def_readonly(BOOST_PP_STRINGIZE(t), BOOST_PP_CAT(I3SimConstants::, t))

// dummy class just used as a namespace.
namespace {
  struct dummy { };
}

void register_I3SimConstants()
{
  class_<dummy>("I3SimConstants")
    .def_readonly("m_e", I3SimConstants::m_e, "electron mass")
    .def_readonly("m_mu", I3SimConstants::m_mu, "muon mass")
    .def_readonly("m_tau", I3SimConstants::m_tau, "tau mass")
    .def_readonly("m_p", I3SimConstants::m_p, "proton mass")
    .def_readonly("m_n", I3SimConstants::m_n, "neutron mass")
    .def_readonly("m_W", I3SimConstants::m_W, "W-boson mass")
    .def_readonly("m_Z", I3SimConstants::m_Z, "Z-boson mass")
    .def_readonly("h", I3SimConstants::h, "Planck Constant")
    .def_readonly("hbar", I3SimConstants::hbar, "reduced Planck Constant")
    .def_readonly("hbarc", I3SimConstants::hbarc, "Unit conversion")
    .def_readonly("G_Fermi", I3SimConstants::G_Fermi, "Fermi constant, / (hbar * c)^3 is assumed")
    .def_readonly("epsilon_0", I3SimConstants::epsilon_0, "Permittivity of free space")
    .def_readonly("mu_0", I3SimConstants::mu_0, "Permeability of free space")
    .def_readonly("alpha", I3SimConstants::alpha, "Fine-structure constant")
    .def_readonly("k_B", I3SimConstants::k_B, "Boltzmann Constant")
    .def_readonly("sinsq_theta_W", I3SimConstants::sinsq_theta_W, "Weak-mixing Angle sin^{2} theta (M_{Z})")
    .def_readonly("cos_theta_W", I3SimConstants::cos_theta_W, "Weak-mixing Angle cos theta_{W}")
    .def_readonly("sin_theta_C", I3SimConstants::sin_theta_C, "Cabibbo Angle sin_theta_C")
    .def_readonly("cos_theta_C", I3SimConstants::cos_theta_C, "Cabibbo Angle cos_theta_C")
    .def( freeze() )
    ;
}

void register_ShowerParameters()
{
    class_<I3SimConstants::ShowerParameters, boost::shared_ptr<I3SimConstants::ShowerParameters> >(
        "ShowerParameters",
        init<I3Particle::ParticleType, double, double>(
            (args("type"),
            args("energy"),
            args("density")=0.9216*(I3Units::g/I3Units::cm3))
        ))
        .def_readonly("a", &I3SimConstants::ShowerParameters::a)
        .def_readonly("b", &I3SimConstants::ShowerParameters::b)
        .def_readonly("emScale", &I3SimConstants::ShowerParameters::emScale)
        .def_readonly("emScaleSigma", &I3SimConstants::ShowerParameters::emScaleSigma)
    ;
}
