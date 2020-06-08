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
# @file FourierToolset.py
# @version $Revision$
# @date $Date$
# @author Ben Jones

import numpy as np
import copy
pi = np.pi

def FourierSeries(x, y, nmodes=-1):
    IsAscending = True
    for i in range(len(x)-1):
        if not x[i] < x[i+1]:
            IsAscending= False  
            break
    if not IsAscending:
        x = x[::-1]
        y = y[::-1]
    L = np.abs(max(x) - min(x))
    a = []
    b = []
    if nmodes == -1 :
        nmodes = int(len(x)/2)
    for n in range(nmodes):
        if n == 0:
            a.append((2/L)*np.trapz(y,x))
            b.append(0)
        else:
            a.append((2/L) * np.trapz(y*np.cos(2*n*pi*x/L),x))
            b.append((2/L) * np.trapz(y*np.sin(2*n*pi*x/L),x))
    a = np.asarray(a)
    b = np.asarray(b)
    A = np.sqrt(a**2 + b**2)
    A[0] = a[0]
    phi = np.arctan2(a,b)    
    f = (A[0]/2)*np.ones(len(x))
    if nmodes == 1 :
        return [x, f, A, phi]
    else:
        for n in range(1,nmodes):
                f = f + A[n] * np.sin((2*n*pi*x/L) + phi[n])
    if not IsAscending:
        x = x[::-1]
        f = f[::-1]
    return [x, f, A, phi]

def PerturbAmplitudes(fs, modes_to_shift, shifts_to_apply, RelativeShift=True):
    x = copy.copy(fs[0])
    L = np.abs(max(x) - min(x))
    A = copy.copy(fs[2])
    phi = copy.copy(fs[3])
    for i in range(len(modes_to_shift)):
        if RelativeShift:
            A[modes_to_shift[i]] = A[modes_to_shift[i]] * (1 + shifts_to_apply[i])
        else:
            A[modes_to_shift[i]] = A[modes_to_shift[i]] + shifts_to_apply[i]    
    f = (A[0]/2)*np.ones(len(x))
    for n in range(1,len(A)):
        f = f + A[n] * np.sin((2*n*pi*x/L) + phi[n])
    return [x, f, A, phi]
    
def PerturbPhases(_fs, modes_to_shift, shifts_to_apply, RelativeShift=False):
    x = copy.copy(_fs[0])
    L = np.abs(max(x) - min(x))
    A = copy.copy(_fs[2])
    phi = copy.copy(_fs[3])
    for i in range(len(modes_to_shift)):
        if RelativeShift:
            phi[modes_to_shift[i]] = phi[modes_to_shift[i]] * (1 + shifts_to_apply[i])
        else:
            phi[modes_to_shift[i]] = phi[modes_to_shift[i]] + shifts_to_apply[i]    
    f = (A[0]/2)*np.ones(len(x))
    for n in range(1,len(A)):
        f = f + A[n] * np.sin((2*n*pi*x/L) + phi[n])
    return [x, f, A, phi]

def Wipe(_var_list):
    for var in _var_list:
        del var
    del _var_list
    return None

###### Simply Loads the Ice Data #####

def LoadIcedata(_icepath):
    _output = None    
    _ice_dep, _ice_abs, _ice_sca, _ice_unk = np.loadtxt(_icepath, unpack=True, usecols=[0,1,2,3])
    _output = [_ice_dep, _ice_abs, _ice_sca, _ice_unk]
    Wipe([ _ice_dep, _ice_abs, _ice_sca, _ice_unk])
    return _output

##########  Loads the ice model given a central model. Additional    ##########
##########  functionality allows for perturbations of the Fourier    ##########

###### Default pertubations are M_plus with relative amplitude shifts and 

def GetIcemodel(ice_data, 
              amp_modes_to_shift = [], 
              phs_modes_to_shift = [],
              amp_shifts = [],
              phs_shifts = [],                         
              Perturb_M_plus = True,
              Relative_Amp_Shifts = True,
              Relative_Phs_Shifts = False):
    _output = None
    if len(amp_modes_to_shift) != len(amp_shifts) or len(phs_modes_to_shift) != len(phs_shifts):
        raise ValueError('Number Of Modes To Shift Does Not Equal Number Of Shifts!!!')
    
    _tru_dep = copy.deepcopy(ice_data[0])
    _tru_sca = copy.deepcopy(ice_data[1])
    _tru_abs = copy.deepcopy(ice_data[2])
    _tru_unk = copy.deepcopy(ice_data[3])
 
    _M_plus = (np.log10(_tru_abs[:-1]*_tru_sca[:-1]))/2
    _M_minus = (np.log10(_tru_abs[:-1]/_tru_sca[:-1]))/2
    
    if Perturb_M_plus:
        _fs  = FourierSeries(_tru_dep[:-1], _M_plus)
    else: 
        _fs  = FourierSeries(_tru_dep[:-1], _M_minus)
    
    _pfs = PerturbAmplitudes( _fs,  
                              amp_modes_to_shift,  
                              amp_shifts,  
                              RelativeShift = Relative_Amp_Shifts)
    _pfs = PerturbPhases(    _pfs,  
                              phs_modes_to_shift,  
                              phs_shifts,  
                              RelativeShift = Relative_Phs_Shifts)       
    if Perturb_M_plus:
        _pfs_abs    = 10**( _pfs[1] + _M_minus)
        _pfs_sca    = 10**( _pfs[1] - _M_minus)
    else: 
        _pfs_abs    = 10**( _M_plus + _pfs[1])
        _pfs_sca    = 10**( _M_plus - _pfs[1])

    _pfs_sca[0]  = _tru_sca[0]
    _pfs_abs[0]  = _tru_abs[0]
    _pfs_sca[-1] = _tru_sca[-2]
    _pfs_abs[-1] = _tru_abs[-2]
    _pfs_sca  = np.append(  _pfs_sca, _tru_sca[-1])    
    _pfs_abs  = np.append(  _pfs_abs, _tru_abs[-1])    
    
    _output = [_tru_dep, _pfs_sca, _pfs_abs, _tru_unk]
    
    Wipe([_tru_dep, _tru_sca, _tru_abs, _M_plus, _M_minus, _fs, _pfs, _pfs_abs, _pfs_sca])
    return _output

###### Returns Amplitude and Phase Parameters of the Ice Decomposition  ###############

def GetFourierParameters(ice_data,
              amp_modes_to_shift = [],
              phs_modes_to_shift = [],
              amp_shifts = [],
              phs_shifts = [],
              Perturb_M_plus = True,
              Relative_Amp_Shifts = True,
              Relative_Phs_Shifts = False):

    _output = None
    if len(amp_modes_to_shift) != len(amp_shifts) or len(phs_modes_to_shift) != len(phs_shifts):
        raise ValueError('Number Of Modes To Shift Does Not Equal Number Of Shifts!!!')

    _tru_dep = copy.deepcopy(ice_data[0])
    _tru_sca = copy.deepcopy(ice_data[1])
    _tru_abs = copy.deepcopy(ice_data[2])
    _tru_unk = copy.deepcopy(ice_data[3])

    _M_plus = (np.log10(_tru_abs[:-1]*_tru_sca[:-1]))/2
    _M_minus = (np.log10(_tru_abs[:-1]/_tru_sca[:-1]))/2

    if Perturb_M_plus:
        _fs  = FourierSeries(_tru_dep[:-1], _M_plus)
    else:
        _fs  = FourierSeries(_tru_dep[:-1], _M_minus)

    _pfs = PerturbAmplitudes( _fs,
                              amp_modes_to_shift,
                              amp_shifts,
                              RelativeShift = Relative_Amp_Shifts)
    _pfs = PerturbPhases(    _pfs,
                              phs_modes_to_shift,
                              phs_shifts,
                              RelativeShift = Relative_Phs_Shifts)
 
    _Amps = _pfs[2]
    _Phss = _pfs[3]

    _output = [_Amps, _Phss]

    Wipe([_tru_dep, _tru_sca, _tru_abs, _M_plus, _M_minus, _fs, _pfs, _Amps, _Phss])
    return _output


########## Load IceModel With Uncertainty Values, Note RelativeValues toggle #######

def LoadUncertainties(_error_path, _ice_data, RelativeValues = False):
    
    _output = None
    
    _err_dep, _err_sca, _err_abs, _bia_sca, _bia_abs  = np.loadtxt( _error_path, unpack=True, 
                                                               usecols=[0,4,5,6,7])
    _tru_dep = copy.deepcopy(_ice_data[0])
    _tru_sca = copy.deepcopy(_ice_data[1])
    _tru_abs = copy.deepcopy(_ice_data[2])
    
    _abs_err = []
    _sca_err = []
    _abs_bia = []
    _sca_bia = []
    
    for depth in _tru_dep:
        if depth in _err_dep:
            i = np.where(_err_dep == depth)[0][0]
            _sca_err.append(_err_sca[i])
            _abs_err.append(_err_abs[i])
            _sca_bia.append(_bia_sca[i])
            _abs_bia.append(_bia_abs[i])
        else:
            _sca_err.append(0)
            _abs_err.append(0)
            _sca_bia.append(0)
            _abs_bia.append(0)

    _abs_err = np.asarray(_abs_err)
    _sca_err = np.asarray(_sca_err)
    _abs_bia = np.asarray(_abs_bia)
    _sca_bia = np.asarray(_sca_bia)

    if not RelativeValues:
        _sca_bia = _tru_sca + ( _tru_sca * _sca_bia)
        _abs_bia = _tru_abs + ( _tru_abs * _abs_bia)
        _sca_err  = _tru_sca  *  _sca_err
        _abs_err  = _tru_abs  *  _abs_err
        
        _output = [_tru_dep, _sca_bia, _sca_err, _abs_bia, _abs_err]
    
    else:
        _sca_bia = np.ones(len(_sca_bia)) + _sca_bia
        _abs_bia = np.ones(len(_abs_bia)) + _abs_bia
        _sca_err  = _sca_err
        _abs_err  = _abs_err
        
        _output = [_tru_dep, _sca_bia, _sca_err, _abs_bia, _abs_err]   

    Wipe([_tru_dep, _tru_sca, _tru_abs, _sca_err, _sca_bia, _abs_err, _abs_bia])
    return _output
 
def MyChi2(_test_model, 
           _error_model, 
           central_model = [], 
           RelativeErrors=False,
           Reduced=True):

    _output = None

    _dep_testcheck =   len(_test_model[0])==len(_error_model[0])
    _sca_testcheck = ( len(_test_model[1])==len(_error_model[1]) ) and ( len(_test_model[1])==len(_error_model[2]) ) 
    _abs_testcheck = ( len(_test_model[2])==len(_error_model[4]) ) and ( len(_test_model[2])==len(_error_model[3]) ) 

    _dep_check = True
    _sca_check = True
    _abs_check = True


    if not (_dep_testcheck and _sca_testcheck and _abs_testcheck):
        raise ValueError('Hypothesis Model and Error Model Have Incompatible Dimensions!')

    _hyp_dep = copy.deepcopy(_test_model[0])
    _hyp_sca = copy.deepcopy(_test_model[1])
    _hyp_abs = copy.deepcopy(_test_model[2])
        
    if RelativeErrors:
        
        _dep_check =   len(central_model[0])==len(_error_model[0])
        _sca_check = ( len(central_model[1])==len(_error_model[1]) ) and ( len(central_model[1])==len(_error_model[2]) ) 
        _abs_check = ( len(central_model[2])==len(_error_model[4]) ) and ( len(central_model[2])==len(_error_model[3]) ) 

        if not (_dep_check and _sca_check and _abs_check):
            raise ValueError('Central Model and Relative Error Model Have Incompatible Dimensions!')
        
        _tru_dep = copy.deepcopy(_error_model[0])
        _tru_sca = copy.deepcopy(central_model[1]) * copy.deepcopy(_error_model[1])
        _err_sca = copy.deepcopy(central_model[1]) * copy.deepcopy(_error_model[2])
        _tru_abs = copy.deepcopy(central_model[2]) * copy.deepcopy(_error_model[3])
        _err_abs = copy.deepcopy(central_model[2]) * copy.deepcopy(_error_model[4])  
    else:
        _tru_dep = copy.deepcopy(_error_model[0])
        _tru_sca = copy.deepcopy(_error_model[1])
        _err_sca = copy.deepcopy(_error_model[2])
        _tru_abs = copy.deepcopy(_error_model[3])
        _err_abs = copy.deepcopy(_error_model[4])
        
    _chi2 = 0
    _nlayers = 0
    for i in range(len(_tru_dep)):
        if _hyp_dep[i] == _tru_dep[i]:
            if _err_abs[i] == 0 or _err_sca[i] == 0:
                _chi2 = _chi2
            else:
                _chi2 = ( _chi2 + 0.5*( ( _hyp_sca[i] - _tru_sca[i])/_err_sca[i] )**2 
                              + 0.5*( ( _hyp_abs[i] - _tru_abs[i])/_err_abs[i] )**2 )
            _nlayers = _nlayers + 1
        else:
            raise ValueError('Incompatible Layers in Hypothesis and Error Models!')
            
    if _nlayers == 0 :
        print("Warning: No Valid Ice Layers Found.")

    if Reduced:
        _output = _chi2/_nlayers
    else:
        _output = [_chi2, _nlayers]

    Wipe([ _dep_testcheck, _sca_testcheck, _abs_testcheck,
               _dep_check,     _sca_check,     _abs_check,          
                 _hyp_dep,       _hyp_sca,       _hyp_abs,
                 _tru_dep,       _tru_sca,       _tru_abs, 
                    _chi2,       _err_sca,       _err_abs,  _nlayers])
    
    return _output
