#!/usr/bin/env python

import sys
from icecube import clsim, icetray, dataclasses
# skip out if I3CLSimServer was not built
try:
    clsim.I3CLSimServer
except AttributeError:
    sys.exit(0)

import time
import tempfile
from multiprocessing import Process
try:
    # PY2
    from Queue import Queue
except ImportError:
    # PY3
    from queue import Queue
from numpy.random import uniform
from numpy.random import seed
seed(0)
from numpy import testing

icetray.logging.set_level('INFO')

def dummy_photon(step):
    photon = clsim.I3CLSimPhoton()
    for attr in 'x', 'y', 'z', 'theta', 'phi', 'time', 'weight', 'id':
        setattr(photon, attr, getattr(step, attr))
    photon.numScatters = 3
    photon.omID = 52
    photon.stringID = 23
    return photon

def dummy_photon_history(photon):
    history = clsim.I3CLSimPhotonHistory()
    for i in range(photon.numScatters):
        history.append(dataclasses.I3Position(i,i+0.5,i+3.14), i)
    return history

class DummyConverter(clsim.I3CLSimStepToPhotonConverter):
    def __init__(self):
        super(DummyConverter, self).__init__()
        self.input_queue = Queue()
    def IsInitialized(self):
        return True
    def GetWorkgroupSize(self):
        return 8
    def GetMaxNumWorkitems(self):
        return 64
    def MorePhotonsAvailable(self):
        return not self.input_queue.empty()
    def EnqueueSteps(self, steps, id):
        # ensure that EnqueueSteps blocks to test backpressure
        time.sleep(0.1)
        self.input_queue.put((steps, id))
    def GetConversionResult(self):
        steps, id = self.input_queue.get()
        icetray.logging.log_debug('{} steps in bunch {}'.format(len(steps), id), unit='clsim.testCLSimServer')
        photons = clsim.I3CLSimPhotonSeries([dummy_photon(step) for step in steps if step.num > 0])
        history = clsim.I3CLSimPhotonHistorySeries(map(dummy_photon_history, photons))
        return clsim.I3CLSimStepToPhotonConverter.ConversionResult_t(id, photons, history)

def test_client(client, num_steps, base=0):
    from threading import Thread
    input_steps = dict()
    def feed():
        for i in range(base,base+num_steps):
            steps = clsim.I3CLSimStepSeries()
            for _ in range(int(uniform(10*num_steps/10.))+1):
                steps.append(clsim.I3CLSimStep())
                steps[-1].pos = dataclasses.I3Position(*uniform(size=3))
                steps[-1].dir = dataclasses.I3Direction(*uniform(size=3))
                steps[-1].time = uniform()
                steps[-1].weight = uniform()
                steps[-1].num = 1
            input_steps[i] = steps
            icetray.logging.log_debug("submitting bunch {}/{} with size {}".format(i+1,num_steps,len(steps)), unit='clsim.testCLSimServer')
            client.EnqueueSteps(steps, i)
        if hasattr(client, 'EnqueueBarrier'):
            client.EnqueueBarrier()
    def drain():
        for i in range(num_steps):
            icetray.logging.log_debug("getting bunch {}/{}".format(i+1,num_steps), unit='clsim.testCLSimServer')
            if hasattr(client, 'GetConversionResultWithBarrierInfo'):
                result, barrier = client.GetConversionResultWithBarrierInfo()
            else:
                result, barrier = client.GetConversionResult(), i == num_steps-1
            icetray.logging.log_debug("got bunch id {}: {} photons (barrier: {})".format(result.identifier, len(result.photons), barrier), unit='clsim.testCLSimServer')
            input = input_steps[result.identifier]
            try:
                assert len(result.photons) == len(input)
            except AssertionError:
                icetray.logging.log_error("{} != {} in bunch {}".format(len(result.photons), len(input), result.identifier), unit='clsim.testCLSimServer')
                raise
            assert len(result.photonHistories) == len(input)
            for step, photon, history in zip(input, result.photons, result.photonHistories):
        
                testing.assert_equal( photon.numScatters, 3 )
                testing.assert_equal( photon.omID, 52 )
                testing.assert_equal( photon.stringID, 23 )
                for attr in 'x', 'y', 'z', 'theta', 'phi', 'time', 'weight':
                    testing.assert_equal( getattr(step, attr), getattr(photon, attr), err_msg='{} not equal'.format(attr))
            
                dummy_history = dummy_photon_history(photon)
                testing.assert_equal( len(dummy_history), len(history) )
                for pos, expected_pos in zip(history, dummy_history):
                    testing.assert_equal( pos, expected_pos )
            if i == num_steps-1:
                assert barrier
            else:
                assert not barrier
    # drain in a thread, since EnqueueSteps() may block
    t = Thread(target=drain)
    t.start()
    feed()
    t.join()
    icetray.logging.log_info("base {} done".format(base), unit='clsim.testCLSimServer')

# icetray.logging.set_level_for_unit('I3CLSimServer', 'TRACE')
# icetray.logging.set_level_for_unit('I3CLSimServer', 'DEBUG')

icetray.logging.set_level_for_unit('I3CLSimClient', 'TRACE')
icetray.logging.set_level_for_unit('clsim.testCLSimServer', 'TRACE')

# First, ensure that the test passes when the converter is called directly
test_client(DummyConverter(), 10)

# Now, call through the server in a separate process
converters = clsim.I3CLSimStepToPhotonConverterSeries([DummyConverter()])
server = clsim.I3CLSimServer('tcp://127.0.0.1:*',converters)
address = server.GetAddress()

def fire_a_few(num_steps=10, base=0):
    # NB: the Python logging bridge deadlocks from secondary threads in Py3
    if sys.version_info.major == 2:
        icetray.logging.BASIC_FORMAT = "{} %(filename)s:%(lineno)s %(levelname)s: %(message)s".format(base)
        icetray.logging.console()
    icetray.logging.set_level_for_unit('clsim.testCLSimServer', 'TRACE')
    icetray.logging.set_level_for_unit('I3CLSimClient', 'TRACE')
    
    icetray.logging.log_debug("client {} connecting to {}".format(base, address), unit='clsim.testCLSimServer')
    client = clsim.I3CLSimClient(address)
    icetray.logging.log_debug("client {} connected".format(base), unit='clsim.testCLSimServer')
    testing.assert_equal( client.workgroupSize, 8 )
    testing.assert_equal( client.maxNumWorkitems, 64 )
    test_client(client, num_steps, base)

procs = [Process(target=fire_a_few, kwargs=dict(num_steps=10, base=10*i)) for i in range(10)]
for p in procs:
    p.start()
for p in procs:
    p.join()
    icetray.logging.log_info("process {} exited with status {}".format(p.pid, p.exitcode))
    assert p.exitcode == 0

icetray.logging.log_info("going to destroy")
del server
icetray.logging.log_info("destroyed")


