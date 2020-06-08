import os, urllib, hashlib

from icecube import icetray, dataio

class I3GCDSanityChecker(icetray.I3Module):
    def __init__(self, context):
        super().__init__(context)

        self.__url = 
        self.AddParameter('GCDFilename', 'Fully resloved filename of GCD file.', None)
        self.AddParameter('Prescale', 'Only check every Nth Q-frame', None)

        self.__q_frame_counter = 0
        self.__validated = False
        self.__stops_called = {'Geometry': 0,
                               'Calbration': 0,
                               'DetectorStatus': 0}
        
    def Configure(self):
        self.gcd_filename = self.GetParameter('GCDFilename')
        self.prescale = self.GetParameter('Prescale')
        
        if not os.path.exists(self.filename):
            icetray.logging.log_fatal("GCD file %s not found." % self.filename)
                
        with dataio.I3File(self.filename) as f:
            self.geometry_frame = None
            self.calibration_frame = None
            self.detector_status_frame = None

            while f.more():
                frame = f.pop_frame()
                if 'I3Geometry' in frame:
                    self.geometry_frame = frame['I3Geometry']
                if 'I3Calibration' in frame:
                    self.calibration_frame = frame['I3Calibration']
                if 'I3DetectorStatus' in frame:
                    self.detector_status_frame = frame['I3DetectorStatus']

        # We'll need to get the checksums from some location
        try:
            url = self.__url + 'checksums'
            response = urllib.urlopen(url)
            self.checksums = pickle.loads(response.read())
        except:
            icetray.logging.log_fatal("Something went wrong checking/loading checksums")
            
    def Geometry(self, frame):
        frame['I3Geometry'] = self.geometry_frame
        self.__stops_called['Geometry'] += 1
        self.PushFrame(frame)
        
    def Calibration(self, frame):
        frame['I3Calibration'] = self.calibration_frame
        self.__stops_called['Calibration'] += 1
        self.PushFrame(frame)
    
    def DetectorStatus(self, frame):
        frame['I3DetectorStatus'] = self.detector_status_frame
        self.__stops_called['DetectorStatus'] += 1
        self.PushFrame(frame)

    def __check_frame(self, frame):
        for frame_key in ['I3Geometry', 'I3Calibration', 'I3DetectorStatus']:
            # TODO: Need to expose something via pybindings; either the underlying
            #       state via pickle (or buffer) or expose a hash (bound to __hash__, perhaps).
            #       There's really no reliable way, that I can see now, to get the
            #       checksum of the memory representation of frame objects.
            m = hashlib.sha256()
            obj = frame[frame_key]
            m.update(obj.buffer) # <--- need pybindings for this.
            if m.hexdigest() != self.checksums[frame_key]:
                icetray.logging.log_fatal('The %s frame has been altered.' % frame_key)

            ### another way might be, if __hash__ is bound
            # checksum = hash(frame[frame_key]) # <-- says nothing about the algo though
            ### ...or...
            # checksum = frame[frame_key].sha256()  # <-- need bindings too
            # if checksum != self.checksums[frame_key]:
            #    icetray.logging.log_fatal('The %s frame has been altered.' % frame_key)
        
    def DAQ(self, frame):
        '''
        Only check the G,C, and D frames at the beginning of the first Q frame.
        '''
        self.__q_frame_counter += 1
        
        if self.__validated or self.__q_frame_counter % self.prescale:
            self.PushFrame(frame)
            return

        self.__check_frame(frame)
        self.PushFrame(frame)
                
    def Finish(self, frame):
        self.__check_frame(frame)

        for stop, ncalls in self.__stops_called.items():
            if ncalls > 1 :
                icetray.logging.log_fatal("%s called %d times." % (stop, ncalls))
