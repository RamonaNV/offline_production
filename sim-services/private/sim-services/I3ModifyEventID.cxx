/**
 * class: I3ModifyEventID.cxx
 * (c) 2004 IceCube Collaboration
 * Version $Id: I3ModifyEventID.cxx,v 1.16 2005/03/17 18:33:43 olivas Exp $
 *
 * Date 08 June 2006
 * @version $Revision: 1.1 $
 * @date $Date: 2005/03/17 18:33:43 $
 * @author Alex Olivas <olivas@icecube.umd.edu>
 *
 */
#include "icetray/I3TrayHeaders.h"
#include "icetray/I3Tray.h"
#include "icetray/I3Module.h"
#include "dataclasses/I3Time.h"
#include "dataclasses/physics/I3EventHeader.h"

/** 
 *@class I3ModifyEventID I3ModifyEventID.h 
 */
class I3ModifyEventID : public I3Module
{
 public:

  I3ModifyEventID(const I3Context& ctx);
  ~I3ModifyEventID();
  
  void Configure();
  void DAQ(I3FramePtr frame);

 private:
  I3ModifyEventID();
  I3ModifyEventID(const I3ModifyEventID&);
  I3ModifyEventID& operator=(const I3ModifyEventID&);

  //Start time of run period
  int year_;
  int64_t daqTime_;
  int mjd_;
  int mjd_s_;
  double mjd_ns_;
  unsigned runNumber_;
  unsigned subRunNumber_;
  unsigned startEID_;

  bool modTime_;
  bool modRunId_;
  bool modEventId_;

  SET_LOGGER("I3ModifyEventID");

};

I3_MODULE(I3ModifyEventID);

const int DEFAULT_YEAR = 2006;
const int64_t DEFAULT_DAQTIME = 0ULL;
const int DEFAULT_MJD_SECONDS(0);
const double DEFAULT_MJD_NANOSECONDS(0.);
const int DEFAULT_RUN_NUMBER(0);
const int DEFAULT_SUBRUN_NUMBER(0);

I3ModifyEventID::I3ModifyEventID(const I3Context& ctx) : 
  I3Module(ctx),
  year_(DEFAULT_YEAR),
  daqTime_(DEFAULT_DAQTIME),
  mjd_(INT_MIN),
  mjd_s_(DEFAULT_MJD_SECONDS),
  mjd_ns_(DEFAULT_MJD_NANOSECONDS),
  runNumber_(DEFAULT_RUN_NUMBER),
  subRunNumber_(DEFAULT_SUBRUN_NUMBER),
  startEID_(0),
  modTime_(0),
  modRunId_(0),
  modEventId_(0)
{ 
  log_debug("Constructor I3ModifyEventID");

  AddParameter("Year", "Year of the run", year_);
  AddParameter("DAQTime", "DAQTime of the run in 1/10 of ns", daqTime_);
  AddParameter("MJD","Modified Julian Date",mjd_);
  AddParameter("MJDSeconds","Number of seconds after the start of the MJD.",mjd_s_);
  AddParameter("MJDNanoSeconds","Number of nanoseconds after the start of the MJD.",mjd_ns_);
  AddParameter("RunNumber", "Run Number", runNumber_);
  AddParameter("SubRunNumber","SubRun Number",subRunNumber_);
  AddParameter("StartEventID","Starting Event ID",startEID_);
  AddParameter("ModifyTime","Flag: Modify event time", modTime_);
  AddParameter("ModifyRunId","Flag: Modify run ID", modRunId_);
  AddParameter("ModifyEventId","Flag: Modify event ID", modEventId_);
  AddOutBox("OutBox");
}

I3ModifyEventID::~I3ModifyEventID(){}

void I3ModifyEventID::Configure()
{
  log_error("This module is deprecated.");
  log_error("If you actually need it please submit a ticket requesting that it not be deprecated.");
  GetParameter("StartEventID", startEID_);
  GetParameter("Year", year_);
  GetParameter("DAQTime", daqTime_);
  GetParameter("MJD",mjd_);
  GetParameter("MJDSeconds",mjd_s_);
  GetParameter("MJDNanoSeconds",mjd_ns_);
  GetParameter("RunNumber", runNumber_);
  GetParameter("SubRunNumber", subRunNumber_);
  GetParameter("ModifyTime", modTime_);
  GetParameter("ModifyRunId", modRunId_);
  GetParameter("ModifyEventId", modEventId_);

  if(mjd_ != INT_MIN && 
     (year_ != DEFAULT_YEAR || daqTime_ != DEFAULT_DAQTIME ))
     log_fatal("Ambiguous settings : Please choose either Mjd or Year and DAQTime.  Not both.");

  if(mjd_ != INT_MIN){
    I3Time t;
    t.SetModJulianTime(mjd_,mjd_s_,mjd_ns_);
    year_ = t.GetUTCYear();
    daqTime_ = t.GetUTCDaqTime();
  }
}

void I3ModifyEventID::DAQ(I3FramePtr frame)
{
  log_debug("DAQ");

  if (frame->Has(I3DefaultName<I3EventHeader>::value()) ) {
          const I3EventHeader& oldHeader = frame->Get<I3EventHeader>();
          I3EventHeaderPtr ehPtr(new I3EventHeader(oldHeader));

          if (modEventId_) {
              ehPtr->SetEventID(startEID_++);
          }
          if (modTime_) {
              I3Time evtTime(year_, daqTime_);
              ehPtr->SetStartTime(evtTime);
          }
          if (modRunId_) {
              ehPtr->SetRunID(runNumber_);
              ehPtr->SetSubRunID(subRunNumber_);
          }
          frame->Delete("I3EventHeader");
          frame->Put(ehPtr);
  } else { 
          log_info("No I3EventHeader to modify.");
  }

  PushFrame(frame,"OutBox");
}//DAQ()
 
