/**
 * Copyright (c) 2013
 * Juan Carlos Diaz-Velez <juancarlos@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3RemoveLargeDT.cxx 137416 2015-09-09 20:15:44Z olivas $
 *
 * @file I3RemoveLargeDT.cxx
 * @version $Revision: 137416 $
 * @date $Date: 2015-09-09 14:15:44 -0600 (Wed, 09 Sep 2015) $
 * @author Juan Carlos Diaz-Velez
 */

#include <boost/foreach.hpp>

#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <simclasses/I3MCPE.h>

// these are utilities we're going to use later
namespace{
  // for time ordering vectors of I3MCPEs
  bool compare(const I3MCPE& lhs, const I3MCPE& rhs){
    if(lhs.time < rhs.time) return true;
    return false;
  }

  // this is mostly for clipping I3MCPEs outside a time range
  struct outside_time_range{
  private:
    double min_time_;
    double max_time_;

  public:
    int n_outliers;

    outside_time_range(double min_time, double max_time) :
      min_time_(min_time),
      max_time_(max_time),
      n_outliers(0)
    {};

    bool operator()(const I3MCPE& pe){
      return ((pe.time < min_time_) || (pe.time > max_time_));
    }
  };
}

/**
 * @brief Removes photo-electron hits that are separated in time by a factor 
 * larger than maxDT/2 from the median time (where maxDT is the maximum size 
 * of the trigger window).
 * 
 * The purpose is to eliminate hits that, while physically related to a triggered
 * event, would never be associated to that event by the DAQ.
 * The long gaps will otherwise be filled by noise and beacon hits in DOMLauncher
 * and will unnecessarily blow up memory consumption.
 * 
 * Clips any PEs later than MaxDeltaT, taken with respect to the earliest
 * PE in the frame.
 */
class I3RemoveLargeDT: public I3Module {
public:
  I3RemoveLargeDT(const I3Context& ctx);
  ~I3RemoveLargeDT(){};

  void Configure();
  void DAQ(I3FramePtr frame);
  void Finish(){};

private:
  std::string inputResponse_;
  std::string outputResponse_;
  double maxdt_;
  bool presorted_;

  SET_LOGGER("I3RemoveLargeDT");
};

I3RemoveLargeDT::I3RemoveLargeDT(const I3Context& ctx) :
  I3Module(ctx),
  inputResponse_("I3MCPESeriesMap"),
  outputResponse_("CleanedI3MCPESeriesMap"),
  maxdt_(100*I3Units::ms),
  presorted_(true)
{
  AddParameter("MaxDeltaT",
              "Largest time span of PEs in an event.",
               maxdt_);
  AddParameter("InputResponse",
               "Name of the input response series",
               inputResponse_);

  AddParameter("OutputResponse",
               "Name of the output response series",
               outputResponse_);

  AddParameter("PreSorted",
               "PEs are already sorted in time",
               presorted_);

  AddOutBox("OutBox");
}

void I3RemoveLargeDT::Configure()
{
    GetParameter("MaxDeltaT", maxdt_);
    GetParameter("InputResponse", inputResponse_);
    GetParameter("OutputResponse", outputResponse_);
    GetParameter("PreSorted", presorted_);
}

void I3RemoveLargeDT::DAQ(I3FramePtr frame)
{
  if(!frame->Has(inputResponse_)){ // push and return
    log_fatal("I3MCPESeriesMap '%s' doesn't exist in the frame.", 
	inputResponse_.c_str());
  }

  I3MCPESeriesMapConstPtr 
	  input = frame->Get<I3MCPESeriesMapConstPtr>(inputResponse_);

  // determine the earliest and latest hits
  // also sorting in-place, which is why we're iterating
  // over the output map, which at this point is just a
  // copy of the input map.
  double earliest_time = std::numeric_limits<double>::max();
  double latest_time = std::numeric_limits<double>::min();

  // we're going to use this to calculate the median
  // ...but only if we need to
  // we could use pe_times to determine the earliest and
  // latest times, but i don't want to have to sort this
  // if i don't have to.
  std::vector<double> pe_times;

  I3MCPESeriesMapPtr output(new I3MCPESeriesMap(*input));
  BOOST_FOREACH(I3MCPESeriesMap::value_type& map_pair, *output){

    if (!presorted_) {
      std::sort(map_pair.second.begin(), map_pair.second.end(), compare);
    }

    earliest_time = std::min(map_pair.second.front().time, earliest_time);
    latest_time = std::max(map_pair.second.back().time, latest_time);

    BOOST_FOREACH(const I3MCPE& pe, map_pair.second){
      pe_times.push_back(pe.time);
    }
  }

  // quickly determine if there's even any heavy lifting to do
  if(latest_time - earliest_time > maxdt_){
    // need to remove anything outside the time range
    // start lifting heavy things...
    std::sort(pe_times.begin(), pe_times.end());

    double median_time = *(pe_times.begin() + pe_times.size()/2);
    double min_time(median_time - maxdt_/2);
    double max_time(median_time + maxdt_/2);
    outside_time_range is_outlier(min_time, max_time);
    BOOST_FOREACH(I3MCPESeriesMap::value_type& map_pair, *output){

      // check to see if the whole series is in range first
      if(!is_outlier(map_pair.second.front()) &&
         !is_outlier(map_pair.second.back())){
        continue; // go to the next DOM
      }

      // nope...gotta clip some
      I3MCPESeries::iterator new_end =
        std::remove_if(
		map_pair.second.begin(), 
		map_pair.second.end(), 
		is_outlier);

      map_pair.second.resize(std::distance(map_pair.second.begin(), new_end));
                
    }
  }

  // In case we are overwriting object
  if (outputResponse_ == inputResponse_) {
    frame->Put("old"+inputResponse_,input);
    frame->Delete(inputResponse_);
  }
  frame->Put(outputResponse_, output);
  PushFrame(frame,"OutBox");
}

I3_MODULE(I3RemoveLargeDT);


