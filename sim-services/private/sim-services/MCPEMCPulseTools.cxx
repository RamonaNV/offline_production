#include <sim-services/MCPEMCPulseTools.hpp>
#include <icetray/I3ConditionalModule.h>

namespace MCHitMerging{

    std::vector<I3ParticleID> findParents(uint32_t pulseIndex, const ParticlePulseIndexMap& idMap){
      std::vector<I3ParticleID> results;
      for(ParticlePulseIndexMap::const_iterator particle=idMap.begin(), particleEnd=idMap.end(); particle!=particleEnd; particle++){
        std::vector<uint32_t>::const_iterator it=std::lower_bound(particle->second.begin(),particle->second.end(),pulseIndex);
        if(it!=particle->second.end() && *it==pulseIndex)
          results.push_back(particle->first);
      }
      return(results);
    }
    
    ParticlePulseIndexMap extractPIDInfoandMerge(std::vector<I3MCPE>& hits, OMKey key){
        ParticlePulseIndexMap particleMap;
        uint32_t counter = 0;
        double lastTime = -std::numeric_limits<double>::max();
        for(I3MCPE& hit : hits){
            I3ParticleID pid(hit.ID.majorID, hit.ID.minorID);
            particleMap[pid].push_back(counter);
            hit.ID.majorID = 0;
            hit.ID.minorID = -1;
            //while we're here, check time ordering
            if(hit.time < lastTime)
                log_fatal_stream("MCPEs on " << key << " are not time ordered: " << lastTime << ", " << hit.time);
            lastTime = hit.time;
            counter++;
        }
        // now we know that the hits are time ordered, and we have the aux structure ready
        timeMergeHits(hits, particleMap);
        log_trace_stream("compressed " << counter << " p.e. on " << key << " to " << hits.size());
        //std::cout << "compressed " << counter << " p.e. on " << domIt->first << " to " << domIt->second.size() << std::endl;
        //unsigned int entries=0;
        //for(const auto& particleRecord : particleMap)
        //    entries += particleRecord.second.size();
        //std::cout << " stored " << entries << " pulse entries from " << particleMap.size() << " particles" << std::endl;

        return(particleMap);
    }
    
    I3ParticleIDMapPtr extractPIDInfo(I3Map<OMKey,std::vector<I3MCPE>>& hits){
        I3ParticleIDMapPtr resultPIDMap(new I3ParticleIDMap());
        for(auto& domRecord : hits){
            const OMKey& dom=domRecord.first;
            ParticlePulseIndexMap& particleMap = (*resultPIDMap)[dom];
            uint32_t counter = 0;
            for(I3MCPE& hit : domRecord.second){
                particleMap[hit.ID].push_back(counter++);
                hit.ID.majorID = 0;
                hit.ID.minorID = -1;
            }
        }
        return(resultPIDMap);
    }
    
    I3ParticleIDMapPtr extractPIDInfoandMerge(I3Map<OMKey,std::vector<I3MCPE>>& hits){
        I3ParticleIDMapPtr resultPIDMap(new I3ParticleIDMap());
        for(auto& domRecord : hits){
            const OMKey& dom=domRecord.first;
            ParticlePulseIndexMap& particleMap = (*resultPIDMap)[dom];
            uint32_t counter = 0;
            double lastTime = -std::numeric_limits<double>::max();
            for(I3MCPE& hit : domRecord.second){
                I3ParticleID pid(hit.ID);
                particleMap[pid].push_back(counter++);
                hit.ID.majorID = 0;
                hit.ID.minorID = -1;
                //while we're here, check time ordering
                if(hit.time<lastTime)
                    log_fatal_stream("MCPEs on " << dom << " are not time ordered: " << lastTime << ", " << hit.time);
                lastTime = hit.time;
            }
            // now we know that the hits are time ordered, and we have the aux structure ready
            timeMergeHits(domRecord.second, particleMap);
            log_trace_stream("compressed " << counter << " p.e. on " << dom << " to " << domRecord.second.size());
            //unsigned int entries=0;
            //for(const auto& particleRecord : particleMap)
            //    entries += particleRecord.second.size();
            //std::cout << " stored " << entries << " pulse entries from " << particleMap.size() << " particles" << std::endl;
        }

        return(resultPIDMap);
    }
}

class I3MCPEMerger : public I3ConditionalModule{
public:
    I3MCPEMerger(const I3Context& context):
    I3ConditionalModule(context),
    inputName("I3MCPESeriesMap")
    {
        AddParameter("Input","Name of the frame object to be compressed",inputName);
        AddOutBox("OutBox");
    }
    
    void Configure(){
        GetParameter("Input",inputName);
    }
    
    void DAQ(boost::shared_ptr<I3Frame> frame){
        if(!frame->Has(inputName)){
            log_warn_stream("Frame does not have \"" << inputName << "\", ignoring");
            PushFrame(frame);
            return;
        }
        if(frame->Has(inputName+"ParticleIDMap")){
            log_fatal_stream("Frame already contains " << (inputName+"ParticleIDMap")
                             << " suggesting that " << inputName << " is already compressed!");
        }
        
        //We have to copy all of the input data. This is gonna hurt.
        boost::shared_ptr<I3Map<OMKey,std::vector<I3MCPE>>> hits=
        boost::make_shared<I3Map<OMKey,std::vector<I3MCPE>>>(
            frame->Get<I3Map<OMKey,std::vector<I3MCPE>>>(inputName)
        );
        //Compute size before
        size_t initialSize=0;
        for(const auto& p : *hits)
            initialSize+=p.second.size()*sizeof(I3MCPE);
        //Do the merge
        I3ParticleIDMapPtr sideTable=MCHitMerging::extractPIDInfoandMerge(*hits);
        //Compute size after
        size_t finalSize=0;
        for(const auto& p : *hits)
            finalSize+=p.second.size()*sizeof(I3MCPE);
        
        for(const auto& p : *sideTable){
            finalSize+=sizeof(ParticlePulseIndexMap)
                +p.second.size()*(sizeof(I3ParticleID)+sizeof(ParticlePulseIndexMap::mapped_type));
            for(const auto& p2 : p.second)
                finalSize+=p2.second.size()*sizeof(ParticlePulseIndexMap::mapped_type::value_type);
        }
        //Note that these size estimates are _in memory_, not serialized, so
        //they also count padding bytes. Size on disk is typically smaller, and
        //compressed files may of course introduce further differences.
        log_trace_stream("Reduced " << initialSize << " bytes of data to " << finalSize);
            
        //Store output, if we actually made an improvement
        if(finalSize<initialSize){
            //inputName is also the output name!
            frame->Delete(inputName);
            frame->Put(inputName,hits);
            frame->Put(inputName+"ParticleIDMap",sideTable);
        }
        
        PushFrame(frame);
    }
private:
    std::string inputName;
    SET_LOGGER("I3MCPEMerger");
};

I3_MODULE(I3MCPEMerger)
