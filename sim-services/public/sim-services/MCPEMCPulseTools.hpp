#ifndef MCPEMCPULSETOOLS_H_INCLUDED
#define MCPEMCPULSETOOLS_H_INCLUDED

#include <simclasses/I3MCPulse.h>
#include <simclasses/I3MCPE.h>
#include <simclasses/I3ParticleIDMap.hpp>

namespace {
    
    struct indexRange{
        std::vector<uint32_t>::iterator cur, end;
        indexRange(std::vector<uint32_t>& container):
        cur(container.begin()),end(container.end()){}
    };
    
    template<typename HitType>
    struct HitTraits{};
    
    template<>
    struct HitTraits<I3MCPE>{
        typedef uint32_t charge_type;
        typedef double time_type;
        static charge_type& getCharge(I3MCPE& h){ return(h.npe); }
        static time_type& getTime(I3MCPE& h){ return(h.time); }
        static const charge_type& getCharge(const I3MCPE& h){ return(h.npe); }
        static const time_type& getTime(const I3MCPE& h){ return(h.time); }
        /*static void printParent(const I3MCPE& h){
            std::cout << '(' << h.major_ID << ',' << h.minor_ID << ')';
        }
        static bool sameParent(const I3MCPE& h1, const I3MCPE& h2){
            return(h1.major_ID==h2.major_ID && h1.minor_ID==h2.minor_ID);
        }*/
    };
    
    template<>
    struct HitTraits<I3MCPulse>{
        typedef float charge_type;
        typedef double time_type;
        static charge_type& getCharge(I3MCPulse& h){ return(h.charge); }
        static time_type& getTime(I3MCPulse& h){ return(h.time); }
        static const charge_type& getCharge(const I3MCPulse& h){ return(h.charge); }
        static const time_type& getTime(const I3MCPulse& h){ return(h.time); }
        static void printParent(const I3MCPulse& h){}
        static bool sameParent(const I3MCPulse& h1, const I3MCPulse& h2){ return(true); }
    };
  
    template<typename T>
    bool earlier(const T& h1, const T& h2){
        return(h1.time < h2.time);
    }
  
    template<typename HitType>
    class AugmentedHit : public HitType{
    public:
        uint32_t index;
        AugmentedHit(const HitType& p, uint32_t i): HitType(p), index(i){}
    };
  
    template<typename HitType>
    struct HitAugmenter{
    private:
        uint32_t index;
    public:
        HitAugmenter():index(0){}
        
        AugmentedHit<HitType> operator()(const HitType& p){
            return(AugmentedHit<HitType>(p,index++));
        }
    };
  
    template<typename HitType>
    struct HitSimplifier{
        HitType operator()(const AugmentedHit<HitType>& p) const{
            return(p);
        }
    };
  
}

namespace MCHitMerging{
    SET_LOGGER("MCHitMerging");

    // Merges together pulses which are within some small time window, compacting the hit series
    //  if it is densely filled with hits.
    // (Here, 'small' is intended to be relative to the resolution of the feature extractor.)
    // \param hits The hit series to be compacted
    // \pre The input hits must be time ordered
    template<typename HitType>
    void timeMergeHits(std::vector<HitType>& hits, ParticlePulseIndexMap& aux){
        //This is the time span in which adjacent hits will be merged
        //The chosen value, 200 picoseconds, is based on the advice of
        //N. Whitehorn with regard to the highest resolution which wavedeform
        //achieves in unfolding pulses
        const double timeWindow = .2; //ns
    
        typedef typename HitTraits<HitType>::charge_type charge_type;
        typedef typename HitTraits<HitType>::time_type time_type;
    
        if(hits.empty())
            return; //easy, nothing to do
    
        std::vector<indexRange> ranges;
        for(std::pair<const I3ParticleID, std::vector<uint32_t>>& pr : aux)
            ranges.push_back(indexRange(pr.second));
    
        uint32_t minOldIndex=0, maxOldIndex=1, newIndex=0, removed=0;
        typename std::vector<HitType>::iterator it1=hits.begin(), it2=hits.begin();
        typename std::vector<HitType> newHits;
        bool madeChanges=false;
        it2++;
        //Iterate over all hits
        while(it1 != hits.end()){
            bool doMerge=false;
            charge_type charge1=HitTraits<HitType>::getCharge(*it1);
            time_type time1=HitTraits<HitType>::getTime(*it1);
            charge_type weightSum=charge1;
            time_type avgTime=time1;
            log_trace_stream("merge base time now " << avgTime);
            //std::cout << "   first index is " << minOldIndex << "; hit parent is ";
            //HitTraits<HitType>::printParent(*it1);
            //std::cout << '\n';
            
            //Scan forward to collect all hits which can be merged
            while(it2 != hits.end() && (it2->time-it1->time) <= timeWindow){
                doMerge=true;
                charge_type charge2=HitTraits<HitType>::getCharge(*it2);
                time_type time2=HitTraits<HitType>::getTime(*it2);
                log_trace_stream("will merge hit " << time2-time1 << " ns later");
                //std::cout << "    merging with hit at index " << maxOldIndex << " whose parent is ";
                //HitTraits<HitType>::printParent(*it2);
                //std::cout << '\n';
                //if(HitTraits<HitType>::sameParent(*it1,*it2))
                //  std::cout << "     DING!\n";
                avgTime=(weightSum*avgTime + charge2*time2)/(weightSum + charge2);
                weightSum+=charge2;
                it2++;
                maxOldIndex++;
            }
            //Merge the set of hits if necessary
            if(doMerge){
                removed+=maxOldIndex-minOldIndex-1;
                if(!madeChanges){
                    //copy all of the unmerged hits that were skipped previously
                    std::copy(hits.begin(),it1,std::back_inserter(newHits));
                    newIndex+=std::distance(hits.begin(),it1);
                    //advance the aux iterators
                    for(indexRange& range : ranges)
                        range.cur=std::lower_bound(range.cur,range.end,newIndex);
                    madeChanges=true;
                }
                log_trace_stream("merged hit time is " << avgTime << ", and weight is " << weightSum);
                newHits.push_back(*it1);
                HitTraits<HitType>::getTime(newHits.back())=avgTime;
                HitTraits<HitType>::getCharge(newHits.back())=weightSum;
                //update aux, replacing all indices in [minOldIndex,maxOldIndex) with newIndex
                for(indexRange& range : ranges){
                    while(range.cur!=range.end && *range.cur<maxOldIndex){
                        *range.cur=newIndex; //overwrite the index with the index of the merged pulse
                        range.cur++; //advance to the next index
                    }
                }
                newIndex++;
            }
            //If we earlier began copying data to a new structure, we must continue doing so
            else if(madeChanges){
                std::copy(it1,it2,std::back_inserter(newHits));
                for(indexRange& range : ranges){
                    while(range.cur!=range.end && *range.cur<maxOldIndex){
                        *range.cur-=removed; //decrease the index by the number of pulses removed
                        range.cur++; //advance to the next index
                    }
                }
                newIndex+=std::distance(it1,it2);
            }
            it1=it2++; //skip to the next possible group
            minOldIndex=maxOldIndex++;
        }
        if(madeChanges){
            hits.swap(newHits);
            //deduplicate the contents of aux
            for(auto& particle : aux){
                std::vector<uint32_t> temp;
                std::unique_copy(particle.second.begin(),particle.second.end(),std::back_inserter(temp));
                particle.second.swap(temp);
            }
        }
    }
    
    //move particle parentage information out into a separate data structure
    //erase particle information from the hits to prevent later confusion
    I3ParticleIDMapPtr extractPIDInfo(I3Map<OMKey,std::vector<I3MCPE>>& hits);

    //move particle parentage information out into a separate data structure
    //erase particle information from the hits to prevent later confusion
    //run the time merging
    boost::shared_ptr<I3ParticleIDMap> extractPIDInfoandMerge(I3Map<OMKey,std::vector<I3MCPE>>& hits);
    
    ParticlePulseIndexMap extractPIDInfoandMerge(std::vector<I3MCPE>& hits, OMKey key);
    
    std::vector<I3ParticleID> findParents(uint32_t pulseIndex, const ParticlePulseIndexMap& idMap);
    
    // Given N pulses derived from M particles, with A being the average number
    // of parent particles per pulse this function should run in (average) time
    // proportional to N*M*log(A*N/M) and use additional space of approximately 
    // N*(sizeof(I3MCPulse)+4) + A*N*4 bytes = N*(sizeof(I3MCPulse) + 4*(A+1)) bytes
    //TODO: This code is horribly cumbersome. Can it be made less so?
    template<typename ForwardIterator>
    void sortMCHits(ForwardIterator begin, ForwardIterator end, ParticlePulseIndexMap& aux){
        log_debug_stream(" Sorting " << std::distance(begin,end) << " pulses");
        //make copies of all of the pulses, augmented with their initial indices
        typedef typename std::iterator_traits<ForwardIterator>::value_type HitType;
        std::vector<AugmentedHit<HitType>> hits;
        hits.reserve(std::distance(begin,end));
        std::transform(begin, end, std::back_inserter(hits), HitAugmenter<HitType>());
		
        //sort the augmented hits
        std::sort(hits.begin(),hits.end(), earlier<AugmentedHit<HitType>>);
		
        //update the particle ID mappings
        uint32_t newIndex=0;
        ParticlePulseIndexMap newAux(aux);
        for(const auto& pulse : hits){
            uint32_t oldIndex=pulse.index;
            //find every instance of the old index, and replace it with the new index
            //(in the copy of aux, so as not to confuse later lookups of old indices)
            ParticlePulseIndexMap::iterator newParticle=newAux.begin();
            for(const auto& particle : aux){
                std::vector<uint32_t>::const_iterator it = std::lower_bound(particle.second.begin(), particle.second.end(), oldIndex);
                if(it!=particle.second.end() && *it==oldIndex)
                    //this should be efficient since particle->second.begin() and it are RandomAccess
                    newParticle->second[std::distance(particle.second.begin(),it)] = newIndex;
                newParticle++;
            }
            newIndex++;
        }
        //replace old mappings with new
        aux.swap(newAux);
        //restore sorted order of mappings
        for(auto& particle : aux)
            std::sort(particle.second.begin(), particle.second.end());
		
        //copy the hits back to the original storage, leaving out the indices
        std::transform(hits.begin(), hits.end(), begin, HitSimplifier<HitType>());
    }


    class MCPEStream {
        
    public:
        bool insert(const I3MCPE &hit) {
            const double timeWindow = .2; //ns
            auto it = std::lower_bound(hits_.begin(), hits_.end(), hit);
            std::pair<double,double> dt(
                it == hits_.begin() ? std::numeric_limits<double>::infinity() : std::abs((it-1)->time - hit.time),
                it == hits_.end() ? std::numeric_limits<double>::infinity() : std::abs(it->time - hit.time)
                );
            bool inserted = true;
            if (dt.first < dt.second && dt.first <= timeWindow) {
                inserted = false;
                --it;
            } else if (dt.second < dt.first && dt.second <= timeWindow) {
                inserted = false;
            }
            
            if (inserted) {
                hits_.insert(it, Bunch(hit));
            } else {
                it->absorb(hit);
            }
            
            return inserted;
        }
        
        std::pair<std::vector<I3MCPE>, ParticlePulseIndexMap> extractMCPEsWithPIDInfo() const
        {
            std::vector<I3MCPE> output;
            ParticlePulseIndexMap particleMap;
            uint32_t counter = 0;
            
            for (const auto &bunch : hits_) {
                if (bunch.parents.size() == 1) {
                    output.emplace_back(bunch.parents.front(),bunch.npe,bunch.time);
                } else {
                    // Erase parentage information to avoid confusion
                    output.emplace_back(0,-1,bunch.npe,bunch.time);
                }
                // Insert parentage info in side map
                for (auto &id : bunch.parents) {
                    particleMap[id].push_back(counter);
                }
                ++counter;
            };
            
            return std::make_pair(std::move(output),std::move(particleMap));
        }
        
        size_t size() const { return hits_.size(); }
    private:
        struct Bunch {
            double time;
            uint32_t npe;
            std::vector<I3ParticleID> parents;
            
            Bunch(const I3MCPE &hit) : time(hit.time), npe(hit.npe), parents(1, hit.ID) {}
            
            void absorb(const I3MCPE &hit) {
                npe += hit.npe;
                auto it = std::lower_bound(parents.begin(), parents.end(), hit.ID);
                if (it == parents.end() || *it < hit.ID)
                    parents.insert(it, hit.ID);
            }
            
            bool operator<(const I3MCPE &hit) const { return time < hit.time; }
        };
        std::vector<Bunch> hits_;
    };
}

#endif
