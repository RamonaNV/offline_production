#include <I3Test.h>
#include <simclasses/I3MCPulse.h>
#include <simclasses/I3MCPE.h>
#include <sim-services/MCPEMCPulseTools.hpp>


TEST_GROUP(HitMerging);

template<typename HitType>
bool isTimeOrdered(const std::vector<HitType>& hits)
{
    if(hits.size()<=1)
        return(true);
    typename std::vector<HitType>::const_iterator last, it, end;
    for(last=hits.begin(), it=hits.begin()+1, end=hits.end(); 
        it!=end; 
        it++,last++)
    {
        if(last->time>it->time)
            return(false);
    }
    return(true);
}

//Check that hit merging works properly:
//Hits with times with .2 ns should be merged
//The merged hit should have a weight equal to the sum of the weights of the original hits,
//and a time given by the weighted average of the times of the original hits.
//Parent particle information should also be correctly updated, in particular when two pulses
//with different parent particles are merged the index of the resulting pulse should appear
//in the index lists of both parents.
TEST(1_MCPEMergingSortingMerging)
{
    std::vector<I3MCPE> hits;
    ParticlePulseIndexMap aux;
    I3ParticleID p1(538,16);
    I3ParticleID p2(619,12);
    uint32_t hitIndex=0;

    //these two hits should be merged
    hits.push_back(I3MCPE(p1,1.0, 0.0));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPE(p1,1.0,0.1));
    aux[p1].push_back(hitIndex++);

    //this hit should not merge with the preceding or following hits
    hits.push_back(I3MCPE(p1,1.0,0.25));
    aux[p1].push_back(hitIndex++);

    //these two hits should be merged
    //the merged hits should be a child of both p1 and p2
    hits.push_back(I3MCPE(p1,1.0,0.5));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPE(p2,1.0,0.6));
    aux[p2].push_back(hitIndex++);

    //neither of these hits should merge with anything
    hits.push_back(I3MCPE(p1,1.0,0.75));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPE(p2,1.0,1.0));
    aux[p2].push_back(hitIndex++);

    MCHitMerging::timeMergeHits(hits, aux);

    ENSURE(hits.size()==5,"Should have 5 pulses after merging");
    ENSURE(aux.size()==2,"Number of parent particles should be unchanged");
    ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
    ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
    ENSURE(aux[p1].size()==4,"Particle 1 should have 4 children");
    ENSURE(aux[p2].size()==2,"Particle 2 should have 2 children");

    uint32_t p1Expected[4]={0,1,2,3};
    uint32_t p2Expected[2]={2,4};

    ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1Expected)==std::make_pair(aux[p1].end(),p1Expected+4));
    ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2Expected)==std::make_pair(aux[p2].end(),p2Expected+2));
    
    hitIndex=hits.size();

    // adding hits that appear to come within
    // the already sorted and merged hits
    I3ParticleID p3(741,13);
    hits.push_back(I3MCPE(p3,1.0,0.35));
    aux[p3].push_back(hitIndex++);
    hits.push_back(I3MCPE(p1,1.0,0.30));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPE(p2,1.0,1.3));
    aux[p2].push_back(hitIndex++);

    MCHitMerging::sortMCHits(hits.begin(), hits.end(), aux);

    // Check again that the state still makes sense
    ENSURE(isTimeOrdered<I3MCPE>(hits));
    ENSURE(hits.size()==8,"Should have 8 pulses after sorting");
    ENSURE(aux.size()==3,"Number of parent particles should be unchanged");
    ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
    ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
    ENSURE(aux[p1].size()==5,"Particle 1 should have 5 children");
    ENSURE(aux[p2].size()==3,"Particle 2 should have 3 children");
    ENSURE(aux[p3].size()==1,"Particle 3 should have 1 child");
    
    uint32_t p1ExpectedSort[5]={0,1,2,4,5};
    uint32_t p2ExpectedSort[3]={4,6,7};
    uint32_t p3ExpectedSort[1]={3};

    ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1ExpectedSort)==std::make_pair(aux[p1].end(),p1ExpectedSort+5));
    ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2ExpectedSort)==std::make_pair(aux[p2].end(),p2ExpectedSort+3));
    ENSURE(std::mismatch(aux[p3].begin(),aux[p3].end(),p3ExpectedSort)==std::make_pair(aux[p3].end(),p3ExpectedSort+1));

    MCHitMerging::timeMergeHits(hits, aux);

    ENSURE(isTimeOrdered<I3MCPE>(hits));
    ENSURE(hits.size()==5,"Should have 6 pulses after sorting");
    ENSURE(aux.size()==3,"Number of parent particles should be unchanged");
    ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
    ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
    ENSURE(aux[p1].size()==3,"Particle 1 should have 3 children");
    ENSURE(aux[p2].size()==3,"Particle 2 should have 3 children");
    ENSURE(aux[p3].size()==1,"Particle 3 should have 1 child");

    uint32_t p1ExpectedSortMerge[3]={0,1,2};
    uint32_t p2ExpectedSortMerge[3]={2,3,4};
    uint32_t p3ExpectedSortMerge[1]={1};

    ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1ExpectedSortMerge)==std::make_pair(aux[p1].end(),p1ExpectedSortMerge+3));
    ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2ExpectedSortMerge)==std::make_pair(aux[p2].end(),p2ExpectedSortMerge+3));
    ENSURE(std::mismatch(aux[p3].begin(),aux[p3].end(),p3ExpectedSortMerge)==std::make_pair(aux[p3].end(),p3ExpectedSortMerge+1));
}

//Check that hit merging works properly:
//Hits with times with .2 ns should be merged
//The merged hit should have a weight equal to the sum of the weights of the original hits,
//and a time given by the time of the first hit encountered.
//Parent particle information should also be correctly updated, in particular when two pulses
//with different parent particles are merged the index of the resulting pulse should appear
//in the index lists of both parents.
TEST(2_MCPEStreaming)
{
    MCHitMerging::MCPEStream stream;
    
    I3ParticleID p1(538,16);
    I3ParticleID p2(619,12);
    
    //these two hits should be merged
    ENSURE(stream.insert(I3MCPE(p1,1.0, 0.0)), "First hit has nothing to merge with");
    ENSURE_EQUAL(stream.size(), size_t(1));
    ENSURE(!stream.insert(I3MCPE(p1,1.0,0.1)), "Second hit should be merged");
    ENSURE_EQUAL(stream.size(), size_t(1));
    
    //this hit should not merge with the preceding or following hits
    ENSURE(stream.insert(I3MCPE(p1,1.0,0.25)), "Hit too far away to merge");
    ENSURE_EQUAL(stream.size(), size_t(2));
    
    //these two hits should be merged
    //the merged hits should be a child of both p1 and p2
    ENSURE(stream.insert(I3MCPE(p1,1.0,0.5)));
    ENSURE(!stream.insert(I3MCPE(p2,1.0,0.6)));
    ENSURE_EQUAL(stream.size(), size_t(3));
    
    //neither of these hits should merge with anything
    ENSURE(stream.insert(I3MCPE(p1,1.0,0.75)));
    ENSURE(stream.insert(I3MCPE(p2,1.0,1.0)));
    
    ENSURE(stream.size()==5,"Should have 5 pulses after merging");
    
    {
        std::vector<I3MCPE> hits;
        ParticlePulseIndexMap aux;
        std::tie(hits, aux) = stream.extractMCPEsWithPIDInfo();
        
        ENSURE(isTimeOrdered<I3MCPE>(hits));
        
        ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
        ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
        ENSURE_EQUAL(aux[p1].size(), 4u,"Particle 1 should have 4 children");
        ENSURE_EQUAL(aux[p2].size(), 2u,"Particle 2 should have 2 children");
        
        uint32_t p1Expected[4]={0,1,2,3};
        uint32_t p2Expected[2]={2,4};
        
        ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1Expected)==std::make_pair(aux[p1].end(),p1Expected+4));
        ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2Expected)==std::make_pair(aux[p2].end(),p2Expected+2));
    }
    
    // adding hits that appear to come within
    // the already sorted and merged hits
    I3ParticleID p3(741,13);
    
    ENSURE(!stream.insert(I3MCPE(p3,1.0,0.35)));
    ENSURE(!stream.insert(I3MCPE(p1,1.0,0.30)));
    ENSURE(stream.insert(I3MCPE(p2,1.0,1.3)));
    
    {
        // Check again that the state still makes sense
        std::vector<I3MCPE> hits;
        ParticlePulseIndexMap aux;
        std::tie(hits, aux) = stream.extractMCPEsWithPIDInfo();
        
        ENSURE(isTimeOrdered<I3MCPE>(hits));
        ENSURE_EQUAL(hits.size(), 6u,"Should have 6 pulses after merging");
        ENSURE(aux.size()==3,"Number of parent particles should be unchanged");
        ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
        ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
        // NB: there are 4 children here instead of 3 because we don't average
        // the hit times in the streaming case
        ENSURE_EQUAL(aux[p1].size(), 4u,"Particle 1 should have 4 children");
        ENSURE_EQUAL(aux[p2].size(), 3u,"Particle 2 should have 3 children");
        ENSURE_EQUAL(aux[p3].size(), 1u,"Particle 3 should have 1 child");

        uint32_t p1ExpectedSort[4]={0,1,2,3};
        uint32_t p2ExpectedSort[3]={2,4,5};
        uint32_t p3ExpectedSort[1]={1};

        ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1ExpectedSort)==std::make_pair(aux[p1].end(),p1ExpectedSort+4));
        ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2ExpectedSort)==std::make_pair(aux[p2].end(),p2ExpectedSort+3));
        ENSURE(std::mismatch(aux[p3].begin(),aux[p3].end(),p3ExpectedSort)==std::make_pair(aux[p3].end(),p3ExpectedSort+1));
    }
}

//Check that hit merging works properly:
//Hits with times with .2 ns should be merged
//The merged hit should have a weight equal to the sum of the weights of the original hits,
//and a time given by the weighted average of the times of the original hits.
//Parent particle information should also be correctly updated, in particular when two pulses
//with different parent particles are merged the index of the resulting pulse should appear
//in the index lists of both parents.
TEST(3_MCPulseMergingSortingMerging)
{
    std::vector<I3MCPulse> hits;
    ParticlePulseIndexMap aux;
    I3ParticleID p1(538,16);
    I3ParticleID p2(619,12);
    uint32_t hitIndex=0;

    //these two hits should be merged
    hits.push_back(I3MCPulse(0.0,1.0));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPulse(0.1,1.0));
    aux[p1].push_back(hitIndex++);

    //this hit should not merge with the preceding or following hits
    hits.push_back(I3MCPulse(0.25,1.0));
    aux[p1].push_back(hitIndex++);

    //these two hits should be merged
    //the merged pulse should be a child of both p1 and p2
    hits.push_back(I3MCPulse(0.5,1.0));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPulse(0.6,1.0));
    aux[p2].push_back(hitIndex++);

    //neither of these hits should merge with anything
    hits.push_back(I3MCPulse(0.75,1.0));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPulse(1.0,1.0));
    aux[p2].push_back(hitIndex++);

    MCHitMerging::timeMergeHits(hits,aux);

    ENSURE(hits.size()==5,"Should have 3 pulses after merging");
    ENSURE(aux.size()==2,"Number of parent particles should be unchanged");
    ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
    ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
    ENSURE(aux[p1].size()==4,"Particle 1 should have 4 children");
    ENSURE(aux[p2].size()==2,"Particle 2 should have 2 children");

    uint32_t p1Expected[4]={0,1,2,3};
    uint32_t p2Expected[2]={2,4};

    ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1Expected)==std::make_pair(aux[p1].end(),p1Expected+4));
    ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2Expected)==std::make_pair(aux[p2].end(),p2Expected+2));
    
    hitIndex=hits.size();
    
    I3ParticleID p3(741,13);
    hits.push_back(I3MCPulse(0.35,1.0));
    aux[p3].push_back(hitIndex++);
    hits.push_back(I3MCPulse(0.30,1.0));
    aux[p1].push_back(hitIndex++);
    hits.push_back(I3MCPulse(1.3,1.0));
    aux[p2].push_back(hitIndex++);
    
    MCHitMerging::sortMCHits(hits.begin(), hits.end(), aux);
    
    // Check again that the state still makes sense
    ENSURE(isTimeOrdered<I3MCPulse>(hits));
    ENSURE(hits.size()==8,"Should have 8 pulses after sorting");
    ENSURE(aux.size()==3,"Number of parent particles should be unchanged");
    ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
    ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
    ENSURE(aux[p1].size()==5,"Particle 1 should have 5 children");
    ENSURE(aux[p2].size()==3,"Particle 2 should have 3 children");
    ENSURE(aux[p3].size()==1,"Particle 3 should have 1 child");
    
    uint32_t p1ExpectedSort[5]={0,1,2,4,5};
    uint32_t p2ExpectedSort[3]={4,6,7};
    uint32_t p3ExpectedSort[1]={3};

    ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1ExpectedSort)==std::make_pair(aux[p1].end(),p1ExpectedSort+5));
    ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2ExpectedSort)==std::make_pair(aux[p2].end(),p2ExpectedSort+3));
    ENSURE(std::mismatch(aux[p3].begin(),aux[p3].end(),p3ExpectedSort)==std::make_pair(aux[p3].end(),p3ExpectedSort+1));

    MCHitMerging::timeMergeHits(hits, aux);

    ENSURE(isTimeOrdered<I3MCPulse>(hits));
    ENSURE(hits.size()==5,"Should have 6 pulses after sorting");
    ENSURE(aux.size()==3,"Number of parent particles should be unchanged");
    ENSURE(aux.find(p1)!=aux.end(),"particle 1 must appear in the output");
    ENSURE(aux.find(p2)!=aux.end(),"particle 2 must appear in the output");
    ENSURE(aux[p1].size()==3,"Particle 1 should have 3 children");
    ENSURE(aux[p2].size()==3,"Particle 2 should have 3 children");
    ENSURE(aux[p3].size()==1,"Particle 3 should have 1 child");

    uint32_t p1ExpectedSortMerge[3]={0,1,2};
    uint32_t p2ExpectedSortMerge[3]={2,3,4};
    uint32_t p3ExpectedSortMerge[1]={1};

    ENSURE(std::mismatch(aux[p1].begin(),aux[p1].end(),p1ExpectedSortMerge)==std::make_pair(aux[p1].end(),p1ExpectedSortMerge+3));
    ENSURE(std::mismatch(aux[p2].begin(),aux[p2].end(),p2ExpectedSortMerge)==std::make_pair(aux[p2].end(),p2ExpectedSortMerge+3));
    ENSURE(std::mismatch(aux[p3].begin(),aux[p3].end(),p3ExpectedSortMerge)==std::make_pair(aux[p3].end(),p3ExpectedSortMerge+1));
}
