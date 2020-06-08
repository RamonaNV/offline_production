#include <I3Test.h>
#include <boost/make_shared.hpp>
#include <icetray/I3Tray.h>
#include <icetray/I3ServiceFactory.h>
#include <simclasses/I3ParticleIDMap.hpp>

#include <sim-services/I3PrimaryPulseMapper.h>

template<typename T>
class param_setter {
	T& owner;
	I3Configuration& config;
	bool final;
public:
	param_setter<T>(T& t, I3Configuration& c):owner(t),config(c),final(true){}
	param_setter<T>(param_setter<T>&& prev):owner(prev.owner),config(prev.config),final(true){
		prev.final=false;
	}
	~param_setter<T>(){
		if(final)
			owner.FinalizeConfiguration();
	}
	template <typename U>
	param_setter& operator()(const std::string& paramName, U value){
		config.Set(paramName, boost::python::object(value));
		return(*this);
	}
};

template<typename ModuleType>
class SingleModuleTestSetup{
private:
	boost::shared_ptr<std::deque<boost::shared_ptr<I3Frame>>> inbox;
	I3Context& context;
	boost::shared_ptr<ModuleType> module;
	struct CollectorModule : public I3Module{
		CollectorModule(I3Context& ctx):I3Module(ctx){}
		boost::shared_ptr<I3Frame> GetFrame(){ return(PopFrame()); }
	};
	boost::shared_ptr<CollectorModule> collector;
public:
	SingleModuleTestSetup(I3Context& context):
	inbox(new std::deque<boost::shared_ptr<I3Frame> >),
	context(context),
	module(new ModuleType(context)),
	collector(new CollectorModule(context))
	{
		module->configuration_.ClassName(I3::name_of(typeid(*module)));
		module->inbox_=inbox;
		module->ConnectOutBox("OutBox",collector);
	}
	
	param_setter<SingleModuleTestSetup<ModuleType>> Configure(){
        return(param_setter<SingleModuleTestSetup<ModuleType> >(*this,module->configuration_));
    }
	void FinalizeConfiguration(){
		module->Configure_();
	}
	boost::shared_ptr<I3Frame> processFrame(boost::shared_ptr<I3Frame> frame){
		if(!module)
			throw std::runtime_error("Test module not configured, or configuration not finalized");
		inbox->push_front(frame);
		module->Process_();
		return(collector->GetFrame());
	}
};

TEST_GROUP(I3PrimaryPulseMapper);

TEST(MapPulses){
	boost::shared_ptr<I3Frame> frame(new I3Frame(I3Frame::DAQ));
	
	boost::shared_ptr<I3MCTree> tree(new I3MCTree);
	I3Particle p1, p2; //primaries
	tree->insert(p2);
	tree->insert(p1);
	I3Particle p1c1, p1c2;
	tree->append_child(p1,p1c1);
	tree->append_child(p1,p1c2);
	I3Particle p1gc1, p1gc2;
	tree->append_child(p1c1,p1gc1);
	tree->append_child(p1c1,p1gc2);
	I3Particle p2c1, p2c2;
	tree->append_child(p2,p2c1);
	tree->append_child(p2,p2c2);
	frame->Put("I3MCTree",tree);
	
	boost::shared_ptr<I3ParticleIDMap> in_map(new I3ParticleIDMap);
	//stuff from the p1 family
	(*in_map)[OMKey(1,1)][p1c1].push_back(0);
	(*in_map)[OMKey(1,1)][p1gc1].push_back(1);
	(*in_map)[OMKey(1,1)][p1c2].push_back(2);
	(*in_map)[OMKey(1,1)][p1gc2].push_back(3);
	//stuff from the p2 family
	(*in_map)[OMKey(1,2)][p2].push_back(0);
	(*in_map)[OMKey(1,2)][p2c1].push_back(1);
	(*in_map)[OMKey(1,2)][p2c2].push_back(2);
	//mixed stuff
	(*in_map)[OMKey(1,3)][p1c1].push_back(0);
	(*in_map)[OMKey(1,3)][p1c2].push_back(0);
	(*in_map)[OMKey(1,3)][p2c1].push_back(0);
	(*in_map)[OMKey(1,3)][p2c2].push_back(0);
	frame->Put("I3MCPulseSeriesMapParticleIDMap",in_map);
	
	I3Context context;
	SingleModuleTestSetup<I3PrimaryPulseMapper> tester(context);
	frame=tester.processFrame(frame);
	
	auto out_map=frame->Get<boost::shared_ptr<const I3ParticleIDMap>>("I3MCPulseSeriesMapPrimaryIDMap");
	
	ENSURE(!!out_map, "Output should exist");
	ENSURE_EQUAL(out_map->size(),3u, "Output should have same entries as input");
	ENSURE(out_map->count(OMKey(1,1)), "Output should have same entries as input");
	ENSURE(out_map->count(OMKey(1,2)), "Output should have same entries as input");
	ENSURE(out_map->count(OMKey(1,3)), "Output should have same entries as input");
	
	{
		const ParticlePulseIndexMap& ppim=out_map->find(OMKey(1,1))->second;
		for(const auto& particle_entry : ppim){
			std::cout << particle_entry.first << ": ";
			for(uint32_t idx : particle_entry.second)
				std::cout << idx << ' ';
			std::cout << std::endl;
		}
		ENSURE_EQUAL(ppim.size(),1u,"Only primary 1 should appear for OM(1,1)");
		ENSURE(ppim.count(p1),"Only primary 1 should appear for OM(1,1)");
		ENSURE_EQUAL(ppim.find(p1)->second.size(),4u,"All pe on OM(1,1) should be associated with primary 1");
		ENSURE(std::is_sorted(ppim.find(p1)->second.begin(),ppim.find(p1)->second.end()),"P.E. indices should be sorted");
	}
	{
		const ParticlePulseIndexMap& ppim=out_map->find(OMKey(1,2))->second;
		for(const auto& particle_entry : ppim){
			std::cout << particle_entry.first << ": ";
			for(uint32_t idx : particle_entry.second)
				std::cout << idx << ' ';
			std::cout << std::endl;
		}
		ENSURE_EQUAL(ppim.size(),1u,"Only primary 2 should appear for OM(1,2)");
		ENSURE(ppim.count(p2),"Only primary 2 should appear for OM(1,2)");
		ENSURE_EQUAL(ppim.find(p2)->second.size(),3u,"All pe on OM(1,2) should be associated with primary 2");
		ENSURE(std::is_sorted(ppim.find(p2)->second.begin(),ppim.find(p2)->second.end()),"P.E. indices should be sorted");
	}
	{
		const ParticlePulseIndexMap& ppim=out_map->find(OMKey(1,3))->second;
		for(const auto& particle_entry : ppim){
			std::cout << particle_entry.first << ": ";
			for(uint32_t idx : particle_entry.second)
				std::cout << idx << ' ';
			std::cout << std::endl;
		}
		ENSURE_EQUAL(ppim.size(),2u,"Both primaries should appear for OM(1,3)");
		ENSURE(ppim.count(p1),"Primary 21should appear for OM(1,3)");
		ENSURE(ppim.count(p2),"Primary 2 should appear for OM(1,3)");
		ENSURE_EQUAL(ppim.find(p1)->second.size(),1u,"The pe on OM(1,3) should be associated with primary 1");
		ENSURE_EQUAL(ppim.find(p2)->second.size(),1u,"The pe on OM(1,3) should be associated with primary 2");
		ENSURE(std::is_sorted(ppim.find(p1)->second.begin(),ppim.find(p1)->second.end()),"P.E. indices should be sorted");
		ENSURE(std::is_sorted(ppim.find(p2)->second.begin(),ppim.find(p2)->second.end()),"P.E. indices should be sorted");
	}
}
