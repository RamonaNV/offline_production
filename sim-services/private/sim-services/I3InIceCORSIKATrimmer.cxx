
#include <icetray/I3ConditionalModule.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/I3Constants.h>

/**
 * @brief Remove muons that have no chance of reaching the detector
 */
class I3InIceCORSIKATrimmer : public I3ConditionalModule {
public:
	I3InIceCORSIKATrimmer(const I3Context &ctx) : I3ConditionalModule(ctx)
	{
		AddOutBox("OutBox");
		AddParameter("MinEnergy", "Minimum surface energy of muons", 273*I3Units::GeV);
		AddParameter("FilterMode", "Drop events with no InIce particles", true);
		AddParameter("RemoveNeutrinos", "Remove neutrinos from the bundle", false);
	}
	
	void Configure()
	{
		double minEnergy;
		GetParameter("MinEnergy", minEnergy);
		minRange_ = muonRange_(minEnergy);
		log_debug_stream("Minimum energy " << minEnergy << " GeV; slant depth " << minRange_ << " meters water-equivalent");
		
		GetParameter("FilterMode", filterMode_);
		GetParameter("RemoveNeutrinos", dropNeutrinos_);
	}
	
	void DAQ(I3FramePtr frame)
	{
		I3MCTreePtr mctree(new I3MCTree(*frame->Get<I3MCTreeConstPtr>("I3MCTree")));
		
		unsigned inIceCount = 0;
		for (I3MCTree::pre_order_iterator p = mctree->begin(); p != mctree->end(); )
			if (Drop(*p)) {
				log_trace_stream("dropped " << p->GetEnergy() << " " << p->GetTypeString());
				p = mctree->erase(p);
			} else {
				if (p->GetLocationType() == I3Particle::InIce)
					inIceCount++;
				p++;
			}
		
		if (filterMode_ && inIceCount == 0)
			return;
		
		// Clean up newly childless parents
		for (I3MCTree::post_order_iterator p = mctree->begin_post(); p != mctree->end_post(); ) {
			if (p->GetLocationType() != I3Particle::InIce && mctree->number_of_children(p) == 0u) {
				p = mctree->erase(p);
			} else {
				p++;
			}
		}
		
		frame->Delete("I3MCTree");
		frame->Put("I3MCTree", mctree);
		PushFrame(frame);
	}
private:
	
	bool Drop(const I3Particle &p)
	{
		if (p.GetLocationType() == I3Particle::InIce) {
			if (p.GetType() == I3Particle::MuMinus ||
			    p.GetType() == I3Particle::MuPlus) {
				double cos_theta = std::cos(p.GetDir().GetZenith());
				return (cos_theta > 0) ?
				    muonRange_(p.GetEnergy()) < minRange_/cos_theta : true;
			} else if (p.IsNeutrino()) {
				return dropNeutrinos_;
			} else {
				return true;
			}
		}
		
		
		return false;
	}
	
	// Calculate the 99.9% range of muons in mwe
	// See: http://arxiv.org/abs/hep-ph/0407075 
	struct Range {
		Range(double a=0.212/1.2, double b=0.251e-3/1.2)
		: a_(a), b_(b) {}
		double a_, b_;
		double operator()(double energy) const
		{
			return std::log(1 + energy*b_/a_)/b_;
		}
	};
	
	double minRange_;
	bool dropNeutrinos_;
	bool filterMode_;
	Range muonRange_;
	
	SET_LOGGER("I3InIceCORSIKATrimmer");
};

I3_MODULE(I3InIceCORSIKATrimmer);
