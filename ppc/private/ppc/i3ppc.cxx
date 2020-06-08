#include <icetray/I3Frame.h>
#include <icetray/I3ConditionalModule.h>
#include <simclasses/I3MMCTrack.h>
#include <simclasses/I3MCPE.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/calibration/I3Calibration.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <phys-services/I3RandomService.h>
#include <clsim/I3CLSimStep.h>
#include <boost/foreach.hpp>

#include "ppc.h"

using namespace std;

class i3ppc : public I3ConditionalModule{
public:
  SET_LOGGER("i3ppc");

  i3ppc(const I3Context &context) : I3ConditionalModule(context),
				    ini(false),
				    verbose(false),
				    save_photons(false),
				    cyl(true),
				    keep_(false),
				    gpu(-1),
				    nph(0),
				    wid(0*I3Units::ns),
                                    efficiency_scaling_(1.),
                                    charge_(1.),
				    infoName_("")
  {
    tau_dnde_vec_.clear();
    AddParameter ("gpu", "GPU to use", gpu);
    AddParameter ("fla", "Flasher position", fla);
    AddParameter ("nph", "Number of photons", nph);
    AddParameter ("wid", "Flasher pulse width", wid);
    AddParameter ("MCTree", "MCTree to use", mct);
    AddParameter ("cyl", "use cylinder (1) or strict +300 m (0) detector volume", cyl);
    AddParameter ("keep", "Keep events that don't produce PEs",keep_);
    AddParameter ("verbose", "Print information messages", verbose);
    AddParameter ("photons", "Save photons that hit DOMs", save_photons);
    AddParameter ("tau_dnde_vec", "vector of pairs of luminescence decay time and dnde, tau in ns, dNdE in gamma/eV", tau_dnde_vec_); // fhl
    AddParameter ("infoName","Name of the ppc info dictionary. Will not be created if set to empty string",infoName_); // fhl
    AddParameter ("efficiency_scaling_factor","Multiplicative factor to modify DOM efficiency",efficiency_scaling_); // jc
    AddParameter ("charge", "Only for exotic particle with non-unitairy charge", charge_); // ward

    AddOutBox("OutBox");

    xppc::start();
    xppc::set_res(0.1f);  // 0.1 ns default resolution for hit binning

    Radius=0, Top=0, Bottom=0, Padding=300;
  }

  ~i3ppc(){
    xppc::stop();
  }

  void Configure(){
    GetParameter ("gpu", gpu);
    GetParameter ("fla", fla);
    GetParameter ("nph", nph);
    GetParameter ("wid", wid); wid/=I3Units::ns;
    GetParameter ("MCTree", mct);
    GetParameter ("cyl", cyl);
    GetParameter ("keep", keep_);
    GetParameter ("verbose", verbose);
    GetParameter ("photons", save_photons);
    GetParameter ("infoName", infoName_); // fhl
    GetParameter ("tau_dnde_vec", tau_dnde_vec_); // fhl
    GetParameter ("efficiency_scaling_factor", efficiency_scaling_); // jc
    GetParameter ("charge", charge_); // Ward
    fb=0, fe=0;

    if(!infoName_.empty()){
      if(verbose) cout<<"Building Info Dict with: "<<endl;
      // create info dict and save it // fhl
      ppcinfo = I3MapStringDoublePtr(new I3MapStringDouble);
      (*ppcinfo)["gpu"] = gpu;
      if(verbose) cout<<"gpu: "<< gpu << endl;
      (*ppcinfo)["nph"] = nph;
      if(verbose) cout<<"nph: "<< nph << endl;
      (*ppcinfo)["verbose"] = verbose;
      if(verbose) cout<<"verbose: "<< verbose << endl;
      (*ppcinfo)["wid"] = wid;
      if(verbose) cout<<"wid: "<< wid << endl;
      (*ppcinfo)["cyl"] = cyl;
      if(verbose) cout<<"cyl: "<< cyl << endl;
      (*ppcinfo)["keep"] = keep_;
      if(verbose) cout<<"keep: "<< keep_ << endl;
      (*ppcinfo)["num_tau_dnde"] = tau_dnde_vec_.size();
      if(verbose) cout<<"efficiency_scaling_factor: "<< efficiency_scaling_ << endl;
      (*ppcinfo)["efficiencyscaling"] = efficiency_scaling_;
      if(verbose) cout<<"charge: "<< charge_ << endl;
      (*ppcinfo)["charge"] = charge_;
      for(unsigned int i=0; i<tau_dnde_vec_.size(); ++i){
        std::stringstream id;
        id << i;
        std::string _tau = "tau_";
        std::string _dnde = "dnde_";
        _tau.append(id.str());
        _dnde.append(id.str());
        (*ppcinfo)[_tau] = tau_dnde_vec_[i].first;
        (*ppcinfo)[_dnde] = tau_dnde_vec_[i].second;
        if(verbose) cout<< (*ppcinfo)[_tau] << endl;
        if(verbose) cout<< (*ppcinfo)[_dnde] << endl;
      }
    }
  }

  void DAQ(I3FramePtr frame){
    log_trace("Entering i3ppc::DAQ()");
    frames.push_back(frame);

    if(!infoName_.empty()){
      frame->Put(infoName_, ppcinfo);
    }

    if(!ini){
      if(verbose) cout<<"Geometry frame"<<endl;
      geo = frame->Get<I3GeometryConstPtr>();
      I3CalibrationConstPtr calib = frame->Get<I3CalibrationConstPtr>();
      I3DetectorStatusConstPtr status = frame->Get<I3DetectorStatusConstPtr>();
      if(verbose){
	if(!geo) cout<<"    --> empty"<<endl;
	else cout<<"    --> contains "<<geo->omgeo.size()<<" entries"<<endl;
      }
      if(geo) if(geo->omgeo.size()>0){
	  for(I3OMGeoMap::const_iterator it=geo->omgeo.begin(); it!=geo->omgeo.end(); ++it){
	    OMKey omkey=it->first;
	    if(it->second.omtype==I3OMGeo::IceCube){
	      xppc::OM om;
	      om.str=omkey.GetString();
	      om.dom=omkey.GetOM();
	      const I3Position& pos = (it->second).position;
	      om.r[0]=pos.GetX()/I3Units::m;
	      om.r[1]=pos.GetY()/I3Units::m;
	      om.r[2]=pos.GetZ()/I3Units::m;
	      xppc::i3oms.push_back(om);

	      if(status){
		map<OMKey, I3DOMStatus>::const_iterator om_stat = status->domStatus.find(omkey);
		if(om_stat!=status->domStatus.end()){
		  double hv=om_stat->second.pmtHV/I3Units::V;
		  if(!isfinite(hv) || hv<0) hv=0;

		  double eff=1.0;
		  double spe_compensation_factor = 1.0;

		  if(calib){
		    map<OMKey, I3DOMCalibration>::const_iterator om_cal = calib->domCal.find(omkey);
		    if(om_cal!=calib->domCal.end()) {
			eff=om_cal->second.GetRelativeDomEff();

			SPEChargeDistribution spe_charge_dist = om_cal->second.GetCombinedSPEChargeDistribution();
			if (!std::isnan(spe_charge_dist.compensation_factor)) {
				spe_compensation_factor = spe_charge_dist.compensation_factor;
			} else {
			  log_debug("OM (%i/%u): spe compensation factor = %g (default (was: NaN))", om.str,om.dom,spe_compensation_factor);
			} 
		    }
		    if(!isfinite(eff) || eff<0) eff=1.0;
		  }

		  xppc::hvs[om]=hv;
		  xppc::rdes[om]=make_pair(eff*spe_compensation_factor*efficiency_scaling_, 0);
		}
	      }

	      double R=sqrt(om.r[0]*om.r[0]+om.r[1]*om.r[1]); if(R>Radius) Radius=R;
	      if(om.r[2]>Top) Top=om.r[2]; else if(om.r[2]<Bottom) Bottom=om.r[2];
	    }
	  }
	}

      {
	xppc::initialize(1.);
	xppc::choose(gpu);
	xppc::ini();

	if(nph>0) xppc::flset(fla.GetString(), fla.GetOM());
      }

      ini=true;
    }

    if(nph>0){
      xppc::sett(0, 0, 1, make_pair(-1, -1), fe);
      xppc::flone((unsigned long long) llroundf(nph));
    }
    else{
      for(I3Frame::const_iterator iter=frame->begin(); iter!=frame->end(); iter++){
	I3MMCTrackListConstPtr mmclist = boost::dynamic_pointer_cast<const I3MMCTrackList>(iter->second);
	if(mmclist){
	  if(verbose) cout<<"  --> This is an I3MMCTrackList: "<<iter->first<<endl;
	  for(I3MMCTrackList::const_iterator it=mmclist->begin(); it!=mmclist->end(); ++it){
	    const I3MMCTrack& m=*it;
	    const I3Particle& p=it->GetI3Particle();
	    i3mmctracks[make_pair(p.GetMinorID(), p.GetMajorID())]=&m;
	  }
	}
      }

      for(I3Frame::const_iterator iter=frame->begin(); iter!=frame->end(); iter++){
	string name = iter->first;
	if(mct.compare("")!=0 && name.compare(mct)!=0) continue;
	I3MCTreeConstPtr mctree = boost::dynamic_pointer_cast<const I3MCTree>(iter->second);
	if(mctree){
	  if(verbose) cout<<"  --> This is an I3MCTree: "<<iter->first<<endl;
	  for(I3MCTree::sibling_const_iterator it=mctree->begin(); it!=mctree->end_sibling(); ++it){
	    if(verbose) cout<<"extracting a particle "<<it->GetType()<<endl;
	    if (it->GetType() == I3Particle::Monopole && it->GetSpeed()/I3Constants::c > 0.995){
	      // leaving this sanity check ungated here as we do not get the MOPO symbol here
	      // and performance should not be impaired by it // fhl
	      log_fatal("Monopole velocities only implemented for v <= 0.995c. ");
	    }
	    pparticle(mctree, it);
	  }
	}
      }

      for(I3Frame::const_iterator iter=frame->begin(); iter!=frame->end(); iter++){
	I3CLSimStepSeriesConstPtr clstep = boost::dynamic_pointer_cast<const I3CLSimStepSeries>(iter->second);
	if(clstep){
	  if(verbose) cout<<"  --> This is an I3CLSimSteperies: "<<iter->first<<endl;
	  for(I3CLSimStepSeries::const_iterator it=clstep->begin(); it!=clstep->end(); ++it){
	    const I3CLSimStep& step = *it;
	    const I3Direction n = *(step.GetDir());
	    float nx=n.GetX(), ny=n.GetY(), nz=n.GetZ();
	    pair<int, unsigned long long> id(step.GetID(), 0ULL);
	    xppc::sett(nx, ny, nz, id, fe);

	    float ll=step.GetLength()/I3Units::m;

	    if(!isfinite(ll) || ll<0) ll=0;

	    float x=step.GetPosX()/I3Units::m;
	    float y=step.GetPosY()/I3Units::m;
	    float z=step.GetPosZ()/I3Units::m;
	    float t0=step.GetTime()/I3Units::ns;

	    unsigned long long num=step.GetNumPhotons();
	    float beta=step.GetBeta();

	    if(isinside(x, y, z, nx, ny, nz, ll)) xppc::addp_clst(x, y, z, t0, num, ll, beta);
	  }
	}
      }

      i3mmctracks.clear();
    }

    fe++;
    pushframe();
  }

  void Finish(){
    if(ini){
      xppc::eout();
      pushframe();
      popframes(fe);
      Flush();
      xppc::fin();
    }
  }

private:

  bool ini, verbose;
  bool save_photons;
  std::map<std::pair<int, unsigned long long>, const I3MMCTrack *> i3mmctracks;

  bool cyl;
  bool keep_;
  I3GeometryConstPtr geo;

  int gpu;
  OMKey fla;
  double nph;
  double wid;
  double efficiency_scaling_; // jc
  double charge_; // ward

  std::string mct;
  int fb, fe;
  std::deque<I3FramePtr> frames;

  std::string infoName_; // name for info dict // fhl
  I3MapStringDoublePtr ppcinfo; 
  I3Vector<std::pair<double,double> > tau_dnde_vec_; // vector of tau,dnde pairs

  class framedata{
  public:
    I3MCPESeriesMapPtr hits;
    std::map< std::pair<int, unsigned long long>, std::vector<I3Particle> > photons;

    framedata(){
      hits=I3MCPESeriesMapPtr(new I3MCPESeriesMap);
    }

    bool empty(){
      return
	photons.empty() && hits->empty();
    }
  };
  std::map<int, framedata> mcpes;

  double Radius, Top, Bottom, Padding;


  // ppc does not process all photons from one IceTray frame but collects photons until it can
  // populate the complete graphics card with photons. The photons are processed,
  // hits are written to the IceTray frames and then pushed back into the IceTray framework // fhl

  // please mind the capitalisation, pushframe pushes to the ppc internal frame stack,
  // PushFrame pushes a frame to the IceTray stack // fhl

  // pops all frames with an number lower than fx from the internal ppc stack to the IceTray stack // fhl
  void popframes(int fx){
    // fb is the internal count number of an event // fhl
    for(; fb<fx; fb++){
      bool empty = true;
      I3FramePtr frame=frames.front();
      map<int, framedata>::iterator mci=mcpes.find(fb);

      if(mci!=mcpes.end()){
	framedata& mc = mci->second;
	empty = mc.empty();

	if(mc.hits) frame->Put("MCPESeriesMap", mc.hits);

	if(save_photons){
	  map<string, I3MCTreePtr> tmp;
	  for(I3Frame::const_iterator iter=frame->begin(); iter!=frame->end(); iter++){
	    string name = iter->first;
	    if(mct.compare("")!=0 && name.compare(mct)!=0) continue;
	    I3MCTreeConstPtr mctree = boost::dynamic_pointer_cast<const I3MCTree>(iter->second);
	    if(mctree){
	      if(verbose) cout<<"  --> This is an I3MCTree: "<<iter->first<<endl;

	      I3MCTreePtr newtree(new I3MCTree(*mctree));
	      for(I3MCTree::post_order_iterator iter = newtree->begin_post(); iter!=newtree->end_post(); ++iter){
		const I3Particle& p = *iter;
		pair<int, unsigned long long> id=make_pair(p.GetMinorID(), p.GetMajorID());
		map< pair<int, unsigned long long>, vector<I3Particle> >::iterator photons = mc.photons.find(id);

		if(photons!=mc.photons.end()){
		  for(vector<I3Particle>::const_iterator j=photons->second.begin(); j!=photons->second.end(); ++j) newtree->append_child(iter, *j);
		  mc.photons.erase(photons);
		}
	      }

	      tmp[name]=newtree;
	    }
	  }

	  for(map<string, I3MCTreePtr>::iterator i=tmp.begin(); i!=tmp.end(); ++i){
	    string name=i->first;
	    frame->Delete(name);
	    frame->Put(name, i->second);
	  }

	  unsigned int num=0;
	  I3MCTreePtr rest(new I3MCTree());

	  for(map< pair<int, unsigned long long>, vector<I3Particle> >::const_iterator i = mc.photons.begin(); i!=mc.photons.end(); ++i){
	    vector<I3Particle> photons = i->second;
	    for(vector<I3Particle>::const_iterator j=photons.begin(); j!=photons.end(); ++j, num++) rest->insert_after(*j);
	  }

	  if(num>0) frame->Put("unattached", rest);
	}

	mcpes.erase(mci);
      }

      if(!empty || keep_){
	log_debug("Pushing the frame");
	PushFrame(frame,"OutBox");
      }
      frames.pop_front();
    }
  }

  void pushframe(){
    int fc=0;

    for(xppc::outz::const_iterator i=xppc::hitz.begin(); i!=xppc::hitz.end(); ++i){
      const xppc::ihit & h = i->first;
      OMKey key(h.omkey.str, h.omkey.dom);

      fc=h.track.frame;
      framedata& mc = mcpes[fc];

      I3MCPE hit(h.track.second, h.track.first);

      const vector<xppc::pout> & p = i->second;
      hit.npe = p.size();

      hit.time = (h.time+(wid>0?wid*xppc::xrnd():0))*I3Units::ns;
      (*(mc.hits))[key].push_back(hit);

      if(save_photons && geo && p.size()>0){
	vector<I3Particle> & photons = mc.photons[h.track];
	const I3Position& pos = geo->omgeo.find(key)->second.position;

	for(vector<xppc::pout>::const_iterator j=p.begin(); j!=p.end(); ++j){
	  I3Particle particle;

	  particle.SetType(I3Particle::CherenkovPhoton);
	  particle.SetLocationType(I3Particle::InIce);
	  particle.SetPos(pos.GetX()+j->r[0]*I3Units::m, pos.GetY()+j->r[1]*I3Units::m, pos.GetZ()+j->r[2]*I3Units::m);
	  particle.SetTime(j->r[3]*I3Units::ns);
	  particle.SetLength(0);
	  particle.SetDir(j->n[0], j->n[1], j->n[2]);
	  particle.SetEnergy((1239.84193/j->n[3])*I3Units::eV);

	  photons.push_back(particle);
	}
      }
    }
    xppc::efin();
    popframes(fc);
  }

  bool isinside(double x, double y, double z){
    // this function checks if a cascade happens close to the geometry of IceCube
    // or wont make any light at all, since then it is not worth processing

    if(cyl){ // use a cylinder along the track # default
      return sqrt(x*x+y*y)<=Radius+Padding && Bottom-Padding<=z && z<=Top+Padding;
    }
    else{ // use a strict 300m detector volume # not default
      for(vector<xppc::OM>::const_iterator i=xppc::i3oms.begin(); i!=xppc::i3oms.end(); ++i){
	double dx=i->r[0]-x;
	double dy=i->r[1]-y;
	double dz=i->r[2]-z;
	if(sqrt(dx*dx+dy*dy+dz*dz)<=Padding) return true;
      }
      return false;
    }
  }

  bool isinside(double x, double y, double z, double nx, double ny, double nz, double l){
    // this function checks if a track passes through the geometry of IceCube
    // or wont make any light at all, since then it is not worth processing
    if(cyl){ // use a cylinder along the track # default
      double R=Radius+Padding, B=Bottom-Padding, T=Top+Padding;
      double a=nx*nx+ny*ny, c=R*R-(x*x+y*y), lo=-1, hi=l+1;

      double h1, h2;
      if(nz>0) h1=(B-z)/nz, h2=(T-z)/nz;
      else if(nz<0) h1=(T-z)/nz, h2=(B-z)/nz;
      else{
	if(z<B || z>T) h1=hi, h2=lo;
	else h1=0, h2=l;
      }

      double r1, r2;
      if(a>0){
	double b=x*nx+y*ny;
	double D=b*b+a*c;
	if(D>=0){
	  D=sqrt(D);
	  r1=(-b-D)/a, r2=(-b+D)/a;
	}
	else r1=hi, r2=lo;
      }
      else{
	if(c<0) r1=hi, r2=lo;
	else r1=0, r2=l;
      }

      double t1=max(r1, h1), t2=min(r2, h2);
      if(t1<0) t1=0; if(t2>l) t2=l;
      return t2>=t1;
    }
    else{ // use a strict 300m detector volume # not default
      for(vector<xppc::OM>::const_iterator i=xppc::i3oms.begin(); i!=xppc::i3oms.end(); ++i){
	double dx=i->r[0]-x;
	double dy=i->r[1]-y;
	double dz=i->r[2]-z;
	double r=dx*nx+dy*ny+dz*nz;
	if(r<0) r=0; else if(r>l) r=l;
	dx-=r*dx, dy-=r*dy, dz-=r*dz;
	if(sqrt(dx*dx+dy*dy+dz*dz)<=Padding) return true;
      }
      return false;
    }
  }

  int iType(const I3Particle& p){
    if(p.GetLocationType()!=I3Particle::InIce) return -4;
    else if(p.GetShape()==I3Particle::Dark) return -3;
    else{
      switch(p.GetType()){
      default:
	return -2;
      case I3Particle::Monopole: // obi
	return -41; // this is the i3particle enum number! // obi
   // Ward (needs implementation in I3Particle)
   // case I3Particle::SMPPlus:
   // case I3Particle::SMPMinus:
   //   return -10;
      case I3Particle::MuPlus:
      case I3Particle::MuMinus:
	return -1;
      case I3Particle::DeltaE:
      case I3Particle::Brems:
      case I3Particle::PairProd:
	return 1;
      case I3Particle::EMinus:
	return 2;
      case I3Particle::EPlus:
	return 3;
      case I3Particle::Gamma:
      case I3Particle::Pi0: // Pi0 decays into two gammas before it has a chance to interact, so pure EM cascades are produced
	return 4;
      case I3Particle::NuclInt:
      case I3Particle::Hadrons:
      case I3Particle::KPlus:
      case I3Particle::KMinus:
      case I3Particle::K0_Short:
      case I3Particle::Eta:
      case I3Particle::Lambda:
      case I3Particle::SigmaPlus:
      case I3Particle::Sigma0:
      case I3Particle::SigmaMinus:
      case I3Particle::Xi0:
      case I3Particle::XiMinus:
      case I3Particle::OmegaMinus:
      case I3Particle::NeutronBar:
      case I3Particle::LambdaBar:
      case I3Particle::SigmaMinusBar:
      case I3Particle::Sigma0Bar:
      case I3Particle::SigmaPlusBar:
      case I3Particle::Xi0Bar:
      case I3Particle::XiPlusBar:
      case I3Particle::OmegaPlusBar:
      case I3Particle::DPlus:
      case I3Particle::DMinus:
      case I3Particle::D0:
      case I3Particle::D0Bar:
      case I3Particle::DsPlus:
      case I3Particle::DsMinusBar:
      case I3Particle::LambdacPlus:
      case I3Particle::WPlus:
      case I3Particle::WMinus:
      case I3Particle::Z0: // listing all hadrons explicitly, no nucleus here 
	return 101;
      case I3Particle::PiPlus:
	return 102;
      case I3Particle::PiMinus:
	return 103;
      case I3Particle::K0_Long:
	return 104;
      case I3Particle::PPlus:
	return 105;
      case I3Particle::Neutron:
	return 106;
      case I3Particle::PMinus:
	return 107;
      }
    }
  }

  void pparticle(I3MCTreeConstPtr tree, I3MCTree::const_iterator sub){
    const I3Particle& p = *sub;
    pair<int, unsigned long long> id=make_pair(p.GetMinorID(), p.GetMajorID());

    double E0=p.GetEnergy()/I3Units::GeV;
    double t0=p.GetTime()/I3Units::ns;
    double ll=p.GetLength()/I3Units::m;

    if(!isfinite(ll) || ll<0) ll=0;

    double x=p.GetX()/I3Units::m;
    double y=p.GetY()/I3Units::m;
    double z=p.GetZ()/I3Units::m;

    const I3Direction& n = p.GetDir();
    double nx=n.GetX(), ny=n.GetY(), nz=n.GetZ();
    xppc::sett(nx, ny, nz, id, fe);

    int ptype = iType(p);

    if(ptype>=0){ // cascade // obi
      if(isfinite(E0)) if(isinside(x, y, z)) xppc::addp(x, y, z, t0, E0, ptype);
    }
    else if (ptype == -41){ // monopoles // obi
      // checks if track segment is near any DOM, use cascade function because track segment is very short < 10m
      if (isinside(x, y, z, nx, ny, nz, ll)){ // use the function for tracks to check if the track passes close to the geometry
	double beta = p.GetSpeed()/I3Constants::c;
	addp_monopole(x, y, z, t0, E0, ll, beta);
      }
    }
    else if(ptype==-1 || ptype==-10){ // track muon // obi // track SMP // Ward
      float scale=ptype==-10?charge_*charge_:1.f;

      // cout << "muon" << endl;
      const double sol=0.299792458;

      map<pair<int, unsigned long long>, const I3MMCTrack *>::const_iterator mmctrack=i3mmctracks.find(id);
      if(mmctrack==i3mmctracks.end()){
	if(isinside(x, y, z, nx, ny, nz, ll)) xppc::addp(x, y, z, t0, E0, (float) ll, scale);
      }
      else{
	const I3MMCTrack * i3mmctrack = mmctrack->second;
	double ti=i3mmctrack->GetTi()/I3Units::ns;
	double Ei=i3mmctrack->GetEi()/I3Units::GeV;
	double tf=i3mmctrack->GetTf()/I3Units::ns;
	double Ef=i3mmctrack->GetEf()/I3Units::GeV;
	if(Ei==0){
	  ti=t0;
	  Ei=E0;
	}
	if(Ef<=0){
	  tf=t0+ll/sol;
	  Ef=0;
	}
	if(Ei>0){
	  double totE=0;
	  for(I3MCTree::sibling_const_iterator it=tree->children(sub); it!=tree->end(sub); ++it){
	    double t=it->GetTime()/I3Units::ns;
	    if(t>=ti && t<=tf) totE+=it->GetEnergy()/I3Units::GeV;
	  }
	  double sl=tf>ti?((Ef-Ei+totE)/(ti-tf)):0;

	  double dr=(ti-p.GetTime()/I3Units::ns)*sol;
	  x+=dr*nx, y+=dr*ny, z+=dr*nz;

	  double t=ti;
	  double l;
	  double E=Ei;
	  double new_t;

	  for(I3MCTree::sibling_const_iterator it=tree->children(sub); it!=tree->end(sub); ++it){
	    new_t=it->GetTime()/I3Units::ns;
	    if(new_t>=ti && new_t<=tf){
	      l=sol*(new_t-t);

	      if(isinside(x, y, z, nx, ny, nz, l)) xppc::addp(x, y, z, t, E, (float) l, scale);

	      x=it->GetX()/I3Units::m;
	      y=it->GetY()/I3Units::m;
	      z=it->GetZ()/I3Units::m;
	      double e=it->GetEnergy()/I3Units::GeV;

	      E-=e+sl*max(new_t-t, 0.);
	      t=new_t;
	    }
	  }
	  l=sol*(tf-t);
	  if(isinside(x, y, z, nx, ny, nz, l)) xppc::addp(x, y, z, t, E, (float) l, scale);
	}
      }
    }

    for(I3MCTree::sibling_const_iterator it=tree->children(sub); it!=tree->end(sub); ++it) pparticle(tree, it);
  }


  double CalculateEnergyLoss(double beta, bool cDensityCorrection){
    // taken from monopole_propagator module

    const double MASS_E = 5.11e-4;  // GeV/c^2
    const double ION_LOSS = 7.5e-8; // GeV
    const double DENSITY=9.17e5;    // g/m^3
    const double Z_OVER_A=0.55509;  // mol/g
    const double CHARGE_E_SQRD=1.439964459e-18; // GeV*m
    // Values given for g=137e/2
    const double QED_CORR=0.406;
    const double BLOCH_CORR=0.248;
    const double pi=3.141592653589793;
    const double NA=6.022140857e23; // mol^-1

    const double COEFF1=(4*pi*NA*DENSITY*Z_OVER_A);     // m^-3
    const double COEFF2=pow(0.5*137.0*CHARGE_E_SQRD,2); // GeV^2 m^2
    const double COEFF3=MASS_E;                         // GeV
    const double COEFF=COEFF1*COEFF2/COEFF3;            // GeV/m

    double densityCorrection = 0.0;
    if(cDensityCorrection){
      // densityCorrection=CalculateDensityCorrection(beta);
      // Constants for Density Correction
      const double A_DC=0.09116;
      const double X0=0.240;
      const double X1=2.8004;
      const double M_DC=3.477;
      const double C_BAR=3.5017;

      double denCorr=0.0;
      double gamma = 1/(sqrt(1-pow(beta,2)));
      double X=log10(beta*gamma);

      if (X>X0&&X<X1){
	denCorr=log(pow(beta*gamma,2))-C_BAR+A_DC*pow((X1-X),M_DC);
      }
      else if (X>X1){
	denCorr=log(pow(beta*gamma,2))-C_BAR;
      }
      densityCorrection= denCorr;
    }

    double dEdX = COEFF*(log((2*MASS_E*pow(beta,2))/(ION_LOSS*(1-pow(beta,2))))+(QED_CORR/2)-0.5-(densityCorrection/2)-BLOCH_CORR);

    return dEdX*1000*1000*1000; // GeV/m -> eV/m
  }

  inline float directCherenkovLightUnnormalizedyield(double beta){
    static const float ppm2=32359.1381; // obi's calculation of the max muon number from Frank Tamm (use this for normalization of obi's values)
    // calculated by Obi, checked by Alex
    // only direct Cherenkov light
    if (beta > 0.758){ // above cherenkov threshold
      // calculation see Thesis Anna Pollmann
      float a=624464828.676, b=1.74108025; // parameters of Frank-Tamm // calculated by Obi, checked by Alex
      return (a * (1-1/( beta*beta * b ))) / ppm2; // Frank-Tamm-Monopole(beta) [photons/meter] / ppm
    }
    return 0;
  }

  inline float indirectCherenkovLightUnnormalizedyield(double beta){
    if (beta >= 0.48) { //indirect Cherenkov threshold
      // Spline Fit, see Thesis Anna Pollmann
      float a=0, b=0, c=0, d=0;
      if (beta < 0.51){ // and beta >= 0.48
	a=1.14630633e-10; b=4.85359497e+01; c=-1.45995816e+00; d=-2.58433114e-01;
      }else if (beta < 0.61){
	a=5.05592842e-03; b=1.47862510e+01; c=-4.20203712e-01; d=-5.30134659e+00;
      }else if (beta < 0.91){
	a=0.70428599; b=5.85816753; c=1.02337763; d=-48.55463142;
      }else if (beta < 0.96){
	a=4.65761093e-04; b=1.36985032e+01; c=4.13888337e-01; d=1.77325577e+02;
      }else if (beta >= 0.96){
	a=2.98686931e-13; b=3.50807192e+01; c=6.12915941e-02; d=4.08217397e+02;
	// if beta > 0.995 c this is a pessimistic estimation of the indirect CHerenkov light since this is increasing heavily after this point
      }else{
	cout << "FATAL: Wrong beta!" << endl; // this does never happen atm // fhl
      }
      return (a * expf((b*beta)+c))+d; // add this to the direct CL if any
    }
    return 0;
  }

  inline float luminescenceLightUnnormalizedyield(double beta, double dnde){
    static const float ppm2=32359.1381; // obi's calculation of the max muon number from Frank Tamm (use this for normalization of obi's values)
    // calculated by Obi, checked by Alex
    double enloss=CalculateEnergyLoss(beta,true); // eV/m
    return enloss*dnde/ppm2; // devided by muon photons per meter
  }

  inline double setParticleProbabilities(double directCherenkov, double indirectCherenkov, double luminescence){
    // function to set the probabilities direct Cherenkov, indirect Cherenkov, and luminescence light
    double extr=directCherenkov + indirectCherenkov;
    return directCherenkov/extr; // indirect = 0 => direct/direct = 1
    // extr+=luminescence;
    // p.p_l = luminescence / extr; // old propability based code
  }

  void addp_monopole(float rx, float ry, float rz, float t, float E, float length, float beta){
    float tau=0, fr=1;

    float directCherenkov, indirectCherenkov, luminescence, extr;
    directCherenkov = directCherenkovLightUnnormalizedyield(beta);
    indirectCherenkov = indirectCherenkovLightUnnormalizedyield(beta);

    luminescence = 0;
    for(I3Vector<std::pair<double,double> >::const_iterator it=tau_dnde_vec_.begin(); it!=tau_dnde_vec_.end(); ++it){
      float tmp = luminescenceLightUnnormalizedyield(beta, it->second);

      // add luminescence photons
      tau = it->first;
      xppc::addp_mopo(rx, ry, rz, t, tmp*length, length, beta, tau, fr);

      luminescence += tmp;
    }

    fr=setParticleProbabilities(directCherenkov, indirectCherenkov, luminescence);
    extr = directCherenkov + indirectCherenkov;

    // add Cherenkov light
    tau = 0; // Cherenkov is emited promptly
    xppc::addp_mopo(rx, ry, rz, t, extr*length, length, beta, tau, fr);
  }

};

I3_MODULE(i3ppc);
