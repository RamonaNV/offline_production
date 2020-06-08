
#include "ppc/I3CLSimStepToPhotonConverterPPC.h"
#include "ppc.h"

boost::weak_ptr<I3CLSimStepToPhotonConverterPPC> I3CLSimStepToPhotonConverterPPC::instance_;

boost::shared_ptr<I3CLSimStepToPhotonConverterPPC> I3CLSimStepToPhotonConverterPPC::GetInstance()
{
    auto shared = instance_.lock();
    if (!shared) {
        shared = boost::shared_ptr<I3CLSimStepToPhotonConverterPPC>(new I3CLSimStepToPhotonConverterPPC);
        instance_ = shared;
    }
    
    return shared;
}

I3CLSimStepToPhotonConverterPPC::I3CLSimStepToPhotonConverterPPC() :
    isInitialized_(false), workgroupSize_(0), maxBunchSize_(0), inbox_(1)
{
    xppc::start();
    xppc::set_res(0.1f);  // 0.1 ns default resolution for hit binning
}

I3CLSimStepToPhotonConverterPPC::~I3CLSimStepToPhotonConverterPPC()
{
    xppc::stop();
}

void I3CLSimStepToPhotonConverterPPC::Initialize()
{
    if (!geometry_)
        log_fatal("Geometry has not yet been configured");
    if (!wavelenthBias_)
        log_fatal("Wavelength bias has not yet been configured!");
    
    if (!xppc::i3oms.empty())
        log_fatal("There can be only one ppc");
    for (size_t i=0; i < geometry_->size(); i++) {
        xppc::OM om;
        om.str=geometry_->GetStringID(i);
        om.dom=geometry_->GetDomID(i);
        om.r[0]=geometry_->GetPosX(i)/I3Units::m;
        om.r[1]=geometry_->GetPosY(i)/I3Units::m;
        om.r[2]=geometry_->GetPosZ(i)/I3Units::m;
        xppc::i3oms.push_back(om);
        
        // Setting PMT high voltage > 0, RDE to 1.0, and om type to 1, and
        // as.dat to 1\n1, prevents xppc::print() from down-sampling photons
        // to account for wavelength and angular efficiency.
        xppc::hvs[om]=1200;
        xppc::rdes[om]=std::make_pair(1.0, 1);
    }

    int gpu(-1);
    xppc::initialize(1.);
    
    xppc::choose(gpu);
    xppc::ini();

    // Harmonize bunch sizes
    maxBunchSize_ = xppc::getMaxBunchSize();
    workgroupSize_ = xppc::getWorkgroupSize();

    isInitialized_ = true;
}

void I3CLSimStepToPhotonConverterPPC::EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier)
{
    inbox_.Put(std::make_pair(steps, identifier));
}

I3CLSimStepToPhotonConverterPPC::ConversionResult_t I3CLSimStepToPhotonConverterPPC::GetConversionResult()
{
    I3CLSimStepSeriesConstPtr steps;
    uint32_t identifier;
    std::tie(steps, identifier) = inbox_.Get();
    
    // sends steps to kernel
    for (const I3CLSimStep &step : *steps) {
        const I3Direction n = *(step.GetDir());
        float nx=n.GetX(), ny=n.GetY(), nz=n.GetZ();
        std::pair<int, unsigned long long> id(static_cast<int>(step.GetID()), 0ULL);
        xppc::sett(nx, ny, nz, id, static_cast<int>(identifier));

        float ll=step.GetLength()/I3Units::m;

        if(!std::isfinite(ll) || ll<0) ll=0;

        float x=step.GetPosX()/I3Units::m;
        float y=step.GetPosY()/I3Units::m;
        float z=step.GetPosZ()/I3Units::m;
        float t0=step.GetTime()/I3Units::ns;

        unsigned long long num=step.GetNumPhotons();
        float beta=step.GetBeta();

        xppc::addp_clst(x, y, z, t0, num, ll, beta);
    }
    // Flush the kernel
    xppc::eout();
    
    ConversionResult_t result(0, boost::make_shared<I3CLSimPhotonSeries>(), NULL);
    
    // collect results
    for (auto &hit : xppc::hitz) {
        if (result.identifier == 0)
            result.identifier = hit.first.track.frame;
        else {
            if ( result.identifier != static_cast<uint32_t>(hit.first.track.frame) )
                log_fatal_stream("Got "<<static_cast<uint32_t>(hit.first.track.frame)<<" after "<<result.photons->size()<<" photons from bunch "<<result.identifier);
            assert( int(result.identifier) == hit.first.track.frame );
        }
        
        for (auto &ppc_photon : hit.second) {
            I3CLSimPhoton photon;
            photon.SetID(hit.first.track.first);
            photon.SetStringID(hit.first.omkey.str);
            photon.SetOMID(hit.first.omkey.dom);
            
            // NB: photon positions are relative to the hit DOM
            photon.SetPosX(ppc_photon.r[0]*I3Units::m);
            photon.SetPosY(ppc_photon.r[1]*I3Units::m);
            photon.SetPosZ(ppc_photon.r[2]*I3Units::m);
            photon.SetTime(ppc_photon.r[3]*I3Units::ns);
            photon.SetDir(ppc_photon.n[0],ppc_photon.n[1],ppc_photon.n[2]);
            photon.SetWavelength(ppc_photon.n[3]*I3Units::nanometer);
            photon.SetWeight(1.0/wavelenthBias_->GetValue(photon.GetWavelength()));
            
            result.photons->push_back(photon);
        }
    }
    // clear output buffer
    xppc::efin();
    
    return result;
}
