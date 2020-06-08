
#ifndef PPC_I3CLSIMSTEPTOPHOTONCONVERTERPPC_H_INCLUDED
#define PPC_I3CLSIMSTEPTOPHOTONCONVERTERPPC_H_INCLUDED

#include "clsim/I3CLSimStepToPhotonConverter.h"
#include "clsim/I3CLSimQueue.h"

class I3CLSimStepToPhotonConverterPPC : public I3CLSimStepToPhotonConverter {
public:
    
   I3CLSimStepToPhotonConverterPPC(const I3CLSimStepToPhotonConverterPPC &) = delete;
   I3CLSimStepToPhotonConverterPPC(const I3CLSimStepToPhotonConverterPPC &&) = delete;
   
   virtual ~I3CLSimStepToPhotonConverterPPC();
   
   /** Get the global shared instance */
   static boost::shared_ptr<I3CLSimStepToPhotonConverterPPC> GetInstance();
   
   /**
    * Sets the wavelength generators. 
    * The first generator (index 0) is assumed to return a Cherenkov
    * spectrum that may have a bias applied to it. This bias factor
    * needs to be set using SetWlenBias().
    * All other generator indices are assumed to be for flasher/laser
    * light generation. During generation, no Cherenkov angle
    * rotation will be applied to those photons with indices >= 1.
    * Will throw if used after the call to Initialize().
    */
   virtual void SetWlenGenerators(const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators) { log_fatal("unimplemented"); };

   /**
    * Sets the wavelength weights. Set this to a constant value
    * of 1 if you do not need biased photon generation.
    * The wavelength spectrum set with SetWlenGenerator()
    * is assumed to have a biassing factor already applied to it.
    * This call sets this factor in order to be able to assign
    * correct weights.
    * Will throw if used after the call to Initialize().
    */
   virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias) { wavelenthBias_=wlenBias; };

   /**
    * Sets the medium properties.
    * Will throw if used after the call to Initialize().
    */
   virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties) { log_fatal("unimplemented"); };

   /**
    * Sets the geometry.
    * Will throw if used after the call to Initialize().
    */
   virtual void SetGeometry(I3CLSimSimpleGeometryConstPtr geometry) { geometry_ = geometry; };
   
   /**
    * Initializes the simulation.
    * Will throw if already initialized.
    */
   virtual void Initialize();

   /**
    * Returns true if initialized.
    * Never throws.
    */
   virtual bool IsInitialized() const { return isInitialized_; };
   
   /**
    * Adds a new I3CLSimStepSeries to the queue.
    * The resulting I3CLSimPhotonSeries can be retrieved from the
    * I3CLSimStepToPhotonConverter after some processing time.
    *
    * Enqueuing a vector after calling EnqueueBarrier 
    * will throw if not all photons have been retrieved.
    *
    * Will throw if not initialized.
    */
   virtual void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier);

   /**
    * Gets the current workgroup size.
    */
   virtual std::size_t GetWorkgroupSize() const { return workgroupSize_; };

   /**
    * Gets the number of parallel work items.
    */
   virtual std::size_t GetMaxNumWorkitems() const { return maxBunchSize_; };

   /**
    * Reports the current queue size. The queue works asynchronously,
    * so this value will probably have changed once you use it.
    *
    * Will throw if not initialized.
    */
   virtual std::size_t QueueSize() const { return inbox_.size(); }; 
   
   /**
    * Returns true if more photons are available.
    * If the return value is false, the current simulation is finished
    * and a new step vector may be set.
    * 
    * Will throw if not initialized.
    */
   virtual bool MorePhotonsAvailable() const  { log_fatal("unimplemented"); };

   /**
    * Returns a bunch of photons stored in a vector<I3CLSimPhoton>.
    *
    * The return value is a pair<uint, vector<I3CLSimPhoton> >.
    * The integer is the same identifier as specified in the call
    * to EnqueueSteps().
    *
    * Might block if no photons are available.
    * 
    * Will throw if not initialized.
    */
   virtual ConversionResult_t GetConversionResult();
   
   virtual std::map<std::string, double> GetStatistics() const { return std::map<std::string, double>(); };

private:
    I3CLSimStepToPhotonConverterPPC();
    static boost::weak_ptr<I3CLSimStepToPhotonConverterPPC> instance_;
    
    bool isInitialized_;
    I3CLSimSimpleGeometryConstPtr geometry_;
    I3CLSimFunctionConstPtr wavelenthBias_;
    
    size_t workgroupSize_, maxBunchSize_;
    
    I3CLSimQueue<std::pair<I3CLSimStepSeriesConstPtr, uint32_t>> inbox_;
};

#endif
