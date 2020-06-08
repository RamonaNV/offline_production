#ifndef I3PROPAGATORSERVICE_H
#define I3PROPAGATORSERVICE_H

#include <map>

#include <icetray/I3PointerTypedefs.h>
#include <icetray/I3Frame.h>
#include <dataclasses/physics/I3Particle.h>
#include <phys-services/I3RandomService.h>

/// @brief Base class for particle propagators
class I3PropagatorService {
  public:
    I3_FORWARD_DECLARATION(DiagnosticMap);
    /// @brief Propagate a particle and return secondaries
    ///
    /// @param[in,out] particle  a particle to be propagated. The propagator
    ///                          may modify its properties, e.g. set its final
    ///                          length.
    /// @param[in,out] frameinfo auxiliary information for the current event.
    ///                          The contents of this map will be added to the
    ///                          I3Frame for the current event after
    ///                          propagation.
    /// @param[in,out] frame     pointer to I3Frame.
    ///                          This is added because NeutrinoPropagator requires
    ///                          to access I3Frame in order to get weights stored
    ///                          by preceding modules. Copying const objects from I3Frame
    ///                          may use computing resources heaviry, so use it carefully.  
    /// @returns                 a vector of secondary particles. These are
    ///                          final states as far as this propagator is
    ///                          concerned; they may be handled by a different
    ///                          propagator, but will never be passed back to
    ///                          this one in a further iteration. 
    virtual std::vector<I3Particle> Propagate(I3Particle& particle, DiagnosticMapPtr frameinfo, I3FramePtr frame) = 0;
    /// @brief Set the random number generator to be used
    virtual void SetRandomNumberGenerator(I3RandomServicePtr random) = 0;

    I3PropagatorService();
    virtual ~I3PropagatorService();
    
    /// @ a read-write replacement for I3Frame
    class DiagnosticMap {
      public:
        bool Has(const std::string &key) const
        {
            return map_.count(key);
        }
        template <typename T>
        T
        Get(const std::string &key,
            typename boost::enable_if<is_shared_ptr<T> >::type * = 0)
        {
            map_type::iterator it = map_.find(key);
            if (it == map_.end())
                return boost::shared_ptr<typename T::element_type>();
            return boost::dynamic_pointer_cast<typename T::element_type>(it->second);
        }
        void
        Put(const std::string &key, I3FrameObjectPtr item)
        {
            map_type::iterator it = map_.find(key);
            if (it != map_.end()) {
                log_fatal_stream("Map already contains " << key);
            }
            
            map_.insert(std::make_pair(key, item));
        }
        
        typedef std::map<std::string, I3FrameObjectPtr> map_type;
        typedef map_type::value_type value_type;
        typedef map_type::iterator iterator;
        
        iterator begin() { return map_.begin(); }
        iterator end() { return map_.end(); }
      private:
        map_type map_;
    };

};

I3_POINTER_TYPEDEFS(I3PropagatorService);

typedef std::map<I3Particle::ParticleType, I3PropagatorServicePtr> I3ParticleTypePropagatorServiceMap;

I3_POINTER_TYPEDEFS(I3ParticleTypePropagatorServiceMap);

#endif //I3PROPAGATORSERVICE_H
