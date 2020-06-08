#ifndef I3_GENERATOR_SERVICE_H
#define I3_GENERATOR_SERVICE_H

#include <icetray/I3PointerTypedefs.h>
#include <icetray/I3Frame.h>
#include <icetray/I3ServiceBase.h>
#include <dataclasses/physics/I3MCTree.h>

/// @brief Base class for MC event generator service
class I3GeneratorService : public I3ServiceBase{
  public:
    virtual I3MCTreePtr GetNextEvent() = 0;
    virtual I3FramePtr GetNextFrame() = 0;

    virtual double GetRate() = 0;

    I3GeneratorService();
    I3GeneratorService(const I3Context& ctx);
    virtual ~I3GeneratorService();
};

I3_POINTER_TYPEDEFS(I3GeneratorService);


#endif // I3_GENERATOR_SERVICE_H
