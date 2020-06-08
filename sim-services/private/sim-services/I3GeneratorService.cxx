#include <sim-services/I3GeneratorService.h>

I3GeneratorService::I3GeneratorService():
I3ServiceBase("I3GeneratorService")
{}

I3GeneratorService::I3GeneratorService(const I3Context& ctx):
I3ServiceBase(ctx)
{}

I3GeneratorService::~I3GeneratorService() 
{}
