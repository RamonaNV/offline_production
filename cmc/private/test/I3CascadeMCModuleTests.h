#ifndef I3CASCADEMCMODULETESTS__H_INCLUDED
#define I3CASCADEMCMODULETESTS__H_INCLUDED

#include<icetray/I3Module.h>

#include<string>

/**
 * Simple class that checks the contents of
 * the MCTrees before and after application of
 * the I3CascadeMCModule
 */
class I3CascadeMCModuleTests : public I3Module
{

  std::string originalTree_;
  std::string cmcTree_;

 public:

  I3CascadeMCModuleTests(const I3Context& context);
  void Configure();
  void DAQ(I3FramePtr frame);

};

#endif
