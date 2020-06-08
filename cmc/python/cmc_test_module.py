from icecube.icetray.I3Test import *
from icecube import icetray

class I3CascadeMCModuleTests(icetray.I3Module):
  def __init__(self, context):
    icetray.I3Module.__init__(self, context)
    self.AddParameter("originaltree","Name of origignal MCTree", "I3MCTree")
    self.AddParameter("cmctree","Name of cmc MCTree", "CMCTree")
    self.AddOutBox("OutBox")

  def Configure(self):
    self.originalTree = self.GetParameter("originaltree")
    self.cmcTree = self.GetParameter("cmctree")

  # checks whether the original MCTree and the one put into the frame
  # by the I3CascadeMCModule are there and the that there are particles
  # from the splitting routine in the the new tree
  # whether the sub-cascade parameters are correct is tested in I3CascadeSplitTests 
  # Test ensures:
  # sub-cascades added at original place
  # sub-cascade enrgies sum up to total energy
  # subcascades have the same direction as the original cascade
  def DAQ(self, frame):
    print(frame)
    original = frame.Get(self.originalTree)
    cmc = frame.Get(self.cmcTree)
  
    # make sure both trees the orifinal and the one after splitting are there
    assert(original)
    assert(cmc)

    # find the cascade in the original tree
    # there's only one cascade, get the depth of it
    # and leave loop
    originalDepth = None
    for p in original:
      if p.is_cascade:
        originalDepth = original.depth(p)
        break

      # check the new tree, find first cascade
      # check whether this is at the same depth (first sub-cascades)
      # check whether it has childs, those are the succeeding sub-cascades
      for p in cmc :
        # check if particle is a cascade
        if p.IsCascade():
          # check that the depth didn't change
          ENSURE_EQUAL(originalDepth, cmc.depth(p), "Sub-cascade list is not at same depth as original cascade")
          # chekc that there are sub-cascades
          ENSURE(cmc.number_of_children(p) > 0, "Sub-cascade list is empty")
          break

    self.PushFrame(frame)
