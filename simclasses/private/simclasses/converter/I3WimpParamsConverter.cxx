/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3WimpParamsConverter.cxx 66215 2010-08-17 16:53:26Z mzoll $
 *
 * @version $Revision: 66215 $
 * @date $LastChangedDate: 2010-08-17 18:53:26 +0200 (Tue, 17 Aug 2010) $
 * @author Marcel Zoll <marcel.zoll@fysik.su.se> $LastChangedBy: mzoll $
 */

#include "I3WimpParamsConverter.h"

I3WimpParamsConverter::I3WimpParamsConverter() : I3ConverterImplementation<I3WimpParams>() {}

I3TableRowDescriptionPtr I3WimpParamsConverter::CreateDescription(const I3WimpParams& wimpparams)
{
  I3TableRowDescriptionPtr desc = I3TableRowDescriptionPtr(new I3TableRowDescription() );

  desc->AddField<WimpSim::SourceType>("source","SourceType","Source of the WIMP");
  desc->AddField<int>("mass","GeV","Mass of the WIMP in GeV");
  desc->AddField<WimpSim::DecayChannel>("channel","DecayChannel","DecayChannel of the WIMP");
  desc->AddField<double>("nu_weight","number","Weight of the primary neutrino");
  desc->AddField<double>("lep_weight","number","Weight of the daughter lepton");
  desc->AddField<double>("had_weight","number","Weight of the daughter hadron");
  desc->AddField<double>("vgen","m3","Generated Volume for this event");
  desc->AddField<double>("time_mjd_double","mjd","Time for this event");
  return desc;
}

size_t I3WimpParamsConverter::FillRows(const I3WimpParams& wimpparams, I3TableRowPtr rows)
{
  rows->Set<WimpSim::SourceType>("source", wimpparams.GetSource()); //DANGER this is enum type py:tableio.I3Datatype(SourceType)
  rows->Set<int>("mass", wimpparams.GetMass());
  rows->Set<WimpSim::DecayChannel>("channel", wimpparams.GetChannel()); //DANGER this is enum type py:tableio.I3Datatype(DecayChannel)
  rows->Set<double>("nu_weight", wimpparams.GetNuWeight());
  rows->Set<double>("lep_weight", wimpparams.GetLepWeight());
  rows->Set<double>("had_weight", wimpparams.GetHadWeight());
  rows->Set<double>("vgen", wimpparams.GetVgen());
  rows->Set<double>("time_mjd_double", wimpparams.GetTime().GetModJulianDayDouble());

  return 1;
}
