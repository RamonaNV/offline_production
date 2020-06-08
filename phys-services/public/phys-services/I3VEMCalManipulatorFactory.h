/**
 * Copyright (C) 2008
 * The IceCube collaboration
 * ID: $Id: I3VEMCalManipulatorFactory.h 136718 2015-08-21 20:23:49Z nega $
 *
 * @file I3VEMCalManipulatorFactory.h
 * @version $Revision: 136718 $
 * @date $Date: 2015-08-21 14:23:49 -0600 (Fri, 21 Aug 2015) $
 * @author $Author: nega $
 */


#ifndef PHYS_SERVICES_VEMCALMANIPULATORFACTORY_H
#define PHYS_SERVICES_VEMCALMANIPULATORFACTORY_H


#include <icetray/I3ServiceFactory.h>
#include <phys-services/I3VEMCalManipulator.h>

class I3VEMCalManipulatorFactory: public I3ServiceFactory
{
public:
    
    I3VEMCalManipulatorFactory(const I3Context& context);
    ~I3VEMCalManipulatorFactory();
    
    bool InstallService(I3Context& services);
    
    void Configure();
    
private:
    
    boost::shared_ptr<I3VEMCalManipulator> vemCalService_;
    
    std::string inCalServiceName_;
    std::string outCalServiceName_;
    
    std::string path_;
    std::string file_;
    bool forceFile_;
    bool useDefaults_;
    
    SET_LOGGER("I3VEMCalManipulatorFactory");
};

#endif 
