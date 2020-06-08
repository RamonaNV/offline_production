/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3CSVTableService.h 94948 2012-11-04 16:21:52Z nwhitehorn $
 *
 * @version $Revision: 94948 $
 * @date $LastChangedDate: 2012-11-04 09:21:52 -0700 (Sun, 04 Nov 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: nwhitehorn $
 */

#ifndef	I3CSVTABLESERVICE_H_INCLUDED
#define I3CSVTABLESERVICE_H_INCLUDED

#include "tableio/I3TableService.h"

class I3CSVTableService : public I3TableService {
    public:
        I3CSVTableService(const std::string& foldername);
        virtual ~I3CSVTableService();

    protected:
        virtual I3TablePtr CreateTable(const std::string& tableName, 
                                       I3TableRowDescriptionConstPtr description);
        virtual void CloseFile();

    private:

        std::string folderName_;

    SET_LOGGER("I3CSVTableService");
};


#endif
