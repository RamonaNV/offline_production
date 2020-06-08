/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3TableTranscriber.h 94949 2012-11-04 16:40:30Z nwhitehorn $
 *
 * @version $Revision: 94949 $
 * @date $LastChangedDate: 2012-11-04 09:40:30 -0700 (Sun, 04 Nov 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: nwhitehorn $
 */

#include <tableio/I3TableService.h>
#include <tableio/I3Table.h>
#include <tableio/converter/I3IndexColumnsGenerator.h>

#ifndef	I3TABLETRANSCRIBER_H_INCLUDED
#define I3TABLETRANSCRIBER_H_INCLUDED

class I3TableTranscriber {
    public:
        I3TableTranscriber(I3TableServicePtr input, I3TableServicePtr output);
        void Execute();
        void Execute(size_t nframes);
        
        void Finish();
        
    protected:
         I3TablePtr ConnectTable(std::string tableName, 
                                 const I3TableRowDescription& description);
         void DisconnectTable(I3TablePtr& table);
    private:
 
        I3TableTranscriber();
        I3TableServicePtr inputService_;
        I3TableServicePtr outputService_;
        // std::map<std::string, I3TablePtr> inputTables_;
	boost::shared_ptr<I3IndexColumnsGenerator> indexer_;
        std::vector<std::pair<I3TablePtr,I3TablePtr> > transcriptions_;
        size_t nEvents_;
    SET_LOGGER("I3TableTranscriber");
};

#endif

