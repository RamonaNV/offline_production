/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3IndexColumnsGenerator.h 136657 2015-08-20 18:50:01Z kkrings $
 *
 * @version $Revision: 136657 $
 * @date $LastChangedDate: 2015-08-20 12:50:01 -0600 (Thu, 20 Aug 2015) $
 * @author Eike Middell <eike.middell@desy.de> Last changed by: $LastChangedBy: kkrings $
 */

#ifndef TABLEIO_I3INDEXCOLUMNSGENERATOR_H_INCLUDED
#define TABLEIO_I3INDEXCOLUMNSGENERATOR_H_INCLUDED

#include "tableio/I3Converter.h"
#include "dataclasses/physics/I3EventHeader.h"
#include <I3/hash_map.h>

class I3IndexColumnsGenerator : public I3ConverterImplementation<I3EventHeader> {
    public:
	I3IndexColumnsGenerator();
	I3IndexColumnsGenerator(const std::vector<std::string> &streams);
	I3IndexColumnsGenerator(I3TableRowDescriptionConstPtr desc);
    private:
	I3TableRowDescriptionPtr CreateDescription(const I3EventHeader& object);
	size_t FillRows(const I3EventHeader& header, I3TableRowPtr rows);
	
	typedef hash_map<std::string,int> stream_map_t;
	typedef std::vector<stream_map_t::key_type> istream_t;
	stream_map_t streams_;
	istream_t istreams_;
    public:
	I3EventHeaderPtr Resurrect(I3TableRowPtr rows);
};

#endif // TABLEIO_I3INDEXCOLUMNSGENERATOR_H_INCLUDED
