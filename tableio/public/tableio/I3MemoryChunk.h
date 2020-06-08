/*
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3MemoryChunk.h 67230 2010-10-05 08:28:09Z kislat $
 *
 * @version $Revision: 67230 $
 * @date $LastChangedDate: 2010-10-05 02:28:09 -0600 (Tue, 05 Oct 2010) $
 * @author Eike Middell <eike.middell@desy.de> Last changed by: $LastChangedBy: kislat $
 */

#ifndef	DATAENTRY_H_INCLUDED
#define DATAENTRY_H_INCLUDED

/**
 * \brief A single field in the in-memory representation of the table.
 *
 * The converters will fill a contigous region of space of which different
 * regions may contain different datatypes. In order to guarantee proper
 * alignment of the fields reserve for every field sizeof(I3MemoryChunk)
 * bytes.
 */
union I3MemoryChunk {
    char char_;
    short short_;
    int int_;
    long long_;
    long long longLong_;
    float float_;
    double double_;
    long double longDouble_;
    bool bool_;
};

/// The amount of memory to be reserved (in the internal storage) for each table field
const size_t I3MEMORYCHUNK_SIZE=sizeof(I3MemoryChunk);

#endif
