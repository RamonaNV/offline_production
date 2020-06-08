/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3PhotonConverter.h 180218 2020-05-20 11:42:16Z jvansanten $
 *
 * @file I3PhotonConverter.h
 * @version $Revision: 180218 $
 * @date $Date: 2020-05-20 05:42:16 -0600 (Wed, 20 May 2020) $
 * @author Claudio Kopper
 */

#include "tableio/I3Converter.h"
#include "tableio/converter/I3MapConverter.h"

#include "simclasses/I3Photon.h"
#include "simclasses/I3CompressedPhoton.h"

class I3PhotonConverter : public I3ConverterImplementation<I3Photon>
{
public:
    typedef booked_type value_type;

    void AddFields(I3TableRowDescriptionPtr, const value_type&  = value_type());
    void FillSingleRow(const value_type&, I3TableRowPtr);

private:
    I3TableRowDescriptionPtr CreateDescription(const I3Photon &photon); 
    std::size_t FillRows(const I3Photon &photon, I3TableRowPtr rows);
};

typedef I3MapModuleKeyVectorConverter<I3PhotonConverter, I3PhotonSeriesMap> I3PhotonSeriesMapConverter;

class I3CompressedPhotonSeriesMapConverter : public I3ConverterImplementation<I3CompressedPhotonSeriesMap>
{
public:
    typedef booked_type value_type;

    void AddFields(I3TableRowDescriptionPtr, const value_type&  = value_type());
    void FillSingleRow(const value_type&, I3TableRowPtr);

private:
    I3TableRowDescriptionPtr CreateDescription(const value_type &example);
    std::size_t GetNumberOfRows(const value_type& value);
    std::size_t FillRows(const value_type &value, I3TableRowPtr rows);
};
