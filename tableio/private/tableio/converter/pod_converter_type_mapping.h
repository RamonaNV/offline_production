/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: pod_converter_type_mapping.h 125932 2014-11-19 20:53:15Z jvansanten $
 *
 * @version $Revision: 125932 $
 * @date $LastChangedDate: 2014-11-19 13:53:15 -0700 (Wed, 19 Nov 2014) $
 * @author Fabian Kislat <fabian.kislat@desy.de> $LastChangedBy: jvansanten $
 */

#ifndef TABLEIO_POD_CONVERTER_TYPE_MAPPING_H_INCLUDED
#define TABLEIO_POD_CONVERTER_TYPE_MAPPING_H_INCLUDED

/// @cond
namespace detail {

  template <typename T>
  struct pod_converter_type_mapping {
    typedef T type;
  };
    
  template <>
  struct pod_converter_type_mapping<signed short> {
    typedef int16_t type;
  };
    
  template <>
  struct pod_converter_type_mapping<signed int> {
    typedef int32_t type;
  };
    
  template <>
  struct pod_converter_type_mapping<signed long> {
    typedef int64_t type;
  };
    
  template <>
  struct pod_converter_type_mapping<signed long long> {
    typedef int64_t type;
  };
    
  template <>
  struct pod_converter_type_mapping<unsigned short> {
    typedef uint16_t type;
  };

  template <>
  struct pod_converter_type_mapping<unsigned int> {
    typedef uint32_t type;
  };
    
  template <>
  struct pod_converter_type_mapping<unsigned long> {
    typedef uint64_t type;
  };
    
  template <>
  struct pod_converter_type_mapping<unsigned long long> {
    typedef uint64_t type;
  };
    
}
/// @endcond

#endif // TABLEIO_POD_CONVERTER_TYPE_MAPPING_H_INCLUDED
