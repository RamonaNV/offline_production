/**
 *    copyright  (C) 2004
 *  the icecube collaboration
 *  @version $Id: copy_if.h 94873 2012-11-01 21:21:27Z claudio.kopper $
 *  @date $Date: 2012-11-01 15:21:27 -0600 (Thu, 01 Nov 2012) $
 */

//
// implementation of copy_if, which, inexplicably, is missing
// from the STL implementation
//
// see Scott Meyers Effective STL p. 156
//
#ifndef COPY_IF_H_INCLUDED
#define COPY_IF_H_INCLUDED

template <typename InputIterator,
          typename OutputIterator,
          typename Predicate>
OutputIterator 
copy_if(InputIterator begin,
        InputIterator end,
        OutputIterator destBegin,
        Predicate p)
{
  while (begin != end)
    {
      if (p(*begin))
        *destBegin++ = *begin;
      ++begin;
    }
  return destBegin;
}

#endif //COPY_IF_H_INCLUDED
