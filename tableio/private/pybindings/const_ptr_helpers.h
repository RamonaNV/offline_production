/**
 * const_ptr_helpers.h
 *
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: const_ptr_helpers.h 136657 2015-08-20 18:50:01Z kkrings $
 *
 * @version $Revision: 136657 $
 * @date $LastChangedDate: 2015-08-20 12:50:01 -0600 (Thu, 20 Aug 2015) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: kkrings $
 */


/* Twiddles to pass const pointers across the language barrier
 *
 * see: http://language-binding.net/pyplusplus/troubleshooting_guide/shared_ptr/shared_ptr.html
 */

#ifndef TABLEIO_CONSTPTRHELPERS_H_INCLUDED
#define TABLEIO_CONSTPTRHELPERS_H_INCLUDED

namespace boost{

    template<class T>
    inline T* get_pointer( boost::shared_ptr<const T> const& p ){
        return const_cast< T* >( p.get() );
    }

}

namespace boost{ namespace python{

    template<class T>
    struct pointee< boost::shared_ptr<T const> >{
        typedef T type;
    };

} } //boost::python

namespace utils{

    template< class T >
    void register_const_ptr(){
        namespace bpl = boost::python;
        // bpl::register_ptr_to_python< boost::shared_ptr< T > >();
        bpl::register_ptr_to_python< boost::shared_ptr< const T > >();
        bpl::implicitly_convertible< boost::shared_ptr< T >, boost::shared_ptr< const T > >();
        // bpl::implicitly_convertible< boost::shared_ptr< const T >, boost::shared_ptr< T > >();
        
    }

}

#endif // TABLEIO_CONSTPTRHELPERS_H_INCLUDED
