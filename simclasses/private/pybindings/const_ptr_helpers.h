#ifndef CONST_PTR_HELPERS_H_INCLUDED
#define CONST_PTR_HELPERS_H_INCLUDED
/**
 * const_ptr_helpers.h
 *
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: const_ptr_helpers.h 136731 2015-08-21 21:52:17Z nega $
 *
 * @version $Revision: 136731 $
 * @date $LastChangedDate: 2015-08-21 16:52:17 -0500 (Fri, 21 Aug 2015) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy$
 */


/* Twiddles to pass const pointers across the language barrier
 *
 * see: http://language-binding.net/pyplusplus/troubleshooting_guide/shared_ptr/shared_ptr.html
 */

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


#endif  // CONST_PTR_HELPERS_H_INCLUDED
