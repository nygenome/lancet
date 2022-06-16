
#ifndef BAMTOOLS_API_EXPORT_H
#define BAMTOOLS_API_EXPORT_H

#ifdef API_STATIC_DEFINE
#  define API_EXPORT
#  define API_NO_EXPORT
#else
#  ifndef API_EXPORT
#    ifdef BamTools_EXPORTS
        /* We are building this library */
#      define API_EXPORT 
#    else
        /* We are using this library */
#      define API_EXPORT 
#    endif
#  endif

#  ifndef API_NO_EXPORT
#    define API_NO_EXPORT 
#  endif
#endif

#ifndef API_DEPRECATED
#  define API_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef API_DEPRECATED_EXPORT
#  define API_DEPRECATED_EXPORT API_EXPORT API_DEPRECATED
#endif

#ifndef API_DEPRECATED_NO_EXPORT
#  define API_DEPRECATED_NO_EXPORT API_NO_EXPORT API_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef API_NO_DEPRECATED
#    define API_NO_DEPRECATED
#  endif
#endif

#endif /* BAMTOOLS_API_EXPORT_H */
