// Copyright (c) 2014-2018 Dr. Colin Hirsch and Daniel Frey
// Please see LICENSE for license or visit https://github.com/taocpp/PEGTL/

#ifndef TAO_PEGTL_INPUT_ERROR_HPP
#define TAO_PEGTL_INPUT_ERROR_HPP

#include <cerrno>
#include <sstream>
#include <stdexcept>
#include <system_error>

#include "config.hpp"

// In PEGTL 3 input_error was changed to std::system_error and
// std::filesystem::filesystem_error (the latter is derived from the former).
// Here we half-backported it - input_error is replaced with system_error.
#if 0
namespace tao
{
   namespace TAO_PEGTL_NAMESPACE
   {
      struct input_error
         : std::runtime_error
      {
         input_error( const std::string& message, const int in_errorno )
            : std::runtime_error( message ),
              errorno( in_errorno )
         {
         }

         int errorno;
      };

   }  // namespace TAO_PEGTL_NAMESPACE

}  // namespace tao
#endif

#define TAO_PEGTL_INTERNAL_UNWRAP( ... ) __VA_ARGS__

#define TAO_PEGTL_THROW_INPUT_ERROR( MESSAGE )                                          \
   do {                                                                                 \
      const int errorno = errno;                                                        \
      std::ostringstream oss;                                                           \
      oss << TAO_PEGTL_INTERNAL_UNWRAP( MESSAGE ); \
      throw std::system_error( errorno, std::system_category(), oss.str() );            \
   } while( false )

#endif
