# ifndef CPPAD_UTILITY_ERROR_HANDLER_HPP
# define CPPAD_UTILITY_ERROR_HANDLER_HPP
/*
  =============================================================================
  adlaplace R-safe shadow of CppAD's error_handler.hpp
  =============================================================================

  WHY THIS FILE EXISTS
    R CMD check flags the shared library for references to std::cerr / exit that
    come from CppAD's default ErrorHandler (see Writing R Extensions,
    "Compiled code"). This header is intentionally placed on the package include
    path *before* the system CppAD headers (see src/Makevars.in PKG_CPPFLAGS:
    -I../src/include appears before @CPPAD_INCLUDE_CPPFLAGS@) so that

      #include <cppad/utility/error_handler.hpp>

    resolves here instead of the Homebrew/system copy.

  WHAT CHANGED VS UPSTREAM
    ErrorHandler::Default uses REprintf + Rf_error instead of std::cerr +
    std::exit. API and the rest of the class match CppAD upstream.

  HOW TO REVERT
    Delete this file (and the sibling README.md in this directory). The build
    will then pick up the system
    <cppad-prefix>/include/cppad/utility/error_handler.hpp again.

  Upstream: CppAD (EPL-2.0 OR GPL-2.0-or-later), Bradley M. Bell.
  =============================================================================
*/

// Do not #include <R.h> here: R's length() macro breaks libc++ / OpenMP headers
// when this file is pulled in early via CppAD. Link against libR only.
// [[noreturn]] must be on the first Rf_error declaration (matches Rinternals).
extern "C" {
   void REprintf(const char *, ...);
   [[noreturn]] void Rf_error(const char *, ...);
}

# include <cppad/configure.hpp>
# include <cppad/local/set_get_in_parallel.hpp>
# include <cassert>

namespace CppAD { // BEGIN CppAD namespace

class ErrorHandler {
   template <class Base>
   friend void parallel_ad(void);
public:
   typedef void (*Handler)
      (bool, int, const char *, const char *, const char *);

   // construct a new handler
   ErrorHandler(Handler handler) : previous( Current() )
   {  if( local::set_get_in_parallel() )
      {  bool known       = true;
         int  line        = __LINE__;
         const char* file = __FILE__;
         const char* exp  = "! local::set_get_in_parallel()";
         const char* msg  =
            "Using ErrorHandler constructor in parallel mode.";
         Call(known, line, file, exp, msg);
      }
      Current() = handler;
   }

   // destructor for an error handler
   ~ErrorHandler(void)
   {  if( local::set_get_in_parallel() )
      {  bool known       = true;
         int  line        = __LINE__;
         const char* file = __FILE__;
         const char* exp  = "! local::set_get_in_parallel()";
         const char* msg  =
            "Using ErrorHandler destructor in parallel mode.";
         Call(known, line, file, exp, msg);
      }
      Current() = previous;
   }

   // report an error
   static void Call(
      bool        known,
      int         line ,
      const char *file ,
      const char *exp  ,
      const char *msg  )
   {  Handler handler = Current();
      handler(known, line, file, exp, msg);
   }

private:
   const Handler previous;

   // R-safe default: no std::cerr / std::exit (CRAN compiled-code check).
   static void Default(
      bool        known,
      int         line ,
      const char *file ,
      const char *exp  ,
      const char *msg  )
   {
      REprintf("%s", CPPAD_PACKAGE_STRING);
      if( known )
         REprintf(" error from a known source:\n");
      else
         REprintf(" error from unknown source\n");
      if( msg != nullptr && msg[0] != '\0' )
         REprintf("%s\n", msg);
      REprintf("Error detected by false result for\n");
      REprintf("    %s\n", exp != nullptr ? exp : "");
      REprintf("at line %d in the file\n", line);
      REprintf("    %s\n", file != nullptr ? file : "");

      // Does not return; replaces assert(false) + std::exit(1).
      Rf_error("CppAD error handler");
   }

   // current error handler
   static Handler &Current(void)
   {  static bool first_call = true;
      static Handler current = Default;
      // CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL
      // code below is like macro above but works when NDEBUG defined
      if( first_call )
      {  if( local::set_get_in_parallel() )
         {  bool known       = false;
            int  line        = __LINE__;
            const char* file = __FILE__;
            const char* exp  = "";
            const char* msg  = "";
            Call(known, line, file, exp, msg);
         }
         first_call = false;
      }
      return current;
   }
};

} // END CppAD namespace

# endif
