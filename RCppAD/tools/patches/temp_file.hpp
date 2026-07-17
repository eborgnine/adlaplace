# ifndef CPPAD_LOCAL_TEMP_FILE_HPP
# define CPPAD_LOCAL_TEMP_FILE_HPP
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later
// SPDX-FileCopyrightText: Bradley M. Bell <bradbell@seanet.com>
// SPDX-FileContributor: 2003-22 Bradley M. Bell
// ----------------------------------------------------------------------------
// RCppAD: provide an inline definition so LinkingTo packages do not need
// libcppad_lib (headers-only packaging). Upstream only declares this symbol.
# include <string>
# include <vector>
# include <cstdio>
# include <cstdlib>
# include <unistd.h>

# include <cppad/configure.hpp>

namespace CppAD { namespace local {
   // Create a unique temporary file path and return its name.
   // Returns "" on failure.
   inline std::string temp_file(void)
   {
      const char* tmp = std::getenv("TMPDIR");
      if( tmp == nullptr )
         tmp = std::getenv("TMP");
      if( tmp == nullptr )
         tmp = "/tmp";
      std::string pattern = std::string(tmp);
      if( !pattern.empty() && pattern.back() != '/' )
         pattern += '/';
      pattern += "cppadXXXXXX";
      std::vector<char> buf(pattern.begin(), pattern.end());
      buf.push_back('\0');
      int fd = mkstemp(buf.data());
      if( fd < 0 )
         return "";
      close(fd);
      return std::string(buf.data());
   }
} }

# endif
