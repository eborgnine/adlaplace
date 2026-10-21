# RCppAD

Headers-only R package shipping [CppAD](https://github.com/coin-or/CppAD)
for use via `LinkingTo: RCppAD`.

## Use from another package

In `DESCRIPTION`:

```
Imports: RCppAD
LinkingTo: RCppAD
```

In C++:

```cpp
#include <cppad/cppad.hpp>
```

R adds `inst/include` from LinkingTo packages to the include path, so
`#include <cppad/cppad.hpp>` resolves to this package.

## Refreshing CppAD

Maintainers (not install-time):

```sh
./tools/update-cppad.sh 20260000.0
# or: ./tools/update-cppad.sh /path/to/cppad-include-parent
```

This copies upstream `include/cppad/`, installs the frozen
`tools/configure.hpp.in`, and reapplies R-safe patches.

## License

This package is licensed under **GPL (>= 2)**. Upstream CppAD is dual-licensed
EPL-2.0 OR GPL-2.0-or-later; this package redistributes it under the GPL
secondary license. See also upstream SPDX headers in `inst/include/cppad/`.
For the EPL option, use CppAD from upstream.
