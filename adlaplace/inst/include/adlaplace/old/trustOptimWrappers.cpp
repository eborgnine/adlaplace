#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Minimal API: query sizes
using plugin_nvars_fn_t = int (*)(SEXP handle);
using plugin_hes_nnz_fn_t = int (*)(SEXP handle);

// CSC Hessian sparsity pattern
// Here: fill p[0..n] and i[0..nnz-1]. Return 0 on success.
using plugin_hes_pattern_csc_fn_t = int (*)(
  SEXP handle,
  int* p, int n_plus1,
  int* i, int nnz
);

// Compute f and g (no allocations)
using plugin_fg_fn_t = int (*)(
  SEXP handle,
  const double* x, int n,
  double* f_out,
  double* g_out
);

// Sparse Hessian values in triplet form (pattern fixed by the plugin)
using plugin_hess_triplet_fn_t = int (*)(
  SEXP handle,
  const double* x, int n,
  const double* w, int m,
  int* ii, int* jj, double* vv, int nnz
);

struct AD_Func_Opt {
  // existing
  std::vector<GroupPack> &tape;
  Eigen::SparseMatrix<int> Htemplate;

  // new: plugin mode
  bool use_plugin = false;
  SEXP plugin_handle = R_NilValue;

  plugin_nvars_fn_t            plugin_nvars = nullptr;
  plugin_hes_nnz_fn_t          plugin_hes_nnz = nullptr;
  plugin_hes_pattern_csc_fn_t  plugin_hes_pattern_csc = nullptr;
  plugin_fg_fn_t               plugin_fg = nullptr;
  plugin_hess_triplet_fn_t     plugin_hess_triplet = nullptr;

  static std::vector<GroupPack>& empty_tape_ref() {
    static std::vector<GroupPack> empty;
    return empty;
  }

  // -------------------------
  // AD pack constructors
  // -------------------------

  // ctor with hessian pattern
  AD_Func_Opt(std::vector<GroupPack> &tape_,
              const std::vector<std::vector<int>> &hessianIJ)
    : tape(tape_), Htemplate() {

    if (hessianIJ.size() < 3) {
      throw std::runtime_error("hessianIJ must have at least 3 components: i, p, x");
    }

    const int n   = hessianIJ[1].empty() ? 0 : (int)hessianIJ[1].size() - 1;
    const int nnz = (int)hessianIJ[0].size();

    Eigen::Map<const Eigen::SparseMatrix<int, Eigen::ColMajor, int>> Hmap(
      n, n, nnz,
      hessianIJ[1].data(), // p
      hessianIJ[0].data(), // i
      hessianIJ[2].data()  // x
    );
    Htemplate = Hmap;
  }

  // -------------------------
  // handle constructors
  // -------------------------

  // Plugin ctor: ask plugin for Hessian CSC pattern and build Htemplate
  AD_Func_Opt(const std::string& pkg, SEXP handle, bool load_hess_pattern)
    : tape(empty_tape_ref()), Htemplate() {
    init_plugin(pkg, handle, load_hess_pattern);
  }

private:
  void init_plugin(const std::string& pkg, SEXP handle, bool want_pattern) {
    use_plugin = true;
    plugin_handle = handle;

    // Look up callables (names must match what newpackage registers!)
    plugin_nvars = (plugin_nvars_fn_t) R_GetCCallable(pkg.c_str(), "mypkg_nvars");
    plugin_hes_nnz = (plugin_hes_nnz_fn_t) R_GetCCallable(pkg.c_str(), "mypkg_hes_nnz");
    plugin_fg = (plugin_fg_fn_t) R_GetCCallable(pkg.c_str(), "mypkg_fg");
    plugin_hess_triplet = (plugin_hess_triplet_fn_t) R_GetCCallable(pkg.c_str(), "mypkg_hess_triplet");

    if (want_pattern) {
      plugin_hes_pattern_csc =
        (plugin_hes_pattern_csc_fn_t) R_GetCCallable(pkg.c_str(), "mypkg_hes_pattern_csc");

      const int n = plugin_nvars(plugin_handle);
      const int nnz = plugin_hes_nnz(plugin_handle);
      if (n < 0 || nnz < 0) Rf_error("plugin returned invalid n/nnz");

      std::vector<int> p(n + 1), i(nnz);

      int rc = plugin_hes_pattern_csc(plugin_handle, p.data(), n + 1, i.data(), nnz);
      if (rc != 0) Rf_error("plugin_hes_pattern_csc failed (rc=%d)", rc);

      // Build an integer sparse matrix template with 1s as values
      std::vector<int> x(nnz, 1);
      Eigen::Map<const Eigen::SparseMatrix<int, Eigen::ColMajor, int>> Hmap(
        n, n, nnz,
        p.data(), i.data(), x.data()
      );
      Htemplate = Hmap;
    }
  }

public:
  // -------------------------
  // trustOptim interface bits
  // -------------------------

  int get_nvars() const {
    if (use_plugin) return plugin_nvars ? plugin_nvars(plugin_handle) : 0;
    if (tape.empty()) return 0;
    return (int)tape[0].fun.Domain();
  }

  int get_nnz() const {
    if (use_plugin) return plugin_hes_nnz ? plugin_hes_nnz(plugin_handle) : 0;
    return (int)Htemplate.nonZeros();
  }

  // f only
  template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX> &x, double &f) {
    if (use_plugin) {
      const int n = (int)x.size();
      std::vector<double> xx(n);
      for (int k = 0; k < n; ++k) xx[k] = x[k];
      f = 0.0;
      // Call plugin fg with g_out = nullptr? (depends on your plugin API)
      // Better: plugin provides a dedicated f-only function. If not, do fg and ignore g.
      std::vector<double> gtmp(n, 0.0);
      int rc = plugin_fg(plugin_handle, xx.data(), n, &f, gtmp.data());
      if (rc != 0) Rf_error("plugin_fg failed (rc=%d)", rc);
      return;
    }

    // existing tape path (your current implementation)
    const std::size_t n = (std::size_t)x.size();
    CppAD::vector<double> xp(n);
    for (std::size_t i = 0; i < n; ++i) xp[i] = x[i];

    f = 0.0;
    #pragma omp parallel
    {
      double fOutHere=0.0;
      #pragma omp for schedule(dynamic,1) nowait
      for(int D=0; D<(int)tape.size(); D++) {
        CppAD::vector<double> y = tape[D].fun.Forward(0, xp);
        fOutHere += y[0];
      }
      #pragma omp critical
      f += fOutHere;
    }
  }

  // f and g
  template <class DerivedX, class DerivedG>
  void get_fdf(const Eigen::MatrixBase<DerivedX> &x, double &f, Eigen::MatrixBase<DerivedG> &g) {
    const int n = (int)x.size();

    if (use_plugin) {
      std::vector<double> xx(n), gg(n);
      for (int k = 0; k < n; ++k) xx[k] = x[k];
      int rc = plugin_fg(plugin_handle, xx.data(), n, &f, gg.data());
      if (rc != 0) Rf_error("plugin_fg failed (rc=%d)", rc);
      for (int k = 0; k < n; ++k) g[k] = gg[k];
      return;
    }

    // existing tape path
    const std::size_t ns = (std::size_t)n;
    std::vector<double> gOut(ns, 0.0);
    CppAD::vector<double> xp(ns);
    for (std::size_t i = 0; i < ns; ++i) xp[i] = x[i];

    grad(xp, f, gOut, tape);
    for (std::size_t i = 0; i < ns; ++i) g[i] = gOut[i];
  }

  // f, g, and H
  template <class DerivedX, class DerivedG>
  void get_fdfh(const Eigen::MatrixBase<DerivedX> &x,
                double &f,
                Eigen::MatrixBase<DerivedG> &g,
                Eigen::SparseMatrix<double> &H) {
    const int n = (int)x.size();

    if (use_plugin) {
      std::vector<double> xx(n), gg(n);
      for (int k = 0; k < n; ++k) xx[k] = x[k];

      int rc = plugin_fg(plugin_handle, xx.data(), n, &f, gg.data());
      if (rc != 0) Rf_error("plugin_fg failed (rc=%d)", rc);
      for (int k = 0; k < n; ++k) g[k] = gg[k];

      // Hessian via triplets
      const int nnz = plugin_hes_nnz(plugin_handle);
      std::vector<int> ii(nnz), jj(nnz);
      std::vector<double> vv(nnz);

      // weights for scalar-output objective
      const double w1 = 1.0;
      rc = plugin_hess_triplet(plugin_handle, xx.data(), n, &w1, 1,
                              ii.data(), jj.data(), vv.data(), nnz);
      if (rc != 0) Rf_error("plugin_hess_triplet failed (rc=%d)", rc);

      // Build Eigen sparse from triplets
      std::vector<Eigen::Triplet<double,int>> trips;
      trips.reserve((size_t)nnz);
      for (int k = 0; k < nnz; ++k)
        trips.emplace_back(ii[k], jj[k], vv[k]);

      H.resize(n, n);
      H.setFromTriplets(trips.begin(), trips.end());
      return;
    }

    // existing tape path
    const std::size_t ns = (std::size_t)n;
    CppAD::vector<double> xp(ns);
    for (std::size_t i = 0; i < ns; ++i) xp[i] = x[i];
    std::vector<double> gOut(ns, 0.0);

    grad(xp, f, gOut, tape);
    get_hess_function(xp, H, tape, Htemplate);

    for (std::size_t i = 0; i < ns; ++i) g[i] = gOut[i];
  }
};
