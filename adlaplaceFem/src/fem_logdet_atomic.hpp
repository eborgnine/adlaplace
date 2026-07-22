#ifndef ADLAPLACEFEM_FEM_LOGDET_ATOMIC_HPP
#define ADLAPLACEFEM_FEM_LOGDET_ATOMIC_HPP

// CppAD atomic op for the FEM precision log-determinant.
//
// The precision is a fixed linear combination of constant Gram matrices,
//   Q(x) = sum_j x_j * M_j    (M_j in {C, G, G2, G3}),
// so log det Q is a function of the m = 3 or 4 scalar coefficients only.
// Derivatives come from the trace identity d log det Q = Tr(Q^{-1} dQ):
//   d/dx_j log det Q  =  Tr(Q^{-1} M_j)
// which needs Q^{-1} only on the sparsity pattern of Q (Takahashi selected
// inverse of the sparse LDL factor). Second order runs the same LDL +
// Takahashi code with a dual-number scalar:
//   d2/dx_i dx_j log det Q = -Tr(Q^{-1} M_i Q^{-1} M_j)
// (the Tr(Q^{-1} d2Q) term vanishes because Q is linear in x).
// Nothing here is ever recorded on the AD tape.

#include <cppad/cppad.hpp>

#include "adlaplace/chol_update_impl.hpp"
#include "adlaplace/takahashi_impl.hpp"

#include <cmath>
#include <cstddef>
#include <deque>
#include <mutex>
#include <utility>
#include <vector>

namespace femlogdet {

// Immutable per-model data, registered at tape-build time and looked up by
// the atomic's call_id during forward/reverse sweeps.
struct Payload {
	std::vector<int> Q_p;
	std::vector<int> Q_i;
	// Symmetry weights per stored (upper-triangle) nonzero: 1 diagonal,
	// 2 off-diagonal.
	std::vector<double> w;
	// Gram value vectors aligned with Q_p/Q_i; one per coefficient.
	std::vector<std::vector<double>> M;
	std::vector<int> perm;
	std::vector<int> perm_inv;
	std::vector<int> L1_p;
	std::vector<int> L1_i;
	std::size_t n = 0;

	std::size_t m() const { return M.size(); }
};

inline bool scalar_invalid(const double& v) { return !std::isfinite(v); }
inline bool scalar_invalid(const adlaplace::chol::Dual& v) {
	return !std::isfinite(v.val);
}

// Assemble Q, factor, and (optionally) compute the gradient of log det Q
// with respect to the m coefficients. With Scalar = Dual the gradient's dot
// components carry the Hessian-direction product.
template <typename Scalar>
Scalar fem_logdet_eval(const Payload& pay, const std::vector<Scalar>& x,
                       std::vector<Scalar>* grad_out) {
	const std::size_t nnz = pay.Q_i.size();
	const std::size_t m = pay.m();

	std::vector<Scalar> Q_x(nnz, Scalar(0));
	for (std::size_t j = 0; j < m; ++j) {
		const std::vector<double>& Mj = pay.M[j];
		for (std::size_t k = 0; k < nnz; ++k) {
			Q_x[k] += x[j] * Scalar(Mj[k]);
		}
	}

	std::vector<Scalar> L_x(pay.L1_i.size(), Scalar(0));
	std::vector<Scalar> D(pay.n, Scalar(0));
	const Scalar log_det = adlaplace::chol::chol_update_csc(
		pay.Q_p, pay.Q_i, Q_x, pay.perm, pay.L1_p, pay.L1_i, L_x, D);

	if (grad_out != nullptr) {
		grad_out->assign(m, Scalar(0));
		if (scalar_invalid(log_det)) {
			const Scalar bad =
				adlaplace::chol::CholLogDetTraits<Scalar>::invalid_result();
			for (std::size_t j = 0; j < m; ++j) {
				(*grad_out)[j] = bad;
			}
			return log_det;
		}
		std::vector<Scalar> sigma;
		adlaplace::chol::takahashi_selected_inv(pay.L1_p, pay.L1_i, L_x, D,
		                                        sigma);
		std::vector<Scalar> S_x;
		adlaplace::chol::selected_inv_scatter(pay.Q_p, pay.Q_i, pay.perm_inv,
		                                      pay.L1_p, pay.L1_i, sigma, S_x);
		for (std::size_t j = 0; j < m; ++j) {
			const std::vector<double>& Mj = pay.M[j];
			Scalar g(0);
			for (std::size_t k = 0; k < nnz; ++k) {
				g += Scalar(pay.w[k] * Mj[k]) * S_x[k];
			}
			(*grad_out)[j] = g;
		}
	}
	return log_det;
}

class atomic_fem_logdet : public CppAD::atomic_four<double> {
public:
	atomic_fem_logdet() : CppAD::atomic_four<double>("fem_logdet") {}

	std::size_t register_payload(Payload&& pay) {
		std::lock_guard<std::mutex> lock(mutex_);
		payloads_.push_back(std::move(pay));
		return payloads_.size() - 1;
	}

private:
	// std::deque keeps element references stable across push_back, so
	// evaluation threads can hold a payload reference while a new model is
	// being taped.
	std::deque<Payload> payloads_;
	mutable std::mutex mutex_;

	const Payload& payload(std::size_t call_id) const {
		std::lock_guard<std::mutex> lock(mutex_);
		return payloads_[call_id];
	}

	bool for_type(size_t call_id,
	              const CppAD::vector<CppAD::ad_type_enum>& type_x,
	              CppAD::vector<CppAD::ad_type_enum>& type_y) override {
		(void)call_id;
		type_y.resize(1);
		CppAD::ad_type_enum t = CppAD::constant_enum;
		for (size_t j = 0; j < type_x.size(); ++j) {
			t = std::max(t, type_x[j]);
		}
		type_y[0] = t;
		return true;
	}

	bool rev_depend(size_t call_id, const CppAD::vector<bool>& ident_zero_x,
	                CppAD::vector<bool>& depend_x,
	                const CppAD::vector<bool>& depend_y) override {
		(void)call_id;
		depend_x.resize(ident_zero_x.size());
		for (size_t j = 0; j < depend_x.size(); ++j) {
			depend_x[j] = depend_y[0];
		}
		return true;
	}

	bool jac_sparsity(
		size_t call_id, bool dependency,
		const CppAD::vector<bool>& ident_zero_x,
		const CppAD::vector<bool>& select_x,
		const CppAD::vector<bool>& select_y,
		CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {
		(void)call_id;
		(void)dependency;
		const size_t n = select_x.size();
		size_t nnz = 0;
		if (select_y[0]) {
			for (size_t j = 0; j < n; ++j) {
				if (select_x[j] && !ident_zero_x[j]) {
					++nnz;
				}
			}
		}
		pattern_out.resize(1, n, nnz);
		size_t k = 0;
		if (select_y[0]) {
			for (size_t j = 0; j < n; ++j) {
				if (select_x[j] && !ident_zero_x[j]) {
					pattern_out.set(k++, 0, j);
				}
			}
		}
		return true;
	}

	bool hes_sparsity(
		size_t call_id, const CppAD::vector<bool>& select_x,
		const CppAD::vector<bool>& select_y,
		CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {
		(void)call_id;
		const size_t n = select_x.size();
		size_t nnz = 0;
		if (select_y[0]) {
			for (size_t i = 0; i < n; ++i) {
				for (size_t j = 0; j < n; ++j) {
					if (select_x[i] && select_x[j]) {
						++nnz;
					}
				}
			}
		}
		pattern_out.resize(n, n, nnz);
		size_t k = 0;
		if (select_y[0]) {
			for (size_t i = 0; i < n; ++i) {
				for (size_t j = 0; j < n; ++j) {
					if (select_x[i] && select_x[j]) {
						pattern_out.set(k++, i, j);
					}
				}
			}
		}
		return true;
	}

	bool forward(size_t call_id, const CppAD::vector<bool>& select_y,
	             size_t order_low, size_t order_up,
	             const CppAD::vector<double>& tx,
	             CppAD::vector<double>& ty) override {
		(void)select_y;
		if (order_up > 1) {
			return false;
		}
		const Payload& pay = payload(call_id);
		const size_t m = pay.m();
		const size_t q = order_up + 1;
		ty.resize(q);

		std::vector<double> x0(m);
		for (size_t j = 0; j < m; ++j) {
			x0[j] = tx[j * q + 0];
		}

		if (order_up == 0) {
			ty[0] = fem_logdet_eval<double>(pay, x0, nullptr);
			return true;
		}

		std::vector<double> grad;
		const double f = fem_logdet_eval<double>(pay, x0, &grad);
		if (order_low <= 0) {
			ty[0] = f;
		}
		double dy = 0.0;
		for (size_t j = 0; j < m; ++j) {
			dy += grad[j] * tx[j * q + 1];
		}
		ty[1] = dy;
		return true;
	}

	bool reverse(size_t call_id, const CppAD::vector<bool>& select_x,
	             size_t order_up, const CppAD::vector<double>& tx,
	             const CppAD::vector<double>& ty,
	             CppAD::vector<double>& px,
	             const CppAD::vector<double>& py) override {
		(void)select_x;
		(void)ty;
		if (order_up > 1) {
			return false;
		}
		const Payload& pay = payload(call_id);
		const size_t m = pay.m();
		const size_t q = order_up + 1;
		px.resize(m * q);
		for (size_t k = 0; k < px.size(); ++k) {
			px[k] = 0.0;
		}

		if (order_up == 0) {
			std::vector<double> x0(m);
			for (size_t j = 0; j < m; ++j) {
				x0[j] = tx[j];
			}
			std::vector<double> grad;
			fem_logdet_eval<double>(pay, x0, &grad);
			for (size_t j = 0; j < m; ++j) {
				px[j] = py[0] * grad[j];
			}
			return true;
		}

		// order_up == 1 (used by sparse_hes: Forward(1) then Reverse(2)).
		// Dual pass seeded with the forward direction x1 gives grad (val)
		// and the Hessian-direction product H * x1 (dot).
		std::vector<adlaplace::chol::Dual> x(m);
		for (size_t j = 0; j < m; ++j) {
			x[j] = adlaplace::chol::Dual(tx[j * q + 0], tx[j * q + 1]);
		}
		std::vector<adlaplace::chol::Dual> grad;
		fem_logdet_eval<adlaplace::chol::Dual>(pay, x, &grad);
		for (size_t j = 0; j < m; ++j) {
			px[j * q + 0] = py[0] * grad[j].val + py[1] * grad[j].dot;
			px[j * q + 1] = py[1] * grad[j].val;
		}
		return true;
	}
};

inline atomic_fem_logdet& fem_logdet_atomic_instance() {
	static atomic_fem_logdet op;
	return op;
}

} // namespace femlogdet

#endif
