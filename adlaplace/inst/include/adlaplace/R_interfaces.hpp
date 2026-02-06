#ifdef R_INTERFACES_HPP
#define R_INTERFACES_HPP



Rcpp::List extractSparsity_h(
	SEXP adPack) {

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];
	return extractSparsity(gp);
}


inline SEXP get_adpack_handle_h(SEXP adPack_xptr, 
	const Rcpp::List map,
	const Rcpp::List& config)
{
  Config configC(config);

  Rcpp::XPtr<std::vector<GroupPack>> xp(adPack_xptr);
  std::vector<GroupPack>& adPack = *xp;

  adlaplace_adpack_handle* h = get_adpack_handle_cpp(adPack, configC);

  // Protect adPack_xptr so ctx->adPack never dangles
  SEXP handle = R_MakeExternalPtr((void*)h, R_NilValue, adPack_xptr);
  R_RegisterCFinalizerEx(handle, handle_finalizer, TRUE);
  Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("adlaplace_handle_ptr"));

  return handle;
}



#endif