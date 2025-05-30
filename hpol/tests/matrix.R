 
Rcpp::sourceCpp("~/Admin/Admin/projects/testcppad/src2/matrix.cpp")

library(Matrix)

# test.sparseCholesky.R <- function() {
suppressMessages(require("Matrix", character.only=TRUE))
data("KNex", package = "Matrix")

rr <- sparseCholesky(KNex)
tinytest::expect_equal(rr[[1]],
             as.vector(solve(crossprod(KNex[[1]]), crossprod(KNex[[1]], KNex[[2]])), mode="numeric"),
             info = "Cholmod solution")




#test.wrapSparse.double.R <- function(){
res <- wrapSparseDouble()
rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
colnames(rr) <- NULL

identical(rr, res)


res <- wrapSparseDoubleColumnMajor()
rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
colnames(rr) <- NULL
tinytest::expect_equal(res, rr, info = "wrap<SparseMatrix<double, Eigen::ColMajor> >")

## test.wrapSparse.int.ColMajor.R <- function(){  ## classes not yet exported from Matrix

##     fx <- cxxfunction( , '

##     Eigen::SparseMatrix<int, Eigen::ColMajor>  mm(9,3);
##     mm.reserve(9);
##     for (int j = 0; j < 3; ++j) {
##         mm.startVec(j);
##         for (int i = 3 * j; i < 3 * (j + 1); ++i)
##             mm.insertBack(i, j) = 1;
##     }
##     mm.finalize();
##     return wrap(mm);
## ' , plugin = "RcppEigen" )

##     #res <- fx()
##     #rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
##     #colnames(rr) <- NULL
##     #expect_equal( res, rr, info = "wrap<SparseMatrix<double, Eigen::ColMajor> >")
##     checkException( fx(), info = "wrap<SparseMatrix<int, Eigen::ColMajor> >" )
## }

#test.wrapSparse.double.RowMajor.R <- function(){
res <- wrapSparseDoubleRowMajor()
rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
colnames(rr) <- NULL
tinytest::expect_equal( res, rr, info = "wrap<SparseMatrix<double, Eigen::RowMajor> >")


## test.wrapSparse.int.RowMajor.R <- function(){

##     fx <- cxxfunction( , '

##     Eigen::SparseMatrix<int, Eigen::RowMajor>  mm(9,3);
##     mm.reserve(9);
##     for (int irow = 0; irow < 9; ++irow) {
##         mm.startVec(irow);
##         mm.insertBack(irow, irow / 3) = 9 - irow;
##     }
##     mm.finalize();
##     return wrap(mm);
## ' , plugin = "RcppEigen" )

##     #res <- fx()
##     #rr <- new( "igRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=9L:1L, Dim=c(9L,3L) )
##     #colnames(rr) <- NULL
##     #expect_equal( res, rr, info = "wrap<SparseMatrix<int, Eigen::RowMajor> >")
##     checkException( fx(), info = "wrap<SparseMatrix<int, Eigen::RowMajor> >" )
## }

#test.asSparse.double.ColMajor.R <- function(){
rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
colnames(rr) <- NULL
res <- asSparseDoubleColumnMajor( rr )
tinytest::expect_equal( res, rr, info = "as<SparseMatrix<double, Eigen::ColMajor> >")


#test.asMappedSparse.double.ColMajor.R <- function(){
rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
colnames(rr) <- NULL
res <- asMappedSparseDoubleColMajor( rr )
tinytest::expect_equal( res, sum(rr), info = "as<Map<SparseMatrix<double, Eigen::ColMajor> > >")


#test.asMappedSparse.deprecated.double.ColMajor.R <- function(){
fx <- asMappedSparseDeprecatedDoubleColMajor
rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
colnames(rr) <- NULL
res <- fx( rr )
tinytest::expect_equal( res, sum(rr), info = "as<MappedSparseMatrix<double, Eigen::ColMajor> >")


#test.asSparse.double.RowMajor.R <- function(){
rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
colnames(rr) <- NULL
res <- asSparseDoubleRowMajor( rr )
tinytest::expect_equal( res, sum(rr), info = "as<SparseMatrix<double, Eigen::RowMajor> >")


#test.asMappedSparse.double.RowMajor.R <- function(){
rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
colnames(rr) <- NULL
res <- asMappedSparseDoubleRowMajor( rr )
tinytest::expect_equal( res, sum(rr), info = "as<Map<SparseMatrix<double, Eigen::RowMajor> > >")


#test.asMappedSparse.deprecated.double.RowMajor.R <- function(){
rr <- new( "dgRMatrix", j=rep(0L:2L, each=3), p=0L:9L, x=as.numeric(9:1), Dim=c(9L,3L) )
colnames(rr) <- NULL
res <- asMappedSparseDeprecatedDoubleRowMajor( rr )
tinytest::expect_equal( res, sum(rr), info = "as<MappedSparseMatrix<double, Eigen::RowMajor> >")



