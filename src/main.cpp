#include <benchmark/benchmark.h>

#include <blaze/Math.h>
#include <blaze/math/Subvector.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/Column.h>
#include <blaze/math/Row.h>
#include <blaze/util/AlignedAllocator.h>

void init( blaze::DynamicVector<double>& v ) {
    for ( size_t m =0; m < v.size(); ++m ) {
        v[m] = blaze::rand<double>(0,10);
    }
}

template <bool SO>
void init( blaze::DynamicMatrix<double, SO>& mat ) {
    for ( size_t m = 0; m < blaze::rows(mat); ++m ) {
        for ( size_t n =0; n < blaze::columns(mat); ++n ) {
            mat(m,n) = blaze::rand<double>(0,10);
        }
    }
}

template <bool SO1, bool SO2, bool SO3>
static void BM_ComputeUgp(benchmark::State& state) {
    blaze::DynamicMatrix<double, SO1> phi_gp(3,3);
    blaze::DynamicMatrix<double, SO2> dofs(state.range(0),3);
    blaze::DynamicMatrix<double, SO3> ugp(state.range(0),3);

    init(phi_gp);
    init(dofs);

    for ( auto _ : state ) {
        ugp = dofs * phi_gp;
    }

}

BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::columnMajor, blaze::columnMajor, blaze::columnMajor) -> Range(1<<9, 1<<13);
BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::columnMajor, blaze::columnMajor, blaze::rowMajor) -> Range(1<<9, 1<<13);
BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::columnMajor, blaze::rowMajor, blaze::columnMajor) -> Range(1<<9, 1<<13);
BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::columnMajor, blaze::rowMajor, blaze::rowMajor) -> Range(1<<9, 1<<13);
BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::rowMajor, blaze::columnMajor, blaze::columnMajor) -> Range(1<<9, 1<<13);
BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::rowMajor, blaze::columnMajor, blaze::rowMajor) -> Range(1<<9, 1<<13);
BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::rowMajor, blaze::rowMajor, blaze::columnMajor) -> Range(1<<9, 1<<13);
BENCHMARK_TEMPLATE(BM_ComputeUgp, blaze::rowMajor, blaze::rowMajor, blaze::rowMajor) -> Range(1<<9, 1<<13);

BENCHMARK_MAIN();