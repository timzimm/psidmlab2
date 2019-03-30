// Unitary split operator. Propagates psi according to
// psi_(n+1) = exp(-i/2 * K) exp(-i * V) exp(-i/2 * K) psi_(n)
// for separable Hamiltionian H = K + V = -1/2*del^2_x + a(t)*V(psi)

class USO_KDK : public Algorithm {
    // Abbreviations
    using cmplx = std::complex<double>;

    // Complex row vector
    using CRV = blaze::DynamicVector<cmplx, blaze::rowVector>;
    // Complex column vector
    using CCV = blaze::DynamicVector<cmplx>;
    // Real column vector
    using RCV = blaze::DynamicVector<double>;
    // Complex dense matrix in row-major order
    using CRM = blaze::DynamicMatrix<cmplx>;
    // Complex sparse diagonal matrix in row-major order
    using CDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<cmplx>>;

    // Data members
    Cosmology cosmo;
    std::unique_ptr<Potential::Algorithm> pot;
    size_t N;
    double L;
    bool firstStep;
    CDM K;
    CDM D;
    RCV wavenumbers;

    fftw_plan forwards;
    fftw_plan backwards;
    fftw_plan forwards_op;
    fftw_plan backwards_op;

    // Buffer to store k representation of psi for KDK
    CRM psis_cached;

    // Transforms each row of matrix_on according to the passed plan and stores
    // the result in matrix_out
    void transform_matrix(const fftw_plan& plan, CRM& matrix_in,
                          CRM& matrix_out);

    // Kick Operator - returns a blaze expression
    auto kick(const CRM& psis_in_k, const double dt, const double weight);

    // Drift Operator - returns a blaze expression
    auto drift(const CRM& psis_in_x, const RCV& V, const double dt,
               const double t, const double weight);

   public:
    USO_KDK(const Parameters& p);
    ~USO_KDK();
    void operator()(SimState& state) override;
};
