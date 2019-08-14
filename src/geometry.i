/* SWIG interface file */

%module geometry
%include "std_vector.i"
%{
    /* Geometry functions */
    extern double growth_rate(const std::vector<double> &, const std::vector<double> &);
    extern double le_angle(const std::vector<double> &, const std::vector<double> &,
                           const std::vector<double> &, const std::vector<double> &);
    extern int check_te_angle(const std::vector<double> &, const std::vector<double> &,
                              const double &, const double &);
    extern std::vector<double> curvature(const std::vector<double> &,
                                         const std::vector<double> &);
    extern std::vector<double> curvature_reversals(const std::vector<double> &,
                                                   const std::vector<double> &,
                                                   const double &);
%}

namespace std {
    %template(vectord) vector<double>;
};

/* Geometry functions */
extern double growth_rate(const std::vector<double> &, const std::vector<double> &);
extern double le_angle(const std::vector<double> &, const std::vector<double> &,
                       const std::vector<double> &, const std::vector<double> &);
extern int check_te_angle(const std::vector<double> &, const std::vector<double> &,
                          const double &, const double &);
extern std::vector<double> curvature(const std::vector<double> &,
                                     const std::vector<double> &);
extern std::vector<double> curvature_reversals(const std::vector<double> &,
                                               const std::vector<double> &,
                                               const double &);
