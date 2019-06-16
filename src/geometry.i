/* SWIG interface file */

%module geometry
%include "std_vector.i"
%{
    /* Geometry functions */
    extern double growth_rate(const std::vector<double> &, const std::vector<double> &);
    extern double le_angle(const std::vector<double> &, const std::vector<double> &,
                           const std::vector<double> &, const std::vector<double> &);
%}

namespace std {
    %template(vectord) vector<double>;
};

/* Geometry functions */
extern double growth_rate(const std::vector<double> &, const std::vector<double> &);
extern double le_angle(const std::vector<double> &, const std::vector<double> &,
                       const std::vector<double> &, const std::vector<double> &);
