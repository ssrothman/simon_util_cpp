#ifndef SIMONTOOLS_TOY_SHOWER_H
#define SIMONTOOLS_TOY_SHOWER_H

#include "jet.h"
#include <string>
#include <queue>
#include <random>
#include <Math/Vector3D.h>
#include <Math/Rotation3D.h>

#include <iostream>
#include <fstream>

namespace simon{
    enum class phiPDF {
        UNIFORM,
        COS2PHI,
    };
    std::ostream& operator<<(std::ostream& os, const phiPDF& pdf){
        switch(pdf){
            case phiPDF::UNIFORM:
                os << "UNIFORM";
                break;
            case phiPDF::COS2PHI:
                os << "COS2PHI";
                break;
        }
        return os;
    }

    enum class zPDF {
        UNIFORM,
        LNX,
        GLUON,
    };
    std::ostream& operator<<(std::ostream& os, const zPDF& pdf){
        switch(pdf){
            case zPDF::UNIFORM:
                os << "UNIFORM";
                break;
            case zPDF::LNX:
                os << "LNX";
                break;
            case zPDF::GLUON:
                os << "GLUON";
                break;
        }
        return os;
    }

    enum class thetaPDF {
        UNIFORM,
        LNX
    };
    std::ostream& operator<<(std::ostream& os, const thetaPDF& pdf){
        switch(pdf){
            case thetaPDF::UNIFORM:
                os << "UNIFORM";
                break;
            case thetaPDF::LNX:
                os << "LNX";
                break;
        }
        return os;
    }

    class LightweightParticle {
    public:
        //no default constructor
        LightweightParticle() = delete;

        //default constructor
        LightweightParticle(const double px, const double py, const double pz,
                            const double dpx, const double dpy, const double dpz) :
            px_(px),
            py_(py),
            pz_(pz),
            pt_(std::sqrt(px*px + py*py)),
            dpx_(dpx),
            dpy_(dpy),
            dpz_(dpz) {}

        LightweightParticle(const ROOT::Math::XYZVector& vec,
                            const ROOT::Math::XYZVector& dvec) :
            px_(vec.x()),
            py_(vec.y()),
            pz_(vec.z()),
            pt_(std::sqrt(px()*px() + py()*py())),
            dpx_(dvec.x()),
            dpy_(dvec.y()),
            dpz_(dvec.z()) {}

        //default destructor
        ~LightweightParticle() = default;

        double pt() const { return pt_; }
        double px() const { return px_; }
        double py() const { return py_; }
        double pz() const { return pz_; }
        double dpx() const { return dpx_; }
        double dpy() const { return dpy_; }
        double dpz() const { return dpz_; }

        bool operator<(const LightweightParticle& other) const {
            return pt() < other.pt();
        }
        bool operator>(const LightweightParticle& other) const {
            return pt() > other.pt();
        }
        bool operator==(const LightweightParticle& other) const {
            return pt() == other.pt();
        }
        bool operator!=(const LightweightParticle& other) const {
            return pt() != other.pt();
        }
        bool operator<=(const LightweightParticle& other) const {
            return pt() <= other.pt();
        }
        bool operator>=(const LightweightParticle& other) const {
            return pt() >= other.pt();
        }

    private:
        //momentum vector
        double px_;
        double py_;
        double pz_;

        //pt, precomputed for sorting
        double pt_;

        //dipole axis
        double dpx_;
        double dpy_;
        double dpz_;
    };

    class ToyShowerer {
    public:
        typedef std::priority_queue<LightweightParticle> particle_queue;

        explicit ToyShowerer(std::string phi_mode,
                             std::string z_mode,
                             std::string theta_mode,
                             double zcut,
                             double theta_min,
                             double theta_max,
                             bool angular_ordered,
                             std::string logfile = "");

        ToyShowerer() = delete;
        ~ToyShowerer() {
            if(logging_) {
                logfile_.close();
            }
        };


        void enable_logging(std::string logfile);

        void shower(const double pt, 
                    const double eta, 
                    const double phi, 
                    const unsigned Npart,
                    jet& result);

        void test_cos2phi(const size_t, const int);
        void test_lnx(const size_t, const int, const double, const double);
        void test_gluon_z(const size_t, const int);
    private:
        double cos2phi(const double phi);
        double cos2phi_cdf(const double phi);
        double sample_cos2phi();
        double sample_phi();

        double lnx_cdf(const double x, const double min, const double max);
        double sample_lnx(double min, double max);

        double gluon_z_pdf(const double z);
        double gluon_z_antiderivative(const double z);
        double gluon_z_cdf(const double z);
        double sample_gluon_z();

        double sample_z();
        double sample_theta();

        void get_rotation_about_axis(const ROOT::Math::XYZVector& axis, 
                                     const double angle, 
                                     ROOT::Math::Rotation3D& result);

        void do_one_splitting(particle_queue& particles,
                double z, double theta, double phi);

        phiPDF phi_mode_;
        zPDF z_mode_;
        thetaPDF theta_mode_;

        double zcut_;
        double gluon_cdf_norm_, gluon_cdf_low_;
        double theta_min_, theta_max_;

        std::default_random_engine rng_;
        std::uniform_real_distribution<double> uniform_dist_;

        std::ofstream logfile_;
        bool logging_;

        bool angular_ordered_;
    };
};

#endif
