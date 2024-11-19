#include "ToyShowerer.h"
#include "util.h"

static void print_queue_contents(const std::priority_queue<LightweightParticle>& particles){
    auto tmp = particles;
    while (!tmp.empty()){
        const LightweightParticle& particle = tmp.top();
        double p_tot = std::sqrt(particle.px()*particle.px() + particle.py()*particle.py() + particle.pz()*particle.pz());
        printf("\tp_tot = %f\n", p_tot);
        tmp.pop();
    }
}

ToyShowerer::ToyShowerer(std::string phi_mode,
                         std::string z_mode,
                         std::string theta_mode,
                         double zcut,
                         double theta_min,
                         double theta_max) :
    zcut_(zcut),
    theta_min_(theta_min),
    theta_max_(theta_max),
    rng_(std::random_device()()),
    uniform_dist_(0.0, 1.0)
{
    if (phi_mode == "UNIFORM"){
        phi_mode_ = phiPDF::UNIFORM;
    } else if (phi_mode == "COS2PHI"){
        phi_mode_ = phiPDF::COS2PHI;
    } else {
        throw std::invalid_argument("Invalid phi mode");
    }

    if (z_mode == "UNIFORM"){
        z_mode_ = zPDF::UNIFORM;
    } else if (z_mode == "LNX"){
        z_mode_ = zPDF::LNX;
    } else if (z_mode == "GLUON"){
        z_mode_ = zPDF::GLUON;
    } else {
        throw std::invalid_argument("Invalid z mode");
    }

    if (theta_mode == "UNIFORM"){
        theta_mode_ = thetaPDF::UNIFORM;
    } else if (theta_mode == "LNX"){
        theta_mode_ = thetaPDF::LNX;
    } else {
        throw std::invalid_argument("Invalid theta mode");
    }

    gluon_cdf_low_ = gluon_z_antiderivative(zcut_);
    gluon_cdf_norm_ = gluon_z_antiderivative(1.0-zcut_) - gluon_z_antiderivative(zcut_);
}

double ToyShowerer::cos2phi(const double phi){
    double tmp = std::cos(phi);
    return tmp*tmp;
}

double ToyShowerer::cos2phi_cdf(const double phi){
    double tmp = 2*phi;
    return (tmp + std::sin(tmp))/(4*M_PI);
}

double ToyShowerer::sample_cos2phi(){
    //use binary search to invert the cdf
    const static double tol = 1e-5;

    const double uniform = uniform_dist_(rng_);

    double low = 0;
    double high = 2*M_PI;
    
    double mid;

    while(high - low > tol){
        mid = (low + high) / 2;
        double cdf = cos2phi_cdf(mid);
        if (cdf < uniform){
            low = mid;
        } else {
            high = mid;
        }
    }
    return mid;
}

void ToyShowerer::test_cos2phi(const size_t N, const int bins){
    printf("Testing cos2phi sampling:\n");
    std::vector<double> samples;
    samples.reserve(N);

    for (size_t i = 0; i < N; ++i){
        samples.push_back(sample_cos2phi());
    }

    std::sort(samples.begin(), samples.end());

    //compare empirical CDF with true CDF 
    for (int i=0; i<bins; ++i){
        double phi = 2*M_PI*(i+1)/bins;
        double cdf = cos2phi_cdf(phi);
        size_t empirical_dist = std::lower_bound(samples.begin(), samples.end(), phi) - samples.begin();
        double empirical_cdf = static_cast<double>(empirical_dist)/N;
        printf("\tphi = %f, cdf = %f, empirical_cdf = %f\n", phi, cdf, empirical_cdf);
    }
    printf("\n");
}

double ToyShowerer::sample_phi_uniform(){
    return 2*M_PI*uniform_dist_(rng_);
}

double ToyShowerer::sample_lnx(double min, double max){
    double uniform = uniform_dist_(rng_);
    double logmin = std::log(min);
    double logmax = std::log(max);
    double logrange = logmax - logmin;

    return std::exp(logmin + uniform*logrange);
}

double ToyShowerer::lnx_cdf(const double x, const double min, const double max){
    double logmin = std::log(min);
    double logmax = std::log(max);
    double logrange = logmax - logmin;

    double logx = std::log(x);
    return (logx - logmin)/logrange;
}

void ToyShowerer::test_lnx(const size_t N, const int bins,
        const double min, const double max){

    printf("Testing lnx sampling:\n");
    
    std::vector<double> samples;
    samples.reserve(N);

    for (size_t i = 0; i < N; ++i){
        samples.push_back(sample_lnx(min, max));
    }
    std::sort(samples.begin(), samples.end());

    for (int i=0; i<bins; ++i){
        double logmin = std::log(min);
        double logrange = std::log(max) - logmin;
        double x = std::exp((i+1)*logrange/bins + logmin);
        double cdf_val = lnx_cdf(x, min, max);
        size_t empirical_dist = std::lower_bound(samples.begin(), samples.end(), x) - samples.begin();
        double empirical_cdf = static_cast<double>(empirical_dist)/N;
        printf("\tx = %f, cdf = %f, empirical_cdf = %f\n", x, cdf_val, empirical_cdf);
    }
    printf("\n");
}

double ToyShowerer::sample_phi(){
    switch (phi_mode_){
        case phiPDF::UNIFORM:
            return sample_phi_uniform();
        case phiPDF::COS2PHI:
            return sample_cos2phi();
        default:
            throw std::invalid_argument("Invalid phi mode");
    }
}

double ToyShowerer::gluon_z_pdf(const double z){
    if (z < zcut_ || z > 1.0-zcut_){
        return 0;
    }
    return square(1 - z + z*z)/(z * (z-1) * gluon_cdf_norm_);
}

double ToyShowerer::gluon_z_antiderivative(const double z){
    return -2*z + (1/2)*z*z - (1/3)*z*z*z - std::log(1-z) + std::log(z);
}

double ToyShowerer::gluon_z_cdf(const double z){
    return (gluon_z_antiderivative(z) - gluon_cdf_low_)/gluon_cdf_norm_;
}

double ToyShowerer::sample_gluon_z(){
    //use binary search to invert the cdf
    const static double tol = 1e-5;

    const double uniform = uniform_dist_(rng_);

    double low = zcut_;
    double high = 1.0 - zcut_;
    
    double mid;

    while(high - low > tol){
        mid = (low + high) / 2;
        double cdf = gluon_z_cdf(mid);
        if (cdf < uniform){
            low = mid;
        } else {
            high = mid;
        }
    }
    return mid;
}

void ToyShowerer::test_gluon_z(const size_t N, const int bins){
    printf("Testing gluon z sampling:\n");

    std::vector<double> samples;
    samples.reserve(N);

    for (size_t i = 0; i < N; ++i){
        samples.push_back(sample_gluon_z());
    }

    std::sort(samples.begin(), samples.end());

    //compare empirical CDF with true CDF 
    for (int i=0; i<bins; ++i){
        double z = zcut_ + (1-2*zcut_)*(i+1)/bins;
        double cdf = gluon_z_cdf(z);
        size_t empirical_dist = std::lower_bound(samples.begin(), samples.end(), z) - samples.begin();
        double empirical_cdf = static_cast<double>(empirical_dist)/N;
        printf("\tz = %f, cdf = %f, empirical_cdf = %f\n", z, cdf, empirical_cdf);
    }
    printf("\n");
}

double ToyShowerer::sample_z(){
    switch (z_mode_){
        case zPDF::UNIFORM:
            return uniform_dist_(rng_) * (1-2*zcut_) + zcut_;
        case zPDF::LNX:
            return sample_lnx(zcut_, 1.0-zcut_);
        case zPDF::GLUON:
            return sample_gluon_z();
        default:
            throw std::invalid_argument("Invalid z mode");
    }
}

double ToyShowerer::sample_theta(){
    switch(theta_mode_){
        case thetaPDF::UNIFORM:
            return uniform_dist_(rng_) * (theta_max_ - theta_min_) + theta_min_;
        case thetaPDF::LNX:
            return sample_lnx(theta_min_, theta_max_);
        default:
            throw std::invalid_argument("Invalid theta mode");
    }
    return sample_theta();
}

void ToyShowerer::get_rotation_about_axis(const ROOT::Math::XYZVector& axis, 
                                          const double angle, 
                                          ROOT::Math::Rotation3D& result){
    //get the rotation matrix for a rotation about an axis
    //by an angle
    //this is a simple application of the Rodrigues formula
    //https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    
    ROOT::Math::XYZVector unit_axis = axis.Unit();
    
    double c = std::cos(angle);
    double s = std::sin(angle);
    double t = 1 - c;

    double x = unit_axis.x();
    double y = unit_axis.y();
    double z = unit_axis.z();

    result.SetComponents(
        c + x*x*t,       x*y*t - z*s,     x*z*t + y*s,
        y*x*t + z*s,     c + y*y*t,       y*z*t - x*s,
        z*x*t - y*s,     z*y*t + x*s,     c + z*z*t
    );
}

void ToyShowerer::do_one_splitting(particle_queue& particles){
    /*
     * The strategy here will be as simple and stupid as I can bear
     * Plan: completely ignore the masses of the particles 
     *              and just track 3-momnta
     *  1. Grab the mother particle from the top of the pT-ordered queue
     *  2. Sample z and theta from some physically motivated distributions
     *  3. Sample phi from a distribution according to the phi_mode_
     *  4. Create two new particles with the appropriate 3-momenta
     *  5. Rotate them as needed
     *  6. Push them back onto the queue
     */
    const LightweightParticle mother = particles.top();
    particles.pop();

    double z = sample_z(); //momentum fraction carried by particle 1
    double theta = sample_theta(); //decay opening angle
    double phi = sample_phi(); //decay azimuthal angle

    /*
     * total momenta of the daughters:
     * p1 = z p
     * but p2 != (1-z) p 
     * Because that doesn't actually conserve momntum 
     * (the traditional zs are for energy)
     *
     * Instead, the math is
     * pvec = p1vec + p2vec
     * pvec^2 = (p1vec + p2vec)^2
     * p^2 = p1^2 + p2^2 + 2 p1 p2 cos(theta)
     * 1 = z^2 + y^2 + 2 z y cos(theta)
     * solving for y:
     * y = - z cos(theta) + sqrt[2 - z^2 (1 - cos(2 theta))]/sqrt[2]
     *
     * NB we do not have z + y = 1
     * Instead z + y >= 1 (not never that much greater)
     */
    double ct = std::cos(theta);
    double y = -z*ct + std::sqrt(2 - z*z*(1 - std::cos(2*theta)))/std::sqrt(2);

    double p = std::sqrt(mother.px()*mother.px() 
            + mother.py()*mother.py() 
            + mother.pz()*mother.pz());
    double p1 = z*p;
    double p2 = y*p;

    /*printf("Running a splitting on mother particle with p_tot = %f\n", p);
    printf("\tz = %f\n", z);
    printf("\ty = %f\n", y);
    printf("\ttheta = %f\n", theta);
    printf("\tphi = %f\n", phi);
    printf("\t\tp1 = %f\n", p1);
    printf("\t\tp2 = %f\n", p2);*/

    /*
     * need to get the angles of the daughters w.r.t. the mother axis
     * the math is:
     *  pvec = p1vec + p2vec
     *  p1vec dot pvec = p1vec dot p1vec + p1vec dot p2vec
     *  p p1 cos (alpha) = p1^2 + p1 p2 cos(theta)
     *  cos(alpha) = (p1 + p2 cos(theta)) / p
     *
     *  I solved this in mathematica to get the below expression
     */
    double cos_alpha = (p1 + p2*ct)/p;
    double alpha = std::acos(cos_alpha);
    //printf("\t\talpha = %f\n", alpha);

    ROOT::Math::XYZVector mothervec(mother.px(), mother.py(), mother.pz());
    ROOT::Math::XYZVector axisvec(mother.dpx(), mother.dpy(), mother.dpz());

    ROOT::Math::Rotation3D rotate_axis;
    get_rotation_about_axis(mothervec, phi, rotate_axis);
    axisvec = rotate_axis*axisvec;

    ROOT::Math::Rotation3D rotate_alpha;
    get_rotation_about_axis(axisvec, -alpha, rotate_alpha);
    ROOT::Math::XYZVector p1vec = rotate_alpha*(z*mothervec);
    //printf("\t\tp1vec.R() = %f\n", p1vec.R());

    ROOT::Math::Rotation3D rotate_beta;
    get_rotation_about_axis(axisvec, theta-alpha, rotate_beta);
    ROOT::Math::XYZVector p2vec = rotate_beta*(y*mothervec);
    //printf("\t\tp2vec.R() = %f\n", p2vec.R());

    //test
    /*double dot = p1vec.Dot(p2vec);
    double cos_angle = dot/(p1vec.R()*p2vec.R());
    //printf("\t\tangle between p1 and p2 = %f\n", std::acos(cos_angle));
    
    double dot2 = p1vec.Dot(mothervec);
    double cos_angle2 = dot2/(p1vec.R()*mothervec.R());
    //printf("\t\tangle between p1 and mother = %f\n", std::acos(cos_angle2));

    double dot3 = p2vec.Dot(mothervec);
    double cos_angle3 = dot3/(p2vec.R()*mothervec.R());
    //printf("\t\tangle between p2 and mother = %f\n", std::acos(cos_angle3));
    */

    particles.emplace(p1vec, axisvec);
    particles.emplace(p2vec, axisvec);
    //printf("Queue contents after splitting:\n");
    //print_queue_contents(particles);
    //printf("\n");
}

void ToyShowerer::shower(const double pt, 
                    const double eta,
                    const double phi,
                    const double mass,
                    const unsigned Npart,
                    jet& result){

    particle_queue particles;

    //random initial dipole axis
    double dtheta = uniform_dist_(rng_) * M_PI;
    double dphi = uniform_dist_(rng_) * 2*M_PI;

    double dx = std::sin(dtheta)*std::cos(dphi);
    double dy = std::sin(dtheta)*std::sin(dphi);
    double dz = std::cos(dtheta);
    
    //pt, eta, phi to cartesian
    double px = pt*std::cos(phi);
    double py = pt*std::sin(phi);
    double pz = pt*std::sinh(eta);

    particles.emplace(px, py, pz, dx, dy, dz);

    while(particles.size() < Npart){
        do_one_splitting(particles);
    }

    while (!particles.empty()){
        const auto& particle = particles.top();
        double pt = particle.pt();
        double phi = std::atan2(particle.py(), particle.px());
        double eta = std::asinh(particle.pz()/pt);
        result.particles.emplace_back(pt, eta, phi);

        result.pt += pt;
        result.sumpt += pt;
        result.rawpt += pt;

        result.eta += pt*eta;
        result.phi += pt*phi;

        particles.pop();
        ++result.nPart;
    }
    result.eta /= result.sumpt;
    result.phi /= result.sumpt;
}
