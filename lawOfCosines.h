#ifndef SIMONTOOLS_LAWOFCOSINES_H
#define SIMONTOOLS_LAWOFCOSINES_H

namespace simon{
    template <typename T>
    inline T getcostheta(T a, T b, T c) {
      return (a*a + b*b - c*c) / (2*a*b);
    }

    template <typename T>
    inline T getC(T a, T b, T costheta) {
      return std::sqrt(a*a + b*b - 2*a*b*costheta);
    }
};

#endif
