#ifndef KEP_TOOLBOX_LAMBERT_EXPOSIN_H
#define KEP_TOOLBOX_LAMBERT_EXPOSIN_H

#include <vector>

namespace kep_toolbox{
  class __KEP_TOOL_VISIBLE lambert_exposin{
  public:
    lambert_exposin(const array3D &r1, const array3D &r2, const double &tof, const double &mu, const int &lw, const int &multi_revs);
    const std::vector<array3D>& get_v1() const;
    const std::vector<array3D>& get_v2() const;
    const std::vector<double>& get_traversal_final_mass(const double &isp, const double &mempty, const double &m) const;
    const std::vector<double>& get_max_thrust(const double &isp, const double &mempty, const double &m) const;
  }
}

#endif
