#ifndef VIFCOPULA_REV_CORE_PRECOMP_DVV_VARI_HPP
#define VIFCOPULA_REV_CORE_PRECOMP_DVV_VARI_HPP


#include <stan/math/rev/core/vari.hpp>
#include <dist/dvv_vari.hpp>

namespace vifcopula {

    // use for single precomputed partials
    class precomp_dvv_vari : public op_dvv_vari {
    protected:
      double db_;
      double dc_;
    public:
      precomp_dvv_vari(double val,
                       double ad,
                       vari* bvi, vari* cvi,
                       double db, double dc)
        : op_vv_vari(val, ad, avi, bvi),
          db_(db),
          dc_(dc) {
      }
      void chain() {
        bvi_->adj_ += adj_ * db_;
        cvi_->adj_ += adj_ * dc_;
      }
    };

  }

#endif
