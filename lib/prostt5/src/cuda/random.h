#pragma once

#ifdef CT2_WITH_CURAND
#include <curand_kernel.h>
#endif

namespace ctranslate2 {
  namespace cuda {

#ifdef CT2_WITH_CURAND
    curandStatePhilox4_32_10_t* get_curand_states(size_t num_states);
#endif

  }
}
