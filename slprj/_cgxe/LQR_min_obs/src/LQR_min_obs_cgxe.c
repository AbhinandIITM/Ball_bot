/* Include files */

#include "LQR_min_obs_cgxe.h"
#include "m_ez6b5k7jBQ3Ohnvv6n8cu.h"

unsigned int cgxe_LQR_min_obs_method_dispatcher(SimStruct* S, int_T method, void*
  data)
{
  if (ssGetChecksum0(S) == 1313378018 &&
      ssGetChecksum1(S) == 717081965 &&
      ssGetChecksum2(S) == 2798384841 &&
      ssGetChecksum3(S) == 1604043053) {
    method_dispatcher_ez6b5k7jBQ3Ohnvv6n8cu(S, method, data);
    return 1;
  }

  return 0;
}
