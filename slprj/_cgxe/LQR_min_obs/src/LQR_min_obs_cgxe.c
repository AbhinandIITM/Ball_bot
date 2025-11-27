/* Include files */

#include "LQR_min_obs_cgxe.h"
#include "m_8EHFMTFjUjqA4ydRIrPDXB.h"

unsigned int cgxe_LQR_min_obs_method_dispatcher(SimStruct* S, int_T method, void*
  data)
{
  if (ssGetChecksum0(S) == 1060106614 &&
      ssGetChecksum1(S) == 3676754321 &&
      ssGetChecksum2(S) == 1083559762 &&
      ssGetChecksum3(S) == 2849408887) {
    method_dispatcher_8EHFMTFjUjqA4ydRIrPDXB(S, method, data);
    return 1;
  }

  return 0;
}
