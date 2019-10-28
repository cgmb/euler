#ifdef AE_FLTDEBUG
#include <fenv.h>
void enable_fpmath_asserts() {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
}
#else
void enable_fpmath_asserts() {
}
#endif
