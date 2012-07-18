#define PROJ_UX(F, I1, I2,  J,  K) f(F, I1,  J,  K) + 0.5* LIM(rx(F, I1, I2,  J,  K)) * (f(F, I2,  J,  K) - f(F, I1,  J,  K))
#define PROJ_UY(F,  I, J1, J2,  K) f(F,  I, J1,  K) + 0.5* LIM(ry(F,  I, J1, J2,  K)) * (f(F,  I, J2,  K) - f(F,  I, J1,  K))
#define PROJ_UZ(F,  I,  J, K1, K2) f(F,  I,  J, K1) + 0.5* LIM(ry(F,  I,  J, K1, K2)) * (f(F,  I,  J, K2) - f(F,  I,  J, K1))


