#include "zp.h"

unsigned int ZP::SIMD_SIZE = 20;

long ZP::_prev_p = 2;

long ZP::_prev_r = 1;

std::function<long(void)> ZP::_getR = ZP::getPrevR;

std::function<long(void)> ZP::_getP = ZP::getPrevP;

