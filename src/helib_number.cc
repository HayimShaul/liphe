#include "helib_keys.h"
#include "helib_number.h"

HelibKeys *HelibNumber::_prev_keys = NULL;
std::function<HelibKeys *(void)> HelibNumber::_getKeys = HelibNumber::getPrevKeys;
