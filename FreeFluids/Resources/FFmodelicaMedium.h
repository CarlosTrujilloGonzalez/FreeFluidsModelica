#ifndef FFMODELICAMEDIUM_H
#define FFMODELICAMEDIUM_H

#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFphysprop.h"
#include "FFeosMix.h"
#include "FFactivity.h"

//Creates a static substances array that is charged each time it is called
void *FF_createSubstanceData(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP);

//Creates a static mixture array that is charged each time it is called
//void *FF_createMixData(const char *name, int numSubs, char *subsNames, const char *resDir, char *eosType, char *cubicMixRule, char *activityModel);
#endif // FFMODELICAMEDIUM_H
