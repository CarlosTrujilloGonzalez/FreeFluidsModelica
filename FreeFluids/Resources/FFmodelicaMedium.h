#ifndef FFMODELICAMEDIUM_H
#define FFMODELICAMEDIUM_H

//Creates a static substances array that is charged each time it is called
void *FF_createSubstanceData(const char *name, const char *resDir, int thermoModel, int refState, double refT, double refP);

//Creates a static mixture array that is charged each time it is called
void *FF_createMixData(const char *name, int numSubs, const char *subsNames[15], const char *resDir, int thermoModel, int eosType, int mixRule, int activityModel, int refCalc, double refT, double refP);
#endif // FFMODELICAMEDIUM_H
