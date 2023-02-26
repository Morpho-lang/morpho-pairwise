#include <stdio.h>
#include <morpho/morpho.h>
#include <morpho/builtin.h>

#define SELINGERNEMATIC_CLASSNAME      "SelingerNematic"

#define SELINGERNEMATIC_ARGS           "SlngrNmtcArgs"
#define SELINGERNEMATIC_ARGS_MSG       "SelingerNematic requires a field as the argument."

#define SELINGERNEMATIC_K11_PROPERTY   "k11"
#define SELINGERNEMATIC_K22_PROPERTY   "k22"
#define SELINGERNEMATIC_K33_PROPERTY   "k33"
#define SELINGERNEMATIC_K24_PROPERTY   "k24"

void selingernematic_initialize(void);