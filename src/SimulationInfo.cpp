#include "InitialPositions.h"
#include <vector>
template <typename T>

struct SimulationStats{
    int jd;
    T initialEnergy;
    T initialMomentum;
    T currentEnergy;
    T currentMomentum;
}

enum numericalIntegrationMethods{
    Euler,
    EulerCromer,
    Verlet,
    Leapfrog,
    Yoshida4thOrder,
    RungeKatta,
    VEFRL,
    PEFRL
}

