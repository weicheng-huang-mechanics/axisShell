#ifndef EXTERNALCONATCTFORCE_H
#define EXTERNALCONATCTFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar);
	~externalContactForce();

	void computeFc();
	void computeJc();

    int ifContact;
    int ifRebound;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double dEnergydD;
    double d2EnergydD2;

    double stiffness;
    double dBar;
};

#endif
