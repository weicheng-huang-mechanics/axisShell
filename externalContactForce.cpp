#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar)
{
	plate = &m_plate;
	stepper = &m_stepper;

    stiffness = m_stiffness;
    dBar = m_dBar;

    ifContact = 0;
    ifRebound = 0;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{
	for(int i = 0; i < plate->nv; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);

		double d = xCurrent(1);

		if (d <= dBar)
		{
			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			stepper->addForce(2 * i + 1, stiffness * dEnergydD);

			ifContact = 1;
		}
	}

	Vector2d xCurrent = plate->getVertex(0);

	if (ifContact == 1 && xCurrent(1) > 0.05)
	{
		ifRebound = 1;
	}
}

void externalContactForce::computeJc()
{
	for(int i = 0; i < plate->nv; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);

		double d = xCurrent(1);

		if (d <= dBar)
		{
			d2EnergydD2 = - 2 * log(d / dBar) - 2 * (d - dBar) / d - 2 * (d - dBar) / d + (d - dBar) * (d - dBar) / (d * d);

			stepper->addJacobian(2 * i + 1, 2 * i + 1, stiffness * d2EnergydD2);
		}
	}
}
