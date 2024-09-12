#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");		

	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");

	thickness = m_inputData.GetScalarOpt("thickness");
	Possion = m_inputData.GetScalarOpt("Possion");

	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");
	nv = m_inputData.GetIntOpt("nv");

	totalCompress = m_inputData.GetScalarOpt("totalCompress");
	inputPressure = m_inputData.GetScalarOpt("inputPressure");
	delta = m_inputData.GetScalarOpt("delta");

	height = m_inputData.GetScalarOpt("height");
	stiffness = m_inputData.GetScalarOpt("stiffness");
	dBar = m_inputData.GetScalarOpt("dBar");

	currentCompress = 0.0;
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDER";
    //name << "_nv_" << nv;
    //name << "_force_" << - gVector(1);
    name << "_height_" << height;
    name << "_delta_" << delta;
    name << "_thickness_" << thickness;
    name << "_dBar_" << dBar;
    name << "_stiffness_" << stiffness / 1e6;
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if ( timeStep % 10 != 0)
	{
		return;
	}

	if (currentTime > 10.0)
	{
		Vector2d xStart = plate->getVertex(0);

		outfile << currentTime - 10.0 << " " << xStart(1) << endl;

		//computeReactionForce();

		//Vector2d xEnd = plate->getVertex(0);

		//outfile << xEnd(1) - 1.0 << " " << reactionForce(plate->ndof - 1) << endl;

		//outfile << currentTime << " " << xEnd(1) << endl;

		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);

			//outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << endl;
		}
	}
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, thickness, Possion, 
		deltaTime, nv, delta, height);

	plateBoundaryCondition();

	plate->setup();

	stepper = new timeStepper(*plate);

	// set up force
	m_inertialForce = new inertialForce(*plate, *stepper);
	m_gravityForce = new externalGravityForce(*plate, *stepper, gVector);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosity);
	m_elasticStretchingForce = new elasticStretchingForce(*plate, *stepper);
	m_elasticBendingForce = new elasticBendingForce(*plate, *stepper);
	m_externalPressureForce = new externalPressureForce(*plate, *stepper);
	m_elasticBoundaryForce = new elasticBoundaryForce(*plate, *stepper);
	m_externalContactForce = new externalContactForce(*plate, *stepper, stiffness, dBar);

	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();
	m_elasticStretchingForce->setFirstJacobian();
	m_elasticBendingForce->setFirstJacobian();

	stepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;
}

void world::plateBoundaryCondition()
{
	Vector2d xStart = plate->getVertex(0);
	plate->setOneBoundaryCondition(xStart(0), 0, 0);

	//Vector2d xStart2 = plate->getVertex(plate->nv_temp_1 - 1);
	//plate->setOneBoundaryCondition(xStart2(1), plate->nv_temp_1 - 1, 1);



	Vector2d xEnd1 = plate->getVertex(plate->nv - 1);
	plate->setOneBoundaryCondition(xEnd1(0), plate->nv - 1, 0);
	
	Vector2d xEnd2 = plate->getVertex(plate->nv - 2);
	plate->setOneBoundaryCondition(xEnd2(0), plate->nv - 2, 0);
}

void world::updateTimeStep()
{
	bool goodSolved = false;

	while (goodSolved == false)
	{
		// Start with a trial solution for our solution x
		plate->updateGuess(); // x = x0 + u * dt

		updateEachStep();

		goodSolved = true;
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
}

void world::updateEachStep()
{
	if (currentTime > 2.0 && currentTime < 2.7)
	{
		m_externalPressureForce->pressure = m_externalPressureForce->pressure + 1000.0 * deltaTime;
	}

	if (currentTime > 2.0)
	{
		deltaTime = 1e-4;
		plate->dt = 1e-4;
	}

	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		m_dampingForce->computeFd();

		//if (currentTime > 20.0)
		{
			m_gravityForce->computeFg();
			m_externalContactForce->computeFc();
		}

		m_elasticStretchingForce->computeFs();
		m_elasticBendingForce->computeFb();
		m_externalPressureForce->computeFp();
		m_elasticBoundaryForce->computeFboundarty();
	
		normf = stepper->GlobalForceVec.norm();

		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		//cout << normf << " ";

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			m_gravityForce->computeJg();
			m_dampingForce->computeJd();

			m_elasticStretchingForce->computeJs();
			m_elasticBendingForce->computeJb();
			m_externalPressureForce->computeJp();
			m_elasticBoundaryForce->computeJboundarty();

			//if (currentTime > 20.0)
			{
				m_externalContactForce->computeJc();
			}
			
			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	//cout << endl;

	if (render)
	{
		cout << "iter " << iter << endl;
	}
}

int world::simulationRunning()
{
	Vector2d uTotal;

	uTotal(0) = 0.0;
	uTotal(1) = 0.0;

	for (int i = 0; i < plate->nv; i++)
	{
		Vector2d uCurrent = plate->getVelocity(i);

		uTotal = uTotal + uCurrent;
	}

	uTotal = uTotal / plate->nv;

	if (uTotal(1) < 0 && m_externalContactForce->ifRebound == 1) 
	{
		//return - 1;
	}
	else 
	{
		//return 1;
	}


	if (currentTime < 35.00)
	{
		return 1;
	}
	else
	{
		return - 1;
	}
}

Vector2d world::getScaledCoordinate(int i, int j)
{
	Vector2d xCurrent;
	
	if (j == 0)
	{
		xCurrent = plate->v_edgeElement[i].x_1 * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = plate->v_edgeElement[i].x_2 * scaleRendering;
	}

	return xCurrent;
}

Vector2d world::getScaledPos(int i)
{
	Vector2d xCurrent;

	xCurrent = plate->getVertex(i) * scaleRendering;

	return xCurrent;
}

int world::getNv()
{
	return plate->nv;
}

int world::numStretchingPair()
{
	return plate->edgeNum;
}

void world::computeReactionForce()
{
	//reactionForce = VectorXd::Zero(plate->ndof);

	//m_elasticStretchingForce->computeFs();
	//m_elasticBendingForce->computeFb();
	//m_externalPressureForce->computeFp();
	//m_elasticBoundaryForce->computeFboundarty();

	//reactionForce = - m_elasticBoundaryForce->totalForce - m_elasticStretchingForce->totalForce - m_elasticBendingForce->totalForce;
}
