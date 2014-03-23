/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  particle.h - Particle for canonical particle swarm optimizer			  *
 *  as defined in :														      *
 *																			  *
 *  J. Kennedy, Some issues and practices for particle swarms, Proc. of the   *
 *  IEEE Swarm Intelligence Symposium, 2007, pp. 162-169.					  *
 *   																		  *
 *  ----------------------------------------------------------------------    *
 *  This program is free software; you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the     *
 *  Free Software Foundation; either version 2 of the License, or (at your    *
 *  option) any later version.                                   			  *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful, but       *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of                *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         *
 *  General Public License for more details.                                  *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License along   *
 *  with this program; if not, write to the Free Software Foundation, Inc.,   *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                 *
 ******************************************************************************/

/*------------------------------- Header files -------------------------------*/
//-- Standard
//-- (none)

//-- Project specific
#include "configuration.h"
#include "objFtn_fam.h"

#ifndef PARTICLE_H
#define PARTICLE_H

using namespace std;

class particle {
public:
	/*------------------------------ Variables -------------------------------*/
	//-- General
	int nDimensions;
	int nObjectives;
	int *lbest;			//- Which particle affects this one

	float alpha, *beta;	//- Inertia weights


	//-- Arrays in the research space
	float  *sMin;		//- Minimal bound
	float  *sMax;		//- Maximal bound

	float  *s;	 		//- Current position (normalized)
	float  *sPrev;		//- Previous position (normalized)

	float **p;			//- Personal best position (normalized)
	artmap  nns, *nnp;	//- Neural network associated with the personal best

	//-- Arrays in the objective space
	result  sFitness; 	//- Current position's performance
	result *pFitness; 	//- Personal best position's performance

	//-- Evaluation object with all the databases
	objFtn_fam optSpace;

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors & destructor
	particle();
	~particle();

	//-- Utility functions
	void initialization(configuration *cfg);
	void randomization(int hard);
	void updateP();

	//-- Set & get real (unnormalized) position
	float  getsReal(int d);		//- Current position
	float  getpReal(int o, int d);		//- Personal best
};
#endif
