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
#include "artmap.h"

#ifndef ENSEMBLE_H
#define ENSEMBLE_H

using namespace std;

class ensemble {
public:
	/*------------------------------ Variables -------------------------------*/
	int      size,				//- Number of members in the ensemble
		    *members,			//- Members of the ensemble
	         prediction,		//- Final prediction
		    *vote;				//- Prediction of all the members

	float   *tieBreakerCl,		//- Tie breaker - Classification rate
			*tieBreakerSz;		//- Tie breaker - Network size

	artmap **nn;				//- Member (i.e. neural networks) addresses
	result **pFitness;			//- Fitness obtained with the neural networks
	dbase   *dbTest,			//- Actual data base used
			 dbInit;			//- For initialization purposes only
	result   performances;		//- Performances of the ensemble

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors & destructor
	ensemble();
	~ensemble();

	//-- Utility functions
	void initialization  (configuration *cfg);
	void reinitialization();
	void assignDbTest    (dbase *db);

	void test            ();
	void save            (char *nameFile, int currentRep, int currentBlock);
};
#endif
