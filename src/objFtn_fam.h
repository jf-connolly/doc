/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  fitness.h - Fitness evaluation object. Fitness is evaluated as the        *
 *  performances after training a classifier and it is possible to evaluate   *
 *  several objective at each iteration.                                      *
 *                                                                            *
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
#include "configuration.h"
#include "dbase.h"
#include "result.h"

#ifndef OBJECTIVEFUNC_H
#define OBJECTIVEFUNC_H

using namespace std;

#define fitEstimation   1
#define fitEvaluation   2

class objFtn_fam {
public:
	/*------------------------------ Variables -------------------------------*/
	//-- General
	int nDimensions,		//- Number of dimensions
		nObjectives,		//- Number of objectives
		nEstimationRep,		//- Number of estimation for fitness estimation
		nBlocks,			//- Number of blocs & current bloc
		nReplications,		//- Number of external replications
		typeSpace,			//- Learn space or test space
	    nAccPatterns;		//- Number of patterns used to TRAIN the networks

	//-- Input and output
	float *solution;		//-- Current position (optimization space)
	result fitness;

	//-- Before/after fitness estimation, validation & fitness neural networks
	artmap nnBfe, nnAfr, nnVal, nnFit;

	//-- Result and result used for validation
	result resultFt, resultVl;

	//-- Data bases
	dbase dbTrn, dbVal, dbFit;

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors & destructor
	objFtn_fam();
	~objFtn_fam();
	void initEstimation( configuration *cfg, int r );
	
	void loadDb  ( configuration *cfg, int time, int r );
	void deleteDb();

	//-- Utility functions
	void estimation();  	//- Fitness estimation

	//-- Training & testing
	void trainV( artmap *bestNN,	 artmap *currentNN,
				 dbase  *dbTrain,	 dbase  *dbValid,
				 result *bestResult, result *valResult );
	void train ( artmap *myNN, dbase *dbTrain );
	void test  ( artmap *myNN, dbase *dbTest, result *myResult );
};
#endif
