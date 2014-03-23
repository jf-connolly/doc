/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  result.h - A result                                                       *
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
//-- (none)

#ifndef RESULT_H
#define RESULT_H

using namespace std;

#define staticImages   1
#define videoSequences 2

class result {
public:
	/*------------------------------ Variables -------------------------------*/
	int   error;			//- Error code

	int	  nClasses,			//- Number of classes
		  sizeDbTest,       //- Total number of test patterns
		  nPatternsTest,	//- Current number of test patterns
		  nPatternsLearned, //- Current number of leaning patterns patterns
		  sizeSwarm, 		//- Size of the swarm (for initialization purposes)

		  sizeEns,			//- Number of ensemble members
		 *members,			//- Member of the ensemble

		 *predictions,		//- Matrix of predictions (for video)
		 *predisemble,		//- Ensemble's predictions
	     *trueClasses,		//- Vector of true classes (for video)
		 *sizeClasses,		//- Size of each classes

		 *membersTime;		//- Convergence time of the ensemble

	//-- Simplified results for fitness estimation
	float nClassOk,			//- Number of correct classifications
		  nClassWg,			//- Number of wrong classifications
		  sizeNn,			//- Size of the neural network
		  convTime,			//- Convergence time
		  clsRate,			//- Classfication rate
		  normCpn;			//- Normalized compression

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors & destructor
	//-- Resultats
	result();
	~result();

	//-- Functions: utilities
	void initialization(int nCls, int nPats, int szSwarm); //- Mem. allocation
	void reinit();										   //- Reinitialization
	void save(char* nameFile, int currentRep);			   //- Save to bin. file
	void copy(result* original);						   //- Copy a result
};
#endif
/*----------------------------------------------------------------------------*/
