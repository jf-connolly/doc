/******************************************************************************
 *  Copyright (C) 2008 by Marcelo Kapp - kapp@livia.etsmtl.ca				  *
 *  Modify by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca			  *
 *                                                                            *
 *  configuration.h - Configuration file with all parameters setting		  *
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
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

//-- Project specific
#include "configuration.h"

/*============================================================================*
 *  Function name			:	configuration::configuration
 *  Originally written by	:	Marcelo Kapp
 *  Modified by				:	Jean-Francois Connolly
 *  Description				:	Constructor with arguments (default)
 *----------------------------------------------------------------------------*
 *  Argument				:   char *nameFile    - File config.txt
 *============================================================================*/
configuration::configuration(char * nameFile)
{
	size_t sizeData = sizeof(char);

	error = 0;

	//-- File opening
	ifstream f(nameFile);
	if(!f)
	{	printf("Error! Can't open the file %s \n", nameFile);
		error = error + 1;
	}

	f.seekg(0,ios::beg);
	//-- Reading the general informations
	f.seekg(160*sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> nameDb;
	f.seekg(40 *sizeData,ios::cur);   f >> nameSc;
	f.seekg(40 *sizeData,ios::cur);   f >> incremental;
	f.seekg(40 *sizeData,ios::cur);   f >> nBlocks;

	f.seekg( 4 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> startingRep;
	f.seekg(40 *sizeData,ios::cur);   f >> nReplications;
	f.seekg(40 *sizeData,ios::cur);   f >> nEstimationRep;

	f.seekg( 4 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> pathDb;
	f.seekg(40 *sizeData,ios::cur);   f >> pathSs;

	//-- Reading the informations for the optimization method
	f.seekg(36 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> nDimensions;
	f.seekg(40 *sizeData,ios::cur);   f >> nObjectives;
	f.seekg(40 *sizeData,ios::cur);   f >> nIterationsMax;
	f.seekg(40 *sizeData,ios::cur);   f >> nIterationsOvr;

	//-- ... rest of the inf for pso
	f.seekg( 4 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> alpha;
	f.seekg(40 *sizeData,ios::cur);   f >> beta;

	f.seekg( 4 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> sizeSwarm;
	f.seekg(40 *sizeData,ios::cur);   f >> sizeHood;
	f.seekg(40 *sizeData,ios::cur);   f >> nSSmax;
	f.seekg(40 *sizeData,ios::cur);   f >> sizeSSmax;
	f.seekg(40 *sizeData,ios::cur);   f >> dSSmin;
	dSSmin = pow(dSSmin,2);
	f.seekg(40 *sizeData,ios::cur);   f >> vFreeMin;
	vFreeMin = pow(vFreeMin,2);

	//-- Memory allocation for the hyperparameters boundaries
	sMin      = new float[nDimensions];
	sMax      = new float[nDimensions];

	f.seekg( 4 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);
		for(int i=0; i<nDimensions; i++) { f >> sMin[i]; }
	f.seekg(40 *sizeData,ios::cur);
		for(int i=0; i<nDimensions; i++) { f >> sMax[i]; }

	//-- ARTMAP informations
	f.seekg(16 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> sizeF2max;
	f.seekg(40 *sizeData,ios::cur);   f >> nEpochMax;

	//-- Archive informations
	f.seekg(16 *sizeData,ios::cur);
	f.seekg(40 *sizeData,ios::cur);   f >> nMemetics;
	f.seekg(40 *sizeData,ios::cur);   f >> sizeMemetics;
	f.seekg(40 *sizeData,ios::cur);   f >> widthMemetics;

	f.close();
}
/*============================================================================*
 *  Function name			:	configuration::configuration
 *  Originally written by	:	Marcelo Kapp
 *  Modified by				:	Jean-Francois Connolly
 *  Description				:	Destructor
 *============================================================================*/
configuration::~configuration()
{
	delete [] sMin;
	delete [] sMax;
}
