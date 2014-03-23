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
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//-- Project specific
#include "result.h"

/*============================================================================*
 *  Function name			:	result::result
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Default constructor
 *----------------------------------------------------------------------------*
 *  Argument				:   int nCls	- Number of classes
 *============================================================================*/
result::result()
{
	error = 300;

	nClasses   		  = 0;
	nPatternsTest     = 0;
	nPatternsLearned  = 0;
	sizeSwarm 		  = 0;

	sizeEns = 0;
    members = (int*)NULL;

    predisemble = (int*)NULL;
	predictions = (int*)NULL;
    trueClasses = (int*)NULL;
    sizeClasses = (int*)NULL;
    membersTime = (int*)NULL;

    nClassOk = 0;
	nClassWg = 0;
	sizeNn   = 0;
	convTime = 0;

	clsRate  = 0;
	normCpn  = 0;
}
/*============================================================================*
 *  Function name			:	result::~result
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Destructor
 *============================================================================*/
result::~result()
{
	error = 300;

	sizeEns			  = 0;
	nClasses   		  = 0;
	nPatternsTest     = 0;
	nPatternsLearned  = 0;
	sizeSwarm 		  = 0;

	//-- Check if it is not a simplified version
	if(members    )   { delete [] members;     }
	if(predisemble)   { delete [] predisemble; }
	if(predictions)   { delete [] predictions; }
	if(trueClasses)   { delete [] trueClasses; }
	if(sizeClasses)   { delete [] sizeClasses; }
	if(membersTime)   { delete [] membersTime; }

	nClassOk = 0;
	nClassWg = 0;
	sizeNn   = 0;
	convTime = 0;

	clsRate  = 0;
	normCpn  = 0;
}
/*============================================================================*
 *  Function name			:	result::initialization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Initialization
 *----------------------------------------------------------------------------*
 *  Argument				:   int nCls	- Number of classes
 *								int nPats	- Number of patterns
 *============================================================================*/
void result::initialization(int nCls, int nPats, int szSwarm)
{
	/*------------------------------ Constants -------------------------------*/
	nClasses       = nCls;
	sizeDbTest     = nPats;
	nPatternsTest  = nPats;
	sizeSwarm      = szSwarm;

	/*--------------------------- Ensemble members ---------------------------*/
	//-- Memory allocation
	members = new int[sizeSwarm];
	if(members == NULL)   { error += 3; return; }

	//-- Initialize
	for(int e=0; e<sizeSwarm; e++)   { members[e] = -1;}

	/*--------------------- Predictions and true classes ---------------------*/
	//-- Memory allocation
	predictions = new int[nPatternsTest*sizeSwarm];
	predisemble = new int[nPatternsTest];
	trueClasses = new int[nPatternsTest];
	if(predictions == NULL || trueClasses == NULL)   { error += 3; return; }

	//-- Initialize
	for(int p=0; p<nPatternsTest*sizeSwarm; p++)   { predictions[p] = -1; }
	for(int p=0; p<nPatternsTest; p++)             { predisemble[p] = -1;
												 	 trueClasses[p] = -1; }

	/*-------------------------- Size of each class --------------------------*/
	//-- Memory allocation
	sizeClasses = new int[sizeSwarm*nClasses];
	if(sizeClasses == NULL)   { error += 3; return; }

	/*-------------------------- Convergence time  ---------------------------*/
	//-- Memory allocation
	membersTime = new int[sizeSwarm];
	if(membersTime == NULL)   { error += 3; return; }

	//-- Initialize
	for(int e=0; e<sizeSwarm; e++)    { membersTime[e] = 0; }
}
/*============================================================================*
 *  Function name			:	result::reinit
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Reinitialize a result
 *============================================================================*/
void result::reinit()
{
	error = 300;

	//-- For simplified results, sizeSwarm & nPatternsTest = 0;
	if(members     == NULL) { sizeSwarm = 0; }
	if(trueClasses == NULL) { nPatternsTest = 0; }

	for(int e=0; e<sizeSwarm;           e++)   { members[e]     = -1; }
	for(int p=0; p<nPatternsTest*sizeSwarm; p++)   { predictions[p] = -1; }
	for(int p=0; p<nPatternsTest;           p++)   { predisemble[p] = -1;
												 trueClasses[p] = -1; }
	for(int k=0; k<nClasses*sizeSwarm;  k++)   { sizeClasses[k] =  0; }
	for(int e=0; e<sizeSwarm;           e++)   { membersTime[e] =  0; }

	sizeEns = 0;

	nClassOk = 0;
	nClassWg = 0;
	sizeNn   = 0;
	convTime = 0;

	nPatternsTest    = 0;
	nPatternsLearned = 0;
	
	clsRate  = 0;
	normCpn  = 0;
}
/*============================================================================*
 *  Function name			:	result::copy
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Copies the result
 *----------------------------------------------------------------------------*
 *  Argument				:   result *original  - Original result
 *============================================================================*/
void result::copy(result* original)
{
	size_t sizeData;

	//-- Check for result existence
	if(original == NULL)   { error += 2; return; }

	//-- Simplified result
	nClassOk = original->nClassOk;
	nClassWg = original->nClassWg;
	sizeNn   = original->sizeNn;
	convTime = original->convTime;

	clsRate  = original->clsRate;
	normCpn  = original->normCpn;

	//-- Complete result
	nClasses 		 = original->nClasses;
	sizeDbTest 		 = original->sizeDbTest;
	nPatternsTest 	 = original->nPatternsTest;
	nPatternsLearned = original->nPatternsLearned;
	sizeSwarm 		 = original->sizeSwarm;

	//-- Ensemble members
	sizeEns = original->sizeEns;
	sizeData =  sizeof(int) *(unsigned)(sizeSwarm);
	memcpy(members, original->members, sizeData);

	//-- Classification for video sequences
	sizeData =  sizeof(int) *(unsigned)(nPatternsTest*sizeSwarm);
	memcpy(predictions, original->predictions, sizeData);

	sizeData =  sizeof(int) *(unsigned)(nPatternsTest);
	memcpy(predisemble, original->predisemble, sizeData);

	sizeData =  sizeof(int) *(unsigned)(nPatternsTest);
	memcpy(trueClasses, original->trueClasses, sizeData);

	//-- Class sizes
	sizeData =  sizeof(int) *(unsigned)(nClasses*sizeSwarm);
	memcpy(sizeClasses, original->sizeClasses, sizeData);

	//-- Ensemble convergence time
	sizeData =  sizeof(int) *(unsigned)(sizeSwarm);
	memcpy(membersTime, original->membersTime, sizeData);

	//-- Check
	if((sizeEns		  && members 	 == NULL) ||
	   (nPatternsTest && predictions == NULL) ||
	   (nPatternsTest && trueClasses == NULL) ||
	   (nClasses 	  && sizeClasses == NULL) )
	{   error += 3; return;   }
}
/*============================================================================*
 *  Function name			:	result::save
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Saves a result. If the file already exists,
 *  							the result is saved at the end of the file.
 *----------------------------------------------------------------------------*
 *  Arguments				:   char nameFile 	- Name of the file
 *============================================================================*/
void result::save(char * nameFile, int currentRep)
{
	//-- File opening
	ofstream f(nameFile, ios::app | ios::binary);
	if(!f)   { printf("Error! Can't open the file %s \n", nameFile);  exit(1); }

	//-- Writing the results
	f.write((char*)&currentRep,		  sizeof(int) );
	f.write((char*)&sizeDbTest,		  sizeof(int) );
	f.write((char*)&nClasses,		  sizeof(int) );
	f.write((char*)&nPatternsTest,    sizeof(int) );
	f.write((char*)&nPatternsLearned, sizeof(int) );
	f.write((char*)&sizeSwarm,		  sizeof(int) );

	f.write((char*)&sizeEns, sizeof(int)            );
	f.write((char*) members, sizeof(int)   *sizeEns );

	f.write((char*)&nClassOk,    sizeof(float)							);
	f.write((char*)&nClassWg,    sizeof(float)							);
	f.write((char*)&clsRate,     sizeof(float)							);
	f.write((char*)&normCpn,     sizeof(float)							);
	f.write((char*) predictions, sizeof(int)   *nPatternsTest * sizeEns );
	f.write((char*) predisemble, sizeof(int)   *nPatternsTest			);
	f.write((char*) trueClasses, sizeof(int)   *nPatternsTest			);
	f.write((char*) sizeClasses, sizeof(int)   *nClasses      * sizeEns );
	f.write((char*) membersTime, sizeof(int)   *sizeEns					);

	//-- Closing
	f.close();
}
/*----------------------------------------------------------------------------*/
