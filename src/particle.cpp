/******************************************************************************
 *   Copyright (C) 2010 by                                                    *
 *   Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca                      *
 *																			  *
 *   particle.cpp - Particle for canonical particle swarm optimizer			  *
 *   as defined in :														  *
 *																			  *
 *   J. Kennedy, Some issues and practices for particle swarms, Proc. of the  *
 *   IEEE Swarm Intelligence Symposium, 2007, pp. 162-169.					  *
 *   																		  *
 *   ----------------------------------------------------------------------   *
 *   This program is free software; you can redistribute it and/or modify it  *
 *   under the terms of the GNU General Public License as published by the    *
 *   Free Software Foundation; either version 2 of the License, or (at your   *
 *   option) any later version.                                   			  *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful, but      *
 *   WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *   General Public License for more details.                                 *
 *                                                                            *
 *   You should have received a copy of the GNU General Public License along  *
 *   with this program; if not, write to the Free Software Foundation, Inc.,  *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                *
 ******************************************************************************/

/*------------------------------- Header files -------------------------------*/
//-- Standard
#include <cstdlib>
#include <stdio.h>

//-- Project specific
#include "particle.h"

using namespace std;

/*============================================================================*
 *  Function name			:	particle::particle
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Constructor without arguments (default)
 *============================================================================*/
particle::particle()
{
	//-- General
	nDimensions = 0;
	nObjectives = 0;
	lbest       = (int*)NULL;

	//-- Inertia
	alpha = 0;
	beta  = (float*)NULL;

	//-- Arrays in the research space
	sMin  = (float*) NULL;
	sMax  = (float*) NULL;
	s     = (float*) NULL;
	sPrev = (float*) NULL;
	p     = (float**)NULL;

	//-- Arrays in the objectives space
	pFitness = (result*)NULL;
}
/*============================================================================*
 *  Function name			:	particle::particle
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Initialization
 *----------------------------------------------------------------------------*
 *  Arguments				:   configuration * cfg  -  Configuration file
 *============================================================================*/
void particle::initialization(configuration *cfg)
{
	//-- General
	nDimensions = cfg->nDimensions;
	nObjectives = cfg->nObjectives;			//- To be modified -> config file
	alpha       = cfg->alpha;
	beta        = new float[nObjectives];
	for(int o=0; o<nObjectives; o++)   { beta[o] = cfg->beta; }

	lbest       = new int[nObjectives];
	for(int o=0; o<nObjectives; o++)   { lbest[o] = -1; }

	//-- Arrays in the research space
	sMin = new float[nDimensions];
	sMax = new float[nDimensions];
	for(int d=0; d<nDimensions; d++)
	{	sMin[d] = cfg->sMin[d];
		sMax[d] = cfg->sMax[d];
	}
    s        = new float [nDimensions];
    sPrev    = new float [nDimensions];
	p        = new float*[nObjectives];
	for(int o=0; o<nObjectives; o++)   { p[o] = new float[nDimensions]; }
    
    char nameFile[100];
	dbase dbTemp;
	sprintf(nameFile, "%s/%s/test.db", cfg->pathDb, cfg->nameDb);
	dbTemp.initialiHeader(nameFile);

	//-- Neural networks (not for theoritical functions)
	nns.initialization( dbTemp.nClasses, dbTemp.nFeatures,
			            cfg->sizeF2max, cfg->nEpochMax);
	nnp = new artmap[nObjectives];
	if(dbTemp.error == 200)
	{	for(int o=0; o<nObjectives; o++)
		{   nnp[o].initialization( dbTemp.nClasses, dbTemp.nFeatures,
								   cfg->sizeF2max, cfg->nEpochMax);   }
	}

	//-- Results in simplified mode (default constructor)!
    pFitness = new result[nObjectives];
}
/*============================================================================*
 *  Function name			:	particle::~particle
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Destructor
 *============================================================================*/
particle::~particle()
{
	//-- General
	nDimensions = 0;
	alpha       = 0;
	delete [] beta;

	//-- Arrays in the research space
	delete [] lbest;
	delete [] sMin;
	delete [] sMax;
	delete [] s;
	delete [] sPrev;
	for(int o=0; o<nObjectives; o++)   { delete [] p[o]; }
	delete [] p;
	delete [] pFitness;
	delete [] nnp;
}
/*============================================================================*
 *  Function name			:	particle::initialization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Random initialization of a particle (can also
 *  							be used for re-initialization).
 *============================================================================*/
void particle::randomization(int hard)
{
	float random;

	//-- Position initialization (pBest = initial position)
	for (int i=0; i<nDimensions; i++)
	{	//-- Current position
		random   = ((float)rand() /(float)RAND_MAX);
     	s[i] = random;
		for(int o=0; o<nObjectives; o++)   { p[o][i] = s[i]; }
		
		//-- Previous position (similiparameter: division factor for sprev)
		random   = ((float)rand() /(float)RAND_MAX);
		sPrev[i] = s[i] + (random-s[i])/4;
	}

	//-- Best artmap and results reinitialization
	//   (NULL for theoritical functions)
	for(int o=0; o<nObjectives; o++)
	{	nnp[o].reinit();
		nnp[o].setHp(getsReal(0),getsReal(1),getsReal(2),getsReal(3));
	}

	//o---------------------- New fitness and networks ----------------------o//
	//-- Learn the new position
	for(int d=0; d<nDimensions; d++)
	{   optSpace.solution[d] = getsReal(d);   }

	//-- Soft - we change position
	if(!hard)
	{   optSpace.estimation();   }
	//-- Hard - we erase every thing and relearn on new data
	else
	{
		optSpace.nnBfe.reinit();
		optSpace.estimation();
	}

	sFitness.copy( &(optSpace.fitness) );
	nns.copy( &(optSpace.nnAfr) );

	//-- Copy to objectives
	for(int o=0; o<nObjectives; o++)
	{
		pFitness[o].copy( &(optSpace.fitness) );
		nnp[o].copy( &(optSpace.nnAfr) );
	}
}
/*============================================================================*
 *  Function name			:	particle::updateP
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Updates cognitive influence if needed
 *----------------------------------------------------------------------------*
 *  Argument				:	artmap* nns  -  Network of the current position
 *============================================================================*/
void particle::updateP()
{
	//-- Performance (check classification rate)
	if( sFitness.clsRate > pFitness[0].clsRate )
	{
		for(int d=0; d<nDimensions; d++) {  p[0][d] = s[d];   }
		pFitness[0].copy(&sFitness);
		nnp[0].copy(&nns);
	}
	
	//-- Compression (check network size)
	if( sFitness.normCpn > pFitness[1].normCpn )
	{
		for (int d =0; d<nDimensions; d++) {  p[1][d] = s[d];   }
		pFitness[1].copy(&sFitness);
		nnp[1].copy(&nns);
	}
}
/*============================================================================*
 *  Function name			:	particle::set & get stuff
 *  Originally written by	:	Marcelo Kapp
 *  Modified by				:	Jean-Francois Connolly
 *  Description				:	Sets and gets particle's current position,
 *  							best position, and their performances.
 *----------------------------------------------------------------------------*
 *  Argument				:	depends
 *----------------------------------------------------------------------------*
 *  Returns					:   depends
 *============================================================================*/
float particle::getsReal(int d)
{   return sMin[d] + s[d]*(sMax[d]-sMin[d]);      }

float particle::getpReal(int o, int d)
{   return sMin[d] + p[o][d]*(sMax[d]-sMin[d]);   }
