/******************************************************************************
 *   Copyright(C) 2008 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca *
 *   based on the work done by Marcelo Kapp									  *
 *   																		  *
 *   optimization.cpp - Optimization method									  *
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

//o------------------------------ Header files ------------------------------o//
//-- Standard
#include <cstdlib>
#include <fstream>
#include <math.h>

//-- Project specific
#include "multiPso.h"

/*============================================================================*
 *  Function name			:	optimization::optimization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Constructor without arguments (default)
 *============================================================================*/
multiPso::multiPso()
{
	//-- Constants
	sprintf(name, "%s", "multi");
	nDimensions = 0;
	size      	= 0;
	nIterations = 0;
	gbest 		= (int*)NULL;

	p = (particle*)NULL;

	//-- For local best topology
	sizeHood     = 0;
	hood         = (int*)   NULL;
	distance     = (float**)NULL;

	//-- For DNPSO
	sizeSSmax = 0;
	nSSmax    = 0;
	vFreeMin  = 0;
	dSSmin    = 0;

	fitness  = (float**)NULL;
	lbests   = (int  **)NULL;
    nSS      = 0;
	members  = (int  **)NULL;
	nMembers = (int  **)NULL;
	freep    = (int  **)NULL;
	nfree    = 0;

	indicator = (float*)NULL;
}
/*============================================================================*
 *  Function name			:	optimization::optimization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Constructor without arguments (default)
 *----------------------------------------------------------------------------*
 *  Argument						nItrs = 0;   bestItr = 0;
	:   Configuration *configPso
 *============================================================================*/
void multiPso::initialization(configuration * cfg)
{
	//-- PSO
	nDimensions = cfg->nDimensions;
	nObjectives = cfg->nObjectives;
	size   		= cfg->sizeSwarm;

	gbest = new int[nObjectives];
	for(int o=0; o<nObjectives; o++)   { gbest[o] = -1; }

	//-- Particles, objective functions, and particle position
	p = new particle[size];
	for(int n=0; n<size; n++)   { p[n].initialization(cfg); }

	//-- Same alpha for everyone, and different beta depending on the group
	for(int n=0; n<size; n++)
	{	p[n].alpha = cfg->alpha;
		for(int o=0; o<nObjectives; o++)   { p[n].beta[o] = cfg->beta; }
	
//		 p[n].beta[1] = 0;  //- for dnpso

		//-- Specialization of the particles
//		if(n < size/3)      { p[n].beta[1] = cfg->beta/2; }
//		if(n >= size*2/3)   { p[n].beta[0] = cfg->beta/2; }
	}

	//-- Archive - performed externally as it needs the replication
	myArchive.initialization(cfg);
	
	//-- Ensemble
	myEnsemble.initialization(cfg);
	int sizeInd = (int)max(size, myArchive.size);
	indicator   = new float[sizeInd];

	//o------------------------------ Subswarms -----------------------------o//
	//-- Local best topology
	sizeHood     = cfg->sizeHood;
	hood         = new int   [sizeHood];
	distance     = new float*[size];
	for(int i=0; i<size; i++)   { distance[i] = new float[size]; }

	//-- Maximal distance between two points
	distanceMax = nDimensions;

	//-- Parameters
	nSSmax     = cfg->nSSmax;
	if(nSSmax == 1) { sprintf(name, "%s", "pso"); }
	sizeSSmax  = cfg->sizeSSmax;
	vFreeMin   = cfg->vFreeMin;
	dSSmin     = cfg->dSSmin;

	//-- Memory allocation
	nSS		 = new int   [nObjectives];
	nfree    = new int   [nObjectives];
	fitness  = new float*[nObjectives];
	lbests   = new int*  [nObjectives];
	nMembers = new int*  [nObjectives];
	members  = new int*  [nObjectives];
	freep    = new int*  [nObjectives];
	for(int i=0; i<nObjectives; i++)
	{	fitness [i]= new float[size];
		lbests  [i]= new int  [size];
		nMembers[i]= new int  [size];
		members [i]= new int  [size*size];
		freep   [i]= new int  [size];
	}

	//-- Initialization
	for(int o=0; o<nObjectives; o++)
	{	for(int i=0; i<size; i++)
		{	lbests  [o][i] = -1;
			nMembers[o][i] = 0;
			fitness [o][i] = 0;
			for(int j=0; j<size; j++) { members[o][i*size +j] = -1; }
			freep   [o][i] = -1;
		}
		nSS  [o] = 0;
		nfree[o] = 0;
	}
											
	for(int o=0; o<nObjectives; o++) {   nSS[o] = 0;   nfree[o] = 0;   }

	//-- Sorting object (Size of the swarm for the multipeak benchmark)
	char nameFile[100];
	dbase dbTest;
	sprintf(nameFile, "%s/%s/test.db", cfg->pathDb, cfg->nameDb);
	dbTest.initialiHeader(nameFile);
	if(dbTest.error > 200)   { sort.initialization(size);        }
	else 					 { sort.initialization(dbTest.nPatterns); }
}
/*============================================================================*
 *  Function name			:	optimization::~optimization
 *  Originally written by	:	Marcelo Kapp
 *  Modified by				:	Jean-Francois Connolly
 *  Description				:	PSO destructor
 *============================================================================*/
multiPso::~multiPso()
{
	//-- Particles
	delete [] p;
	delete [] gbest;

	//-- Subswarms
	delete [] nSS;
	delete [] nfree;
	for(int o=0; o<nObjectives; o++)
	{
		delete [] fitness [o];
		delete [] lbests  [o];
		delete [] nMembers[o];
		delete [] members [o];
		delete [] freep   [o];
	}
	delete [] fitness;
	delete [] lbests;
	delete [] nMembers;
	delete [] members;
	delete [] freep;

	//-- For local best topology
	sizeHood = 0;
	delete [] hood;
	delete [] distance;

	//- Ensemble Diversity
	delete [] indicator;

	//-- Variables - Swarm
	nDimensions = 0;
	size   = 0;
	gbest       = 0;

	//-- Variables - DNPSO
	sizeSSmax = 0;
	nSSmax    = 0;
	dSSmin = 0;
	vFreeMin  = 0;
}
/*============================================================================*
 *  Function name			:	optimization::defineSubSwarm
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Defines the subswarms, i.e. the masters,
 *  							followers and free particles (algorithm of
 *  							Figure 1 in Nickabadi and al., "DNPSO",
 *  							p. 2642, 2008)
 *============================================================================*/
void multiPso::defineSubSwarms()
{
  	//-- Distance between particles
	for(int i=0; i<size; i++)
	{   for(int j=0; j<i; j++)
		{	distance[i][j] = euclidean2(&p[i],&p[j]);
			distance[j][i] = distance[i][j];
		}
		distance[i][i] = 0;
	}
	//o------------------------ Assign fitness values -----------------------o//
	//-- Assign for the descending sort!
	for(int n=0; n<size; n++)
	{
		fitness[0][n] = 0 - p[n].pFitness[0].clsRate;
		fitness[1][n] = 0 - p[n].pFitness[1].normCpn;
	}
	
	//o------------------------- Define the masters -------------------------o//
	for(int o=0; o<nObjectives; o++)
	{	//-- Initialized without any masters
		for(int n=0; n<size; n++)
	    {   lbests[o][n] = -1;   p[n].lbest[o] = -1;   }

	  	//-- Distance^2 between particles *********************************
		for(int i=0; i<size; i++)
		{   for(int j=0; j<i; j++)
			{	distance[i][j] = euclidean2p(i,j,o);
				distance[j][i] = distance[i][j];
			}
			distance[i][i] = 0;
		}

		//-- Define the masters (particles who are their own local best)
		//-- Finds the local best and assign if necessary
		nSS[o] = 0;
		for(int n=0; n<size; n++)
		{
			findlBest(n,o);
			if( p[n].lbest[o] == n ) { lbests[o][nSS[o]] = n;   nSS[o]++; }
		}

		//-- Merge the masters if to close
		for(int i=0; i<nSS[o]; i++)
		{
			for(int j=i+1; j<nSS[o]; j++)
			{
				//-- If they are to close ... sort!
				if(distance[ lbests[o][i] ][ lbests[o][j] ] < dSSmin)
				{	int first  = o,			        //- 0 -> 0, and 1 -> 1
						second = (int)pow(o-1,2);   //- 0 -> 1, and 1 -> 0

					sort.size          = 2;
					sort.order     [0] = 0;
					sort.fitness   [0] = fitness[first ][ lbests[o][i] ];
					sort.tieBreaker[0] = fitness[second][ lbests[o][i] ];

					sort.order     [1] = 1;
					sort.fitness   [1] = fitness[first ][ lbests[o][j] ];
					sort.tieBreaker[1] = fitness[second][ lbests[o][j] ];

					sort.heapSort();

					//-- redefine the master depending on the sort.
					int better, worst;
					if(sort.order[0])	{ better = j;  worst = i; }
					else				{ better = i;  worst = j; }

					p[ lbests[o][worst] ].lbest[o] = lbests[o][better];

					//-- ... eliminate it from the list!
					for(int k=worst; k<nSS[o]; k++)
					{   lbests[o][k] = lbests[o][k+1];   }
					nSS[o]--;
				}
			}
		}
	}
	//o------------------------ Define the followers ------------------------o//
	for(int o=0; o<nObjectives; o++)
	{
		//-- Empty at first, and ...
		for(int n=0; n<size; n++)
		{	nMembers[o][n] = 0;
			for(int m=0; m<size; m++)  { members[o][n*size +m] = -1; }
		}

		//-- ... add the masters (the +0)
		for(int m=0; m<nSS[o]; m++)
		{   members [o][lbests[o][m]*size +0] = lbests[o][m];
		    nMembers[o][lbests[o][m]        ] = 1;
		}

		for(int n=0; n<size; n++)
		{	//-- If the local best is not a master is is redefined
			//   by its local best until a master is reached
			if(p[n].lbest[o] != n)
			{	int m = p[n].lbest[o];
				//-- ****possible_infinit_loop****
				while(p[m].lbest[o] != m)   { m = p[m].lbest[o]; }

				//-- Add a follower
				members [o][ m*size +nMembers[o][m] ] = n;
				nMembers[o][m]++;
			}
		}
	}
	//o--------- Adjust the number of followers in each sub swarms ----------o//
	//o----------------------- (if there are to many) -----------------------o//
	for(int o=0; o<nObjectives; o++)
	{
		for(int m=0; m<nSS[o]; m++)
		{	int clbest = lbests[o][m];
			while(nMembers[o][clbest] > sizeSSmax)
			{
				//-- Find the farthest ...
				float worstDist =  0;
				int   worst     = -1;
				for(int k=1; k<nMembers[o][clbest]; k++)
				{   int testedMember = members[o][clbest*size +k];
					if(distance[clbest][testedMember] > worstDist)
					{ worstDist = distance[clbest][testedMember];   worst = k; }
				}
				//-- and eliminate it!
				for(int k=worst; k<nMembers[o][clbest]; k++)
				{	int addr = clbest*size +k;
					members[o][addr] = members[o][addr+1];
				}
				nMembers[o][clbest]--;
			}
		}
	}
	//o---------------------- Define the free particles ---------------------o//
	for(int o=0; o<nObjectives; o++)
	{
		//-- Initialization
		for(int i=0; i<size; i++)   { freep[o][i] = -1; }
		nfree[o] = 0;

		for(int i=0; i<size; i++)
		{
			//-- Search amongst the subswarms, and ...
			int flagfree = 1;
			for(int j=0; j<nSS[o]; j++)
			{	//-- ... amongst the members
				int clbest = lbests[o][j];
				for(int k=0; k<nMembers[o][clbest]; k++)
				{	if(i == members[o][ clbest*size +k ])
					{	flagfree = 0;   }
				}
			}
			//-- If free
			if(flagfree)   { freep[o][nfree[o]] = i;   nfree[o]++; }
		}
	}
	//o----------------------- Redefine the lbest ... -----------------------o//
	for(int o=0; o<nObjectives; o++)
	{
		for(int j=0; j<nSS[o]; j++)
		{	//-- ... of the masters , and ...
			p[ lbests[o][j] ].lbest[o] = lbests[o][j];

			//-- ... the followers.
			for(int k=0; k<nMembers[o][lbests[o][j]]; k++)
			{	int cfollower = members[o][ lbests[o][j]*size +k ];
				p[ cfollower ].lbest[o] = lbests[o][j];
			}
		}

		//-- Free particle don't have any lbest!
		for(int i=0; i<nfree[o]; i++)   { p[freep[o][i]].lbest[o] = -1; }
	}
}

/*============================================================================*
 *  Function name			:	optimization::moveSwarm
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Moves the particle swarm as defined in
 *  							algorithm of Figure 2 in Nickabadi and al.,
 *  							"DNPSO", p. 2642, 2008 (1 iteration).
 *============================================================================*/
void multiPso::movePopulation()
{
	for(int n=0; n<size; n++)
	{
		for(int o=0; o<nObjectives; o++)
		{
			//o------------------------- Influences -------------------------o//
			//-- Cognitive influence (i.e. personal best)
			float betaCg = p[n].beta[o] / (nObjectives*2);

			//-- Social influence - only if has a local best
			//   (i.e. != free particle)
			float betaSc = 0;
			if(p[n].lbest[o] > -1)  { betaSc = p[n].beta[o] / (nObjectives*2); }

			//o----------- Check minimal distance with local best -----------o//
			//-- Test for distance and better personal best
			//   (does not apply for free particles and local bests)
			float testDistance = 0;  //- Automatic pass for free particles
			int   lb           = p[n].lbest[o];
			if(lb>-1 && lb!=n)
			{	//-- If to pbest is to near form lbest ...
				testDistance = euclidean2p(n, lb, o);
				if(testDistance < dSSmin )
				{	//-- ... erase the particle's memory, and
					p[n].nnp[o].reinit();

					//-- ... remove all influences form that objective
					betaCg = 0;   betaSc = 0;
				}
			}

			//o--------------------- Moving a particle ----------------------o//
			for(int d=0; d<nDimensions; d++)
			{
				float currPos  = p[   n           ].s    [d],
					  prevPos  = p[   n           ].sPrev[d],
					  pbestPos = p[   n           ].p[o][d],
					  lbestPos, r;

				//-- Local best (if exist, i.e. not a free particle)
				if(p[n].lbest[o] > -1)  { lbestPos = p[p[n].lbest[o]].p[o][d]; }

				//o-- Move 0/3 - new position -------------------------------o//
				float newPos = currPos;

				//o-- Move 1/3 - inertia ------------------------------------o//
				newPos += p[n].alpha *(currPos  - prevPos);

				//o-- Move 2/3 - Cognitive influence (i.e. personal best) ---o//
				r = ((float)rand() / ( (float)(RAND_MAX)+(float)(1) ));
				newPos += r *betaCg *(pbestPos - currPos);

				//o-- Move 3/3 - Social influence (i.e. local best) ---------o//
				//-- If it is a free particle -> remove influence
				if(p[n].lbest[o] > -1)
				{	r = ((float)rand() / ( (float)(RAND_MAX)+(float)(1) ));
					newPos += r *betaSc *(lbestPos - currPos);
				}

				//-- Set the new position
				p[n].s    [d] = newPos;
				p[n].sPrev[d] = currPos;
			}
		}

		//-- For free particles (class. rate only), check for reinit.
		//   with the velocity
		if(p[n].lbest[0] < 0 && p[n].lbest[1] < 0)
		{
			//-- Distance between old and new position)
			float dTest = 0;
			for(int d=0; d<nDimensions; d++)
			{   dTest += pow( p[n].s[d]-p[n].sPrev[d], 2 );   }

			//-- Test
			if(dTest < vFreeMin)      { p[n].randomization(0); }
		}
	}
}
/*============================================================================*
 *  Function name			:	optimization::updatePopulation
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Move randomly free particles and update each 
 * 								particle current & best position fitness 
 * 								(necessary for dynamic environment)
 *============================================================================*/
void multiPso::reevaluateSwarm()
{
	//o------------------------- Parallelized update ------------------------o//
 	#pragma omp parallel for default(shared) schedule(dynamic)
	for(int n=0; n<size; n++)
	{
		//-- Define pointers for parallelization (it will be copied)
		particle *pPtr = &p[n];

		//o------------------ Move randomly free particles ------------------o//
		if(pPtr->lbest[0] < 0 && pPtr->lbest[1] < 0)
		{	p[n].randomization(0);   }
		
    	//o---------------------- The current position ----------------------o//
		for(int d=0; d<nDimensions; d++)
		{   pPtr->optSpace.solution[d] = pPtr->getsReal(d);   } //- Position
		pPtr->optSpace.estimation();           //- Fitness est. * Bottle neck *
		pPtr->sFitness.copy( &(pPtr->optSpace.fitness) );       //- New fitness
		pPtr->nns.copy( &(pPtr->optSpace.nnAfr) );      		//- New NN
		
    	//o----------------------- The pbest solution -----------------------o//
		for(int o=0; o<nObjectives; o++)
		{
			for(int d=0; d<nDimensions; d++)
			{ pPtr->optSpace.solution[d] = pPtr->getpReal(o,d); } //- Position
			pPtr->optSpace.estimation();        //- Fitness est. * Bottle neck *
			pPtr->pFitness[o].copy( &(pPtr->optSpace.fitness) );  //- Fitness
			pPtr->nnp[o].copy( &(pPtr->optSpace.nnAfr) );		  //- New NN
		}
	}
}
/*============================================================================*
 *  Function name			:	optimization::updatePopulation
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Update population's fitness
 *============================================================================*/
void multiPso::updatePopulation()
{
	//-- Parallelized update
	#pragma omp parallel for default(shared) schedule(dynamic)
	for(int n=0; n<size; n++)
	{
		//-- Define pointers for parallelization
		particle *pPtr = &p[n];

		//-- Check if the particle is within boundaries
		//   (if not, fitness is not evaluated)
		int flagOpt = 1;
		for(int d=0; d<nDimensions; d++)
		{   if( pPtr->s[d] > 1  ||  pPtr->s[d] < 0 ) { flagOpt = 0; }   }

		if(flagOpt)
		{	//--  the new solution (unnormalized position)
			for(int d=0; d<nDimensions; d++)
			{	float newSolution = pPtr->getsReal(d);
				pPtr->optSpace.solution[d] = newSolution;
			}
						
			//-- Estimate and update the fitness  ***** Bottle neck *****
			pPtr->optSpace.estimation();
			pPtr->sFitness.copy( &(pPtr->optSpace.fitness) );
			pPtr->nns.copy( &(pPtr->optSpace.nnAfr) );      //- New NN

			//-- Update personal best for both objectives
			pPtr->updateP();
		}
		else { pPtr->sFitness.reinit(); }
	}
}
/*============================================================================*
 *  Function name			:	multiPso::archiving
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Update the archive
 *============================================================================*/
void multiPso::archiving()
{
//	printf("      archiving\n");
	for(int n=0; n<size; n++)
	{
//		printf("      %d - s(%1.2f, %1.2f, %1.2f, %1.2f) f(%1.2f, %1.2f <- %d) - %d",
//		         n, p[n].s[0], p[n].s[1], p[n].s[2], p[n].s[3],
//				 1-p[n].sFitness.clsRate,  p[n].sFitness.sizeNn,
//			     p[n].nns.sizeF2, p[n].nns.nPatternsLearned);

		//-- Check if the network exist
		//     -> if not, the particle was outside boundaries
		int flagArc = 1;
		for(int d=0; d<nDimensions; d++)
		{   if( p[n].s[d] > 1  ||  p[n].s[d] < 0 ) { flagArc = 0; }   }
	
		if( flagArc )
		{//	printf(" - check");
			int reinitialized = myArchive.checkDominance(&p[n]);
		  
			//-- Reinitialize with the archive occupancy
			if(reinitialized)   { p[n].randomization(0); }
//			printf(" - done\n");
		}
//		else
//		{   printf("\n");   }
	}
}
/*============================================================================*
 *  Function name			:	multiPso::updateArchiveFit
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Update the archive
 *============================================================================*/
void multiPso::updateArchiveFitness()
{
	//-- Db test for the archive update (we get theone from the particles)
	dbase *dbTest = &p[0].optSpace.dbFit;
	
	#pragma omp parallel for default(shared) schedule(dynamic)
	for(int m=0; m<myArchive.size; m++)
	{
		//-- Define pointers for parallelization
		particle *pPtr   = & myArchive.p[m];

		//-- If the solution exist
		if(myArchive.filled[m])
		{
			artmap *nnTest     = & pPtr->nns;
			result *resultTest = & pPtr->optSpace.fitness;

			pPtr->optSpace.test( nnTest, dbTest, resultTest );	//- Update 
			pPtr->sFitness.copy( resultTest );					//- Copy
		}
	}

	//-- Update the archive accordingly
	myArchive.update();
}
/*============================================================================*
 *  Function name			:	multiPso::eucledean2
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Eucledean distance ^2 between two particles
 *  							and two particle's personal best positions
 *----------------------------------------------------------------------------*
 *  Arguments				:   particles pi & pj - particles
 *  Returns					:	float     dist    - (euclidean distance)^2
 *============================================================================*/
//o-------------------------- Two particle postions -------------------------o//
float multiPso::euclidean2(particle *pi, particle *pj)
{
	float dist = 0;

	for(int d=0; d<nDimensions; d++) { dist += pow(pi->s[d]-pj->s[d],2); }
	return dist;
}
//o----------------- Two particle's personal best positions -----------------o//
float multiPso::euclidean2p(int i, int j, int obj)
{
	float dist = 0;

	for(int d=0; d<nDimensions; d++)
	{	dist += pow(p[i].p[obj][d] - p[j].p[obj][d], 2);   }
	return dist;
}
/*============================================================================*
 *  Function name			:	optimization::findgBest
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Find the global best
 *============================================================================*/
void multiPso::findgbest()
{
	sort.size = size;

	//-- Global best for each objective
	for(int o=0; o<nObjectives; o++)
	{	//-- We want a descending class. rate and descending compression
		for(int n=0; n<size; n++)  { fitness[o][n] = 0; }
		
		for(int n=0; n<size; n++)
		{   fitness[0][n] = 0 - p[n].pFitness[o].clsRate;
			fitness[1][n] = 0 - p[n].pFitness[o].normCpn;
		}
		
		int first  =  o,				//- 0 -> 0, and 1 -> 1
			second = (int)pow(o-1,2);	//- 0 -> 1, and 1 -> 0

		for(int n=0; n<size; n++)
		{	sort.order     [n] =                       n ;
			sort.fitness   [n] = fitness[first ][n];
			sort.tieBreaker[n] = fitness[second][n];
		}

		//-- Sorts and the global best is the first element
		sort.heapSort();
		gbest[o] = sort.order[0];
	}
}
/*============================================================================*
 *  Function name			:	optimization::findLbest
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Find the local best of current particles with
 *  							two consecutive sorts
 *----------------------------------------------------------------------------*
 *  Arguments				:   int current   - Current particle
 *============================================================================*/
void multiPso::findlBest(int current, int type)
{
	//o----------------------- Finds the neighborhood -----------------------o//
	//-- Sort by distances
	sort.size = size;
	for(int n=0; n<size; n++)
	{	//-- We want an ascending sort, so we leave things as is
		sort.fitness   [n] = distance[current][n];
		sort.tieBreaker[n] = 0;
	}
	sort.heapSort();

	//o------------------- Define the parameter to sort by ------------------o//
	sort.size = sizeHood;
	//-- Sorts the neighborhood by performance (for two objectives) ...
	int first  = type,			        //- 0 -> 0, and 1 -> 1
		second = (int)pow(type-1,2);	//- 0 -> 1, and 1 -> 0
	for(int n=0; n<sizeHood; n++)
	{	//-- We want a descending class. rate and ascending sizes
		hood           [n] =                  sort.order[n]  ;
		sort.fitness   [n] = fitness[first ][sort.order[n]];
		sort.tieBreaker[n] = fitness[second][sort.order[n]];
	}

	//o----------- Finds the local best (or global best for pso) ------------o//
	sort.heapSort();
	p[current].lbest[type] = hood[sort.order[0]];

	//-- If not a master -> check for mutual local best references
	//   for the ****possible_infinit_loop****
	if(p[current].lbest[type] != current)
	{	int m = p[current].lbest[type];
		while( (m != -1) && (p[m].lbest[type] != m))
		{	//-- In case of infinit loop, make the particle a master.
			m = p[m].lbest[type];
			if(m == current)   { p[current].lbest[type] = current;  break; }
		}
	}
}
/*============================================================================*
 *  Function name			:	optimization::findEnsembleMembers
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Creates the ensemble - Best of each memetic.
 *============================================================================*/
void multiPso::setensembleMemetic()
{
	//-- Reinitialize the sorting order and the ensemble
	myEnsemble.reinitialization();
	//-- Add members coresponding the the best ones of each memetic
	for(int m=0; m<myArchive.nMemetics; m++)
	{
		if(myArchive.nMembers[m])
		{
			int addr = myArchive.sizeMemetics * m + myArchive.nMembers[m]-1;

			myEnsemble.members [myEnsemble.size] =              addr;
			myEnsemble.pFitness[myEnsemble.size] = &myArchive.p[addr].sFitness;
			myEnsemble.nn      [myEnsemble.size] = &myArchive.p[addr].nns     ;
			myEnsemble.size++;
		}
	}
}
/*============================================================================*
 *  Function name			:	optimization::findEnsembleMembers
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Creates the ensemble - Best archive members of 
 *  							each memetic + diversity greedy search.
 *============================================================================*/
void multiPso::setensembleGreedyArchive()
{
	int sizeArc = myArchive.nFilled;

	//-- Reinitialize the sorting order and the ensemble
	myEnsemble.reinitialization();

	//o--------------------- Compute initial diversity ----------------------o//
	//-- Initialize the indicators
	for(int e=0; e<sizeArc; e++)   { indicator[e] = 0; }

	//-- Set the initial ensemble made with the archive
	setensembleMemetic();

	//--  Initial ensemble diversity (average distances between p best),
	//    or zero if the ensemble contains only one network.
	float diversity = 0;
	if(myEnsemble.size > 1)
	{	for(int e1=0; e1<myEnsemble.size-1; e1++)
	  	{	for(int e2=e1+1; e2<myEnsemble.size; e2++)
		  	{	int i = myEnsemble.members[e1], j = myEnsemble.members[e2];

				float div = 0;	
				for(int d=0; d<nDimensions; d++)
				{   div += pow(myArchive.p[i].s[d] - myArchive.p[j].s[d], 2); }
				diversity += sqrt(div);
		  }
	  }
	  diversity = diversity * 2 / ( myEnsemble.size * (myEnsemble.size-1) );
	}
	float bestDiversity = diversity;

	//o------------- te = Time for the greedy algo. (ulas 2009) -------------o//
	for(int te=myEnsemble.size; te<sizeArc; te++)
	{	int better = 0, newMember = 0, sizeNew = te+1;

		//-- Check for all particles in the archive (i.e. networks)
		for(int e=0; e<sizeArc; e++)
		{
			//-- Don't check if network is already present
			//     -> diversity won't increase if so! 
			if(myArchive.filled[e])
			{	//-- Assign the member and network
				myEnsemble.members [te] = e;

				//o------------------- Compute diversity --------------------o//
				diversity = 0;
				for(int e1=0; e1<sizeNew-1; e1++)
				{	for(int e2=e1+1; e2<sizeNew; e2++)
					{	int i = myEnsemble.members[e1],
							j = myEnsemble.members[e2];
						float div = 0;	
						for(int d=0; d<nDimensions; d++)
						{   div += pow(myArchive.p[i].s[d] -
						               myArchive.p[j].s[d], 2); }
						diversity += sqrt(div);
					}
				}
				diversity = diversity * 2 / ( sizeNew * (sizeNew-1) );

				//-- Does the new member bring something new
				if(diversity > bestDiversity)
				{	//-- If so, assign stuff
					better        = 1;
					newMember	  = e;
					bestDiversity = diversity;
				}//-- For validation
			}
		}

		//-- If there is no change, we erase the last tested member and get out.
		if(better)
		{	myEnsemble.members [te] =                newMember            ;
			myEnsemble.nn	   [te] = &( myArchive.p[newMember].nns      );
			myEnsemble.pFitness[te] = &( myArchive.p[newMember].sFitness );
			myEnsemble.size         = sizeNew;
		}
		else
		{   myEnsemble.members [myEnsemble.size]  = -1;
			myEnsemble.nn      [myEnsemble.size]  = (artmap*)NULL;
			myEnsemble.pFitness[myEnsemble.size]  = (result*)NULL;
			break;
		}
	}
}
/*============================================================================*
 *  Function name			:	optimization::findEnsembleMembers
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Creates the ensemble - All archive members.
 *============================================================================*/
void multiPso::setensembleArchive()
{
	//-- Reinitialize the sorting order and the ensemble
	myEnsemble.reinitialization();
	
	//-- Add members coresponding the the best ones of each memetic
	for(int m=0; m<myArchive.size; m++)
	{
		if(myArchive.filled[m])
		{	myEnsemble.members [myEnsemble.size] =               m;
			myEnsemble.pFitness[myEnsemble.size] = & myArchive.p[m].sFitness;
			myEnsemble.nn      [myEnsemble.size] = & myArchive.p[m].nns     ;
			myEnsemble.size++;
		}
	}
}
/*============================================================================*
 *  Function name			:	optimization::findEnsembleMembers
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Creates the ensemble - lbest only
 *============================================================================*/
void multiPso::setensembleLocalBests()
{
	//-- Reinitialize the sorting order and the ensemble
	findgbest();
	myEnsemble.reinitialization();

	//-- Add local best - performances
	myEnsemble.size = nSS[0];
	for(int e=0; e<nSS[0]; e++)
	{	//-- Assign the member ... for objective 0
		myEnsemble.members [e] =       lbests[0][e];
		myEnsemble.pFitness[e] = &( p[ lbests[0][e] ].pFitness[0] );
		myEnsemble.nn      [e] = &( p[ lbests[0][e] ].nnp     [0] );
	}
}
/*============================================================================*
 *  Function name			:	optimization::findEnsembleMembers
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Creates the ensemble - the good way in the
 *  							search space.
 *============================================================================*/
void multiPso::setensembleGreedy()
{
	//-- Reinitialize the sorting order and the ensemble
	findgbest();
	myEnsemble.reinitialization();

	//o--------------------- Compute initial diversity ----------------------o//
	//-- Initialize the indicators
	for(int e=0; e<size; e++)   { indicator[e] = 0; }

	//-- Set the initial ensemble made of the local best (both objectives)
	setensembleLocalBests();
	
	//--  Initial ensemble diversity (average distances between p best),
	//    or zero if the ensemble contains only one network.
	float diversity = 0;
	if(myEnsemble.size > 1)
	{
	  for(int e1=0; e1<myEnsemble.size-1; e1++)
	  {	for(int e2=e1+1; e2<myEnsemble.size; e2++)
		  {	int i = myEnsemble.members[e1], j = myEnsemble.members[e2];
			  for(int d=0; d<nDimensions; d++)
			  {   diversity += pow(p[i].p[0][d] - p[j].p[0][d], 2);   }
		  }
	  }
	  diversity = diversity * 2 / ( myEnsemble.size * (myEnsemble.size-1) );
	}
	float bestDiversity = diversity;

	//o------------- te = Time for the greedy algo. (ulas 2009) -------------o//
	for(int te=myEnsemble.size; te<size; te++)
	{	int better = 0, newMember = 0, sizeNew = te+1;

		//-- Check for all particles (i.e. networks)
		for(int e=0; e<size; e++)
		{
			int flagCheck = 1;
			for(int check=0; check<myEnsemble.size; check++)
			{   if(e == myEnsemble.members[check])   { flagCheck = 0; }   }

			if(flagCheck)
			{	//-- Assign the member and network
				//   (1) We hack the order vector and (2) 1st = gbest
				myEnsemble.members [te] = e;

				//o------------------- Compute diversity --------------------o//
				diversity = 0;
				for(int e1=0; e1<sizeNew-1; e1++)
				{	for(int e2=e1+1; e2<sizeNew; e2++)
					{	int i = myEnsemble.members[e1],
							j = myEnsemble.members[e2];
						float div = 0;
						for(int d=0; d<nDimensions; d++)
						{   div += pow(p[i].p[d] - p[j].p[d], 2);   }
						diversity += sqrt(div);
					}
				}
				diversity = diversity * 2 / ( sizeNew * (sizeNew-1) );

				//-- Does the new member bring something new
				if(diversity > bestDiversity)
				{	//-- If so, assign stuff
					better        = 1;
					newMember	  = e;
					bestDiversity = diversity;
				}//-- For validation
			}
		}

		//-- If there is no change, we erase the last tested member and get out.
		if(better)
		{	myEnsemble.members [te] =      newMember            ;
			myEnsemble.nn	   [te] = &( p[newMember].nnp     [0] );
			myEnsemble.pFitness[te] = &( p[newMember].pFitness[0] );
			myEnsemble.size         = sizeNew;
		}
		else
		{   myEnsemble.members [myEnsemble.size]  = -1;
			myEnsemble.nn      [myEnsemble.size]  = (artmap*)NULL;
			myEnsemble.pFitness[myEnsemble.size]  = (result*)NULL;
			break;
		}
	}
}
/*============================================================================*
 *  Function name			:	optimization::setensembleTotal
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Creates the ensemble with all the swarm
 *                              (aimed at performance)
 *============================================================================*/
void multiPso::setensembleSwarm()
{
	//-- Reinitialize the ensemble
	myEnsemble.reinitialization();

	myEnsemble.size = size;
	for(int e=0; e<size; e++)
	{	myEnsemble.members [e] =      e            ;
		myEnsemble.pFitness[e] = &( p[e].pFitness[0] );
		myEnsemble.nn      [e] = &( p[e].nnp     [0] );
	}
}
/*============================================================================*
 *  Function name			:	optimization::setensembleGbest
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Creates an ensemble with only the gbest
 *                              network (aimed at performance).
 *============================================================================*/
void multiPso::setensembleGbest()
{
	//-- Reinitialize the ensemble
	myEnsemble.reinitialization();

	myEnsemble.size        = 1;
	myEnsemble.members [0] =      gbest[0]               ;
	myEnsemble.pFitness[0] = &( p[gbest[0]].pFitness[0] );
	myEnsemble.nn      [0] = &( p[gbest[0]].nnp     [0] );
}
/*============================================================================*
 *  Function name			:	optimization::setobjFtn
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	s and s the gbest neural network and
 *  							copies it all other particles for the next
 *  							optimization process
 *----------------------------------------------------------------------------*
 *  Arguments & returns		:   depends
 *============================================================================*/
void multiPso::setobjFtn()
{
	//--  the initial conditions for the next opt. process (i.e. block)
	for(int n=0; n<size; n++)
	{   p[n].optSpace.nnBfe.copy( &(p[n].nns) );   }
}
/*============================================================================*
 *  Function name			:	optimization::savePopulation
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Saves the population
 *----------------------------------------------------------------------------*
 *  Arguments				:   char *nameFile
 *  							int   currentRep  - Current replication
 *============================================================================*/
void multiPso::save(char *nameFile, int currentRep)
{
	//-- File opening
	ofstream f(nameFile, ios::app | ios::binary);
	if(!f)   { printf("Error! Can't open the file %s \n", nameFile);  exit(1); }

	//-- General informations
	f.write((char*)&currentRep,  sizeof(int));
	f.write((char*)&size,   sizeof(int));
	f.write((char*)&nDimensions, sizeof(int));
	f.write((char*)&nObjectives, sizeof(int));
	f.write((char*)&nIterations, sizeof(int));

	for(int o=0; o<nObjectives; o++)
	{   f.write((char*)&(gbest[o]), sizeof(int));   }

	//-- For each particle
	for(int i=0; i<size; i++)
	{
		//-- Position, previous position (for velocity), & best position
		f.write((char*) p[i].s    , sizeof(float) * (unsigned)nDimensions);
		f.write((char*) p[i].sPrev, sizeof(float) * (unsigned)nDimensions);
		
		for(int o=0; o<nObjectives; o++)
		{ f.write((char*) p[i].p[o], sizeof(float) * (unsigned)nDimensions); }

		//-- Local best (to identify the subswarm))
		f.write((char*) p[i].lbest, sizeof(int) * (unsigned)nObjectives);

		//-- Performances (current and best position)
		f.write((char*)&(p[i].sFitness.clsRate), sizeof(float));
		f.write((char*)&(p[i].sFitness.sizeNn ), sizeof(float));
		f.write((char*)&(p[i].sFitness.normCpn), sizeof(float));

		for(int o=0; o<nObjectives; o++)
		{	f.write((char*)&(p[i].pFitness[o].clsRate  ), sizeof(float));
			f.write((char*)&(p[i].pFitness[o].sizeNn   ), sizeof(float));
			f.write((char*)&(p[i].pFitness[o].normCpn  ), sizeof(float));   }
	}

	//-- File closing
	f.close();
}
