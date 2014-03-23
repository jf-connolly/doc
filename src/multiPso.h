/******************************************************************************
 *   Copyright(C) 2008 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca *
 *   based on the work done by Marcelo Kapp									  *
 *																			  *
 *   pso.h - Canonical particle swarm optimizer as defined in :				  *
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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU         *
 *   General Public License for more details.                                 *
 *                                                                            *
 *   You should have received a copy of the GNU General Public License along  *
 *   with this program; if not, write to the Free Software Foundation, Inc.,  *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                *
 ******************************************************************************/
/*------------------------------- Header files -------------------------------*/
//-- Standard
//-- (none)

//-- Project specific
#include "archive.h"
#include "configuration.h"
#include "ensemble.h"
#include "particle.h"
#include "sorting.h"

#ifndef PSO_H
#define PSO_H

using namespace std;

#define argMax(val_a,val_b,a,b)   ( ((val_a>val_b)||(val_a==val_b))?a:b )
#define argMin(val_a,val_b,a,b)   ( ((val_a<val_b)||(val_a==val_b))?a:b )
#define performance               0
#define compression               1

/*--------------------------------- Class PSO --------------------------------*/
class multiPso
{
public:
	/*------------------------------ Variables -------------------------------*/
	char name[10];			//- Algorithm name
	int time;           	//- Time (in the sense of learning blocks)
	int size,				//- Total number of particles in the swarm.
		nDimensions, 		//- Number of Dimensions of the Swarm
		nObjectives, 		//- Number of Dimensions of the Swarm

	   *gbest,				//- Global best (i.e. THE solution!)
		nIterations;		//- Global best (i.e. THE solution!)

	particle *p;			//- Particles

	artmap nnOptimal;		//- Best neural network after an iteration

	//-- For local best topology
	int     sizeHood,		//- Size of the neighborhhood
	      * hood;			//- Neighborhood of a particle
	float   distanceMax,	//- Maximal distance (for initialization)
	      **distance; 		//- Distance between particles

	//-- For DNPSO
	int   sizeSSmax,		//- Maximal size of subswarms
	      nSSmax;			//- Maximal number of subswarms
	float vFreeMin,			//- Minimal velocities of free particles
	      dSSmin;			//- Minimal distance between two subswarm

	float **fitness;		//- Fitness values for the objectives
	int **lbests,			//- d'hu
	     *nSS,				//- Number of subswarms
	    **freep,	  		//- Free particles
	     *nfree,			//- Number of free particles
	    **members,			//- Subswarms members
	    **nMembers;			//- Number of members

	//-- Ensemble
	ensemble myEnsemble;    //- The ensemble itself
	float   *indicator;		//- Diversity indicator (values)

	//-- Sorting object
	sorting sort;

	//-- Sorting object
	archive myArchive;

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors and destructor
	multiPso();
	void initialization(configuration *cfg);
	~multiPso();

	//-- PSO algorithm
	void defineSubSwarms();
	void movePopulation();
	void updatePopulation();
	void reevaluateSwarm();
    void archiving();
    void updateArchiveFitness();

	//-- Save
	void save(char *nameFile, int currentRep);

	//-- Ensemble selection
	void  setensembleMemetic      ();
	void  setensembleGreedyArchive();
	void  setensembleArchive      ();
	void  setensembleLocalBests   ();
	void  setensembleGreedy       ();
	void  setensembleSwarm        ();
	void  setensembleGbest        ();

	//-- Utilities
	float euclidean2 (particle *pi, particle *pj);
	float euclidean2p(int i, int j, int obj);
	void  setobjFtn();
	void  findgbest();
	void  findlBest(int current, int type);
};
#endif
