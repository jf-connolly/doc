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
//-- (none)

//-- Project specific
//-- (none)

#ifndef configuration_H
#define configuration_H

using namespace std;

class configuration {
public:
	/*------------------------------ Variables -------------------------------*/
	//-- General
	int  error;

	char nameDb[10],		//-- Name of the data base used
		 nameSc[4];			//-- Name of the learning scenario
	int  incremental,		//-- Incremental learning (vs. batch)?
		 nBlocks;			//-- Number of learning blocks

	int	startingRep,		//-- Starting replication
		nReplications,		//-- Number of replications
		nEstimationRep;		//-- Number of internal replications

	char pathDb[200],		//-- Path to the data bases
		 pathSs[200];		//-- Path to the saved stuff


	//-- Optimization method (pso)
	int nDimensions, 		//-- Number of Dimensions of the Swarm
		nObjectives,
		nIterationsMax,
		nIterationsOvr;

	float alpha, beta;		//-- Inertia

	int   sizeSwarm,		//-- Total number of particles in the swarm.
		  sizeHood,       	//-- Size of the neighborhood
		  sizeSSmax,    	//-- Maximal size of the subswarms
		  nSSmax;			//-- Maximal number of subswarms
	float dSSmin,		//-- Minimal distance between two subswarm
		  vFreeMin;			//-- Minimal velocities of free particles


	float *sMin, *sMax;

	//-- ARTMAP
	int sizeF2max,		//-- Max. nb. of allocated neurons in memory
		nEpochMax;		//-- Maximal nb. of training epochs

	//-- Archive
	int nMemetics,		//- Number of memetic regions
	    sizeMemetics,	//- Size of each regions
	    widthMemetics;	//- Size of each regions

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors & destructor
	configuration(char * nameFile);
	~configuration();
};
#endif
