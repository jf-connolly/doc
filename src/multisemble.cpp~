/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  multisemble.cpp - Incremental learning with ensemble of FAM               *
 *                                                                            *
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

//o------------------------------ Header files ------------------------------o//
//-- Standard
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>


//-- Project specific
#include "multiPso.h"			//- Optimization method
#include "configuration.h"		//- Configuration

//o--------------------------------------------------------------------------o//
//o---------------------------------- Main ----------------------------------o//

using namespace std;

int main(int argc, char *argv[])
{
	char nameArtmap[150], namePso[150], nameArc[150], nameResult[150],
	     nameOp[10], nameLn[5], typeHp[20];

	//-- Configuration parameters
	if (argc < 2)
	{	sprintf(nameArtmap, "config.txt");
	    printf("/-- Usage: %s <%s>\n", argv[0], nameArtmap);
	}
	else
	{	sprintf(nameArtmap, "%s", argv[1]);
		printf("/-- Usage: %s %s\n", argv[0], nameArtmap);
}

//o----------------------------- Initialization -----------------------------o//
    configuration cfg(nameArtmap);
    multiPso myPso;
    dbase    dbTest;

    //-- Fixes stuff for incremental and batch learning with PSO
    int reset;
    if( cfg.incremental ) { reset = 0;
							sprintf(nameLn, "%s", "Inc");
							sprintf(nameOp, "%s", myPso.name); }
	else				  { reset = 1;
							sprintf(nameLn, "%s", "Bth");      }

    int nTest = 4;

	int nItrs    = 0,
		bestItr  = 0,
		nItrsMax = cfg.nIterationsMax,
		nItrsOvr = cfg.nIterationsOvr;

	int nBlocks       = cfg.nBlocks,
		startingRep   = cfg.startingRep,
	    nReplications = cfg.nReplications,
	    size     = cfg.sizeSwarm,
	    nObjectives   = cfg.nObjectives;

    float *bestFit = new float[nObjectives];
    int   *nStar   = new int  [nObjectives];

	//-- Times
	time_t tStart, tEnd;//, tInitialization, tEnd;

	//o------------------------------- Start --------------------------------o//
	if(cfg.error) { printf("Never started (wrong cfg file)\n");  exit(1); }

	printf("\n/---------------------- Start ----------------------/\n");
	//-- Name of the pso


	//-- Initializing time & random number generator
	tStart = time(NULL);
//	srand( 0 );
 	srand( (unsigned)tStart );

	printf("/-- Configuration file:                      %s\n"
		   "/-- Data base and learning mode:             %s & %s\n"
		   "/-- Nb. blocks, nb. rep., nb estimation rep: %d, %d, %d\n"
		   "/-- Nb. de particules:                       %d\n/\n",
		   nameArtmap, cfg.nameDb, nameLn, cfg.nBlocks,
		   cfg.nReplications, cfg.nEstimationRep, cfg.sizeSwarm);

	//o------------------------ Delete exiting files ------------------------o//
	for(int typeEs=0; typeEs<nTest; typeEs++)
	{
		switch(typeEs)
		{
			case 0 : sprintf(typeHp, "%s", "EoMA+d");		break;
			case 1 : sprintf(typeHp, "%s", "EoMA");			break;
			case 2 : sprintf(typeHp, "%s", "EoMA+t");		break;
			case 3 : sprintf(typeHp, "%s", "Hdnc");			break;
			case 4 : sprintf(typeHp, "%s", "EoLB_p");		break;
			case 5 : sprintf(typeHp, "%s", "EoLB+d");		break;
			case 6 : sprintf(typeHp, "%s", "EoFAM+t");		break;
		}

		//-- Results for the ensemble
		sprintf(nameResult, "%s/%s_%s%s%s%s_%dBlocks_wd%d.result", cfg.pathSs,
				cfg.nameDb, cfg.nameSc, nameLn, typeHp, nameOp, cfg.nBlocks,
				cfg.widthMemetics);
		remove(nameResult);

		//-- Ensemble of networks
		sprintf(nameArtmap, "%s/%s_%s%s%s%s_%dBlocks_wd%d.artmap", cfg.pathSs,
				cfg.nameDb, cfg.nameSc, nameLn, typeHp, nameOp, cfg.nBlocks,
				cfg.widthMemetics);
		remove(nameArtmap);
	}

	//-- Populations of solutions
		sprintf(namePso, "%s/%s_%s%s%s%s_%dBlocks_wd%d.pso", cfg.pathSs,
				cfg.nameDb, cfg.nameSc, nameLn, typeHp, nameOp, cfg.nBlocks,
				cfg.widthMemetics);
	remove(namePso);

	//-- Archives
	sprintf(nameArc, "%s/%s_%s%s%s%s_%dBlocks_wd%d.arc", cfg.pathSs,
			cfg.nameDb, cfg.nameSc, nameLn, typeHp, nameOp, cfg.nBlocks,
			cfg.widthMemetics);
	remove(nameArc);


	//-- Load the test data base
	//-- Populations of solutions
// 	printf("dbTest\n");
	sprintf(nameResult, "%s/%s/test.db", cfg.pathDb, cfg.nameDb);
	dbTest.initialization(nameResult);

	//-- Initialization of the whole shebang (pso, archive, and ensemble)
// 	printf("pos init\n");
	myPso.initialization(&cfg);

	//o------------------------ Optimization process ------------------------o//
	//-- For all blocs and replications
	for(int r = startingRep-1; r < startingRep-1 +nReplications; r++)
	{	printf("/\n/-- Replication %d\n", r+1);

		//-- Initialize class activity for the ensemble test process
// 		printf("active classes\n");
		for(int k=0; k<dbTest.nClasses; k++)  {	dbTest.activeClasses[k] = 0; }
	
		//-- Objective space initialization - Swarm and archive
// 		printf("obj space init\n");
		for(int n=0; n<size; n++)
		{   myPso.p[n].optSpace.initEstimation(&cfg, r);   }

		//o-------------------- For each learning blocks --------------------o//
		//-- Incremental -> start from one, batch -> start from nBlocks
		int tStart = 0;
		if(!cfg.incremental)   { tStart = nBlocks-1; }
		for(int t=tStart; t<nBlocks; t++)
		{	
			//-- Set the time
			printf("/      Bloc %d", t+1);
			myPso.time = t;

			//-- Load data bases - Swarm and archive
			for(int n=0; n<size; n++)  { myPso.p[n].optSpace.loadDb(&cfg,t,r); }

			//-- Redefine the objective function with the pbest neural networks
			//   and reevaluated the fitness of the pbest position
			myPso.setobjFtn();
			myPso.reevaluateSwarm();
			myPso.updateArchiveFitness();
			int maxDomaine = myPso.p[0].optSpace.nAccPatterns;
			myPso.myArchive.defineBoundaries(maxDomaine);

			//-- Define the subswarms for a first time
			myPso.defineSubSwarms();

			//-- Until maximal number of iterations are reached or gbest
			//   did not got better for a maximal "over" number of iterations
			nItrs = 0;   bestItr = 0;
			while( nItrs < nItrsMax  &&  (nItrs-bestItr) < nItrsOvr )
			{
				//-- Optimization - move, update, subswarms, & adjust pbest
				myPso.movePopulation();
 				myPso.updatePopulation();
				myPso.defineSubSwarms();
				myPso.archiving();
				nItrs++;
				myPso.nIterations = nItrs;
				
				//-- Find and assign the global best
				myPso.findgbest();
				nStar[0] = myPso.gbest[0];
				nStar[1] = myPso.gbest[1];

				//-- Check if it is the best iteration. If the case, update
				//   the best fitness and the best iteration
				if( myPso.p[ nStar[0] ].pFitness[0].clsRate > bestFit[0] )
				{   bestFit[0] = myPso.p[ nStar[0] ].pFitness[0].clsRate;   }
				if(	myPso.p[ nStar[1] ].pFitness[1].normCpn > bestFit[1] )
				{	bestFit[1] = myPso.p[ nStar[1] ].pFitness[1].normCpn;   }

				if( myPso.p[ nStar[0] ].pFitness[0].clsRate > bestFit[0] ||
					myPso.p[ nStar[1] ].pFitness[1].normCpn > bestFit[1]    )
				{	bestItr = nItrs;   }

				//-- Save the swarm and archive
				myPso.save(namePso, r);
				myPso.myArchive.save(nameArc, r);
			}
			printf(" (nb. iterations: %d, archive size: %d)\n", nItrs,
			       myPso.myArchive.nFilled);

			//-- Update class activity with the class(es) in the current block
			for(int k=0; k<dbTest.nClasses; k++)
			{   if(myPso.p[0].optSpace.dbTrn.sizeClasses[k])
				{ dbTest.activeClasses[k] = 1; }
			}

			//o----------------- Test for all ensemble type -----------------o//
			for(int typeEs=0; typeEs<nTest; typeEs++)
			{
				//-- Type of ensemble
				switch(typeEs)
				{
					case 0 : sprintf(typeHp, "%s", "EoMA+d");
							 myPso.setensembleGreedyArchive();		break;
					case 1 : sprintf(typeHp, "%s", "EoMA");
							 myPso.setensembleMemetic();			break;
					case 2 : sprintf(typeHp, "%s", "EoMA+t");
							 myPso.setensembleArchive();			break;
					case 3 : sprintf(typeHp, "%s", "Hdnc");
							 myPso.setensembleGbest();				break;
					case 4 : sprintf(typeHp, "%s", "EoLB_p");
							 myPso.setensembleLocalBests(); 		break;
					case 5 : sprintf(typeHp, "%s", "EoLB+d");
							 myPso.setensembleGreedy();				break;
					case 6 : sprintf(typeHp, "%s", "EoFAM+t");
							 myPso.setensembleSwarm();				break;
				}

				//-- Names for the result and ensemble
				sprintf(nameResult, "%s/%s_%s%s%s%s_%dBlocks_wd%d.result",
				cfg.pathSs, cfg.nameDb, cfg.nameSc, nameLn, typeHp, nameOp,
				cfg.nBlocks, cfg.widthMemetics);

				sprintf(nameArtmap, "%s/%s_%s%s%s%s_%dBlocks_wd%d.artmap",
				cfg.pathSs,	cfg.nameDb, cfg.nameSc, nameLn, typeHp, nameOp,
				cfg.nBlocks, cfg.widthMemetics);
				
				//-- Set the optimal neural network, set it to all other
				//   particles for the next iteration and save it
//				myPso.myEnsemble.save(nameArtmap, r, t);

				//-- Test & save the performance of the gbest network
				//   and the ensemble
				myPso.myEnsemble.assignDbTest(&dbTest);
				myPso.myEnsemble.test();
				
				myPso.myEnsemble.performances.nPatternsLearned = 
										    myPso.p[0].optSpace.nAccPatterns;
				myPso.myEnsemble.performances.save(nameResult, r);

				//o----------------- Print results on screen ----------------o//
				printf("/         %s: [ ", typeHp);

				for(int e=0; e<myPso.myEnsemble.size; e++)
				{   printf("%d (%d) ", myPso.myEnsemble.members[e],
				                       myPso.myEnsemble.nn[e]->sizeF2);   }

				printf("],\n/              ");
				printf("class. rate: %1.2f, sizeEns: %d, sizeNn: %1.0f\n",
					   myPso.myEnsemble.performances.clsRate,
					   myPso.myEnsemble.size,
					   myPso.myEnsemble.performances.sizeNn);
				//o-------------- end - Print results on screen -------------o//
				
			}//** end : for all types of ensemble

			//-- Delete the data bases (calls the destructor MANUALLY),
			//   except for the last block of the last replication
			if(t != nBlocks-1 || r != startingRep-1 + nReplications-1)
			{
				for(int n=0; n<size; n++)
				{	myPso.p[n].optSpace.deleteDb();   }
			}

			//-- If we are doing batch learning, reinitialize neural networks,
			//   and the particle positions
			printf("/      Reset (maybe)\n");
			if( reset )
			{	for(int n=0; n<size; n++)
				{   myPso.p[n].optSpace.nnBfe.reinit();
					myPso.p[n].optSpace.nnAfr.reinit();
					myPso.p[n].randomization(0);
				}
			}
			printf("/      End\n");
		}//** end : for t Blocks

		//-- After each replication: we reinitialize everything
		printf("/   Reinitialize the swarm");
		for(int n=0; n<size; n++)
		{   myPso.p[n].optSpace.nnBfe.reinit();
			myPso.p[n].optSpace.nnAfr.reinit();
			myPso.p[n].randomization(0);
		}
	
		//-- Archive reinitialization
		printf(" & the archive");
		myPso.myArchive.reinit();
		printf(" - done\n");
    }

	//o----------------------- We delete every thing ------------------------o//
	tEnd = time(NULL);
	printf("\n");
	printf("Start time   : %s", ctime(&tStart));
	printf("End time     : %s", ctime(&tEnd));
	double diff = difftime(tEnd, tStart);
	int hours   = (int)(  diff/(60*60) );
	int minutes = (int)( (diff-hours*60*60)/60 );
    int seconds = (int)( (diff-hours*60*60-minutes*60) );
	printf("Process time : %*d:%*d:%*d\n", 2, hours, 2, minutes, 2, seconds );
	printf("/----------------------- End ------------------------/\n");
	return 0;
}
//o---------------------------------- Fin -----------------------------------o//
//o--------------------------------------------------------------------------o//
