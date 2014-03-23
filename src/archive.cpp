/******************************************************************************
 *   Copyright(C) 2008 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca *
 *																			  *
 *   archive.cpp - Archive for multiobjective optimization					  *
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
#include "archive.h"

/*============================================================================*
 *  Function name			:	archive::archive
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Constructor without arguments (default)
 *============================================================================*/
archive::archive()
{
	sizeMemetics  = 0;
	widthMemetics = 0;
	nDimensions   = 0;
	nMemetics     = 0;
	size          = 0;
	upperBound    = 0;

	nFilled       = 0;
  
	filled     = (int*)     NULL;
	filledOld  = (int*)     NULL;
	boundaries = (int*)     NULL;
	nMembers   = (int*)     NULL;
	members    = (int*)     NULL;
	p          = (particle*)NULL;
	o          = (particle*)NULL;
}
 /*===========================================================================*
 *  Function name			:	archive::initialization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Constructor without arguments (default)
 *----------------------------------------------------------------------------*
 *  Argument				:   Configuration *configPso
 *============================================================================*/
void archive::initialization(configuration *cfg)
{
	dbug = 0;
	
	//-- Parameters
	sizeMemetics  = cfg->sizeMemetics;
	widthMemetics = cfg->widthMemetics;
	nDimensions   = cfg->nDimensions;
	upperBound    = cfg->sizeF2max;
	
	//-- Other parameters
	if(cfg->sizeF2max)   { size = sizeMemetics*(upperBound/widthMemetics +1); }
	
	else				 { size = sizeMemetics;                                }
	boundaries = new int[size];
	
	//-- Which particles are in the archive
	nFilled   = 0;
	filled    = new int[size];
	filledOld = new int[size];
	for(int m=0; m<size; m++) { filled[m] = 0;   filledOld[m] = 0; }

	//-- Particles - current (& objective function) and old
	p = new particle[size];
	o = new particle[size];
	for(int n=0; n<size; n++)   { p[n].initialization(cfg);
								  o[n].initialization(cfg); }

	//-- Tested existing members
	members  = new int[size];
	nMembers = new int[size];
	for(int m=0; m<size; m++)   { nMembers[m] = 0; }
}
/*============================================================================*
 *  Function name			:	archive::~archive
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	PSO destructor
 *============================================================================*/
archive::~archive()
{
	delete [] filled;
	delete [] filledOld;
	delete [] nMembers;
	delete [] o;
	delete [] p;
	delete [] boundaries;

	sizeMemetics = 0;
	nMemetics    = 0;
	size         = 0;
	nDimensions  = 0;
	upperBound   = 0;
}
/*============================================================================*
 *  Function name			:	archive::reinit
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Reinitialization
 *============================================================================*/
void archive::reinit()
{	for(int m=size-1; m>-1; m--)   { if(filled[m])  { remove(m); } }   }
/*============================================================================*
 *  Function name			:	archive::reinit
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Reinitialization
 *============================================================================*/
void archive::defineBoundaries(int maxDomain)
{
	
	maxDomain = upperBound;
	
	//-- Reinitialize the boundaries
	for(int m =0; m<nMemetics; m++)   { boundaries[m] = 0; }
	nMemetics = 0;
	
	//-- Redefine the boundaries (upper bound = 1000)
	int position = 0;
	for(int m =0; m<9999; m++)
	{	position += widthMemetics;
		if (position <= maxDomain)   { boundaries[m] = position; 
										nMemetics     = m+1;	  }
		else						  { break;					  }
	}
	
	//-- Adjust the last boundaries up to infinity!
	boundaries[nMemetics] = 99999;
	nMemetics++;
	
 	if(dbug)
 	{
		printf("Nb. of memetics: %d - boundaries:", nMemetics);
		for(int m=0; m<nMemetics; m++)   { printf(" %d", boundaries[m]); }
		printf("\n");
 	}
}
/*============================================================================*
 *  Function name			:	archive::checkDominance
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Check if the solution is to be added in the 
 *  							archive.
 *----------------------------------------------------------------------------*
 *  Argument				:   particle *test        - Tested particle
 *----------------------------------------------------------------------------*
 *  Returns					:   int toBeReinitialized - Indicates what to do
 *============================================================================*/
int archive::checkDominance(particle *test)
{
	int nObjectives = test->nObjectives,
		toBeReinitialized = 0;
	float err, errTest, sze, szeTest;

	//o--------------- Definition of the objective values ---------------o//
	//o--------- Classification rate and normalized compression ---------o//
	//-- Personal best
	if(test->nns.sizeF2)
	{	errTest = 1 - test->sFitness.clsRate;
		szeTest =     test->sFitness.sizeNn;   }
	else
	{	return 0;							   }

	//o--- Check for which memetic the solution belongs (compression) ---o//
	int cMemetic = whichMemetic(szeTest);

	//-- Boundaries addresses
	int addrStart =  cMemetic * sizeMemetics,
		addrStop  = addrStart + nMembers[cMemetic];

	int flagNonDominated = 0;
	if(errTest || szeTest)
	{
		if(dbug)   { printf("  mem %d (%d to %d)\n",
							cMemetic, addrStart, addrStop-1); }

		//o-------------------- Check for dominance ---------------------o//
		flagNonDominated = 1;	//- Dominant by default
		for(int m = addrStart; m < addrStop; m++)
		{
			err = 1 - p[m].sFitness.clsRate;
			sze =     p[m].sFitness.sizeNn;

			if(dbug)  { printf("    Check dom%d - tested: (%1.3f, %1.3f)"
							   "-vs-(%1.3f, %1.3f) -> ",
							   m, errTest, szeTest, err, sze); }

			//o------------------- Check if dominated -------------------o//
			//-- Same
			if( errTest == err   &&   szeTest == sze )
			{   flagNonDominated = 0;
				if(dbug)   { printf(" same\n"); }
				break;  }
			//-- Dominated
			if( errTest >= err   &&   szeTest >= sze )
			{   flagNonDominated = 0;
				if(dbug)   { printf(" dominated\n"); }
				break;  }
			//-- Old solution is dominated -> we remove and reverify
			else if(errTest <= err   &&   szeTest <= sze)
			{
				if(dbug)   { printf(" remove dominated\n"); }
				remove(m);   addrStop--;   m--;  }
			else
			{	if(dbug)   { printf(" equivalent\n"); }   }
		}
	}

	//o------------------- Check if it must be added --------------------o//
	int flagAdd = 0;
	if(flagNonDominated)
	{
		//o---------- If the memetic if full - Its cuboid time! ---------o//
		if( nMembers[cMemetic] == sizeMemetics )
		{
			//-- Same solution - personal best -vs current position
			if(szeTest == p[addrStart].sFitness.normCpn &&
			   szeTest == p[addrStop-1].sFitness.normCpn)
			{   flagAdd = 0;   }
			//-- Extreme left - replace left solution
			else if( szeTest < p[addrStart].sFitness.sizeNn  )
			{ remove(addrStart);      flagAdd = 1;   }
			else if( szeTest > p[addrStop-1].sFitness.sizeNn )
			{   remove(addrStop-1);   flagAdd = 1;   }
			//-- Define the tested fitness - compare with the rest
			else
			{
				//o------ Crowding distance for the added particle ------o//
				//-- We search for the nearest solution
				float minDist   = 9999, szeNear;
				int   nearest   = -1;
				for(int m = addrStart; m < addrStop; m++)
				{
					float dist = pow(errTest - p[m].sFitness.clsRate, 2) +
								 pow(szeTest - p[m].sFitness.normCpn, 2);

					if(dist < minDist)
					{	minDist = dist;   nearest = m;
						szeNear = p[m].sFitness.normCpn;
					}
				}

				//-- The so called cuboid area is a freaking L1 distance!!!!
				int n;
				if(szeTest < szeNear)   { n = nearest-1; }
				else                    { n = nearest;   }

				float cubTest =
				  pow(p[n].sFitness.clsRate - p[n+1].sFitness.clsRate,2) +
				  pow(p[n].sFitness.normCpn  - p[n+1].sFitness.normCpn, 2);

				//o----------- Crowding distance for the rest -----------o//
				float cubWorst = 9999;
				int worst = -1;
				for(int m = addrStart+1; m < addrStop-1; m++)
				{
					float cuboid = pow(p[m-1].sFitness.clsRate -
									   p[m+1].sFitness.clsRate, 2) +
								   pow(p[m-1].sFitness.normCpn  -
									   p[m+1].sFitness.normCpn, 2);

					if(cuboid<cubWorst)   { cubWorst = cuboid;  worst = m; }
				}

				//-- Make place for the new (with better crowding)!
				if(cubWorst < cubTest)
				{
					if(dbug)
					{ printf("    Remove diversity (%1.3f, %1.3f)\n",
					  1-p[worst].sFitness.clsRate,
					  p[worst].sFitness.sizeNn); }
					remove(worst);   flagAdd = 1;
				}
				//-- Tie breaker - max of the min technique
				else if(worst == nearest)
				{
					float minLeft, minRght, minTest, minWorst;

					minLeft = pow(errTest - p[worst-1].sFitness.clsRate,2) +
							  pow(szeTest - p[worst-1].sFitness.normCpn,2);

					minRght = pow(errTest - p[worst+1].sFitness.clsRate,2) +
								pow(szeTest - p[worst+1].sFitness.normCpn,2);

					minTest = min(minLeft,minRght);

					minLeft = pow(p[worst  ].sFitness.clsRate -
								  p[worst-1].sFitness.clsRate,2) +
							  pow(p[worst  ].sFitness.normCpn -
								  p[worst-1].sFitness.normCpn,2);

					minRght = pow(p[worst ].sFitness.clsRate -
								  p[worst+1].sFitness.clsRate,2) + 
							  pow(p[worst  ].sFitness.normCpn -
								  p[worst+1].sFitness.normCpn,2);

					minWorst = min(minLeft, minRght);

					if(minTest > minWorst)
					{
						if(dbug)   { printf("    Remove diversity ",
											"tie-breaker (%1.2f, %1.2f)\n",
											1-p[worst].sFitness.clsRate,
											p[worst].sFitness.sizeNn); }
						remove(worst);   flagAdd = 1;
					}

					//-- If the memetic is full and and the new solution
					//   does not increase diversity -> reinitialize
					//   (for the current position only)
					if( !flagAdd )
					{ toBeReinitialized = 1; 
					  if(dbug)   { printf("to be reinit\n"); }
					}
				}
			}
		}
		//-- If the memetic is NOT full - add automatically
		else   { flagAdd = 1; }
	}

	//o---------------------- Add the new solution ----------------------o//
	if(flagAdd)   { add(test, cMemetic); }

	if(dbug)   { for(int m = 0 ; m < size; m++)
				 { if(filled[m])
				   { printf("    %d - %1.2f, %1.2f <- %d\n", m,
						1-p[m].sFitness.clsRate, p[m].sFitness.sizeNn,
						p[m].nns.sizeF2);                              }
				 }
			   }

	return toBeReinitialized;
}
/*============================================================================*
 *  Function name			:	archive::whichMemetic
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Check in which memetic is the tested
 *  							compression
 *----------------------------------------------------------------------------*
 *  Argument				:   float testedSze - D'hu!
 *----------------------------------------------------------------------------*
 *  Returns					:	int memetic    - D'hu! (bis)
 *============================================================================*/
int archive::whichMemetic(float testedSze)
{
	int cMemetic = (int)floor( testedSze / (float)widthMemetics );
	if(cMemetic > nMemetics)   { cMemetic = nMemetics-1; }
	return cMemetic;
}
/*============================================================================*
 *  Function name			:	optimization::updatePopulation
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Add particles in an orderly fashion
 *----------------------------------------------------------------------------*
 *  Argument				:   particle *test - Tested particle
 *  							int cMemetic   - Current memetic
 *  							int o          - Objective
 *============================================================================*/
void archive::update()
{
	//-- Copy the current archive
	for(int m = 0; m < size; m++)
	{
		//-- Reinitialize the old
		o[m].nns.reinit();
		o[m].sFitness.reinit();
		for(int d=0; d<nDimensions; d++)   { o[m].s[d] = 0; }

		//-- Copy in the old if necessary
		filledOld[m] = filled[m];
		if(filled[m])
		{	o[m].nns.copy( &(p[m].nns ) );
			o[m].sFitness.copy( &(p[m].sFitness) );
			for(int d=0; d<nDimensions; d++)   { o[m].s[d] = p[m].s[d]; }
		}
	}

	//-- Remove every thing in the current archive
	for(int m=size-1; m>-1; m--)   { remove(m); }

	//-- Update the current archive accordingly
	if(dbug)   { printf("Updated solutions\n");
				 for(int m = 0; m < size; m++)
				 { if(filled[m])
				   { printf(" m:%d - err. rate:%1.2f, size: %1.0f - %d\n",
						    m, 1-o[m].sFitness.clsRate, o[m].sFitness.sizeNn,
			         filledOld[m]);  }
				 }
			   }

	for(int m = 0; m < size; m++)
	{	if(filledOld[m])
		{
			if(dbug)   { printf("Upd %d - \n",m); }
			checkDominance( &o[m] );
		}
	}

/*	if(dbug)
	{	printf("\n/*-- Update             /*-- After\n");
		for(int m = 0 ; m < size; m++)
		{	printf("    %d - %1.4f, %1.4f - %d    %1.4f, %1.4f - %d\n", m,
					1-o[m].sFitness.clsRate, o[m].sFitness.sizeNn, filledOld[m],
					1-p[m].sFitness.clsRate, p[m].sFitness.sizeNn, filled[m]);
		}
	}
*/
}
/*============================================================================*
 *  Function name			:	optimization::updatePopulation
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Add particles in an orderly fashion
 *----------------------------------------------------------------------------*
 *  Argument				:   particle *test - Tested particle
 *  							int cMemetic   - Current memetic
 *  							int o          - Objective
 *============================================================================*/
void archive::add(particle *test, int cMemetic)
{
	//-- Boundaries addresses
	int addrStart = cMemetic  * sizeMemetics,
		addrStop  = addrStart + nMembers[cMemetic];

	float errTest = 1 - test->sFitness.clsRate,
	      szeTest =     test->sFitness.sizeNn;

	//o---------------------- Find the insertion point ----------------------o//
	int addrAdd = addrStart;
	for(int m = addrStart; m < addrStop; m++)
	{	//-- If a smaller network size, we add here, else we increment.
		if(szeTest < p[m].sFitness.sizeNn )   { break;     }
		else 								  { addrAdd++; }
	}

	if(dbug)   { printf("    +++ @ %d, cMem %d, nMembers %d, "
				 "nFilled %d (%1.2f, %1.2f)\n", addrAdd, cMemetic,
				 nMembers[cMemetic], nFilled, errTest, szeTest); }

	//o---------------------- Push the other solutions ----------------------o//
	for(int m = addrStop; addrAdd < m; m--)
	{
		if(dbug)   { printf("    pushing (%1.2f, %1.2f ) from %d to %d \n",
							errTest, szeTest, m-1, m); }
		filled[m] = filled[m-1];
		p[m].nns.     copy( &(p[m-1].nns     ) );
		p[m].sFitness.copy( &(p[m-1].sFitness) );
		for(int d=0; d<nDimensions; d++)   { p[m].s[d] = p[m-1].s[d]; }
	}
	//o--------------------------- Add the new one --------------------------o//
	nFilled++;
	nMembers[cMemetic]++;
	filled[addrAdd] = 1;
	
	p[addrAdd].nns.copy( &(test->nns) );
	p[addrAdd].sFitness.copy( &(test->sFitness) );
	for(int d=0; d<nDimensions; d++)   { p[addrAdd].s[d] = test->s[d]; }
	
	if(dbug)   { printf("    added @ %d, cMem %d, nMembers %d, "
				 "nFilled %d (%1.2f, %d)\n", addrAdd, cMemetic,
				 nMembers[cMemetic], nFilled,
				 1-p[addrAdd].sFitness.clsRate, p[addrAdd].nns.sizeF2); }
}
/*============================================================================*
 *  Function name			:	optimization::updatePopulation
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Update each particle current & best position
 *  							fitness (necessary for dynamic environment)
 *----------------------------------------------------------------------------*
 *  Argument				:   int addr     - Removed address
 *  							int cMemetic - Current memetic
 *============================================================================*/
void archive::remove(int addr)
{
	//-- Find the memetic
	int cMemetic = (int)floor( (float)addr / (float)sizeMemetics );

	//-- Erase only existing networks
	if( filled[addr] )
	{
		nMembers[cMemetic]--;
		nFilled--;
		filled[addr] = 0;
		p[addr].nns.reinit();
		p[addr].sFitness.reinit();
		for(int d=0; d<nDimensions; d++)   { p[addr].s[d] = -1; }
	}

	//o----- If removed WAS not last -> pull to fill in the blank space -----o//
	int last = cMemetic *sizeMemetics + nMembers[cMemetic];

	if(dbug)   { printf("    --- @ addr %d (last %d), cMem %d, nMembers %d, "
						"nFilled %d\n", addr, last, cMemetic,
					    nMembers[cMemetic], nFilled); }

	if(addr != last)
	{	for(int m = addr; m < last; m++)
		{
			if(dbug)   { printf("    pulling from %d to %d (%1.2f, %1.2f)\n",
								m+1, m, 1-p[m+1].sFitness.clsRate,
								p[m+1].sFitness.sizeNn); }
			filled[m] = filled[m+1];
			p[m].nns.     copy( &(p[m+1].nns     ) );
			p[m].sFitness.copy( &(p[m+1].sFitness) );
			for(int d=0; d<nDimensions; d++)   { p[m].s[d] = p[m+1].s[d]; }
		}

		//-- Erase the last one
		filled[last] = 0;
		p[last].nns.reinit();
		p[last].sFitness.reinit();
		for(int d=0; d<nDimensions; d++)   { p[last].s[d] = -1; }
	}
}
/*============================================================================*
 *  Function name			:	optimization::savePopulation
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Saves the population
 *----------------------------------------------------------------------------*
 *  Arguments				:   char *nameFile
 *  							int   currentRep  - Current replication
 *============================================================================*/
void archive::save(char *nameFile, int currentRep)
{
	//-- File opening
	ofstream f(nameFile, ios::app | ios::binary);
	if(!f)   { printf("Error! Can't open the file %s \n", nameFile);  exit(1); }

	//-- General informations
	f.write((char*)&currentRep,    sizeof(int));
	f.write((char*)&size,          sizeof(int));
	f.write((char*)&sizeMemetics,  sizeof(int));
	f.write((char*)&widthMemetics, sizeof(int));
	f.write((char*)&nMemetics,     sizeof(int));
	f.write((char*)&nDimensions,   sizeof(int));

	f.write((char*) boundaries,  sizeof(int)   * (unsigned)nMemetics);
	f.write((char*)&nFilled,     sizeof(int)                        );
	f.write((char*) filled,      sizeof(int)   * (unsigned)size     );
	f.write((char*) nMembers,    sizeof(int)   * (unsigned)nMemetics);

	for(int m=0; m<size; m++)
	{	if(filled[m])
		{	//-- Position, fitness for the two objectives
			f.write((char*)  p[m].s, sizeof(float) * (unsigned)nDimensions);
			f.write((char*)&(p[m].sFitness.clsRate), sizeof(float));
			f.write((char*)&(p[m].sFitness.sizeNn ), sizeof(float));
		}
	}

	//-- File closing
	f.close();
}
