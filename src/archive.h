/******************************************************************************
 *   Copyright(C) 2008 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca *
 *																			  *
 *   archive.h - Archive for multiobjective optimization					  *
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
#include "particle.h"
#include "sorting.h"

#ifndef ARCHIVE_H
#define ARCHIVE_H

using namespace std;

/*--------------------------------- Class PSO --------------------------------*/
class archive
{
public:
	/*------------------------------ Variables -------------------------------*/
	int size,			//- Maximal size of the archive (current = nFilled)
		sizeMemetics,	//- Size of each memetic
		widthMemetics,	//- Width of each memetic
		nMemetics,		//- Number of categories
		upperBound,
		nDimensions;	//- Number of dimensions

	int nFilled,		//- Number of particles in the archive
	   *filled,			//- Indicate if slot is filled
	   *filledOld,		//- Indicate if slot is filled (old)
	   *nMembers,		//- Number of members for each memetic
	   *members;		//- Active members of the archive

	particle *p, *o;	//- Current and old archive

	int *boundaries;	//- Fitness of the tested member
	
	int dbug;

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors and destructor
	archive();
	void initialization(configuration *cfg);
	~archive();

	//-- Other
	void reinit();
	void defineBoundaries(int time);
	int  checkDominance(particle *tested);
	void update();
	int  whichMemetic(float testedCpn);
	void add(particle *test, int cMemetic);
	void remove(int addr);
	void save(char *nameFile, int currentRep);
};
#endif
