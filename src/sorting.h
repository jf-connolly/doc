/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  sort.h - Ascending sort with a tie breaker.									  *
 *     																		  *
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

#ifndef SORT_H
#define SORT_H

using namespace std;

class sorting {
public:
	/*------------------------------ Variables -------------------------------*/
	int	   size,			//- Number of elements to sort
		  *order;			//- Ascending order

	float *fitness,			//- Values to sort
	      *tieBreaker;		//- D'hu!

	/*---------------------------- Functions used ----------------------------*/
	//-- Constructors & destructor
	//-- Resultats
	sorting();
	~sorting();

	//-- Functions: utilities
	void initialization(int sz);
	void heapSort();

	/*--------------------- Functions to do the heapsort ----------------------*/
	void doHeap();
	void siftDown(int start, int end);
	void swap    (int *a, int *b);
};
#endif
/*----------------------------------------------------------------------------*/
