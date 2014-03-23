/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  dbase.h - Database                                                        *
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

#ifndef DBASE_H
#define DBASE_H

using namespace std;

class dbase {
public:
	/*------------------------------ Variables -------------------------------*/
	int		nPatterns;		//- Number of patterns in the DB
	int		nFeatures;		//- Number of features in the DB
	int		nClasses;		//- Number of classes in the DB

	int	   *sizeClasses;	//- Number of patterns per class
	int    *addrClasses;	//- Addresses of the classes
	int    *activeClasses;  //- Is that class part of the system

	int	   *labels;			//- Label (or tags)
	float  *data;			//- Data - Vector of size 2*nFeatures*nPatterns

	int		error;          //- Error code

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors & destructor
	dbase();
	~dbase();
	void initialization(char* nameFile);
	void initialiHeader(char* nameFile);

	//-- Utilities
	float* getPatternAddr(int p);
};
#endif
