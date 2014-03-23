/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  dbase.c - Database                                                        *
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
#include <stdio.h>

//-- Project specific
#include "dbase.h"

/*============================================================================*
 *  Function name			:	dbase::dbase
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Default constructor
 *============================================================================*/
dbase::dbase()
{
	error       = 200;

	nPatterns   = 0;
	nFeatures   = 0;
	nClasses    = 0;

	sizeClasses   = (int*)NULL;
	addrClasses   = (int*)NULL;
	activeClasses = (int*)NULL;

	labels		= (int*)  NULL;
	data		= (float*)NULL;
}
/*============================================================================*
 *  Function name			:	dbase::initialization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Default constructor
 *----------------------------------------------------------------------------*
 *  Argument				:   char nameFile   - File containing the database
 *============================================================================*/
void dbase::initialization(char* nameFile)
{
	float *doubPtr;
	FILE*	f;
	size_t  read;

	//-- Binary file opening for reading (rb)
	f = fopen(nameFile,"rb");

	//-- File existence
	if(f == NULL)	{ error += 1; return; }

	doubPtr = data;
	/*-------------------------- Reading the header --------------------------*/
	//-- Reading the number of patterns
	read = fread(&nPatterns, sizeof(int), 1, f);
	if(read != 1) {	error +=2; return; }

	//-- Reading the number of features
	read = fread(&nFeatures, sizeof(int), 1, f);
	if(read != 1) {	error += 2; return; }

	//-- Reading the number of classes
	read = fread(&nClasses, sizeof(int), 1, f);
	if(read != 1) {	error += 2; return; }

	/*------------------------- Reading the patterns -------------------------*/
	//-- Memory allocation: (float) x 2 x (n Pattern) x (m features)
	data = new float[2 *nFeatures *nPatterns];
	if(data == NULL)   { error += 3; return; }

	//-- Complement coding
	doubPtr = data;
	for(int p=0; p<nPatterns; p++)
	{
		//-- Reading a pattern
		read = fread(doubPtr, sizeof(float), nFeatures, f);
		if(read != (unsigned)nFeatures ) { error += 4; return; }

		//-- Positioning to complement
		doubPtr += nFeatures;

		//-- Complement
		for(int i=0; i<nFeatures; i++)
		{	*doubPtr = 1 - *(doubPtr - nFeatures);
			doubPtr ++; //-- We reposition each time
		}
	}

	/*-------------------------- Reading the labels --------------------------*/
	//-- Memory allocation: (int) x (n Pattern)
	labels = new int[nPatterns];
	if(labels == NULL)   { error += 3; return; }

	//-- Lecture des etiquettes
	read = fread(labels, sizeof(int), nPatterns, f);
	if(read !=	(unsigned)nPatterns )   { error += 4; return; }

	/*-------------------- Number of patterns per class ----------------------*/
	//-- Memory allocation: (int) x (n classes)
	sizeClasses = new int[nClasses];
	if(sizeClasses == NULL)   { error += 3; return; }

	//-- Initialization and assignement
	for(int k=0; k<nClasses;  k++)     { sizeClasses[k]=0; }
	for(int p=0; p<nPatterns; p++)
	{ sizeClasses[labels[p]] = sizeClasses[labels[p]] + 1; }

	/*-------------- Addresses of the beginning of each classes --------------*/
	//-- Memory allocation: (int) x (n classes)
	addrClasses = new int[nClasses];
	if(addrClasses == NULL)   { error += 3; return; }

	//-- Initialization and assignement
	int addrTemp = 0;
	for(int k=0; k<nClasses;  k++)
	{	addrClasses[k] = addrTemp;   addrTemp += sizeClasses[k];   }

	//-- Initial class activity
	activeClasses = new int[nClasses];
	if(activeClasses == NULL)	{ error += 3; return; }

	for(int k=0; k<nClasses;  k++)   { activeClasses[k] = 0; }

	/*--------------------------------- End ----------------------------------*/
	//- Closing the file
	if( fclose(f) )   { error += 6; return; }
}
/*============================================================================*
 *  Function name			:	dbase::initialiHeader
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Default constructor
 *----------------------------------------------------------------------------*
 *  Argument				:   char nameFile   - File containing the database
 *============================================================================*/
void dbase::initialiHeader(char* nameFile)
{
	FILE*	f;
	size_t  read;

	//-- Binary file opening for reading (rb)
	f = fopen(nameFile,"rb");

	//-- File existence
	if(f == NULL)	{ error += 1; return; }

	/*-------------------------- Reading the header --------------------------*/
	//-- Reading the number of patterns
	read = fread(&nPatterns, sizeof(int), 1, f);
	if(read != 1) {	error +=2; return; }

	//-- Reading the number of features
	read = fread(&nFeatures, sizeof(int), 1, f);
	if(read != 1) {	error += 2; return; }

	//-- Reading the number of classes
	read = fread(&nClasses, sizeof(int), 1, f);
	if(read != 1) {	error += 2; return; }

	/*--------------------------------- End ----------------------------------*/
	//- Closing the file
	if( fclose(f) )		{ error += 6; return; }
}
/*============================================================================*
 *  Function name			:	dbase::~dbase
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Destructor
 *============================================================================*/
dbase::~dbase()
{
	//-- Freeing memory
	if( data != NULL )
	{
		delete [] sizeClasses;
		delete [] addrClasses;
		delete [] activeClasses;

		delete [] data;
		delete [] labels;
	}

	//-- Re-initializing constants
	nClasses	= 0;
	nFeatures	= 0;
	nPatterns	= 0;
	error       = 0;
}
/*============================================================================*
 *  Function name			:	dbase::getPatternAddr
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Gets a a pattern's address
 *----------------------------------------------------------------------------*
 *  Argument				:   int    p
 *  Returns					:   float  the address of a pattern
 *============================================================================*/
float* dbase::getPatternAddr(int p) { return data + (2 *nFeatures *p); }
/*----------------------------------------------------------------------------*/
