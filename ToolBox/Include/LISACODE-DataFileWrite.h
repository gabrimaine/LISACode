/*
 *  LISACODE-DataFileWrite.h
 *  LC20_xcode
 *
 *  Created by Antoine Petiteau on 15/02/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \ingroup ToolBox 
 * \defgroup LCDataFileWrite_  Class LCDataFileWrite
 * (See class #LCDataFileWrite for a detailed description)
 * \{
 */

#ifndef __LCDATAFILEWRITE_H
#define __LCDATAFILEWRITE_H

#include "LISACODE-Tools.h"
#include "LISACODE-Serie2.h"

class LCDataFileWrite 
{
protected:
	
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	
	/** \brief Number of records */
	int NRec;
	
	/** \brief List of current value */
	double ** RecData;
	
	/** \brief Name of record. Each name have a maximal size of 128 characters */
	char ** RecName;
	
	/** \brief Output file */ 
	std::ofstream fOut; 
	
	/** \brief Name of output file */
	char fNOut[2048];
	
	/** \brief Type of output file. Coding :
	 *	 ASCII  : ASCII
	 *	 BINARY : Binary with included header (3 first lines of the file)
	 *	 XML    : Binary without 3 lines header
	 */
	TypeFile typeOutFile;
	
	
	/** \brief If true, write reference */
	bool WriteRef;
	
	/** \brief Initial time (to write in the header and computing time if needed) */ 
	double t0;
	
	/** \brief Time step (to write in the header and computing time if needed) */
	double dt;
	
	
	/** \brief Number of data (to write in the header and computing time if needed) */
	int NDatExpect;
	
	/** \brief Number of data already write in file. Note : -1 if the file has not been opened */
	int iDatWrite; 
	
public:
	/********** Constructor **********/
	
	/** \brief Default constructor */
	LCDataFileWrite();
	
	/** \brief Constructor with name of output file */
	LCDataFileWrite(LCTools * MT_n, const char * fNOut_n, TypeFile typeOutFile_n);
	
	/** \brief Destructor */
	~LCDataFileWrite();
	
	
	/**********  Access methods  **********/
	void sett0(double t0_n){t0=t0_n;};
	void setdt(double dt_n){dt=dt_n;};
	void setNDatExpect(int NDatExpect_n){NDatExpect=NDatExpect_n;};
	
	/**********  Initialisation methods  **********/
	
	/** \brief Initialization of everything at NULL or equivalent 
	 *	@param[in] CleanMem  If true clean all memory
	 */
	void initNULL(bool CleanMem);
	
	/**  \brief Add one record 
	 *	RETURN : Adress of the record location
	 */
	double * AddRecord(char * RecNameAdd);
	
	/** \brief Initialization : last step before using during the running */
	void init(std::ostream * GlobalHeader, int indent);
	
	
	/**********  Running methods  **********/
	
	
	/** \brief Add one record : record the data in the file */
	void RecordData();
	
	/** \brief Add one record : record the data in the file 
	 *	vRef : IN : Value of the reference (ex: time)
	 */
	void RecordData(double vRef);
	
	
	/**********  Other methods  **********/
	/**  \brief Write the XML header bloc corresponding to the file in the output with a NIndent tabulation at the beginning of each line */
	void WriteXMLHeader(std::ostream * fHeader, int NIndent);
	
	
	/** \brief Display info */
	void DispInfo(char * BTab);
	
};

#endif //__LCDATAFILEWRITE_H

/**\}*/

// end of LISACODE-DataFileWrite.h