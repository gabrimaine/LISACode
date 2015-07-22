/*
 *  LISACODE-DataFileRead.h
 *  LC20
 *
 *  Created by Antoine Petiteau on 15/02/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */


/** \ingroup ToolBox 
 * \defgroup LCDataFileRead_  Class LCDataFileRead
 * (See class #LCDataFileRead for a detailed description)
 * \{
 */

#ifndef __LCDATAFILEREAD_H
#define __LCDATAFILEREAD_H

#include "LISACODE-Tools.h"
#include "LISACODE-Serie2.h"


class LCDataFileRead 
{
protected:
	/** \brief Pointer on toolbox class */
	LCTools * MT;
	
	/** \brief Name of input file */
	char fNIn[2048];
	
	/** \brief File type : 
	 * - ASCII  : ASCII,
	 * - BINARY : BINARY,
	 * - XML : XML (binary file + XML header in another file). 
	 * If the header of the 'xml' file is the noisexmlbloc passed by config, the type can be ASCII, BINARY or XML.
	 */
	TypeFile FileType;
	
	/*! \brief File encoding :
	 * - BIGENDIAN 
	 * - LITTLEENDIAN 
	 */
	TypeEncoding FileEncoding; 
	
	/*! \brief Number of data : number of lines */
	long NData;
	
	/*! \brief Number of records : number of columns ( -1 if the time is present) */
	int NRec;
		
	/*! \brief Reference step */
	double dx;

	/*! \brief Initial reference value */
	double x0;
	
	/*! \brief Final reference value */
	double xend;
	
	/*! \brief True if we don't want to store the first record (usually because it's the reference) */
	bool SkipFirstRec;
	
	/*! \brief Contain header words (only read with ASCII file) */
	std::vector<char*> Header; 
	
	/*! \brief Position of the beginning of the binary part in the binary file  */
	long PosStartBin;
	
	/*! \brief List of series which contain the data */
	LCSerie2 ** Data;
	
	
public:
	/********** Constructor **********/
	
	/** \brief Default constructor */
	LCDataFileRead();
	
	/** \brief Standard constructor */
	LCDataFileRead(LCTools * MT_n);
	
	/** \brief Standard constructor */
	LCDataFileRead(LCTools * MT_n, char * NewFileName);
	
	/** \brief Destructor */
	~LCDataFileRead();
	
	/** \brief Initialization of everything at NULL or equivalent 
	 *	@param[in] CleanMem  If true clean all memory
	 */
	void initNULL(bool CleanMem);
	
	
	/**********  Access methods  **********/
	double getx0(){return(x0);};
	double getdx(){return(dx);};
	double getxend(){return(xend);};
	int getNRec() {return(NRec);};
	int getNDat() {return(NData);};
	int ReadFirstCol() { return(!SkipFirstRec); }; 
	
	/*! \brief Return the pointer on signal iSig */
	LCSerie2 * getPtrData(int iSig) { return(Data[iSig]); };
	
	
	/*! \brief Copy the word of the header in list of words 
	 *	@param[out]	Words		List of words
	 *	@param[out]	NWords		Number of words
	 *	@param[in]	NCharWord	Number of character in a word
	 */
	void getHeader(char ** & Words, int & NWords, int NCharWord);
    
    /*! \brief Copy the file name
     *	@param[out]	FileNameCopy	Copy of file name
     */
    void getFileName(char * FileNameCopy) {strcpy(FileNameCopy, fNIn); };
	
	/** \brief Set the name of the input file */
	void setFile(char * NewFileName);
	
	/**********  Configuration methods  **********/
	
	
	/*! \brief Configuration from the xml bloc : <XSIL Name="..." Type="Noise"> ... </XSIL> */
	void config(ezxml_t noisexmlbloc);
	
	
	/**********  Linking and initalization methods  **********/	
		
	/** \brief Initialization : last step before using during the running */
	void init();
	
	/** \brief Initialization but just return pointer on data and associated size */
	void init(double ** & DataExt, int & NRecExt, int & NDataExt);
	
	
	/**********  Other methods  **********/
	
	/*! \brief Analyzing file's parameters */
	void FileAnalysis();
	
	/*! \brief Analyzing ASCII file's parameters */
	void FileAnalysisASCII();
	
	/*! \brief Analysis Binary file's parameters */
	void FileAnalysisBinary();
	
	/*! \brief Analysis XML file's parameters */
	void FileAnalysisXML();
	
	/*! \brief Read data only for the specified record */
	void ReadDataSeries(double ** ExtData);
	
	/*! \brief Read data from ASCII file (only for the specified record) */
	void ReadDataSeriesASCII(double ** ExtData);
	
	/*! \brief Read data from Binary file (only for the specified record) */
	void ReadDataSeriesBinary(double ** ExtData);
	
	/*! \brief Display file's parameters  */
	void ControlDisplay();

	/**********  Access data methods  **********/
	
	/** \brief Return the corresponding data  
	 *	@param[in]	iS				Index of signal (or record)
	 *	@param[in]	t				Time
	 *	@param[in]	InterpType		Type of interpolation
	 *	@param[in]	InterpUtilValue Value used for interpolation (for Lagrange interpolation it's the order)
	 *	@return Corresponding value
	 */
	double gData(int iS, double t, INTERP InterpType, double InterpUtilValue);
	
	/** \brief Return the corresponding data
	 *	@param[in]	iS				Index of signal (or record)
	 *	@param[in]	iT				Index of time
	 *	@return Corresponding value
	 */
	double gDataBin(int iS, int iT);
	
	/** \brief Return the corresponding refernce 
	 *	@param[in]	iS				Index of signal (or record)
	 *	@param[in]	iT				Index of time
	 *	@return Corresponding reference
	 */
	double gRefBin(int iT){ return(x0+iT*dx);};
	
	/*! \brief Compute raw derivative using current bin and bin at -NDbin */
	double DerivRawCur(int iS, int NDbin);
	
	/*! \brief Compute backward derivative at second order of precision (use the 3 last bins)  */
	double DerivBackOrder2Cur(int iS);
	
	/*! \brief Compute raw derivative using x and x - Dx points,   */
	double DerivRawSpe(int iS, double x, double Dx, INTERP InterpType, double InterpUtilValue);
	
	/*! \brief Compute backward derivative et second order of precision (use x, x-Dx and x-2Dx points)  */
	double DerivBackOrder2Spe(int iS, double x, double Dx, INTERP InterpType, double InterpUtilValue);
	
	
};

#endif //__LCDATAFILEREAD_H

/**\}*/

// end of LISACODE-DataFileRead.h