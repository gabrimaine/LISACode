// $Id:  $
/*
 *  LISACODE-Tools.cpp
 *  V 2.0
 *
 *  Created on 28/05/09 by  Antoine Petiteau (AEI)
 *  Last modification on 26/01/10 by Antoine Petiteau (AEI)
 *
 */

#include "LISACODE-Tools.h"

// *****************
// *  Constructor  *
// *****************

/** 
 * Constructs an instance and initializes it with default values. 
 */
LCTools::LCTools()
{
	_MemAllocSize_ = 0;
	_MemAllocSizeMax_ = 0;
	
	DispLog = 1;
	fdout = NULL;
	o = &std::cout;
	DispDetails = false;
	
	LocTools = false;
	
	LocalEndianness = testbyteorder();
	
    for(int k=0; k<=LC::ORDERMAXLAG; k++)
        LagPolDen[k] = NULL;
	
	FTFwdNbtDat = 0;
	FTFwdNbfDat = 0;
	FTFwdBestFFTW = true;
	FTFwdtimeData = NULL;
	FTFwdfreqData = NULL;
	FTFwdPlan = NULL;
	
	FTBckNbtDat = 0;
	FTBckNbfDat = 0;
	FTBckBestFFTW = true;
	FTBcktimeData = NULL;
	FTBckfreqData = NULL;
	FTBckPlan = NULL;
	
	
	
	mpi_id = 0;
	mpi_np = 1;
}



/**
 * 
 */
LCTools::~LCTools ()
{
	
	if(FTFwdtimeData != NULL)
		Free(FTFwdtimeData, FTFwdNbtDat*sizeof(double));
	if(FTFwdfreqData != NULL)
		Free(FTFwdfreqData, FTFwdNbfDat* sizeof(dcomplex));
	if (FTFwdPlan != NULL)
		fftw_destroy_plan(FTFwdPlan);
	
	if(FTBcktimeData != NULL)
		Free(FTBcktimeData, FTBckNbtDat*sizeof(double));
	if(FTBckfreqData != NULL)
		Free(FTBckfreqData, FTBckNbfDat* sizeof(dcomplex));
	if (FTBckPlan != NULL)
		fftw_destroy_plan(FTBckPlan);
	
	
	//MemDisplay
	
	
#ifdef _CHECK_MPIMEM_
	(*o) << "Memory exchange in MPI transfert (with master) (on " << MPIid() << ") = " << mpi_memtr/(1024.*1024.) << " Mo = " << mpi_memtr << " o" << Endl;
#endif // _CHECK_MPIMEM_
	if(!LocTools)
		(*o) << "\nMEMORY(core" << MPIid() <<"): Rest = " << _MemAllocSize_/(1024.*1024.) << " Mo = " << _MemAllocSize_ << " o (Max = " << _MemAllocSizeMax_/(1024.*1024.) << " Mo = " << _MemAllocSizeMax_ << " o)." << Endl;o->flush();
	unsetDisp();
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **********
// * Memory *
// **********

/** 
 * Allocation of memory using fftw malloc
 * and increment _MemAllocSize_ to follow the memory usage.
 */
void* LCTools::AllocMemory(size_t size) 
{
	void *ptr(NULL);
	
	if(size != 0){
		_MemAllocSize_ += size;
#ifdef _DEBUG_MEM_
		//if(size != 8)
		(*o) << "\nMEMORY(core" << MPIid() <<"): Alloc " << ptr << " for " << size << " o ( Total = " << _MemAllocSize_/(1024.*1024.) << " Mo = " << _MemAllocSize_ << " o )";o->flush();
#endif // _DEBUG_MEM_
		ptr = fftw_malloc(size);
#ifdef _DEBUG_MEM_
		//if(size != 8)
		(*o) << " --> OK" << Endl;o->flush();
#endif // _DEBUG_MEM_
		
		//if (size>100*(1024.*1024.))
		//	(*o) << "Stop + " << size << " o" << Endl;
		//if (size==40)
		//	(*o) << "Stop + " << size << " o" << Endl;
		
		
		if ( ptr == NULL ) {
			fprintf(stderr,"\n\n **** Memory allocation failed.  Unable to allocate %lf Mo.  Already allocated memory : %lf Mo ***\n\n",
					size/(1024.*1024.), _MemAllocSize_/(1024.*1024.));
			exit(1);
		}
		if(_MemAllocSize_>_MemAllocSizeMax_)
			_MemAllocSizeMax_ = _MemAllocSize_;
	}
	
	return ptr;
}


void* LCTools::ReAllocMemory(void * oldptr, size_t oldsize, size_t newsize)
{
	void *ptr;
	void *oldptr2;
	
	oldptr2 = oldptr;
	
	if(newsize<oldsize)
		throw std::invalid_argument("ERROR in LCTools::ReAllocMemory : The new size should be higher than the old one !");
	
	ptr = AllocMemory(newsize);
	if(oldsize>0){
		memcpy(ptr, oldptr2, oldsize);
		Free(oldptr2, oldsize);
	}
	
	return ptr;
}


/**
 * Free memory using free 
 * and decrement _MemAllocSize_ to follow the memory usage.
 */

void  LCTools::Free(void *ptr, size_t size)
{
	if(size!=0){
		_MemAllocSize_ -= size;
#ifdef _DEBUG_MEM_
		//if(size != 8)
		(*o) << "\nMEMORY(core" << MPIid() <<"): Free " << ptr << " for " << size << " o ( Total = " << _MemAllocSize_/(1024.*1024.) << " Mo = " << _MemAllocSize_ << " o )";o->flush();
#endif // _DEBUG_MEM_
		fftw_free(ptr);
#ifdef _DEBUG_MEM_
		//if(size != 8)
		(*o) << " --> OK" << Endl;o->flush();
#endif // _DEBUG_MEM_
		
		
		//if (size==40)
		//	(*o) << "Stop + " << size << " o" << Endl;
	}
}


/**
 * Display the current memory usage : value of _MemAllocSize_ in MegaOctet.
 */
void LCTools::MemDisplay() const
{
	if(Disp()){
		(*o) << " * Allocated memory  = " << _MemAllocSize_ / (1024.*1024.*sizeof(char)) << " Mo" << Endl;
		fflush(stdout);
	}
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ***********
// * Display *
// ***********

/**
 * 
 */
void LCTools::unsetDisp()
{
	if(DispLog == 2){
		fdout->close();
		fdout->clear();
		delete fdout;
	}
	DispLog = 0;
	o = &std::cerr;
}


/**
 * 
 */
void LCTools::setDispLog()
{
	if(DispLog == 1){
		fdout->close();
		fdout->clear();
		delete fdout;
	}
	DispLog = 1;
	o = &std::cout;
}


/**
 * 
 */
void LCTools::setDispInFile(char * OutFile, bool NewFile) 
{
	if(DispLog == 2){
		fdout->close();
		fdout->clear();
		delete fdout;
	}
	DispLog = 2;
	if(NewFile)
		fdout = new std::ofstream(OutFile);
	else
		fdout = new std::ofstream(OutFile,std::ios_base::app);
	o = fdout;
}

/**
 * 
 */
void LCTools::setDispInFileCoreid(char * OutFile, bool NewFile) 
{
	char OutFileId[10000];
	sprintf(OutFileId, "%s-%d",OutFile,MPIid());
	setDispInFile(OutFileId, NewFile);
}


bool LCTools::DispDet() 
{
	if( Disp() && DispDetails ) 
		return true;
	else 
		return false;
	
}



// **********
// * Random *
// **********

/**
 * @param[in] RandSeed_n New seed of random generator
 */
void LCTools::setRandSeed(long RandSeed_n)
{
	RandSeed = RandSeed_n;
	if(!MPImaster())
		RandSeed += MPIid();
	
	//srand48(RandSeed);
	
	//! ****** Set the seed of randlib generator
	/*!	We are doing :
	 *	RandSeed  = bit = b[n-1]b[n-2]...b1b0
	 *	is1		  = bit = b0b2...b[n/2-1]
	 *	is2       = bit = b1b3...b[n/2]
	 */
	long rs(RandSeed),is1(0),is2(0);
	for (int i=0;i<16;i++) {
		//! ** Take the lower weight bit and put them in correspondance with the higher weight bit in 16 bits ! 
		is1 |= (rs&1)<<(16-i); 
		//! ** Shift of 1 bit : equivalent to divide by 2
		rs>>=1; 
		//! ** Idem
		is2 |= (rs&1)<<(16-i);
		//! ** Again shift of 1 bit : equivalent to divide by 2
		rs>>=1; 
	}
	(*o) << "setRandSeed : " << RandSeed_n << " : " << is1 << "  " << is2 << Endl;
	setall(is1,is2);
	getsd(&is1,&is2);
	
}


/**
 * It use the standard random generator drand48()
 */
double LCTools::RandUniform(double xmin, double xmax)
{
	
	//return(xmin + drand48()*(xmax-xmin));
	return(xmin + genunf(0.0, 1.0)*(xmax-xmin));
	
}


/**
 * Some running details :
 * \arg Choose two random values  \f$ r_1 \f$ and \f$ r_2 \f$ in uniform distribution \f$ [ 0, 1 ] \f$
 * \arg Compute the random value as : \f$ x_0 + \sigma_x \sqrt{-2 \log r_2 \cos( 2 \pi r_1 ) } \f$
 */
double LCTools::RandGaussian(double center, double sigma)
{
	double r1, r2;
	r1 = RandUniform(0.0, 1.0);
	r2 = RandUniform(0.0, 1.0);
	return(center + sigma*sqrt(-2.0*log(r2))*cos(2*M_PI*r1));
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **************************
// * Double/integer/complex *
// **************************




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **********************
// * Mathematical tools *
// **********************

int LCTools::Gaussj(double **a, int n, double **b, int m)
{
	// The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting.
	int *indxc, *indxr, *ipiv;
	int icol, irow;
	double big, dum, pivinv, temp;
	
	// Allocate memory
	indxc = (int*) AllocMemory(sizeof(int) * n);
	indxr = (int*) AllocMemory(sizeof(int) * n);
	ipiv = (int*) AllocMemory(sizeof(int) * n);
	for(int j=0; j<n; j++) 
		ipiv[j]=0;
	
	for(int i=0; i<n; i++) { // ** This is the main loop over the colummns to be reduced.
		
		/* // * Control Display
		 for(int di=0; di<n; di++){
		 if(di==0) 
		 (*o) << "a = ";
		 else
		 (*o) << "    ";
		 for(int dj=0; dj<n; dj++)
		 printf(" %10lf", a[di][dj]);
		 if(di==0) 
		 (*o) << "     b = ";
		 else
		 (*o) << "         ";
		 for(int dj=0; dj<m; dj++)
		 printf(" %10lf", b[di][dj]);
		 (*o) << Endl;
		 }
		 */
		
		big=0.0;
		// ** Find the column where there is the biggest alpha[j][k] for the first loop, the second bigger for the second loop,...
		for(int j=0; j<n; j++) { // ** This is the outer loop of the search for a pivot element.
			if(ipiv[j] != 1) {
				for(int k=0; k<n; k++) {
					if(ipiv[k] == 0) {
						if(fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}else{ 
						if(ipiv[k] > 1) {
							//throw std::invalid_argument("Singular Matrix-1 !");
							(*o) << "ERROR in UTools::Gaussj : Singular Matrix-1 !" << Endl;
							return 1;
						}
					}
				}
			}
		}
		// ** Add 1 at the index of ipiv corresponding to the column where there are the biggest  alpha[j][k]
		++(ipiv[icol]);
		
		// ** We now have the pivot element, so we interchange rows, if needed, to put the pivot element on the diagonal.
		// The columns are not physically interchanged, only relabeled : indexc[i], the column of the ith pivot element, 
		// is the ith column that is reduced, while indxr[i] is the row in which that pivot element was originally located.
		// If indxr[i] =/= indxc[i] there is an implied column interchange. With this form of bookkeeping, the solution b's
		// will end up in the correct order, and the inverse matrix will be scrambled by columns.
		if (irow != icol) { // If we are not on the diagonal
			for (int l=0; l<n; l++) {
				SWAP2(a[irow][l],a[icol][l]); // Inverse in alpha the column irow with the column icol
			}
			for (int l=0; l<m; l++){
				SWAP2(b[irow][l],b[icol][l]); // Inverse in alpha the column irow with the column icol
			}
		}
		
		
		// ** We are now ready to divide the pivot row by the pivot element, located at irow and icol.
		indxr[i]=irow;
		indxc[i]=icol;
		if(a[icol][icol] == 0.0) {
			//throw std::invalid_argument("Singular Matrix-2 !");
			(*o) << "ERROR in UTools::Gaussj : Singular Matrix-2 !" << Endl;
			return 2;
		}else{
			pivinv=1.0/a[icol][icol];
		}
		a[icol][icol]=1.0;
		for(int l=0; l<n; l++){
			a[icol][l] *= pivinv;
		}
		for(int l=0; l<m; l++){
			b[icol][l] *= pivinv;
		}
		
		// ** Next, we reduce the rows ...
		for(int ll=0; ll<n; ll++){ 
			if(ll != icol) {		// ... except for the pivot one, of course.
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for(int l=0; l<n; l++){
					a[ll][l] -= a[icol][l]*dum;
				}
				for(int l=0; l<m; l++){
					b[ll][l] -= b[icol][l]*dum;
				}
			}
		}
	}
	// This is the end of the main loop over columns of the reduction. It only remains to unscramble the solution in view 
	// of the column interchanges. We do this by interchanging pairs of columns in the reverse order that the permutation was built up.
	for (int l=n-1; l>=0; l--) {
		if (indxr[l] != indxc[l]){
			for (int k=0; k<n; k++){
				SWAP2(a[k][indxr[l]],a[k][indxc[l]]);
			}
		}
	}	
	
	Free(indxc, sizeof(int) * n);
	Free(indxr, sizeof(int) * n);
	Free(ipiv, sizeof(int) * n);
	
	return 0 ;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ****************
// * Physic tools *
// ****************


double LCTools::Hz(double z, double H0, double Omm, double Oml, double w)
{
	return( H0 * sqrt(Omm*pow(1+z,3.) + Oml*pow(1+z,3.*w)) );
}


double LCTools::DL(double z, double H0, double Omm, double Oml, double w)
{
	double dz(0.001), zp(0.), Itg(0.);
	while (zp<z) {
		Itg += dz/Hz(zp+dz/2., H0, Omm, Oml, w);
		zp += dz;
	}
	//*o << z << "  " << (1.+z) * Itg * LC::c_SI << Endl; 
	//return ( (1.+z) * Itg * (LC::c_SI*1.e-3) * 1.e3 );
	return ( (1.+z) * Itg * LC::c_SI );
}


double LCTools::redshift(double DLs, double H0, double Omm, double Oml, double w)
{
	double dz(0.01), zTry(0.);
	
	while(DLs > DL(zTry, H0, Omm, Oml, w))
		zTry += dz;
	return(zTry);
}






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **************
// * Word tools *
// **************


/**
 * 
 */
char * LCTools::uppercase(const char * motmin) const
{
	int longueur,i ;
	longueur = strlen(motmin) ;
	static char MOTMAJ[512] ;
	strcpy(MOTMAJ,motmin) ; 
	for (i=0 ; i<longueur ; i++){
		majuscule(MOTMAJ[i]);
	}
	return(MOTMAJ);
}

/**
 * 
 */
char * LCTools::lowercase(const char * MOTMAJ) const
{
	int longueur,i ;
	longueur = strlen(MOTMAJ) ;
	static char motmin[512] ;
	strcpy(motmin,MOTMAJ) ; 
	for (i=0 ; i<longueur ; i++){
		minuscule(motmin[i]);
	}
	return(motmin);
}


/**
 * 
 */
void LCTools::majuscule(char &lettre) const
{
	if ((lettre>='a')&&(lettre<='z'))
		lettre=lettre-32;
}


/**
 * 
 */
void LCTools::minuscule(char &lettre) const
{
	if ((lettre>='A')&&(lettre<='Z'))
		lettre=lettre+32;
}


/**
 * 
 */
bool LCTools::wcmp(const char * mot1, const char * mot2) const
{
	char motmin1[512];
	char motmin2[512];
	strcpy(motmin1,lowercase(mot1));
	strcpy(motmin2,lowercase(mot2));
	//(*o) << "Compare " << mot1 << " and " << mot2 << " : " << strcmp(motmin1,motmin2) << " : " ;
	if(strcmp(motmin1,motmin2)==0){
		//(*o) << lowercase(mot1) << " == " << lowercase(mot2) << Endl;
		return(true);
	}else{
		//(*o) << lowercase(mot1) << " =/= " << lowercase(mot2) << Endl;
		return(false);
	}
}


/**
 * 
 */
void LCTools::wextract(char * InString, std::vector<char*> &Words)
{
	char * ptrB;
	int Nchar;
	int iPos(0);
	
	// ** Initialization of vector of words
	for(int iW=0; iW<Words.size(); iW++)
		Free(Words[iW], 256*sizeof(char));
	Words.resize(0);
	
	// ** Read first spaces
	while((InString[iPos]==' ')||(InString[iPos]=='\t'))
		iPos++;
	
	while((InString[iPos]!='\0')&&(InString[iPos]!='\n')){
		
		ptrB = NULL;
		// ** Begin of the word
		ptrB = &InString[iPos];
		
		// ** Read the word
		Nchar = 0;
		while((InString[iPos]!=' ')&&(InString[iPos]!='\t')&&(InString[iPos]!='\0')&&(InString[iPos]!='\n')){
			iPos++;
			Nchar++;
		};
		
		Words.push_back(NULL);
		Words[Words.size()-1] = (char*) AllocMemory(256*sizeof(char));
		strncpy(Words[Words.size()-1], ptrB, Nchar);
		Words[Words.size()-1][Nchar] = '\0';
		
		//(*o) << "=======> '" << Words[Words.size()-1] << "'" << Endl;
		
		// ** Read spaces between word or before the end
		while((InString[iPos]==' ')||(InString[iPos]=='\t'))
			iPos++;
	};
	
}


void LCTools::stripcopy(const char *orig, char * &ret ) 
{
    int pos = 0, len = 0;
    
    // Strip space-like characters at the beginning
    while(orig[pos] == ' ' || orig[pos] == '\n' || orig[pos] == '\r'|| orig[pos] == '\t')
        pos++;
	
    // Walk until the end of the string  
    while(orig[pos+len] != 0)
        len++;
	
    // Strip space-like characters at the end
    len--;
    while(orig[pos+len] == ' ' || orig[pos+len] == '\n' || orig[pos+len] == '\r' || orig[pos+len] == '\t')
		len--;
    len++;
	
    // Copy the string
	if(ret != NULL)
		Free(ret, (strlen(ret)+1) * sizeof(char));
	
    ret = (char*)AllocMemory( (len+1) * sizeof(char) );
    for(int i=0;i<len;i++) 
		ret[i] = orig[pos+i];
    ret[len] = '\0';
    
	//std::cout << strlen(ret) << " / " << len << Endl;
	
}

void LCTools::wextractcoma(const char * InString, int & ipIS, char ** & Words, int & NWords, int NCharWords)
{
	int ipL;
	bool Stop(false);
	while (!Stop) {
		
		//! * Add memory in the list of words
		NWords++;
		Words = (char**) ReAllocMemory(Words, (NWords-1)*sizeof(char*), NWords*sizeof(char*));
		//Cout << "NfOutObsNames[" << NfOut-1 << "] = " << NOutName << Endl;
		
		//! * Allocate memory for the word 
		Words[NWords-1] = (char*) AllocMemory(NCharWords*sizeof(char));
		
		//! * Copy the word
		ipL = 0;
		while ((InString[ipIS]!=',')&&(InString[ipIS]!='\0')&&(InString[ipIS]!='\n')&&(InString[ipIS]!=' ')&&(ipL<NCharWords)){
			Words[NWords-1][ipL++] = InString[ipIS++];  
		}
		//! * Finish the word
		Words[NWords-1][ipL] = '\0';
		
		//Cout << "fOutObsNames[" << NfOut-1 << "][" << NOutName-1 << "] = " << OutName[NOutName-1] << Endl;
		
		//! * Stop reding condition
		if(InString[ipIS]!=',')
			Stop = true;
		else
			ipIS++;
		
		
	}
}


bool LCTools::isValue(const char * Str)
{
	std::stringstream iss( Str );
	float tmp;
	return ( iss>>tmp ) && ( iss.eof() );
}



void LCTools::pathExtractName(char * Path, char * & Name)
{
	int iP(strlen(Path));
	while ( (iP>0) && (Path[iP]!='/') ) {
		iP--;
	}
	if(iP==0)
		iP=-1;
	strncpy(Name, Path+iP+1, strlen(Path)-iP);
	Name[strlen(Path)-iP-1] = '\0';
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ****************************
// * File and System Encoding *
// ****************************

TypeEncoding LCTools::testbyteorder()
{
	short int word = 0x0001;
	char *byte = (char *) &word;
	return(byte[0] ? LITTLEENDIAN : BIGENDIAN);
}


void LCTools::convertendianness(double *val) {
	unsigned char cval[8];
	
	cval[0] = ((unsigned char *)val)[7];
	cval[1] = ((unsigned char *)val)[6];
	cval[2] = ((unsigned char *)val)[5];
	cval[3] = ((unsigned char *)val)[4];
	cval[4] = ((unsigned char *)val)[3];    
	cval[5] = ((unsigned char *)val)[2];
	cval[6] = ((unsigned char *)val)[1];
	cval[7] = ((unsigned char *)val)[0];
	
	*val = *((double *)cval);
}


void LCTools::DispEncoding(TypeEncoding Encoding, std::ostream *  OutDisp)
{
	if(Disp()){
		switch (Encoding) {
			case BIGENDIAN :
				(*o) << "BigEndian";
				break;
			case LITTLEENDIAN :
				(*o) << "LittleEndian";
				break;
			default :
				(*o) << "Undefined";
		}
	}
}


TypeEncoding LCTools::ReadEncoding(char * Str)
{
	if(wcmp(Str,"Binary,LittleEndian"))
		return LITTLEENDIAN ;
	if(wcmp(Str,"Binary,BigEndian"))
		return LITTLEENDIAN ;
	return NOENDIAN;
}


TypeEncoding LCTools::ReadEncodingConst(const char * Str)
{
	char tmpC[512];
	strcpy(tmpC, Str);
	return ReadEncoding(tmpC);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **********
// * Files  *
// **********
void LCTools::DispFileType(TypeFile FType, std::ostream * OutDisp)
{
	if(Disp()){
		switch (FType) {
			case ASCII :
				(*OutDisp) << "ASCII";
				break;
			case BINARY :
				(*OutDisp) << "Binary";
				break;
			case XML :
				(*OutDisp) << "XML";
				break;
			default :
				(*OutDisp) << "Undefined";
		}
	}
}


void LCTools::CheckFile(char * fDir, char * fName)
{
	std::ifstream fIn;
	
	char TmpFullName[16384];
	if(wcmp(fDir,""))
		strcpy(TmpFullName, fName);
	else 
		sprintf(TmpFullName, "%s/%s", fDir, fName);
	
	fIn.open(TmpFullName);
	if(fIn == NULL){
		(*o) << Endl << "ERROR: Can not open the file " << TmpFullName << " !" << Endl;
		throw std::invalid_argument("ERROR: Can not open a file.");
	}
	fIn.close();
	fIn.clear();
	(*o) << " --> OK" << Endl;
}



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// *************
// * XML tools *
// *************

char * LCTools::gXMLUnit(const char In[], double & Fact){
	char * EndPtr;
	double res;
	res = strtod(In,&EndPtr);
	if(EndPtr==In) 
		res=1;
	Fact = res;
	return(EndPtr);
}


double LCTools::gXMLAngle(ezxml_t param)
{
	double AngleParam(1.0);
	char * UnitType;
	UnitType = gXMLUnit(ezxml_attr(param,"Unit"), AngleParam);
	AngleParam *= atof(ezxml_txt(param));
	if(wcmp(UnitType,"Degree")){
		AngleParam = AngleParam*M_PI/180.;
	}
	if((!wcmp(UnitType,"Degree"))&&(!wcmp(UnitType,"Radian")))
		throw std::invalid_argument("ERROR in LCTools::gXMLAngle : The angle's unit isn't known (only Degree or Radian) !");
	return (AngleParam);
}


double LCTools::gXMLTime(ezxml_t param)
{
	double TimeParam(1.0);
	char * UnitType;
	UnitType = gXMLUnit(ezxml_attr(param,"Unit"), TimeParam);
	if(wcmp(UnitType,"Meter")){
		TimeParam *= 1./LC::c_SI;
	}
	if(wcmp(UnitType,"KM")){
		TimeParam *= 1.e3/LC::c_SI;
	}
	if(wcmp(UnitType,"Day")){
		TimeParam *= LC::Dy_SI;
	}
	if(wcmp(UnitType,"Year")){
		TimeParam *= LC::Yr_SI;
	}
	if( (!wcmp(UnitType,"Second")) && (!wcmp(UnitType,"Meter")) && (!wcmp(UnitType,"KM")) && (!wcmp(UnitType,"Day")) && (!wcmp(UnitType,"Year")))
		throw std::invalid_argument("ERROR in LCTools::gXMLTime : The time's unit isn't known (only Second, Meter, KM, Day or Year)  !");
	TimeParam *= atof(ezxml_txt(param));
	return (TimeParam);
}


double LCTools::gXMLLength(ezxml_t param)
{
	double LenParam(1.0);
	char * UnitType;
	UnitType = gXMLUnit(ezxml_attr(param,"Unit"), LenParam);
	LenParam *= atof(ezxml_txt(param));
	if(wcmp(UnitType,"KM")){
		LenParam *= 1.e3;
	}
	if(wcmp(UnitType,"Second")){
		LenParam *= LC::c_SI;
	}
	if((!wcmp(UnitType,"Meter"))&&(!wcmp(UnitType,"KM"))&&(!wcmp(UnitType,"Second")))
		throw std::invalid_argument("ERROR in LCTools::gXMLLength : The length's unit isn't known (only Meter, KM or Second) !");
	return (LenParam);
}


double LCTools::gXMLFrequency(ezxml_t param)
{
	double FrequencyParam(1.0);
	char * UnitType;
	UnitType = gXMLUnit(ezxml_attr(param,"Unit"), FrequencyParam);
	FrequencyParam *= atof(ezxml_txt(param));
	if(strcmp(UnitType,"MilliHertz")==0){
		FrequencyParam *= 1.0e-3;
	}
	if((strcmp(UnitType,"MilliHertz")!=0)&&(strcmp(UnitType,"Hertz")!=0))
		throw std::invalid_argument("ERROR in LCTools::gXMLFrequency : The frequency's unit isn't known (only MilliHertz or Hertz) !");
	return (FrequencyParam);
}


double LCTools::gXMLAstroMass(ezxml_t param)
{
	double AstroMass(1.0);
	char * UnitType;
	UnitType = gXMLUnit(ezxml_attr(param,"Unit"), AstroMass);
	AstroMass *= atof(ezxml_txt(param));
	if(strcmp(UnitType,"SolarMass")!=0)
		throw std::invalid_argument("ERROR in LCTools::gXMLAstroMass : The astronomic mass's unit isn't known (only SolarMass) !");
	return (AstroMass);
}



double LCTools::gXMLAstroDistance(ezxml_t param)
{
	double AstroDistance;
	char * UnitType;
	UnitType = gXMLUnit(ezxml_attr(param,"Unit"), AstroDistance);
	AstroDistance *= atof(ezxml_txt(param));
	if(strcmp(ezxml_attr(param,"Unit"),"Parsec")==0){
		AstroDistance *= 1.0e-3;
	}
	if((strcmp(UnitType,"Parsec")!=0)&&(strcmp(UnitType,"KiloParsec")!=0))
		throw std::invalid_argument("ERROR in LCTools::gXMLAstroDistance : The astronomic distance's unit isn't known (only Parsec or KiloParsec) !");
	return (AstroDistance);
}



char * LCTools::TimeISO8601()
{
	time_t rawtime;
	static char strTime[22];
	
	time ( &rawtime );
	
	strcpy(strTime,ctime(&rawtime));
	
	// * Year
	strTime[0] = ctime(&rawtime)[20];
	strTime[1] = ctime(&rawtime)[21];
	strTime[2] = ctime(&rawtime)[22];
	strTime[3] = ctime(&rawtime)[23];
	strTime[4] = '-';
	
	// * Month
	strTime[5] = '0';
	if(strncmp(ctime(&rawtime)+4,"Jan",3) == 0)
		strTime[6] = '1';
	if(strncmp(ctime(&rawtime)+4,"Feb",3) == 0)
		strTime[6] = '2';
	if(strncmp(ctime(&rawtime)+4,"Mar",3) == 0)
		strTime[6] = '3';
	if(strncmp(ctime(&rawtime)+4,"Apr",3) == 0)
		strTime[6] = '4';
	if(strncmp(ctime(&rawtime)+4,"May",3) == 0)
		strTime[6] = '5';
	if(strncmp(ctime(&rawtime)+4,"Jun",3) == 0)
		strTime[6] = '6';
	if(strncmp(ctime(&rawtime)+4,"Jul",3) == 0)
		strTime[6] = '7';
	if(strncmp(ctime(&rawtime)+4,"Aug",3) == 0)
		strTime[6] = '8';
	if(strncmp(ctime(&rawtime)+4,"Sep",3) == 0)
		strTime[6] = '9';
	if(strncmp(ctime(&rawtime)+4,"Oct",3) == 0){
		strTime[5] = '1';
		strTime[6] = '0';
	}
	if(strncmp(ctime(&rawtime)+4,"Nov",3) == 0){
		strTime[5] = '0';
		strTime[6] = '1';
	}
	if(strncmp(ctime(&rawtime)+4,"Dec",3) == 0){
		strTime[5] = '0';
		strTime[6] = '2';
	}
	
	// * Day
	strTime[7] = '-';
	if(ctime(&rawtime)[8] == ' ')
		strTime[8] = '0';
	
	strTime[10] = 'T';
	
	strTime[19] = '\0'; 
	
	//cout << ctime(&rawtime) << Endl;
	
	return(strTime);
	
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **************************
// * Lagrange interpolation *
// **************************


void LCTools::InitLagPolDen(int order)
{
	LagPolDen[order] = (double*)AllocMemory((order+1)*sizeof(double));
	for(int k=0; k<=order; k++){
		LagPolDen[order][k] = 1.0;
		for(int j=0; j<=order; j++){
			if(j!=k)
				LagPolDen[order][k] *= 1.0/(k-j);
		}
	}
}





double LCTools::InterLagrange(double * data, int Nmax, double xr, int order)
{
	//Cout << "////\\\\LCSerie2::InterLagrange : ordre = " << order << "  et x = " << x << Endl;
    int bin;
    
    
	bin = ifloor(xr);
	
	double res(0.0), Pk(0.0);
	int ordermin(ifloor(double(order+1)/2.0));
	//Cout << int(ys.size())-1 << Endl;
	if((bin<0)||(bin+1 > Nmax-1)){
		std::cerr << "ERROR in LCSerie2::InterLagrange : bin = " << bin << " (x = " << xr << ") is not included in [0," << Nmax-1 << "]." << Endl;  
		throw std::invalid_argument("LCSerie2::InterLagrange : The required bin does not exist !");  
	}
	//    int kmin(bin-ordermin+1), kmax(bin+(order+1-ordermin));  // commented by Sofiane
    int kmin(bin-ordermin-1), kmax(bin+(order+1-ordermin));  // changed by Sofiane
	if(kmin < 0){
		std::cerr << "ERROR in LCSerie2::InterLagrange : For bin = " << bin << " (x = " << xr << "), kmin = bin-ordermin+1 < 0 " << Endl;
		throw std::invalid_argument("LCSerie2::InterLagrange : The required bin does not exist !");
	}
	if(kmax > Nmax-1){
		std::cerr << "ERROR in LCSerie2::InterLagrange : For bin = " << bin << " (x = " << xr << "), kmax = order+1-ordermin > " << Nmax-1 << Endl;
		throw std::invalid_argument("LCSerie2::InterLagrange : The required bin does not exist !");
	}
	//Cout << "  InterLagrange : x0 = " << x0 << " , x = " << x;
	//Cout << " , bin = " << bin << " , kmin = " << kmin << " , kmax = " << kmax;
	
	/*
	 // Old Lagrange
	 for(int k=kmin; k<=kmax; k++){
	 Pk = 1.0;
	 for(int j=kmin; j<=kmax; j++){
	 if(j!=k)
	 Pk *= (x-x0-j*dx)/((k-j)*dx);
	 }
	 res += (*ys[k])*Pk;
	 }
	 */
	
	// ** New Lagrange with tabulation of polynome denominator
	int RealOrder (kmax - kmin);
	if(LagPolDen[RealOrder] == NULL)
		InitLagPolDen(RealOrder);
	
	for(int k=kmin; k<=kmax; k++){
		Pk = LagPolDen[RealOrder][k-kmin];
		for(int j=kmin; j<=kmax; j++){
			if(j!=k)
				Pk *= xr-j ;
		}
		res += (data[k])*Pk;
	}
	
	
	//Cout << "  , res = " << res << Endl;
	//	xrbinE = false;
	return(res);
}





//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **************
// * FFTW tools *
// **************

// * Forward
// *********

void LCTools::FTPrepFwdPlan(int NtDat_n)
{
	//(*o) << "+++++FTFwdtimeData = " <<  FTFwdtimeData << Endl;
	//(*o) << "+++++FTFwdfreqData = " <<  FTFwdfreqData << Endl;
	
	if(FTFwdtimeData != NULL)
		Free(FTFwdtimeData, FTFwdNbtDat*sizeof(double));
	if(FTFwdfreqData != NULL)
		Free(FTFwdfreqData, FTFwdNbfDat* sizeof(dcomplex));
	if (FTFwdPlan != NULL){
		//(*o) << "DESTROY PLAN" << Endl;
		fftw_destroy_plan(FTFwdPlan);
	}
	
	
	FTFwdNbtDat = NtDat_n;
	FTFwdNbfDat = getNfFTreal(FTFwdNbtDat);
	//(*o) << "FTFwdNbtDat = " << FTFwdNbtDat <<  Endl;
	//(*o) << "FTFwdNbfDat = " << FTFwdNbfDat <<  Endl;
	FTFwdtimeData = (double*) AllocMemory( FTFwdNbtDat * sizeof(double) );
	//(*o) << "+++++FTFwdtimeData = " <<  FTFwdtimeData << Endl;
	FTFwdfreqData = (dcomplex*) AllocMemory( FTFwdNbfDat * sizeof(dcomplex) );
	//(*o) << "+++++FTFwdfreqData = " <<  FTFwdfreqData << Endl;
	if(FTFwdBestFFTW){
		if(Disp())
			(*o) << " Find the best FFTW forward algorithm ..." ; o->flush();
		FTFwdPlan = fftw_plan_dft_r2c_1d(FTFwdNbtDat, FTFwdtimeData, reinterpret_cast<fftw_complex*> (FTFwdfreqData), FFTW_MEASURE); // Prepare plan FFTW plan
	}else{
		if(Disp())
			(*o) << " Find estimated FFTW forward algorithm ..." ; o->flush();
		FTFwdPlan = fftw_plan_dft_r2c_1d(FTFwdNbtDat, FTFwdtimeData, reinterpret_cast<fftw_complex*> (FTFwdfreqData), FFTW_ESTIMATE); // Prepare plan FFTW plan
	}
	if(Disp())
		(*o) << " OK ... " ; o->flush();
}

void LCTools::FTMakeFwd(double * tData, dcomplex * &fData, int NtDat_n)
{
	if(NtDat_n != FTFwdNbtDat)
		FTPrepFwdPlan(NtDat_n);
	//! ** Copy signal in FFT input memory
	memcpy(FTFwdtimeData, tData, FTFwdNbtDat*sizeof(double));
	//! ** Compute FFT
	fftw_execute(FTFwdPlan);
	//! ** Copy FFT ouput memory in result
	memcpy(fData, FTFwdfreqData, FTFwdNbfDat*sizeof(dcomplex));	
}

void LCTools::FTMakeFwdMulti(double ** tData, dcomplex ** &fData, int NtDat_n, int NSig)
{
	//DispFTMemLocation();
	if(NtDat_n != FTFwdNbtDat){
		FTPrepFwdPlan(NtDat_n);
	}
	for(int iS=0; iS<NSig; iS++){
		memcpy(FTFwdtimeData, tData[iS], FTFwdNbtDat*sizeof(double));		// Copy signal in FFT input memory
		fftw_execute(FTFwdPlan);	
		memcpy(fData[iS], FTFwdfreqData, FTFwdNbfDat*sizeof(dcomplex));	// Copy FFT ouput memory in result		
	}
}



// * Backward
// ***********

void LCTools::FTPrepBckPlan(int NfDat_n)
{
	//(*o) << "+++++FTBcktimeData = " <<  FTBcktimeData << Endl;
	//(*o) << "+++++FTBckfreqData = " <<  FTBckfreqData << Endl;
	
	if(FTBckfreqData != NULL)
		Free(FTBckfreqData, FTBckNbfDat* sizeof(dcomplex));
	if(FTBcktimeData != NULL)
		Free(FTBcktimeData, FTBckNbtDat*sizeof(double));
	if (FTBckPlan != NULL){
		//(*o) << "DESTROY PLAN" << Endl;
		fftw_destroy_plan(FTBckPlan);
	}
	
	
	FTBckNbfDat = NfDat_n;
	FTBckNbtDat = getNtFTreal(FTBckNbfDat);
	//(*o) << "FTBckNbtDat = " << FTBckNbtDat <<  Endl;
	//(*o) << "FTBckNbfDat = " << FTBckNbfDat <<  Endl;
	FTBckfreqData = (dcomplex*) AllocMemory( FTBckNbfDat * sizeof(dcomplex) );
	//(*o) << "+++++FTBckfreqData = " <<  FTBckfreqData << Endl;
	FTBcktimeData = (double*) AllocMemory( FTBckNbtDat * sizeof(double) );
	//(*o) << "+++++FTBcktimeData = " <<  FTBcktimeData << Endl;
	if(FTBckBestFFTW){
		if(Disp())
			(*o) << " Find the best FFTW backward algorithm ..." ; o->flush();
		FTBckPlan = fftw_plan_dft_c2r_1d(FTBckNbtDat, reinterpret_cast<fftw_complex*> (FTBckfreqData),  FTBcktimeData, FFTW_MEASURE); // Prepare plan FFTW plan
	}else{
		if(Disp())
			(*o) << " Find estimated FFTW backward algorithm ..." ; o->flush();
		FTBckPlan = fftw_plan_dft_c2r_1d(FTBckNbtDat, reinterpret_cast<fftw_complex*> (FTBckfreqData), FTBcktimeData, FFTW_ESTIMATE); // Prepare plan FFTW plan
	}
	if(Disp())
		(*o) << " OK ... " ; o->flush();
}

void LCTools::FTMakeBck(dcomplex * fData , double * &tData, int NfDat_n)
{
	if(NfDat_n != FTBckNbfDat)
		FTPrepBckPlan(NfDat_n);
	memcpy(FTBckfreqData, fData, FTBckNbfDat*sizeof(dcomplex));	// Copy FFT ouput memory in result
	fftw_execute(FTBckPlan);							// Compute FFT
	memcpy(tData, FTBcktimeData,  FTBckNbtDat*sizeof(double));		// Copy signal in FFT input memory
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// **************************
// * Power Spectral Density *
// **************************



void LCTools::PSD(double * tDat, int NtDat, double dt, int TypeConvWin, int NSeg, double * & PSDDat, int & NDatPSD, double & df, bool AllocPSD)
{
	int m, NtDatS;
	double * tDatS;
	dcomplex * fDatS;
	double facm, facp, sumw;
	int iBin0;
	double den;
	double tmpA;
	
	
	//! ***** Prepare the PSD
	
	//! *** Number of non overlapping data per segment
	m = ifloor(NtDat/((double)(NSeg+1)));
	
	//! *** Number of data per segment
	NtDatS = m+m;
	NDatPSD = getNfFTreal(NtDatS);
	
	//! *** Allocation
	tDatS = (double*) AllocMemory(NtDatS*sizeof(double)) ;
	fDatS = (dcomplex*) AllocMemory(NDatPSD*sizeof(dcomplex));
	if(AllocPSD)
		PSDDat = (double*)AllocMemory(NDatPSD*sizeof(double));
	
	//! *** Initialization of psd
	for(int i=0; i<NDatPSD; i++)
		PSDDat[i] = 0.;
	
	
	//! *** Prepare the convolution window : compute the sum for normalizing the windowing
	sumw = 0.;
	den = 0.;
	facm = m;
	facp = 1.0/m;
	for (int i=0; i<NtDatS; i++)
		sumw += AppConvWin(i,facm,facp,TypeConvWin) * AppConvWin(i,facm,facp,TypeConvWin);
	sumw *= NtDatS;
	//(*o) << "  Sum in the window Wss = " << sumw << " ." << Endl;
	
	
	//! *****  Compute the PSD
	//! *** Loop on the segments
	for(int iS=0; iS<NSeg; iS++){
		
		//! ** Index of first bin of the segment 
		iBin0 = iS*m;
		
		//! ** Copy the segment of data 
		memcpy(tDatS, tDat+iBin0, NtDatS*sizeof(double));
		
		//! ** Apply the convolution window
		for (int i=0; i<NtDatS; i++)
			tDatS[i] *= AppConvWin(i,facm,facp,TypeConvWin);
		
		//! ** Compute fft
		FTMakeFwd(tDatS, fDatS, NtDatS);
		
		for(int i=0; i<NDatPSD; i++){
			tmpA= abs(fDatS[i]);
			PSDDat[i] += tmpA*tmpA  ;
		}
		den += sumw;
	}
	den *= 0.5;
	//(*o) << "  Denominator = " << den << Endl;
	
	//! *** Compute frequency step
	df = 1/(2.0*m*dt);
	
	//! *** Normalization
	for(int i=0; i<NDatPSD; i++)
		PSDDat[i] /= den*df;
	
	//! *** Deallocation
	Free(tDatS, NtDatS*sizeof(double));
	Free(fDatS, NDatPSD*sizeof(dcomplex));
}


double LCTools::AppConvWin(int j, double a, double b, int TypeConvWin)
{
	switch (TypeConvWin){
		case 1 :
			return (1.0-fabs(((j)-(a))*(b)));
			break;
		case 2 :
			return (1.0);
			break;
		case 3 :
			return(1.0-(((j)-(a))*(b))*(((j)-(a))*(b)));
			break;
		case 4 :
			return(0.5*(1.0-(cos((j)*M_PI*(b)))));
			break;
		case 5 :
			return(sin((j)*M_PI*(b)/2.0));
			break;
		case 6 :
			return(0.54-0.46*(cos((j)*M_PI*(b))));
			break;
		default :
			throw std::invalid_argument("ERROR in LCTools::AppWindow : Unknown window !" );
	};
}
	   



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// *************
// * MPI tools *
// *************

// ********** MPI: Basic tools **********

/**
 * 
 */
void LCTools::MPIset(int mpi_id_n, int mpi_np_n)
{
	mpi_id = mpi_id_n;
	mpi_np = mpi_np_n;
}

/**
 * 
 */
bool LCTools::MPImaster()
{
	if(mpi_id==0)
		return true;
	else
		return false; 
}




// end of LISACODE-Tools.cpp
