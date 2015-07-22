/*
 *  LISACODE-GWGalBin.cpp
 *  LC20
 *
 *  Created by Antoine Petiteau on 06/05/11.
 *  Copyright 2011 Max-Planck-Institut für Gravitationsphysik - Albert Einstein Institute. All rights reserved.
 *
 */

#include "LISACODE-GWCosmicString.h"


// *********************
// ***  Constructor  ***
// *********************

LCGWCosmicString::LCGWCosmicString()
: LCGW()
{
	initNULL(false);
}


LCGWCosmicString::LCGWCosmicString(LCTools * MT_n)
: LCGW(MT_n)
{
	initNULL(false);
}


LCGWCosmicString::~LCGWCosmicString()
{
	initNULL(true);
}




// ********************************************
// ***        Required methods              ***
// ***  (to be present in derived classes)  ***
// ********************************************

void LCGWCosmicString::initNULL(bool CleanMem)
{
    //Cout << "LCGWCosmicString::initNULL ..." << Endl;
    
	initNULLBase(CleanMem);
	
    
	if (CleanMem){
        if(hp_time!=NULL)
            MT->Free(hp_time,N*sizeof(double));
        if(hp_freq!=NULL)
            MT->Free(hp_freq,Nf*sizeof(dcomplex));
        
        if(amplitudeInit!=NULL){
            MT->Free(amplitudeInit,Nb*sizeof(double));
            MT->Free(centralTime,Nb*sizeof(double));
            MT->Free(cos2Pol,Nb*sizeof(double));
            MT->Free(msin2Pol,Nb*sizeof(double));
        }
	}
    
	NParams = 12;
	
    
	
	//Gmu = -1.;
	//DL = -1.;
	//ParentSize = -1.;
	//alpha=-1;
	flow = 1.e-3;
	//tcosmo=0.;
    
	fhigh=-1;
    
	// !!!!!!!!! on peut peut etre definir le dt dirrecment ici pour avoir une bonne précision  ?
	
    
	
    
	fmin = 1.e-4;
	fmax = 1.e4;
	N=1;
	df=1;
    
    
	
	
	hp_time = NULL;
	hp_freq = NULL;
    
    Nb=0;
	amplitudeInit=NULL;
	centralTime=NULL;
	cos2Pol=NULL;
	msin2Pol=NULL;
    PolarizationDelta=M_PI; 
    
    strcpy(FileNameBurstList,"None");
    

	q=-1;
    strcpy(BurstType,"None");
    
	//Cout<<"initNULL"<< " " << NParams << Endl;
    
    
	/*ecc = 0.;
     DPeriod = 0.;
     
     
     DontChangeDFreq = false;
     
     hp0 = 0.;
     hc0 = 0.;
     om = 1.e-3;
     Dom = 0.;*/
    
	//Cout << "LCGWCosmicString::initNULL --> OK" << Endl;
	
}


// ***************************
// *  Configuration methods  *
// ***************************


void LCGWCosmicString::config(ezxml_t xmlbloc)
{
    double TMPamplitudeInit(0.), TMPcentralTime(0.);
    
    
    Cout << "Configuration of LCGWCosmicString : Read parameters : ";
    //! *** Read sky position, polarization and name
	configBase(xmlbloc);
    
	ezxml_t param;
	for(param = ezxml_child(xmlbloc,"Param"); param; param = param->next){
        
        ParameterReaded = false;
        
		//=======================================
        
		/*if(MT->wcmp(ezxml_attr(param,"Name"),"MassL"))
         setParam(3, atof((*param).txt));*/
        
        
        if(MT->wcmp(ezxml_attr(param,"Name"),"FileBusrtList")||(MT->wcmp(ezxml_attr(param,"Name"),"ListBusrtFile"))){
            char * tmpF(NULL);
            MT->stripcopy((*param).txt, tmpF); 
            strcpy(FileNameBurstList, tmpF);
            if(tmpF != NULL)
                MT->Free(tmpF, (strlen(tmpF)+1) * sizeof(char));
        }
        
        if(MT->wcmp(ezxml_attr(param,"Name"),"amplitudeInit"))
             TMPamplitudeInit = atof((*param).txt);
        
		if(MT->wcmp(ezxml_attr(param,"Name"),"centralTime"))
            TMPcentralTime = MT->gXMLTime(param);
        
        if(MT->wcmp(ezxml_attr(param,"Name"),"PolarisationDelta"))
            PolarizationDelta = MT->gXMLAngle(param);
        
		if(MT->wcmp(ezxml_attr(param,"Name"),"q"))
            setParam(5, atof((*param).txt));
        
		if(MT->wcmp(ezxml_attr(param,"Name"),"flow"))
            setParam(6, atof((*param).txt));
		
		if(MT->wcmp(ezxml_attr(param,"Name"),"BurstType")){
            char * BurstType(NULL);
            MT->stripcopy((*param).txt, BurstType);
            if(MT->wcmp(BurstType,"cusp")){
                setParam(5,1./3.);
            }else{
                if(MT->wcmp(BurstType,"kink")){
                    setParam(5,2./3.);
                }else{
                    Cout << "ERROR in LCGWCosmicString::config : The burst type >" << (*param).txt << "< is unknown !" << Endl;
                    throw std::invalid_argument("ERROR in LCGWCosmicString::config : The burst type is unknown !");
                    break;
                }
            }
            if(BurstType != NULL)
                MT->Free(BurstType, (strlen(BurstType)+1) * sizeof(char));
            
		}
        
        
        
		/*if(MT->wcmp(ezxml_attr(param,"Name"),"Gamma"))
         setParam(5, atof((*param).txt));*/
        
		/*if(MT->wcmp(ezxml_attr(param,"Name"),"epsilon"))
         setParam(7 ,atof((*param).txt));*/
        
        
		//if(MT->wcmp(ezxml_attr(param,"Name"),"fhigh"))
        //    setParam(7, atof((*param).txt));
        
		//Cout <<param<< Endl ;
		
		//=======================================
        
        //if(ParameterReaded)
        //    Cout << " " << ezxml_attr(param,"Name") ; 
		
	}
    
    if(MT->wcmp(FileNameBurstList,"None")){
        if(amplitudeInit!=NULL){
            MT->Free(amplitudeInit,Nb*sizeof(double));
            MT->Free(centralTime,Nb*sizeof(double));
            MT->Free(cos2Pol,Nb*sizeof(double));
            MT->Free(msin2Pol,Nb*sizeof(double));
        }
        Nb=1;
        amplitudeInit = (double*)MT->AllocMemory(Nb*sizeof(double));
        centralTime   = (double*)MT->AllocMemory(Nb*sizeof(double));
        cos2Pol       = (double*)MT->AllocMemory(Nb*sizeof(double));
        msin2Pol      = (double*)MT->AllocMemory(Nb*sizeof(double));
        amplitudeInit[0] = TMPamplitudeInit;
        centralTime[0] = TMPcentralTime;
        cos2Pol[0]  = cos(2*Polarization);
        msin2Pol[0] = -sin(2*Polarization);
    }
    
	//Cout << Endl;
}


void LCGWCosmicString::config(int iParam, double ParamVal)
{
	
}


// ********************
// *  Access methods  *
// ********************


double LCGWCosmicString::getParam(int iP)
{

    if((Nb==0)&&((iP==2)||(iP==3)||(iP==4))){
        Cout << "ERROR: LCGWCosmicString::getParam : there is no burst in the list." << Endl; 
        throw std::invalid_argument("ERROR: LCGWCosmicString::getParam : there is no burst in the list.");
    }
    
    if((Nb>1)&&((iP==2)||(iP==3)||(iP==4))){
        Cout << "WARNING: LCGWCosmicString::getParam : there is more than one burst in the list." << Endl; 
        std::cerr << "WARNING: LCGWCosmicString::getParam : there is more than one burst in the list." << Endl;
    }
    
	switch (iP) {
		case 0:
			return(Beta);
			break;
            
		case 1:
			return(Lambda);
			break;
            
		case 2:
			return(acos(cos2Pol[0])/2.);
			break;
            
		case 3:
			return(amplitudeInit[0]);
			break;
            
		case 4:
			return(centralTime[0]);
			break;
            
        case 5:
            return(q);
            break;
            
        case 6:
            return(flow);
            break;
            
        
        //case 7:
        //    return(fhigh);
        //    break;
        
			
		default:
			Cout << "ERROR in LCGWCosmicString::getParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWCosmicString::getParam : The parameter is unknown !");
			break;
	}
	//Cout <<"getparam"<< " " <<  Gmu<<  Endl;
	return(0.);
}


void LCGWCosmicString::setParam(int iP, double Param_n)
{
    ParameterReaded = true;
    
    if((Nb==0)&&((iP==2)||(iP==3)||(iP==4))){
        Cout << "ERROR: LCGWCosmicString::setParam : there is no burst in the list." << Endl; 
        throw std::invalid_argument("ERROR: LCGWCosmicString::getParam : there is no burst in the list.");
    }
    
    if((Nb>1)&&((iP==2)||(iP==3)||(iP==4))){
        Cout << "WARNING: LCGWCosmicString::setParam : there is more than one burst in the list." << Endl; 
        std::cerr << "WARNING: LCGWCosmicString::getParam : there is more than one burst in the list." << Endl;
    }
    
    switch (iP) {
        case 0:
			Beta = Param_n;
			break;
		case 1:
			Lambda = Param_n;
			break;
            
        case 2:
            cos2Pol[0]  = cos(2*Param_n);
            msin2Pol[0] = -sin(2*Param_n);
            break;
            
		case 3:
			amplitudeInit[0] = Param_n;
			break;
            
		case 4:
			centralTime[0] = Param_n;
			break;
            
        case 5:
            q=Param_n;
            if(MT->deq(q,1./3.))
                strcpy(BurstType,"Cusp");
            if(MT->deq(q,2./3.))
                strcpy(BurstType,"Kink");
            //q=4./3.;
            break;
            
        case 6:
			flow = Param_n;
			break;
            
		//case 7:
        //    fhigh= Param_n;
        //    break;
            
		default:
			Cout << "ERROR in LCGWCosmicString::setParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWCosmicString::setParam : The parameter is unknown !");
			break;
            
	}
	//Cout <<"setparam"<< " " << Param_n <<  Endl;	
    
}


void LCGWCosmicString::getRange(int iP, double &Pmin, double &Pmax)
{
	switch (iP) {
		case 0:
			Pmin = -M_PI/2.0;
			Pmax = M_PI/2.0;
			break;
		case 1:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
		case 2:
			Pmin = 0.0;
			Pmax = 2.0*M_PI;
			break;
            
            
        case 3:
			Pmin = 1e-30;
			Pmax = 1e-15;
			break;
			/*case 9: 
             Pmin = 5e3;
             Pmax = 5e4;
             break;
             //en al avec l=alpha*t=epsilon*gamma*G*mu*t*1/c**2 ; où : epsilon =1, gamma = 50, tcosmo=4.3e17s => 1e4 al*/
            
		case 4:
			Pmin = 0;
			Pmax = 2*31556926;//2 years : tObsTot
			break;
		case 5:
            Pmin = 1./3;
            Pmax = 2./3;
            
		case 6:
			Pmin = 5e-4;
			Pmax = 5e-5;
            
			/*case 7:
             Pmin = 1e-12; //electroweak phase transition 1e-35 pour grande unification
             Pmax =1e-6; //quark-hadron transition 
             Cout <<"getRange"<< " " << Pmin << " " << Pmax << " " <<"tcosmo"<< Endl;
             
             
             /*case 8: 
             Pmin = 0.;
             Pmax = 1.;*/
            
		case 7: 
            Pmin = 0;
            Pmax = 1e10;
            
            
			
            
		default:
			Cout << "ERROR in LCGWCosmicString::setParam : The parameter " << iP << " is unknown !" << Endl;
			throw std::invalid_argument("ERROR in LCGWCosmicString::setParam : The parameter is unknown !");
            
            break;
	}
    
}


double LCGWCosmicString::getDelta(int iP)
{
	double tmpMin, tmpMax, DRes;
	getRange(iP, tmpMin, tmpMax);
	
	switch (iP){
		case 0:
			DRes = 2.e-6;
			break;
		case 1:
			DRes = 2.e-6;
			break;
		case 2:
			DRes = 1.e-5;
			break;
        	
            
		default:
			DRes = (tmpMax-tmpMin)*5.0e-8;
	};
	
	//Cout << iP << " --> Delta = " << DRes<< Endl;
	//Cout <<"getDelta"<< Endl;
	return (DRes);
}


void LCGWCosmicString::setSpecialParam(int iPS, double SpecParam_n)
{
	
}


double LCGWCosmicString::getSpecialParam(int iPS)
{
	return(0.);
}


void LCGWCosmicString::setTimeInfo(double t0_n, double dt_n, double TObs_n, double tAskMin_n, double tDOrbMax_n,  double tMaxDiff)
{	
    t0Real = t0_n; 
	dtReal = dt_n;
	TObsReal = TObs_n;
	tAskMin = tAskMin_n;
	tDOrbMax = tDOrbMax_n;
	Cout <<"setTimeInfo"<< " " <<TObsReal<<  Endl;
}


void LCGWCosmicString::RandParam(int iP)
{
	double Pmin(0.0);
	double Pmax(1.0);
	
	getRange(iP, Pmin, Pmax);
	setParam(iP,MT->RandUniform(Pmin, Pmax));
	
	//! ** Specifical case : non-uniform
	if(iP == 0)
		setParam(0, M_PI/2.0 - acos(MT->RandUniform(-1.0, 1.0)));
	//if(iP == 3)
	//setParam(3, acos(MT->RandUniform(-1.0, 1.0)));
	
	/*if(iP == 4)
     DontChangeAmp = false;
     if(iP == 6)
     DontChangeDFreq = false;*/
}

// ***************************************
// * Linking and initialization methods  *
// ***************************************

// !!!!!!!!!!!!!!!!!!!!!! a modif : mais je capte pas trop pour le moment

int LCGWCosmicString::init()
{
	initBase();
	Cout <<"Init GW CosmicString ..."<< Endl;
#ifdef _DEBUG_GW_    
	char fNCheck[512];
	int iRCheck(MT->ifloor(MT->RandUniform(0, 10000.)));
	sprintf(fNCheck,"CheckGWCosmicString_f_%d.txt",iRCheck);
	std::cerr << "DEBUG:GW:LCGWCosmicString : File = " << fNCheck << Endl;
	DEBUGfCheck = new std::ofstream(fNCheck);
    (*DEBUGfCheck) << "#f, hp(f), |hp(f)|, arg(hp(f))" << Endl;
#endif
    
	MT->unsetBestFTBck();
    
    
    
    
    
    
    if(!(MT->wcmp(FileNameBurstList,"None"))){
        std::ifstream fIn;
        char junk[128];
        char TMPType[128];
        int iBurst(0);
        double TMPamp(0.), TMPcent(0.), TMPpolar(0.);
        
        
        //! **** Free memory previously allocated
        if(amplitudeInit!=NULL){
            MT->Free(amplitudeInit,Nb*sizeof(double));
            MT->Free(centralTime,Nb*sizeof(double));
            MT->Free(cos2Pol,Nb*sizeof(double));
            MT->Free(msin2Pol,Nb*sizeof(double));
        }
        
        Cout << "Load list of bursts from " << FileNameBurstList << " reading the " << BurstType << " only ..." << Endl;
        //! *** Count the number of busrts with type (Cusp/Kink) corresponding to the one of GW 
        fIn.open(FileNameBurstList);
        if(fIn == NULL){
            Cout << " ERROR: LCGWCosmicString::init : Can not open the file " << FileNameBurstList << " ." << Endl;
            throw std::invalid_argument(" ERROR: LCGWCosmicString::init : Can not open the file.");
        }
		if(fIn.peek() == '#')
            fIn.ignore(16384,'\n');
        Nb=0;
        while(!fIn.eof()){
            fIn >> junk >> junk >> junk >> TMPType;
            if(MT->wcmp(TMPType, BurstType))
                Nb++;
        }
        fIn.close();
        fIn.clear();
        
        //! *** Allocate memory
        amplitudeInit = (double*)MT->AllocMemory(Nb*sizeof(double));
        centralTime   = (double*)MT->AllocMemory(Nb*sizeof(double));
        cos2Pol  = (double*)MT->AllocMemory(Nb*sizeof(double));
        msin2Pol = (double*)MT->AllocMemory(Nb*sizeof(double));
        
        //! *** Load data
        fIn.open(FileNameBurstList);
        if(fIn.peek() == '#')
            fIn.ignore(16384,'\n');
        iBurst=0;
        while(!fIn.eof()){
            fIn >> TMPcent >> TMPamp >> junk >> TMPType;
            if(MT->wcmp(TMPType, BurstType)){
                amplitudeInit[iBurst] = TMPamp;
                centralTime[iBurst] = TMPcent*LC::Yr_SI;
                TMPpolar = Polarization + MT->RandUniform(0., PolarizationDelta);
                cos2Pol[iBurst]  = cos(2*TMPpolar);
                msin2Pol[iBurst] = -sin(2*TMPpolar);
                iBurst++;
            }
        }
        fIn.close();
        fIn.clear();
    }
    
    //! *** Desactivate the computation of polarisation with standard GW
    setParam(2, 0.);
    ComputePolarization = false;
	
	/*! *** Compute some fixed quantities :
	 * \f[ h{+,0} = A (1+ \cos^2 \iota)  \f]
	 * \f[ h{\times,0} = - 2 A (\cos \iota) \f]
	 * \f[ \omega = 2 \pi f \f]
	 * \f[ \dot{\omega} = \pi \dot{f} \f]
	 */
	
	
	/*hp0 = Amp * ( 1. + cos(inc)*cos(inc) ) ;
     hc0 = -2. * Amp * cos(inc) ;
     om = 2. * M_PI * Freq;
     Dom = M_PI * DFreq * 1.e-16;
     
     FreqMin = 1.e-6;
     FreqMax = 2.*Freq;
     */
    
	//------------------------------------
    
    int Nmin;
	double hp_temp;
	dcomplex hphigh;
	double temp1;
	double cc;
    double tmp, tmp2;
    double fftwNorm;
    double freqTemp;
    bool  highFlag;
    
    //! From LISACode simulation
	fmin=1/TObsReal;
    fmax=1/(2.*dtReal);
    
	fminSim=flow/20.;
    //fminSim=fmin;
    //fminSim=1e-6;
	TObsSim=1/fminSim;

    dtSim=dtReal/2.;
	Nmin=TObsSim/dtSim;
    N=2;
    while (N<Nmin) {
        N *= 2;
    }
    Nf=MT->getNfFTreal(N);
    dtSim=TObsSim/N;
    fmaxSim=1./(2.*dtSim);
    fhigh=fmaxSim/10.; 
    //fhigh = 0.1;
    
	
	df= (fmaxSim-fminSim)/Nf;
	Cout << "TObSim N dt Nf df = " <<TObsSim<< " " <<N<< " " <<dtSim<< " " <<Nf<< " " <<df<< Endl ;
    Cout <<"Sim: fmin fmax = " << fminSim << " " << fmaxSim << Endl;
    Cout <<"flow fhigh = " << flow << " " << fhigh <<  Endl;
    
    
    
	//ParentSize=epsilon*Gamma*MassL*LC::G_SI/(LC::c_SI*LC::c_SI)*tcosmo; //tcosmo est le temps cosmique de création voir constraint on scomic strings (la taille de la corde est définie comme l=a*t => a franction of the horizon at the time of formation 
	
	//amplitudeInit=(LC::G_SI*MassL*pow(ParentSize,2./3))/DL;
	
	//getNtFTreal
    
    
	hp_time = (double*)MT->AllocMemory(N*sizeof(double)); //h+ en fonction du temsp : tableau N
	hp_freq = (dcomplex*)MT->AllocMemory(Nf*sizeof(dcomplex)); 
    
	freqTemp = fminSim; //on initialise
	highFlag = true;

    //! *** FFTW normalisation
    fftwNorm = 1./(dtSim*N); 
    Cout << "FFTW normalisation = " << fftwNorm << Endl;
    
    for(int i=0; i<=Nf; i++)
    {
        
        
        
        
        hp_freq[i] = fftwNorm*pow(freqTemp,-q);
        
        //! ** Low frequency reduction
        tmp=flow/freqTemp;
        tmp2=1./(1.+tmp*tmp);
        hp_freq[i] *= tmp2*tmp2*tmp2*tmp2;

        
        
        //! ** High freuqncy reduction
        
        if  (freqTemp>=fhigh) {
            
            if (highFlag) {
                //temp1=1+(flow/freqTemp)*(flow/freqTemp); // >>> AP : Pas utilise ?
                //hp_temp=amplitudeInit * pow(freqTemp,-q);
                //hp_freq[i]=hp_temp * dcomplex( cos(-M_PI*freqTemp*centralTime) , sin(-M_PI*freqTemp*centralTime) ); 
                
                hphigh=hp_freq[i];
                cc=0;
                highFlag=false;
                
            }
            
            //hp_freq[i]=hphigh*cos((M_PI/2)*((freqTemp-fhigh)/(fmaxSim-fhigh)));
            //hp_freq[i]=hphigh*(1/(exp(cc))); // >>> AP : Justification mathematique (fonction) ?
            
            cc=cc+1;
            //Cout <<hp_freq[i]<< " " <<cc<< Endl ;
            
            //hp_freq[i] = hphigh*exp(1-freqTemp/fhigh);
            
            hp_freq[i] *= exp(1-freqTemp/fhigh);
            //hp_freq[i] *= cos((M_PI/2)*(freqTemp-fhigh)/(fmaxSim-fhigh));
        }
        
        
        hp_freq[i] *= dcomplex( cos(-M_PI*freqTemp*TObsSim) , sin(-M_PI*freqTemp*TObsSim) );
        
        
        freqTemp+=df; //on augmente de df en fréquence  
	    //Cout <<hp_freq[i]<< " " <<freqTemp<< " " <<" hp_freq[i] freqTemp"<< Endl ;
#ifdef _DEBUG_GW_
        DEBUGfCheck->precision(15);
	    (*DEBUGfCheck) << freqTemp << " " <<hp_freq[i]<< " " << sqrt(hp_freq[i].real()*hp_freq[i].real() + hp_freq[i].imag()*hp_freq[i].imag()) << " " << atan2(hp_freq[i].imag(),hp_freq[i].real()) << Endl;
        
#endif
    }
	
	//Cout <<"A="<< " " <<amplitudeInit<< Endl ;
    
    MT->FTMakeBck(hp_freq, hp_time, Nf);// on fait la tf de h+(freqTemp) et on ajoute la valeur à h+(t)
    
	MT->Free(hp_freq,Nf*sizeof(dcomplex)); 
    hp_freq = NULL;
	
	/*for (int i=0; i<=N; i++)
     {
     
     Cout <<"hp en fonction du temps (tf de hp_f)"<< " " <<hp_time[i]<< " " <<"i=" << " " <<i<< " " <<"/"<< " " << N << Endl;
     }*/
#ifdef _DEBUG_GW_
	DEBUGfCheck->close();
    DEBUGfCheck->clear();
    sprintf(fNCheck,"CheckGWCosmicString_t_%d.txt",iRCheck);
	std::cerr << "DEBUG:GW:LCGWCosmicString : t : File = " << fNCheck << Endl;
	DEBUGfCheck = new std::ofstream(fNCheck);
    for(int i=0; i<N;i++){
        DEBUGfCheck->precision(15);
        (*DEBUGfCheck) << i*dtSim << " " << hp_time[i] << Endl;
    }
    DEBUGfCheck->close();
#endif
	/*	for  (int i(centralTime-10) ; i<=centralTime+10 ; i ++)
     {
     Cout <<i<< " " <<hp_time[i]<< Endl ;
     }*/
    return 0;
    
}
// *********************
// *  Running methods  *
// *********************


//Calculation on h+(t) based on h+(f)
/*
 double LCGWCosmicString::hp_time(dcomplex fmin,dcomplex df, int N)
 {
 
 dcomplex 
 hp_time=(double*)MT->AllocMemory(N*sizeof(double)); //h+ en fonction du temsp : tableau N
 hp_freq=(dcomplex*)MT->AllocMemory(N*sizeof(double)); 
 
 dcomplex freqTemp=fmin; //on initialise
 
 for(int i=0; i<=N; i++)
 {
 //dcomplex hp_freqTemp=amplitudeInit*pow(freqTemp,-4/3)*pow((pow((1+(fLow/freqTemp)),2)),-4); //on calcul h+(f) avec f=freqTemp
 
 hp_freq[i]=amplitudeInit*pow(freqTemp,-4/3)*pow((pow((1+(fLow/freqTemp)),2)),-4);
 
 freqTemp+=df; //on augmente de df en fréquence     
 
 }
 hp_time=FTMakeBck(h+_freq, hp_time, N);// on fait la tf de h+(freqTemp) et on ajoute la valeur à h+(t)
 
 return hp_time;
 }
 */
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!  Interpolation de lagrange  !!!!!!!!!!!!!!!!!!!!!!!!!!!


void LCGWCosmicString::Computehpc(double t)
{
	/*! *** Compute the component of the waveform (in source frame) :
	 * \f[ h^S_{+}(t) = h{+,0} \cos [ \omega t + \dot{\omega} t^2 ] \f]
	 * \f[ h^S_{\times}(t) = h{\times,0} \sin [ \omega t + \dot{\omega} t^2 ] \f]
	 */
    //Cout <<"t="<< " " <<t<< Endl ;
    int ib(0);
    double timeIndex(0.), hSp(0.);
    
    hBpLast=0.;
    hBcLast=0.;
    for(int iBurst=0; iBurst<Nb; iBurst++){
        timeIndex = (t-t0Real+TObsSim/2.-centralTime[iBurst])/dtSim; //NOTE : timeIndex porte mal son nom, ce n'est pas un entier
        ib = MT->ifloor(timeIndex);
        //Cout <<t<< " " <<timeIndex<< Endl ;
        if ((timeIndex>0)&&(timeIndex<N-1)){

            //Cout <<t<< " " <<timeIndex<< Endl ;
            //if ( ( (t-t0Real)<(centralTime-TObsSim/2) ) || ( (t-t0Real)>(centralTime+TObsSim/2) ) ) { // pour simuler que là où il y a le signal.
            //    hBpLast=0;
            //Cout <<"2"<< Endl ;
            //}
            //else {
            //Cout <<"3"<< Endl ;
            
            //hBpLast=((hp_time[MT->ifloor(timeIndex)]-hp_time[MT->iceil(timeIndex)])/(MT->ifloor(timeIndex)-MT->iceil(timeIndex)))*timeIndex + ((MT->ifloor(timeIndex)*hp_time[MT->iceil(timeIndex)]-MT->iceil(timeIndex)*hp_time[MT->ifloor(timeIndex)]) / (MT->ifloor(timeIndex)-MT->iceil(timeIndex))); /* linear interpolation f(x)=((ya-yb)/(xa-xb))*x + (xa*yb -xb*ya)/(xa-xb) where f(xa)=ya and f(xb)=yb*/
            hSp = amplitudeInit[iBurst]*(((hp_time[ib]-hp_time[ib+1]))*timeIndex + ((ib*hp_time[ib+1]-(ib+1)*hp_time[ib]))); /* linear interpolation f(x)=((ya-yb)/(xa-xb))*x + (xa*yb -xb*ya)/(xa-xb) where f(xa)=ya and f(xb)=yb*/
            
            hBpLast += cos2Pol[iBurst]  * hSp;     
            hBcLast += msin2Pol[iBurst] * hSp;  
            
            /* if (t<=centralTime-10) { //&& t<=centralTime+10) {
             Cout << t << " " <<centralTime-10<< " " << hBpLast << Endl ; 
             }*/
            /*																       //hBpLast=(hp_time[MT->ifloor(timeIndex)]+hp_time[MT->iceil(timeIndex)])/2; //We take the lower and higher int of timeIndex and return the mean of hp_time at both index.
             //Cout <<"Computehpc, TimeIndex"<< " " << timeIndex <<  Endl;*/
            //}
        }
    }
    
    
    
    //hBcLast = hc0 * sin( om*t + Dom*t*t + Phi0);
	
	
	/*
     if((t>2.1460e7)&&(t<2.1461e7)){
     if(DEBUGfCheck==NULL){
     DEBUGfCheck = new std::ofstream("CheckData.txt");
     DEBUGfCheck->precision(15);
     }
     (*DEBUGfCheck) << t << " " << hBpLast << " " << cos( om*t + Dom*t*t + Phi0) << " " << om*t + Dom*t*t << " " << (om*t + Dom*t*t + Phi0)  << Endl;
     }
     /
     }
     
     */
}

// *******************
// *  Other methods  *
// *******************

void LCGWCosmicString::DispInfo(char * BTab)
{
	if(MT->Disp()){
		DispInfoBase(BTab);
	}
    
	//Cout << BTab << "\t- Polarization   = " << polarization << " rad" << Endl;
    Cout << BTab << "\t- q                 = " << q << Endl;
    Cout << BTab << "\t- Low frequency     = " << flow << " Hz" << Endl;
    //Cout << BTab << "\t- High frequency = " << fhigh << " Hz" << Endl;
    Cout << BTab << "\t- Burst type        = " << BurstType << Endl ;
    Cout << BTab << "\t- Bursts' list file = " << FileNameBurstList << Endl ;
    Cout << BTab << "\t- " << Nb << " bursts : Index Amplitude CentralTime Polarization cos(2Pol)" << Endl;
    for (int iBurst=0; iBurst<Nb; iBurst++)
        Cout << BTab << "\t" << iBurst << " " << amplitudeInit[iBurst] << " " << centralTime[iBurst] << " " << acos(cos2Pol[iBurst])/2. << " " << cos2Pol[iBurst] << Endl; 
    
	
}


void LCGWCosmicString::DispAllParam(std::ostream * out)
{
	for(int iP=0; iP<NParams; iP++)
		(*out) << " " << getParam(iP);
}


void LCGWCosmicString::DispAllParamName(std::ostream * out)
{
	(*out) << "Bet Lam Pol Gmu DL alpha flow tcosmo";
}




// end of LISACODE-GWSpinBBH2.cpp
