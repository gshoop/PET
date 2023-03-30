#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <list>
#include <cstdlib>
#include <sstream>
#include <random>
#include <iomanip>
#include <limits>
using namespace std;

#define PI 3.14159265

struct debuginfo
{
	int totclus;
	int clussamep;
	int totgroup;
	int clusright;
};

struct debuginfo debugclus[2]={{0,0,0,0}, {0,0,0,0}};


struct hits
{
	int eventid;
	int particleid;
	int numcompton;
	int numrayleigh;
	int vid1;
	int vid2;
	int vid3;
	double time;
	float energy;
	double blurtime;
	float blurenergy;
	float x;
	float y;
	float z;
	float bx;
	float by;
	float bz;
	float srcx;
	float srcy;
	float srcz;
};

struct lor
{
	struct hits hits1;
	struct hits hits2;
	int Misaln;
	int Incorrseq;
	int IncorrseqCS;	// only for cross strip detectors.
	int Multiclus;
	int failtype;
	float dda;
	float totale1;
	float totale2;
};

struct csortingpara
{
	int checksecdiff;
	int minsecdiff;
	int totnumsec;
	int dbinfo;
};

vector<string> explode( string s, char c); 
vector<hits> BinCrystal(vector<hits>  group);
lor * csorting (vector<hits> group, lor * plorf,int numcolumn, int numcolumnm, int numrowr, int numrowm, float doires, float Elow, float Ehigh, float Ethresh, int * failnum, string blurtype,int * pgroupwarn, int * pnumoverl, int maxclus, int maxevent, int minsecdiff, csortingpara otherpara);
float Etheta(float energy);	
float Ptheta(float pos1[], float pos2[], float pos3[], int * poverlap);
int PanelBlur(float pos[], float * pblur, int crystal[], int numcr, int numcm, int numrowr, int numrowm, float doires, string blurtype);
int checkclus(vector<hits> group, vector<vector<int>> ematrix, int kk);
hits newCrossStrip(hits event1, hits event2, int numcolumnm);
vector<vector<int>> EnergyCheck(const vector<vector<int>>& fmatrix, const vector<hits>& group, const float &Elow, const float &Ehigh, int &numclusafter, int &maxnumevent,const vector<int> &emindex, vector<int> &nemindex);


int main(int argc,char* argv[]){
	vector<string> vpara;
	float Elimlo, Elimhi, CTSthresh, Ethresh, Eblur_perc, Tblur, doires;
    int numcolumnr=0, numrowr=0, numcolumnm=0, numrowm=0;
	int maxclus=2, maxevent=2, binevent=1, minsecdiff=1;
	csortingpara otherpara;

	string line;
	stringstream ss;
	ifstream config ("config.txt");

	if (config.is_open())
	{
		while ( getline (config,line) )
		{
			vpara=explode(line,' ');
			if(vpara[0]=="Elimlo") {ss<<vpara[2];ss>>Elimlo;ss.clear();} else
			if(vpara[0]=="Elimhi") {ss<<vpara[2];ss>>Elimhi;ss.clear();} else
			if(vpara[0]=="CTSthresh") {ss<<vpara[2];ss>>CTSthresh;ss.clear();} else
			if(vpara[0]=="Ethresh") {ss<<vpara[2];ss>>Ethresh;ss.clear();} else
			if(vpara[0]=="Eblur") {ss<<vpara[2];ss>>Eblur_perc;ss.clear();} else
			if(vpara[0]=="Tblur") {ss<<vpara[2];ss>>Tblur;ss.clear();} else
			if(vpara[0]=="DOIres") {ss<<vpara[2];ss>>doires;ss.clear();} else
			if(vpara[0]=="numcolumnInRsector") {ss<<vpara[2];ss>>numcolumnr;ss.clear();} else
			if(vpara[0]=="numrowInRsector") {ss<<vpara[2];ss>>numrowr;ss.clear();} else
			if(vpara[0]=="numcolumnInModule") {ss<<vpara[2];ss>>numcolumnm;ss.clear();} else
			if(vpara[0]=="numrowInModule") {ss<<vpara[2];ss>>numrowm;ss.clear();} else
			if(vpara[0]=="MaximumNumberClusters") {ss<<vpara[2];ss>>maxclus;ss.clear();} else
			if(vpara[0]=="MaximumNumberEventsInACluster") {ss<<vpara[2];ss>>maxevent;ss.clear();} else
			if(vpara[0]=="BinEvents") {ss<<vpara[2];ss>>binevent;ss.clear();} else
			if(vpara[0]=="CheckSecDiff") {ss<<vpara[2];ss>>otherpara.checksecdiff;ss.clear();} else
			if(vpara[0]=="MinSecDiff") {ss<<vpara[2];ss>>otherpara.minsecdiff;ss.clear();} else
			if(vpara[0]=="TotNumSec") {ss<<vpara[2];ss>>otherpara.totnumsec;ss.clear();} else
			if(vpara[0]=="Debug") {ss<<vpara[2];ss>>otherpara.dbinfo;ss.clear();}	

		}
		config.close();
	}
	else cout << "Unable to open file"<<endl;

	cout<<"-------------------------------------------"<<endl;
	cout<<"Input parameters:"<<endl;
	cout<<"E lower threshold: "<<Elimlo<<" keV"<<endl;
	cout<<"E upper threshold: "<<Elimhi<<" keV"<<endl;
	cout<<"Coincidence time window: "<<CTSthresh<<" s"<<endl;
	cout<<"Lower energy threshold: "<<Ethresh<<" keV"<<endl;
	cout<<"Energy resolution: "<<Eblur_perc<<"%"<<endl;
	cout<<"Time resolution: "<<Tblur<<" s"<<endl;
	cout<<"DOI resolution: "<<doires<<" mm"<<endl;
	cout<<"Number of columns in rsector: "<<numcolumnr<< endl;
	cout<<"Number of rows in rsector: "<<numrowr<< endl;
	cout<<"Number of columns in module: "<<numcolumnm << endl;
	cout<<"Number of rows in module: "<<numrowm << endl;
	cout<<"Maximum number of clusters: "<<maxclus << endl;
	cout<<"Maximum number of events in a cluster: "<<maxevent << endl;
	cout<<"Bin event according to crystal (0:no, 1:yes)? "<<binevent<<endl;
	cout<<"Check minimum rsector difference during coincidence sorting (0:no, 1:yes)? "<<otherpara.checksecdiff<<endl;
	if(otherpara.checksecdiff==1)
	{
		cout<<"Minimum rsector difference: "<<otherpara.minsecdiff<<endl;
		cout<<"Total number of rsectors: "<<otherpara.totnumsec<<endl;
		
	}
	else otherpara.minsecdiff=1;
	cout<<"-------------------------------------------"<<endl;

	string failtype[8] = {"less than 2Node","less than 2Node after Ethresh","less than Ewin", "more than Ewin","More than Max Num of clusters","More than Max Num of events in a cluster","All clusters in the group are less than MinSecDiff (excluding less than 2node case)", "At least one clusters have two events with same anode or same cathod"};
	int failnum[8] = {0,0,0,0,0,0,0,0};

	string filein=argv[1];
	string ifreorg=argv[2];
	string blurtype=argv[3];
	string ifSrcPos=argv[4];
	ifstream fin (filein);
	ofstream fout (filein.append(".lst"));
	vector<hits> group;
	lor lorf;
	lor * plorf;
	plorf=&lorf;
	hits event;
	hits * pevent;
	pevent = & event;
	int tR=0,tS=0,tSP=0,tM=0,tI=0,tICS=0, tT=0,tngroup=0,tngroupaccept=0;
	int R,S,SPhoton,M,I,ICS,T;
	vector<string> l;
	string inttype;
	float energy, esigma, eblur;
	double time, tsigma, tblur, trigger;
	unsigned seede = 100, seedt = 200;
	default_random_engine generatore (seede), generatort (seedt);
	//normal_distribution<double> distributione (energy,esigma), distributiont (time,tsigma);

	int groupwarn=0;
	int * pgroupwarn;
	pgroupwarn= &groupwarn;

	//number of group accepted that is influenced by position overlap among events.
	int numoverl=0;
	int * pnumoverl = &numoverl;

	vector<hits> newgroup;


	if(!fin.is_open() || !fout.is_open())
	{
		cout<<"Error: input or output files are not opened!!"<<endl;
		abort();
	}

	while ( getline (fin,line) )
	{
		R=0;
		S=0;
		SPhoton=0;
		M=0;
		I=0;
		T=0;
		l.clear();
		l=explode(line,' ');
		if(ifreorg=="reorg") {l.insert(l.begin()+4,"0");l.erase(l.begin()+10);}
		
		//cout<<l.at(5)<<" "<<l.at(6)<<" "<< l.at(7)<<endl;

		inttype=l.at(22);
		if(inttype=="ElectronIonisation" || inttype=="Bremsstrahlung" || inttype=="annihil" || inttype=="MultipleScattering") continue;

		energy=stof(l.at(11))*1000.;
		esigma=Eblur_perc / 100. * 511. * sqrt(energy/511.) / 2.355;
		normal_distribution<double> distributione (energy,esigma);
		if(energy > 0.){
			eblur=distributione(generatore);
			if(eblur<0.) eblur=0.;
		}
		else if(energy==0.) eblur=0.;
		else cout<<line<<endl;

		time=stod(l.at(10));
		tsigma=Tblur/2.355;
		normal_distribution<double> distributiont (time,tsigma);
		tblur=distributiont(generatort);


		pevent->eventid=stoi(l.at(1));
		pevent->particleid=stoi(l.at(17));
		pevent->numcompton=stoi(l.at(20));
		pevent->numrayleigh=stoi(l.at(21));
		pevent->vid1=stoi(l.at(5));	//panelID
		pevent->vid2=stoi(l.at(6));	//moduleID, in CZT, it's crystal ID
		pevent->vid3=stoi(l.at(7));	//anode and cathod cellID
		pevent->time=time;
		pevent->energy=energy;
		pevent->blurtime=tblur;
		pevent->blurenergy=eblur;
		pevent->x=stof(l.at(13));
		pevent->y=stof(l.at(14));
		pevent->z=stof(l.at(15));

		if(ifSrcPos == "norm")	//For including position of source particle
		{
			pevent->srcx=stof(l.at(25));
			pevent->srcy=stof(l.at(26));
			pevent->srcz=stof(l.at(27));
		}

		//if(tngroup<5) {cout<<pevent->x<<endl;}

		if(group.size()==0) {
			group.push_back(event);
			trigger=tblur;
		}
		else {
			if(abs(tblur-trigger)<CTSthresh) group.push_back(event);
			else {
				tngroup+=1;
				plorf->failtype=0;

				if(binevent==1)
				{
					newgroup = BinCrystal(group);
					plorf=csorting(newgroup,plorf,numcolumnr,numcolumnm, numrowr, numrowm, doires, Elimlo, Elimhi, Ethresh, failnum, blurtype, pgroupwarn, pnumoverl, maxclus, maxevent, minsecdiff, otherpara);
				}
				else
				{
					plorf=csorting(group,plorf,numcolumnr,numcolumnm, numrowr, numrowm, doires, Elimlo, Elimhi, Ethresh, failnum, blurtype, pgroupwarn, pnumoverl, maxclus, maxevent, minsecdiff, otherpara);
				}

			//	//weird things about the group. need printing out.
			//	if(groupwarn==2) 
			//	{
			//		for(int gi=0; gi<group.size(); gi++)
			//		{
			//			cout<<group.at(gi).eventid<<endl;
			//		
			//		}
			//		cout<<"End of one group"<<endl;
			//		abort();
			//	}

				if(plorf->failtype==0)
				{
					tngroupaccept+=1;
					if(plorf->hits1.eventid != plorf->hits2.eventid) {R=1; tR+=1;}
					else if(plorf->hits1.numcompton != 0 || plorf->hits2.numcompton != 0 || plorf->hits1.numrayleigh != 0 || plorf->hits2.numrayleigh != 0) {S=1; tS+=1;}
					else if(plorf->hits1.particleid == plorf->hits2.particleid) {SPhoton=1; tSP+=1;}
					else if(plorf->Misaln==1) {M=1; tM+=1;}
					else 
					{
						T=1;
						tT+=1;
						if(plorf->Incorrseq==1) {I=1; tI+=1;}
						if(plorf->IncorrseqCS==1) tICS+=1;
					}

					//write to file
					if(tngroupaccept==1)
					{

						fout<<"EventID1(1) "<<"BlurredPosX1(2) "<<"BlurredPosY1(3) "<<"BlurredPosZ1(4) "<<"AccuratePosX1(5) "<<"AccuratePosY1(6) "<<"AccuratePosZ1(7) "<<"Rsector1(8) "<<"Module1(9) "<<"SubModule1(10) "<<"TotalEnergy1(11) "<<"EventID2(12) "<<"BlurredPosX2(13) "<<"BlurredPosY2(14) "<<"BlurredPosZ2(15) "<<"AccuratePosX2(16) "<<"AccuratePosY2(17) "<<"AccuratePosZ2(18) "<<"Rsector2(19) "<<"Module2(20) "<<"SubModule2(21) "<<"TotalEnergy2(22) "<<"TrueIndicator(23) "<<"RandomIndicator(24) "<<"ScatterIndicator(25) "<<"SamePhotonIndicator(26) "<<"MisalignedIndicator(27) "<<"IncorrectIndicator(28) "<<"DDA(29) "<<"SourcePosX(30) "<<"SourcePosY(31) "<<"SourcePosZ(32) "<<"TimeStamp(33)"<<endl;
					}

					fout << plorf->hits1.eventid <<" "<< plorf->hits1.bx <<" "<< plorf->hits1.by <<" "<< plorf->hits1.bz <<" "<< plorf->hits1.x <<" "<< plorf->hits1.y <<" "<< plorf->hits1.z <<" "<< plorf->hits1.vid1 <<" "<< plorf->hits1.vid2 <<" "<< plorf->hits1.vid3 <<" "<< plorf->totale1 <<" "<< plorf->hits2.eventid <<" "<< plorf->hits2.bx <<" "<< plorf->hits2.by <<" "<< plorf->hits2.bz <<" "<< plorf->hits2.x <<" "<< plorf->hits2.y <<" "<< plorf->hits2.z <<" "<< plorf->hits2.vid1 <<" "<< plorf->hits2.vid2 <<" "<< plorf->hits2.vid3 <<" "<< plorf->totale2 <<" "<< T <<" "<< R <<" "<< S <<" "<< SPhoton <<" "<< M <<" "<< I <<" "<< plorf->dda <<" ";
					
					if(ifSrcPos == "norm") fout<< plorf->hits1.srcx <<" "<< plorf->hits1.srcy <<" "<< plorf->hits1.srcz<<" "<< plorf->hits1.blurtime<<endl;
					else fout<< plorf->hits1.blurtime<<endl;
						 
				}

				group.clear();
				group.push_back(event);
				trigger=tblur;
			}
		}


		//cout<<setprecision(20)<<tblur<<' '<<eblur<<' '<<stof(l.at(26))<<' '<<stof(l.at(25))<<endl;
		//cout<<l.at(10)<<' '<<l.at(11)<<' '<<setprecision(20)<<time<<' '<<energy<<endl;
		//cout<<l.at(10)<<' '<<l.at(11)<<' '<<setprecision(20)<<time<<' '<<energy<<' '<<tblur<<' '<<eblur<<' '<<l.at(25)<<' '<<l.at(26)<<endl;

	}
	
	cout<<"----------------------------------------"<<endl;
	cout<<tngroup<<" total number of group found based on time window"<<endl;
	cout<<tngroupaccept <<" total number of group accepted. "<< ((double)tngroupaccept)/tngroup * 100. <<"%" <<endl;
	cout<<" "<<endl;
	cout<<"less than 2Node: "<< failnum[0] <<endl;
	cout<<"less than 2Node after Ethresh: "<< failnum[1] <<endl;
	//cout<<"less than Ewin: "<< failnum[2]<<endl;
	//cout<<"more than Ewin: "<< failnum[3] <<endl;
	cout<<"Out of energy window: "<< failnum[2]<<endl;
	cout<<"More than "<< maxclus<<" clusters: "<< failnum[4] <<endl;
	cout<<"More than "<<maxevent<<" events in a cluster: " << failnum[5] <<endl;
	cout<<"All clusters in the group are less than MinSecDiff (excluding less than 2node case): "<<failnum[6]<<endl;
	cout<< failtype[7] << ": " << failnum[7]<<endl;
	cout<<"-----------------------------------------"<<endl;
	cout<<"Number of true coincidences: "<< tT << endl;
	cout<<"Number of random coincidences: " << tR <<endl;
	cout<<"Number of scattering coincidences: "<< tS <<endl;
	cout<<"Number of coincidences generated by the same photon: "<< tSP << endl;
	cout<<"Number of misaligned coincidences (caused by event with energy below Ethresh): " << tM << endl;
	cout<<"Number of true coincidences with incorrected sequences: " << tI << endl;
	cout<<"Incorrected sequences for cross strip: " << tICS << endl;
	cout<<"------------------------------------------"<<endl;
	if(otherpara.dbinfo==1)
	{
		cout<<"Debug information:"<<endl;
		cout<<"Number of accepted coincidences affected by position overlap among events: "<< numoverl << endl;
		cout<<"Number of clusters: "<< debugclus[0].totclus << ". Number of clusters with all same particleID and eventID: "<< debugclus[0].clussamep << ". "<< debugclus[0].clussamep/(double)debugclus[0].totclus * 100. <<"% (Unaccepted groups are also considered)" <<endl;
		cout<<"Number of groups: " <<debugclus[0].totgroup << ". Number of groups with right clustering: " << debugclus[0].clusright<< ". " << (double)debugclus[0].clusright/debugclus[0].totgroup * 100. <<"% (Unaccepted groups are also considered)"<<endl;

		cout<<"Number of clusters: "<< debugclus[1].totclus << ". Number of clusters with all same particleID and eventID: "<< debugclus[1].clussamep << ". "<< debugclus[1].clussamep/(double)debugclus[1].totclus * 100. <<"% (Only accepted groups are considered)" <<endl;
		         cout<<"Number of groups: " <<debugclus[1].totgroup << ". Number of groups with right clustering: " << debugclus[1].clusright<< ". " << (double)debugclus[1].clusright/debugclus[1].totgroup * 100. <<"% (Only accepted groups are considered)"<<endl;
		cout<<"------------------------------------------"<<endl;
	}

	fin.close();
	fout.close();

	return 0;
}

//break string down to string vector based on space 
vector<string> explode(string s, char c)
{
	string buff="";
	vector<string> v;
	char n;
	
	for(unsigned i=0; i<s.length(); ++i)
	{
		n=s.at(i);
		if(n != c) buff+=n; else
		if(n == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if(buff != "") v.push_back(buff);
	
	return v;
}

//bin events into single according to crystal number.
vector<hits> BinCrystal(vector<hits>  group)
{
	vector<hits> newgroup;
	newgroup.clear();
	vector<int> p,m,c,relay;
	vector<vector<int>> ematrix;
	int indcnew = 1;

    p.push_back(group.at(0).vid1);
    m.push_back(group.at(0).vid2);
    c.push_back(group.at(0).vid3);
    relay.assign(1,0);
    ematrix.push_back(relay);

    for(int i=1;i<group.size();i++)
    {
        indcnew=1;
        for(int j=0;j<p.size();j++)
        {
            if(group.at(i).vid1==p.at(j) && group.at(i).vid2==m.at(j) && group.at(i).vid3==c.at(j)) 
            {
                indcnew=0;
                ematrix.at(j).push_back(i);
                break;
            }
        }
        if(indcnew==1)
        {
            p.push_back(group.at(i).vid1);
            m.push_back(group.at(i).vid2);
            c.push_back(group.at(i).vid3);
            relay.clear();
            relay.assign(1,i);
            ematrix.push_back(relay);

        }
    }

	int num1=-1, num2=-1;
	float px=0., py=0., pz=0.; 
	for(int i=0; i<ematrix.size(); i++)
	{
		num1= ematrix.at(i).at(0);
		px= group.at(num1).blurenergy * group.at(num1).x;
		py= group.at(num1).blurenergy * group.at(num1).y;
		pz= group.at(num1).blurenergy * group.at(num1).z;

		for(int j=1; j<ematrix.at(i).size(); j++)
		{
			num2= ematrix.at(i).at(j);
			group.at(num1).energy += group.at(num2).energy;
			group.at(num1).blurenergy += group.at(num2).blurenergy;
			px += group.at(num2).blurenergy * group.at(num2).x;
			py += group.at(num2).blurenergy * group.at(num2).y;
			pz += group.at(num2).blurenergy * group.at(num2).z;
		}

		group.at(num1).x = px / group.at(num1).blurenergy;
		group.at(num1).y = py / group.at(num1).blurenergy;
		group.at(num1).z = pz / group.at(num1).blurenergy;
		newgroup.push_back(group.at(num1));
	}


	return newgroup;
}

//coincidence sorting function
lor * csorting (vector<hits> group, lor * plorf,int numcolumn, int numcolumnm, int numrowr, int numrowm, float doires, float Elow, float Ehigh, float Ethresh, int * failnum, string blurtype, int * pgroupwarn, int * pnumoverl, int maxclus, int maxevent, int minsecdiff, csortingpara otherpara)
{
	vector<int> p,m,mr,mc,relay;
	vector<vector<int>> ematrix;
	int numnewevent;
	int indcnew;
	int numpanel=0, numpanelafter=0;
	vector<int> pnool, pnoolafter;
	int indnewpanel;
	//int p=-1,m=-1,mr=-1,mc=-1;
	
	//check overlap of position in different events
	int overl=0;
	int * poverlap = & overl;


	p.push_back(group.at(0).vid1);
   	m.push_back(group.at(0).vid2); 
	mr.push_back(group.at(0).vid2/numcolumn);
   	mc.push_back(group.at(0).vid2%numcolumn);
	numnewevent=1;
	relay.assign(1,0);
	ematrix.push_back(relay);

	//p, m, mr, mc: panel number, module number, module row, module column
	//numnewevent: number of clusters
	//neweventindex: index of first event of a cluster in the group
	//numev: number of events in each cluster
	//indcnew: =1 if the event does not belong to any existing cluster, thus assigning a new cluster
	//ematrix: each vector inside it is a cluster showing all indexes in the group
	

	// Clustering process
	for(int i=1;i<group.size();i++)
	{
		indcnew=1;
		for(int j=0;j<p.size();j++)
		{
			// The number within if is only for CZT system.
			if(group.at(i).vid1==p.at(j) && abs(group.at(i).vid2/numcolumn - mr.at(j)) < 11 && abs(group.at(i).vid2%numcolumn - mc.at(j)) < 2) 
			{
				indcnew=0;
				ematrix.at(j).push_back(i);
				
				// Add on 10/05/2018. If the event is within a cluster, then updata the mr and mc
				mr.at(j) = (mr.at(j) * (ematrix.at(j).size() - 1) + group.at(i).vid2/numcolumn) / ematrix.at(j).size();
				mc.at(j) = (mc.at(j) * (ematrix.at(j).size() - 1) + group.at(i).vid2%numcolumn) / ematrix.at(j).size();
				break;
			}
		}
		if(indcnew==1)
		{
		    p.push_back(group.at(i).vid1);
			m.push_back(group.at(i).vid2); 
			mr.push_back(group.at(i).vid2/numcolumn);
			mc.push_back(group.at(i).vid2%numcolumn);
			numnewevent+=1;
			relay.clear();
			relay.assign(1,i);
			ematrix.push_back(relay);
	
		}
	}


	//check if clustering is correct. First, inside same cluster, eventid and photonid should be same. Second: between clusters, if eventid is same, then photonid should be different.
	if(otherpara.dbinfo==1)
	{
		checkclus(group, ematrix, 0);

//		int eventid=-1;
//		int pid=-1;
//		int sameclus= 1;
//		int diffclus = 1;
//
//		for(int i=0;i<ematrix.size(); i++)
//		{
//			debugclus.totclus += 1;
//			eventid=-1;
//			pid=-1;
//			sameclus=1;
//			for(int j=0; j<ematrix.at(i).size(); j++)
//			{
//				if(eventid==-1) {eventid=group.at(ematrix.at(i).at(j)).eventid; pid=group.at(ematrix.at(i).at(j)).particleid;}
//				else
//				{
//					if(group.at(ematrix.at(i).at(j)).eventid != eventid || group.at(ematrix.at(i).at(j)).particleid != pid)
//					{
//						sameclus=0;
//						break;
//					}
//				}
//			}
//			if(sameclus==1)
//			{
//				debugclus.clussamep += 1;
//			}
//		}
//
//		debugclus.totgroup += 1;
//		for(int i=0;i<ematrix.size()-1; i++)
//		{
//			for(int j=0;j<ematrix.at(i).size(); j++)
//			{
//				for(int ii=i+1; ii< ematrix.size(); ii++)
//				{
//					for(int jj=0; jj<ematrix.at(ii).size(); jj++)
//					{
//						
//						if(group.at(ematrix.at(i).at(j)).eventid == group.at(ematrix.at(ii).at(jj)).eventid && group.at(ematrix.at(i).at(j)).particleid == group.at(ematrix.at(ii).at(jj)).particleid)
//						{
//							diffclus=0;
//							goto groupcluscheck;
//
//						}
//					}
//				}
//			}
//		}
//groupcluscheck:
//		if(diffclus==1) debugclus.clusright+=1;

	}

//	//test if there is overlap in ematrix
//	for(int i=0; i<ematrix.size(); i++)
//	{
//		for(int j=0; j<ematrix.at(i).size(); j++)
//		{
//			for(int k1=0; k1<ematrix.size(); k1++)
//			{
//				for(int k2=0; k2<ematrix.at(k1).size(); k2++)
//				{
//					if(k1!=i || k2!=j)
//					{
//						if(group.at(ematrix.at(i).at(j)).x==group.at(ematrix.at(k1).at(k2)).x && group.at(ematrix.at(i).at(j)).y==group.at(ematrix.at(k1).at(k2)).y && group.at(ematrix.at(i).at(j)).z==group.at(ematrix.at(k1).at(k2)).z ) {cout<<"Warning: overlap in cluster!!"<<endl; *pgroupwarn+=1; return plorf;}
//					}
//				}
//			}
//		}
//	}
	
	//fmatrix: matrix of event with energy larger than Ethresh
	//numclusafter: number of clusters after checking energy
	//maxnumevent: maximum number of event in a cluster after checking energy
	int numclusafter=0;
	int maxnumevent=0, numevent=0;
	vector<vector<int>> fmatrix;
	vector<int> frelay;
	vector<int> emindex; //the index of each vector in fmatrix in the original ematrix.
	for(int i=0; i< ematrix.size();i++)
	{
		numevent=0;
		frelay.clear();
		for(int j=0; j<ematrix.at(i).size();j++)
		{
			if (group.at(ematrix.at(i).at(j)).blurenergy > Ethresh) {numevent+=1;frelay.push_back(ematrix.at(i).at(j));}
		}
		if(numevent> 0) {numclusafter+=1; emindex.push_back(i);fmatrix.push_back(frelay);}
		if(numevent> maxnumevent) maxnumevent=numevent;
	}

	//check number of panels before Ethresh
	numpanel = 1;
	pnool.clear();
	pnool.push_back(group[ematrix[0][0]].vid1);
	for(int k=1; k < ematrix.size(); k++){
		indnewpanel = 1;
		int panelnum = group.at(ematrix.at(k).at(0)).vid1;
		for(int j = 0; j < pnool.size(); j++) if(panelnum == pnool[j]) {indnewpanel = 0; break;}
		if(indnewpanel == 1) {
			numpanel += 1;
			pnool.push_back(panelnum);
		}
	}

	//check number of panels after Ethresh
    numpanelafter = 1;
	pnoolafter.clear();
	if(fmatrix.size() > 0){
	pnoolafter.push_back(group[fmatrix[0][0]].vid1);
    for(int k=1; k < fmatrix.size(); k++){
        indnewpanel = 1;
		int panelnum = group.at(fmatrix.at(k).at(0)).vid1;
        for(int j = 0; j < pnoolafter.size(); j++) if(panelnum == pnoolafter[j]) {indnewpanel = 0; break;}
        if(indnewpanel == 1) {
			numpanelafter += 1;
			pnoolafter.push_back(panelnum);
		}
    }
	}
	else numpanelafter = 0;

    if(numpanel<2) {plorf->failtype=1; *failnum+=1; return plorf;}
	if(numpanel>1 && numpanelafter<2) {plorf->failtype=2; *(failnum+1)+=1; return plorf;}

	//check total energy of a cluster. only keep clusters within energy range
	vector<int> nemindex;
	fmatrix = EnergyCheck(fmatrix, group, Elow, Ehigh, numclusafter, maxnumevent, emindex, nemindex);
	emindex = nemindex;
    numpanelafter = 1;
    pnoolafter.clear();
    if(fmatrix.size() > 0){
    pnoolafter.push_back(group[fmatrix[0][0]].vid1);
    for(int k=1; k < fmatrix.size(); k++){
        indnewpanel = 1;
        int panelnum = group.at(fmatrix.at(k).at(0)).vid1;
        for(int j = 0; j < pnoolafter.size(); j++) if(panelnum == pnoolafter[j]) {indnewpanel = 0; break;}
        if(indnewpanel == 1) {
            numpanelafter += 1;
            pnoolafter.push_back(panelnum);
        }
    }
    }
    else numpanelafter = 0;  
	if(numpanelafter<2) {plorf->failtype=3; *(failnum+2)+=1; return plorf;}	//out of energy window	

	//check if there are too many  clusters.
	if(numclusafter>maxclus) {plorf->failtype=5;*(failnum+4)+=1; return plorf;}

	//check if there are too many events in a cluster.
	if(maxnumevent>maxevent) {plorf->failtype=6; *(failnum+5)+=1; return plorf;}


	//compton kinematic to find two clusters, and first event in each cluster
	//if(numnewevent>2) plorf->Multiclus=1;
	if(fmatrix.size() > 2) plorf->Multiclus = 1;	// Change on 10/05/1028
	else plorf->Multiclus=0; // number of clusters more than 2.(not considering Ethresh)

	float ddamin=1000.,dda,dda1,dda2,thetae,thetap;
	int clus1=-1, clus2=-1, elem1=-1, elem2=-1; //cluster number and element number for the lor that results in minimum dda
	float pos1[3], pos2[3], pos3[3]; //sequence: lor point 1, lor point 2 first interaction, lor point 2 second interaction
	int numpp=0;
    vector<int>	penum;

	//first check if there is a pair of PP
	int clusind[2]; //event index of PP in the group
    float pos[3], blurpos[3];
	int crystalind[3]; //for position blurring
//	for(int i=0; i<fmatrix.size(); i++)
//	{
//		if(fmatrix.at(i).size()==1)
//		{
//			if(panelnum==-1) 
//			{
//				numpp+=1;
//				panelnum= group.at(fmatrix.at(i).at(0)).vid1;
//				clus1=i;
//				elem1=0;
//			}
//			
//			else
//			{
//				if(panelnum!= group.at(fmatrix.at(i).at(0)).vid1)
//				{
//					numpp+=1;
//					clus2=i;
//					elem2=0;
//					break;
//				}
//			}
//		}
//	}

	int ifCrossStrip = 0;   //indicator for new hits.


	for(int i=0; i<fmatrix.size(); i++)
	{
		if(fmatrix.at(i).size()==1) penum.push_back(i);
	}
	int secdiff;
	if(penum.size()>1)
	{
		for(int i=0; i<penum.size()-1;i++)
		{
			for(int j= 1 + i; j<penum.size();j++)
			{
				secdiff=abs(group.at(fmatrix.at(penum.at(i)).at(0)).vid1 - group.at(fmatrix.at(penum.at(j)).at(0)).vid1);

				if((otherpara.checksecdiff==0 && secdiff >=1) || (otherpara.checksecdiff==1 && secdiff >= otherpara.minsecdiff && secdiff <= (otherpara.totnumsec - otherpara.minsecdiff)))
				{
					clus1=penum.at(i);
					clus2=penum.at(j);
					elem1=0;
					elem2=0;
					numpp=2;
					ddamin=0.;
					goto blurpp;
				}	
			}
		}
	}

	//cout<<numpp<<" "<<clus1<<" "<<clus2<<endl;
	//if there are two clusters with only 1 event. Then only blur the position of the two event and directly go to energy check.
blurpp:
	if(numpp>1) {
		clusind[0]=fmatrix.at(clus1).at(0);
		clusind[1]=fmatrix.at(clus2).at(0);
		for(int i=0; i<2; i++)
		{
			pos[0]=group.at(clusind[i]).x;
       		pos[1]=group.at(clusind[i]).y;
       		pos[2]=group.at(clusind[i]).z;
			crystalind[0]=group.at(clusind[i]).vid1;
            crystalind[1]=group.at(clusind[i]).vid2;
            crystalind[2]=group.at(clusind[i]).vid3;

       		PanelBlur(pos, blurpos,crystalind, numcolumn, numcolumnm,numrowr, numrowm, doires, blurtype);
       		group.at(clusind[i]).bx=blurpos[0];
       		group.at(clusind[i]).by=blurpos[1];
       		group.at(clusind[i]).bz=blurpos[2];
       		if(sqrt(pow(blurpos[0]-pos[0], 2) + pow(blurpos[1]-pos[1], 2) + pow(blurpos[2]-pos[2], 2))>2.6) cout<<crystalind[0]<<" "<<crystalind[1]<<" "<<crystalind[2]<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<blurpos[0]<<" "<<blurpos[1]<<" "<<blurpos[2]<<endl;
	
		}
		goto energycheck;
	}
	

	// Cannot find a pair of PP. Deal with MIPE. First check if there is 2A/1C or 1A/2C.
	for(int i=0; i< fmatrix.size(); i++){
		for(int j1 = 0; j1 < fmatrix.at(i).size() - 1; j1++){
			for(int j2 = j1 + 1; j2 < fmatrix.at(i).size(); j2++){
				hits element1 = group.at(fmatrix.at(i).at(j1));
				hits element2 = group.at(fmatrix.at(i).at(j2));
				if(element1.vid1 == element2.vid1 && element1.vid2 == element2.vid2 && (element1.vid3%numcolumnm == element2.vid3%numcolumnm || element1.vid3/numcolumnm == element2.vid3/numcolumnm)) {
					plorf->failtype=8; 
					*(failnum+7) += 1;
				   	return plorf;
				}
			}
		}
	}

	//blur position of all events in the group in order to perform compton kinematic calculation
	for(int i=0; i<group.size(); i++)
	{
		pos[0]=group.at(i).x;
		pos[1]=group.at(i).y;
		pos[2]=group.at(i).z;
        crystalind[0]=group.at(i).vid1;
        crystalind[1]=group.at(i).vid2;
        crystalind[2]=group.at(i).vid3;
		PanelBlur(pos, blurpos, crystalind, numcolumn, numcolumnm, numrowr, numrowm, doires, blurtype);
		group.at(i).bx=blurpos[0];
		group.at(i).by=blurpos[1];
		group.at(i).bz=blurpos[2];

		if(sqrt(pow(blurpos[0]-pos[0], 2) + pow(blurpos[1]-pos[1], 2) + pow(blurpos[2]-pos[2], 2))>2.6) cout<<crystalind[0]<<" "<<crystalind[1]<<" "<<crystalind[2]<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<blurpos[0]<<" "<<blurpos[1]<<" "<<blurpos[2]<<endl;
	}

	//cout<<pos[0]<<" "<<blurpos[0]<<endl;

	for(int i=0; i<fmatrix.size() -1 ; i++)
	{
		for(int j=1+i; j<fmatrix.size();j++)
		{
			secdiff = abs(group.at(fmatrix.at(i).at(0)).vid1 - group.at(fmatrix.at(j).at(0)).vid1);
			//if(secdiff < otherpara.minsecdiff || secdiff > (otherpara.totnumsec - otherpara.minsecdiff)) continue;
			if((otherpara.checksecdiff==0 && secdiff >=1) || (otherpara.checksecdiff==1 && secdiff >= otherpara.minsecdiff && secdiff <= (otherpara.totnumsec - otherpara.minsecdiff)))
			{
			//	if(fmatrix.at(i).size()==1 && fmatrix.at(j).size()==1)
			//	{
			//		ddamin=0.;
			//		clus1=i;
			//		clus2=j;
			//		elem1=0;
			//		elem2=0;
			//		goto energycheck;
			//	}
				if(fmatrix.at(i).size()==1 && fmatrix.at(j).size()==1) cout<<"Something seems wrong, as PP event has been pre-filtered!!"<<endl;

				else if(fmatrix.at(i).size()==1)
				{
					for(int kt1=0; kt1<fmatrix.at(j).size(); kt1++)		//kt1, kt2 are indexes of assumed first event and second event in a cluster.
					{
						for(int kt2=0; kt2<fmatrix.at(j).size(); kt2++)
						{
							if(kt1==kt2) continue;
							thetae=Etheta(group.at(fmatrix.at(j).at(kt1)).blurenergy);
							pos1[0]=group.at(fmatrix.at(i).at(0)).bx;
							pos1[1]=group.at(fmatrix.at(i).at(0)).by;
							pos1[2]=group.at(fmatrix.at(i).at(0)).bz;
							pos2[0]=group.at(fmatrix.at(j).at(kt1)).bx;
							pos2[1]=group.at(fmatrix.at(j).at(kt1)).by;
							pos2[2]=group.at(fmatrix.at(j).at(kt1)).bz;
							pos3[0]=group.at(fmatrix.at(j).at(kt2)).bx;
							pos3[1]=group.at(fmatrix.at(j).at(kt2)).by;
							pos3[2]=group.at(fmatrix.at(j).at(kt2)).bz;
							thetap=Ptheta(pos1, pos2, pos3, poverlap);
							dda = abs(thetae - thetap);
							if(dda < ddamin)
							{
								clus1=i;
								clus2=j;
								elem1=0;
								elem2=kt1;
								ddamin=dda;
                                ifCrossStrip = 0;
							}

							if(group[fmatrix[j][kt1]].vid1 == group[fmatrix[j][kt2]].vid1 && group[fmatrix[j][kt1]].vid2 == group[fmatrix[j][kt2]].vid2){
							thetae=Etheta(group.at(fmatrix.at(j).at(kt1)).blurenergy);
							pos1[0]=group.at(fmatrix.at(i).at(0)).bx;
							pos1[1]=group.at(fmatrix.at(i).at(0)).by;
							pos1[2]=group.at(fmatrix.at(i).at(0)).bz;
							pos2[0]=group.at(fmatrix.at(j).at(kt2)).bx;
							pos2[1]=group.at(fmatrix.at(j).at(kt1)).by;
							pos2[2]=group.at(fmatrix.at(j).at(kt2)).bz;
							pos3[0]=group.at(fmatrix.at(j).at(kt1)).bx;
							pos3[1]=group.at(fmatrix.at(j).at(kt2)).by;
							pos3[2]=group.at(fmatrix.at(j).at(kt1)).bz;
							thetap=Ptheta(pos1, pos2, pos3, poverlap);
							dda = abs(thetae - thetap);
							if(dda < ddamin)
							{
								clus1=i;
								clus2=j;
								elem1=0;
								elem2=kt1;
								ddamin=dda;
                                ifCrossStrip = 2;
                                //crossStripAdd = createNew(group.at(fmatrix.at(j).at(kt1)), group.at(fmatrix.at(j).at(kt2)));

							}
							}

						}
					}
				}

				else if(fmatrix.at(j).size()==1)
				{
                    for(int kt1=0; kt1<fmatrix.at(i).size(); kt1++)     //kt1, kt2 are indexes of assumed first event and second event in a cluster.
                    {
                        for(int kt2=0; kt2<fmatrix.at(i).size(); kt2++)
                        {
							if(kt1==kt2) continue;
                            thetae=Etheta(group.at(fmatrix.at(i).at(kt1)).blurenergy);
                            pos1[0]=group.at(fmatrix.at(j).at(0)).bx;
                            pos1[1]=group.at(fmatrix.at(j).at(0)).by;
                            pos1[2]=group.at(fmatrix.at(j).at(0)).bz;
                            pos2[0]=group.at(fmatrix.at(i).at(kt1)).bx;
                            pos2[1]=group.at(fmatrix.at(i).at(kt1)).by;
                            pos2[2]=group.at(fmatrix.at(i).at(kt1)).bz;
                            pos3[0]=group.at(fmatrix.at(i).at(kt2)).bx;
                            pos3[1]=group.at(fmatrix.at(i).at(kt2)).by;
                            pos3[2]=group.at(fmatrix.at(i).at(kt2)).bz;
                            thetap=Ptheta(pos1, pos2, pos3, poverlap);
                            dda = abs(thetae - thetap);
                            if(dda < ddamin)
                            {
                                clus1=i;
                                clus2=j;
                                elem1=kt1;
                                elem2=0;
								ddamin=dda;
                                ifCrossStrip = 0;
                            }

							if(group[fmatrix[i][kt1]].vid1 == group[fmatrix[i][kt2]].vid1 && group[fmatrix[i][kt1]].vid2 == group[fmatrix[i][kt2]].vid2){
                            thetae=Etheta(group.at(fmatrix.at(i).at(kt1)).blurenergy);
                            pos1[0]=group.at(fmatrix.at(j).at(0)).bx;
                            pos1[1]=group.at(fmatrix.at(j).at(0)).by;
                            pos1[2]=group.at(fmatrix.at(j).at(0)).bz;
                            pos2[0]=group.at(fmatrix.at(i).at(kt2)).bx;
                            pos2[1]=group.at(fmatrix.at(i).at(kt1)).by;
                            pos2[2]=group.at(fmatrix.at(i).at(kt2)).bz;
                            pos3[0]=group.at(fmatrix.at(i).at(kt1)).bx;
                            pos3[1]=group.at(fmatrix.at(i).at(kt2)).by;
                            pos3[2]=group.at(fmatrix.at(i).at(kt1)).bz;
                            thetap=Ptheta(pos1, pos2, pos3, poverlap);
                            dda = abs(thetae - thetap);
                            if(dda < ddamin)
                            {
                                clus1=i;
                                clus2=j;
                                elem1=kt1;
                                elem2=0;
								ddamin=dda;
                                ifCrossStrip = 1;
                            }
							}

                        }                                                                                                                                                 
                    }
	
				}

				else
				{
				
					for(int ko1=0; ko1<fmatrix.at(i).size(); ko1++)
					{
						for(int ko2=0; ko2<fmatrix.at(i).size();ko2++)
						{
							if(ko1==ko2) continue;
							for(int kt1=0; kt1<fmatrix.at(j).size(); kt1++)
							{
								for(int kt2=0; kt2<fmatrix.at(j).size(); kt2++)
								{
									if(kt1 == kt2) continue;

									bool sC1 = (group[fmatrix[i][ko1]].vid1 == group[fmatrix[i][ko2]].vid1 && group[fmatrix[i][ko1]].vid2 == group[fmatrix[i][ko2]].vid2);
									bool sC2 = (group[fmatrix[j][kt1]].vid1 == group[fmatrix[j][kt2]].vid1 && group[fmatrix[j][kt1]].vid2 == group[fmatrix[j][kt2]].vid2);
									//first calculate dda on ko side
									thetae=Etheta(group.at(fmatrix.at(i).at(ko1)).blurenergy);
                           			pos1[0]=group.at(fmatrix.at(j).at(kt1)).bx;
                           			pos1[1]=group.at(fmatrix.at(j).at(kt1)).by;
                           			pos1[2]=group.at(fmatrix.at(j).at(kt1)).bz;
                           			pos2[0]=group.at(fmatrix.at(i).at(ko1)).bx;
                           			pos2[1]=group.at(fmatrix.at(i).at(ko1)).by;
                           			pos2[2]=group.at(fmatrix.at(i).at(ko1)).bz;
                           			pos3[0]=group.at(fmatrix.at(i).at(ko2)).bx;
                           			pos3[1]=group.at(fmatrix.at(i).at(ko2)).by;
                           			pos3[2]=group.at(fmatrix.at(i).at(ko2)).bz;
                           			thetap=Ptheta(pos1, pos2, pos3, poverlap);
                           			dda1 = abs(thetae - thetap);

                                    //second calculate dda on kt side
									thetae=Etheta(group.at(fmatrix.at(j).at(kt1)).blurenergy);
                                    pos1[0]=group.at(fmatrix.at(i).at(ko1)).bx;
                                    pos1[1]=group.at(fmatrix.at(i).at(ko1)).by;
                                    pos1[2]=group.at(fmatrix.at(i).at(ko1)).bz;
                                    pos2[0]=group.at(fmatrix.at(j).at(kt1)).bx;
                                    pos2[1]=group.at(fmatrix.at(j).at(kt1)).by;
                                    pos2[2]=group.at(fmatrix.at(j).at(kt1)).bz;                                                             
                                    pos3[0]=group.at(fmatrix.at(j).at(kt2)).bx;
                                    pos3[1]=group.at(fmatrix.at(j).at(kt2)).by;
                                    pos3[2]=group.at(fmatrix.at(j).at(kt2)).bz;
                                    thetap=Ptheta(pos1, pos2, pos3, poverlap);
                                    dda2 = abs(thetae - thetap);

									//calculate average
									dda= (dda1 + dda2) / 2.;
									if(dda< ddamin)
									{
										clus1=i;
										clus2=j;
										elem1=ko1;
										elem2=kt1;
										ddamin=dda;
                                        ifCrossStrip = 0;
									}

									if(sC1){
									//first calculate dda on ko side. change ko coordinate
									thetae=Etheta(group.at(fmatrix.at(i).at(ko1)).blurenergy);
                           			pos1[0]=group.at(fmatrix.at(j).at(kt1)).bx;
                           			pos1[1]=group.at(fmatrix.at(j).at(kt1)).by;
                           			pos1[2]=group.at(fmatrix.at(j).at(kt1)).bz;
                           			pos2[0]=group.at(fmatrix.at(i).at(ko2)).bx;
                           			pos2[1]=group.at(fmatrix.at(i).at(ko1)).by;
                           			pos2[2]=group.at(fmatrix.at(i).at(ko2)).bz;
                           			pos3[0]=group.at(fmatrix.at(i).at(ko1)).bx;
                           			pos3[1]=group.at(fmatrix.at(i).at(ko2)).by;
                           			pos3[2]=group.at(fmatrix.at(i).at(ko1)).bz;
                           			thetap=Ptheta(pos1, pos2, pos3, poverlap);
                           			dda1 = abs(thetae - thetap);

                                    //second calculate dda on kt side
									thetae=Etheta(group.at(fmatrix.at(j).at(kt1)).blurenergy);
                                    pos1[0]=group.at(fmatrix.at(i).at(ko2)).bx;
                                    pos1[1]=group.at(fmatrix.at(i).at(ko1)).by;
                                    pos1[2]=group.at(fmatrix.at(i).at(ko2)).bz;
                                    pos2[0]=group.at(fmatrix.at(j).at(kt1)).bx;
                                    pos2[1]=group.at(fmatrix.at(j).at(kt1)).by;
                                    pos2[2]=group.at(fmatrix.at(j).at(kt1)).bz;                                                             
                                    pos3[0]=group.at(fmatrix.at(j).at(kt2)).bx;
                                    pos3[1]=group.at(fmatrix.at(j).at(kt2)).by;
                                    pos3[2]=group.at(fmatrix.at(j).at(kt2)).bz;
                                    thetap=Ptheta(pos1, pos2, pos3, poverlap);
                                    dda2 = abs(thetae - thetap);

									//calculate average
									dda= (dda1 + dda2) / 2.;
									if(dda< ddamin)
									{
										clus1=i;
										clus2=j;
										elem1=ko1;
										elem2=kt1;
										ddamin=dda;
                                        ifCrossStrip = 41;
									}
									}

									if(sC2){
									//first calculate dda on ko side. change kt coordinate
									thetae=Etheta(group.at(fmatrix.at(i).at(ko1)).blurenergy);
                           			pos1[0]=group.at(fmatrix.at(j).at(kt2)).bx;
                           			pos1[1]=group.at(fmatrix.at(j).at(kt1)).by;
                           			pos1[2]=group.at(fmatrix.at(j).at(kt2)).bz;
                           			pos2[0]=group.at(fmatrix.at(i).at(ko1)).bx;
                           			pos2[1]=group.at(fmatrix.at(i).at(ko1)).by;
                           			pos2[2]=group.at(fmatrix.at(i).at(ko1)).bz;
                           			pos3[0]=group.at(fmatrix.at(i).at(ko2)).bx;
                           			pos3[1]=group.at(fmatrix.at(i).at(ko2)).by;
                           			pos3[2]=group.at(fmatrix.at(i).at(ko2)).bz;
                           			thetap=Ptheta(pos1, pos2, pos3, poverlap);
                           			dda1 = abs(thetae - thetap);

                                    //second calculate dda on kt side
									thetae=Etheta(group.at(fmatrix.at(j).at(kt1)).blurenergy);
                                    pos1[0]=group.at(fmatrix.at(i).at(ko2)).bx;
                                    pos1[1]=group.at(fmatrix.at(i).at(ko1)).by;
                                    pos1[2]=group.at(fmatrix.at(i).at(ko2)).bz;
                                    pos2[0]=group.at(fmatrix.at(j).at(kt2)).bx;
                                    pos2[1]=group.at(fmatrix.at(j).at(kt1)).by;
                                    pos2[2]=group.at(fmatrix.at(j).at(kt2)).bz;                                                             
                                    pos3[0]=group.at(fmatrix.at(j).at(kt1)).bx;
                                    pos3[1]=group.at(fmatrix.at(j).at(kt2)).by;
                                    pos3[2]=group.at(fmatrix.at(j).at(kt1)).bz;
                                    thetap=Ptheta(pos1, pos2, pos3, poverlap);
                                    dda2 = abs(thetae - thetap);

									//calculate average
									dda= (dda1 + dda2) / 2.;
									if(dda< ddamin)
									{
										clus1=i;
										clus2=j;
										elem1=ko1;
										elem2=kt1;
										ddamin=dda;
                                        ifCrossStrip = 42;
									}
									}

									if(sC1 && sC2){
									//first calculate dda on ko side. change kt and ko coordinate
									thetae=Etheta(group.at(fmatrix.at(i).at(ko1)).blurenergy);
                           			pos1[0]=group.at(fmatrix.at(j).at(kt2)).bx;
                           			pos1[1]=group.at(fmatrix.at(j).at(kt1)).by;
                           			pos1[2]=group.at(fmatrix.at(j).at(kt2)).bz;
                           			pos2[0]=group.at(fmatrix.at(i).at(ko2)).bx;
                           			pos2[1]=group.at(fmatrix.at(i).at(ko1)).by;
                           			pos2[2]=group.at(fmatrix.at(i).at(ko2)).bz;
                           			pos3[0]=group.at(fmatrix.at(i).at(ko1)).bx;
                           			pos3[1]=group.at(fmatrix.at(i).at(ko2)).by;
                           			pos3[2]=group.at(fmatrix.at(i).at(ko1)).bz;
                           			thetap=Ptheta(pos1, pos2, pos3, poverlap);
                           			dda1 = abs(thetae - thetap);

                                    //second calculate dda on kt side
									thetae=Etheta(group.at(fmatrix.at(j).at(kt1)).blurenergy);
                                    pos1[0]=group.at(fmatrix.at(i).at(ko2)).bx;
                                    pos1[1]=group.at(fmatrix.at(i).at(ko1)).by;
                                    pos1[2]=group.at(fmatrix.at(i).at(ko2)).bz;
                                    pos2[0]=group.at(fmatrix.at(j).at(kt2)).bx;
                                    pos2[1]=group.at(fmatrix.at(j).at(kt1)).by;
                                    pos2[2]=group.at(fmatrix.at(j).at(kt2)).bz;                                                             
                                    pos3[0]=group.at(fmatrix.at(j).at(kt1)).bx;
                                    pos3[1]=group.at(fmatrix.at(j).at(kt2)).by;
                                    pos3[2]=group.at(fmatrix.at(j).at(kt1)).bz;
                                    thetap=Ptheta(pos1, pos2, pos3, poverlap);
                                    dda2 = abs(thetae - thetap);

									//calculate average
									dda= (dda1 + dda2) / 2.;
									if(dda< ddamin)
									{
										clus1=i;
										clus2=j;
										elem1=ko1;
										elem2=kt1;
										ddamin=dda;
                                        ifCrossStrip = 44;
									}
									}

								}
							}
						}
					}
				}


			}
		}
	}

	if(ddamin == 1000.) {plorf->failtype=7; *(failnum+6)+=1; return plorf;}

energycheck:
	//int low1=0,low2=0,high1=0,high2=0;
	float totale1=0., totale2=0.;

	for(int i=0; i< fmatrix.at(clus1).size(); i++) totale1+=group.at(fmatrix.at(clus1).at(i)).blurenergy;

	for(int i=0; i< fmatrix.at(clus2).size(); i++) totale2+=group.at(fmatrix.at(clus2).at(i)).blurenergy;

	//if(totale1 < Elow || totale2 < Elow) {plorf->failtype=3; *(failnum+2)+=1; return plorf;}

	//if(totale1 > Ehigh || totale2 > Ehigh) {plorf->failtype=4; *(failnum+3)+=1; return plorf;}

	
	//check if overlap is existing
	if(overl>0) *pnumoverl+=1;

	//check if clustering is correct. First, inside same cluster, eventid and photonid should be same. Second: between clusters, if eventid is same, then photonid should be different.
	if(otherpara.dbinfo==1)
	{
		checkclus(group, ematrix, 1);
	}


	//check if misaligned
	if(ematrix.at(emindex.at(clus1)).at(0) != fmatrix.at(clus1).at(0) || ematrix.at(emindex.at(clus2)).at(0) != fmatrix.at(clus2).at(0)) plorf->Misaln=1;
	else plorf->Misaln=0;

	//check if first event is incorrectly identified
	if(plorf->Misaln==0) 
	{
		if(elem1 != 0 || elem2 != 0 || ifCrossStrip != 0) plorf->Incorrseq=1;
		else plorf->Incorrseq=0;
		if(ifCrossStrip != 0) plorf->IncorrseqCS = 1;
		else plorf->IncorrseqCS = 0;
	}

	//if not failed, fill plorf
    //hits crossStripAdd; // New hits element due to position ambiguity in cross strip detector.
    if(ifCrossStrip == 0){
	    plorf->hits1=group.at(fmatrix.at(clus1).at(elem1));
	    plorf->hits2=group.at(fmatrix.at(clus2).at(elem2));
    }
    else if(ifCrossStrip == 1 || ifCrossStrip == 41){
        int num1 = fmatrix.at(clus1).at(elem1);
        int num2 = fmatrix.at(clus1).at((elem1 +1) % 2);
        plorf->hits1=newCrossStrip(group.at(num1), group.at(num2), numcolumnm);
    	plorf->hits2=group.at(fmatrix.at(clus2).at(elem2));
    }
    else if(ifCrossStrip == 2 || ifCrossStrip == 42){
        int num1 = fmatrix.at(clus2).at(elem2);
        int num2 = fmatrix.at(clus2).at((elem2 +1) % 2);
        plorf->hits2=newCrossStrip(group.at(num1), group.at(num2), numcolumnm);
    	plorf->hits1=group.at(fmatrix.at(clus1).at(elem1));
    }
    else if(ifCrossStrip == 44){
        int num1 = fmatrix.at(clus1).at(elem1);
        int num2 = fmatrix.at(clus1).at((elem1 +1) % 2);
        plorf->hits1=newCrossStrip(group.at(num1), group.at(num2), numcolumnm);
    	
        num1 = fmatrix.at(clus2).at(elem2);
        num2 = fmatrix.at(clus2).at((elem2 +1) % 2);
        plorf->hits2=newCrossStrip(group.at(num1), group.at(num2), numcolumnm);
    }
    else cout<< "Unknown indicator for assigning events to LOR!!"<<endl;

	plorf->failtype=0;
	plorf->dda=ddamin;
	plorf->totale1=totale1;
	plorf->totale2=totale2;

	return plorf;

}

float Etheta(float energy)
{
	float e1, thetae;
	e1=511. - energy;
	if(e1 < 171.) e1=171.;
	thetae=acos(1. - 511. * (1./e1 - 1./511.))* 180.0 / PI;
	return thetae;
}

float Ptheta(float pos1[], float pos2[], float pos3[], int * poverlap)
{
	float thetap, c1=0., c2=0., c3=0.;
	float dpos1[3], dpos2[3];
	for(int i=0; i<3; i++)
	{
		dpos1[i]=pos2[i] - pos1[i];
		dpos2[i]=pos3[i] - pos2[i];
		c1 += dpos1[i] * dpos2[i];
		c2 += pow(dpos1[i],2);
		c3 += pow(dpos2[i], 2);
	}
	if(c3==0 || c2==0 ) 
	{
		//cout<<"Warning: In compton kinematics, at least two points are the same."<<endl;
		*poverlap += 1;
		return 0.;
	}
	thetap=acos(c1 / sqrt(c2 * c3))* 180.0 / PI;
	return thetap;
}

//blur the position of event. Assign a random location inside a crystal to the location of event. Consider DOI. Distribution of location inside crystal is Uniform distribution. Some numbers need modifying based on a specific system.
int PanelBlur(float pos[], float * pblur, int crystal[], int numcr, int numcm, int numrowr, int numrowm, float doires, string blurtype)
{
    float xdm, ydm, zdm,xds, yds, zds, xd, yd, zd, xmea, ymea, zmea,xreal, yreal, zreal, xdreal, xdmea;
	int nrse, nmod, nsub;
	int mdm[2], sdm[2];
	float theta;
	float numcrhalf, numcmhalf;
	static unsigned seed = 300;
	static default_random_engine generator (seed);
	uniform_real_distribution<double> distribution(-1.0,1.0);

	xreal=pos[0];
	yreal=pos[1];
	zreal=pos[2];
	nrse=crystal[0];
	nmod=crystal[1];
	nsub=crystal[2];
	mdm[0]=nmod/numcr;	// z
	mdm[1]=nmod % numcr;	// y
    sdm[0]=nsub/numcm;	// y
	sdm[1]=nsub % numcm;	// x
	numcrhalf=(numcr-1)/2.;
	numcmhalf=(numcm-1)/2.;
	float numrowrhalf = (numrowr-1)/2.0;
	float numrowmhalf = (numrowm-1)/2.0;
	//cout<<nrse<<endl;

	if(blurtype=="bin2")
	{
         xdm = 0.0;
		 ydm = (mdm[1] - numcrhalf) * 40.1;
		 zdm = (mdm[0] - numrowrhalf) * 5.1;
         yds = (sdm[0] - numrowmhalf) * 1.0;
         xds = (sdm[1] - numcmhalf) * 5.0;

         if(nrse==0)
		 {
			xmea = 120.0 + xds;
			ymea = ydm + yds;
			zmea = zdm - 0.5 * 5.0 + (trunc((zreal-zdm + 0.5 * 5.0)/doires)+0.5)*doires;
		 }
		 else if(nrse==1)
		 {
			xmea = -(120.0 + xds);
			ymea = -(ydm + yds);
			zmea = zdm - 0.5 * 5.0 + (trunc((zreal-zdm + 0.5 * 5.0)/doires)+0.5)*doires;

		 }
	}


	else if(blurtype=="bin31")
	{
        ydm=mdm[1] * 9 - numcrhalf * 9;
        zdm=mdm[0] * 9 - numcrhalf * 9;
        yds=sdm[1] * 2.2 - numcmhalf * 2.2;
        zds=sdm[0] * 2.2 - numcmhalf * 2.2;
        yd=ydm + yds + distribution(generator);
        zd=zdm + zds + distribution(generator);
		//cout<<yd<<" "<<zd<<" "<<endl;

        if (nrse==0)
		{
            xmea=(trunc((xreal - 100.)/doires)+0.5)*doires +100.+ doires/2. * distribution(generator);
            ymea=yd;
            zmea=zd;
			//cout<<xmea<<" "<<ymea<<" "<<zmea<<endl;
		}
		else if(nrse==1)
		{
            theta=28. * PI / 180.;
            xdreal=-(xreal + 61.03)*sin(theta)+(yreal - 114.78)*cos(theta);
            xdmea=(trunc((xdreal+10.)/doires)+0.5)*doires -10.+ doires/2. * distribution(generator);
			//cout<<xdmea<<endl;
            xmea=-yd * cos(theta)-xdmea * sin(theta) - 61.03;
            ymea=-yd * sin(theta) +xdmea * cos(theta) + 114.78;
            zmea=zd;
		}
		else if(nrse==2)
		{
            theta=28. * PI / 180.;
            xdreal=-(xreal + 61.03)*sin(theta)-(yreal + 114.78)*cos(theta);
            xdmea=(trunc((xdreal+10.)/doires)+0.5)*doires -10.+ doires/2. * distribution(generator);
            //cout<<xdmea<<endl;
			xmea=yd * cos(theta) -xdmea * sin(theta)  - 61.03;
            ymea=-yd * sin(theta) -xdmea * cos(theta) - 114.78;
            zmea=zd;
		}
	}


	else if(blurtype=="bin32")
	{
        ydm=mdm[1] * 9 - 11.5 * 9;
        zdm=mdm[0] * 9 - 11.5 * 9;
        yds=sdm[1] * 2.2 - 1.5 * 2.2;
        zds=sdm[0] * 2.2 - 1.5 * 2.2;
        yd=ydm + yds + distribution(generator);
        zd=zdm + zds + distribution(generator);
        if(nrse==0)
		{
            xmea=(trunc((xreal-100.)/doires)+0.5)*doires +100.+ doires/2. * distribution(generator);
            ymea=yd;
            zmea=zd;
		}
		else if(nrse==1)
		{
            xmea=-yd - 48.;
            ymea=(trunc((yreal-120.)/doires)+0.5)*doires +120.+ doires/2. * distribution(generator);
            zmea=zd;
		}
		else if(nrse==2)
		{
            xmea=yd - 48.;
            ymea=(trunc((-yreal-120.)/doires)+0.5)*(-doires)-120.+ doires/2. * distribution(generator);
            zmea=zd;
		}
	}


	else if(blurtype=="bin4")
	{
        ydm=mdm[1] * 9 - 11.5 * 9;
        zdm=mdm[0] * 9 - 11.5 * 9;
        yds=sdm[1] * 2.2 - 1.5 * 2.2;
        zds=sdm[0] * 2.2 - 1.5 * 2.2;
        yd=ydm + yds + distribution(generator);
        zd=zdm + zds + distribution(generator);
        if(nrse==0)
		{
            xmea=(trunc((xreal-100.)/doires)+0.5)*doires +100.+ doires/2. * distribution(generator);
            ymea=yd;
            zmea=zd;
		}
		else if(nrse==1)
		{
            xmea=-yd - 8.;
            ymea=(trunc((yreal-120.)/doires)+0.5)*doires +120.+ doires/2. * distribution(generator);
            zmea=zd;
		}
		else if(nrse==2)
		{
            xmea=yd - 8.;
            ymea=(trunc((-yreal-120.)/doires)+0.5)*(-doires)-120.+ doires/2. * distribution(generator);
            zmea=zd;
		}
        else if(nrse==3)
		{
            xmea=(trunc((xreal+155.)/doires)+0.5)*doires - 155.+ doires/2. * distribution(generator);
            ymea=yd;
            zmea=zd;
		}
	}

	
	else
	{	
		cout<<"Warning: Unknown blurring type!! Will use precise position as blurred position"<<endl;
		xmea=xreal;
		ymea=yreal;
		zmea=zreal;
	}

	*pblur= xmea;
	*(pblur+ 1) = ymea;
	*(pblur +2 ) = zmea;
	//cout<<xmea<<" "<<ymea<<" "<<zmea<<endl;
	return 0;
}

int checkclus(vector<hits> group, vector<vector<int>> ematrix, int kk)
{
    int eventid=-1;
    int pid=-1;
    int sameclus= 1;
    int diffclus = 1;

    for(int i=0;i<ematrix.size(); i++)
    {
        debugclus[kk].totclus += 1;
        eventid=-1;
        pid=-1;
        sameclus=1;
        for(int j=0; j<ematrix.at(i).size(); j++)
        {
            if(eventid==-1) {eventid=group.at(ematrix.at(i).at(j)).eventid; pid=group.at(ematrix.at(i).at(j)).particleid;}
            else
            {
                if(group.at(ematrix.at(i).at(j)).eventid != eventid || group.at(ematrix.at(i).at(j)).particleid != pid)
                {
                    sameclus=0;
                    break;
                }
            }
        }
        if(sameclus==1)
        {
            debugclus[kk].clussamep += 1;
        }
    }

    debugclus[kk].totgroup += 1;
    for(int i=0;i<ematrix.size()-1; i++)
    {
        for(int j=0;j<ematrix.at(i).size(); j++)
        {
            for(int ii=i+1; ii< ematrix.size(); ii++)
            {
                for(int jj=0; jj<ematrix.at(ii).size(); jj++)
                {
                    
                    if(group.at(ematrix.at(i).at(j)).eventid == group.at(ematrix.at(ii).at(jj)).eventid && group.at(ematrix.at(i).at(j)).particleid == group.at(ematrix.at(ii).at(jj)).particleid)
                    {
                        diffclus=0;
                        goto groupcluscheck;

                    }
                }
            }                                                                                                                                                         
        }
    }
groupcluscheck:
    if(diffclus==1) debugclus[kk].clusright+=1;

	return 0;
}

hits newCrossStrip(hits event1, hits event2, int numcolumnm){
	hits newEvent;
	newEvent = event1;
	newEvent.x = event2.x;
	newEvent.z = event2.z;
	newEvent.bx = event2.bx;
	newEvent.bz = event2.bz;
	newEvent.vid3 = (event2.vid3)%numcolumnm + ((event1.vid3)/numcolumnm) * numcolumnm;
	
	return newEvent;
}

vector<vector<int>> EnergyCheck(const vector<vector<int>>& fmatrix, const vector<hits>& group, const float &Elow, const float &Ehigh, int &numclusafter, int &maxnumevent,const vector<int> &emindex, vector<int> &nemindex){
	vector<vector<int>> withinEW;
	numclusafter = 0;
	maxnumevent = 0;
	nemindex.clear();

	for(int i=0; i< fmatrix.size(); i++){
		float totalE = 0.0f;
		int numevent = 0;
		for(int j = 0; j< fmatrix[i].size(); j++) {totalE += group[fmatrix[i][j]].blurenergy; numevent += 1;}
		if(totalE > Elow && totalE < Ehigh) {
			withinEW.push_back(fmatrix[i]);
			numclusafter += 1;
			if(maxnumevent < numevent) maxnumevent = numevent;
			nemindex.push_back(emindex[i]);
		}
	}
	return withinEW;
}


