#ifndef XCOM_H
#define XCOM_H 1

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

using namespace std;

#define ME 1500
#define MEA 600
#define MEB 108
#define AVOG 0.60221367
#define EPAIR1 1.022007E+06
#define EPAIR2 2.044014E+06
#define CROSS_SECTION_TYPES 7

const vector<string> OutputNames{
    "PhotoelectricEdge","PhotonEnergy","CoherentScattering",
    "IncoherentScattering","PhotoelectricAbsorption","NuclearFieldPairProduction",
    "ElectronFieldPairProduction","TotalWithCoherent","TotalWithoutCoherent"
    };
//  ATWTS: atomic mass
const double ATWTS[100] = {
    1.00794,       4.002602,      6.941,         9.012182,
    10.811,        12.011,        14.00674,      15.9994,
    18.9984032,    20.1797,       22.989768,     24.3050,
    26.981539,     28.0855,       30.973762,     32.066,
    35.4527,       39.948,        39.0983,       40.078,
    44.955910,     47.88,         50.9415,       51.9961,
    54.93805,      55.847,        58.93320,      58.69,
    63.546,        65.39,         69.723,        72.61,
    74.92159,      78.96,         79.904,        83.80,
    85.4678,       87.62,         88.90585,      91.224,
    92.90638,      95.94,         97.9072,      101.07,
    102.9055,      106.42,        107.8682,      112.411,
    114.82,        118.710,       121.75,        127.60,
    126.90447,     131.29,        132.90543,     137.327,
    138.9055,      140.115,       140.90765,     144.24,
    144.9127,      150.36,        151.965,       157.25,
    158.92534,     162.50,        164.93032,     167.26,
    168.93421,     173.04,        174.967,       178.49,
    180.9479,      183.85,        186.207,       190.2,
    192.22,        195.08,        196.96654,     200.59,
    204.3833,      207.2,         208.98037,     208.9824,
    209.9871,      222.0176,      223.0197,      226.0254,
    227.0278,      232.0381,      231.03588,     238.0289,
    237.0482,      239.0522,      243.0614,      247.0703,
    247.0703,      251.0796,      252.083,       257.0951
};

// periodic table

const  map<string, int> PeriodicTable = {
{"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"Ne",10},
{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},{"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},
{"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},{"Cu",29},{"Zn",30},
{"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},{"Rb",37},{"Sr",38},{"Y",39},{"Zr",40},
{"Nb",41},{"Mo",42},{"Tc",43},{"Ru",44},{"Rh",45},{"Pd",46},{"Ag",47},{"Cd",48},{"In",49},{"Sn",50},
{"Sb",51},{"Te",52},{"I",53},{"Xe",54},{"Cs",55},{"Ba",56},{"La",57},{"Ce",58},{"Pr",59},{"Nd",60},
{"Pm",61},{"Sm",62},{"Eu",63},{"Gd",64},{"Tb",65},{"Dy",66},{"Ho",67},{"Er",68},{"Tm",69},{"Yb",70},
{"Lu",71},{"Hf",72},{"Ta",73},{"W",74},{"Re",75},{"Os",76},{"Ir",77},{"Pt",78},{"Au",79},{"Hg",80},
{"Tl",81},{"Pb",82},{"Bi",83},{"Po",84},{"At",85},{"Rn",86},{"Fr",87},{"Ra",88},{"Ac",89},{"Th",90},
{"Pa",91},{"U",92},{"Np",93},{"Pu",94},{"Am",95},{"Cm",96},{"Bk",97},{"Cf",98},{"Es",99},{"Fm",100},
{"Md",101},{"No",102},{"Lr",103},{"Rf",104},{"Db",105},{"Sg",106},{"Bh",107},{"Hs",108},{"Mt",109},{"Ds",110},
{"Rg",111},{"Cn",112},{"Nh",113},{"Fl",114},{"Mc",115},{"Lv",116},{"Ts",117},{"Og",118}
};

// read MDATX3 file par
typedef struct MDATX3{
	int IDG[14],KMX[14];
	float EDGEN[14];
	string ADG[14];
	float ENG[14][35],PHC[14][35];
	float E[MEB],SCATCO[MEB],SCATIN[MEB];
	float PHOT[MEB],PAIRAT[MEB],PAIREL[MEB];
	float ATWT;
	int IZ,MAXEDG,MAXE,LAX;
}MDATX3;

// read MDATX3 datebase
int ReadMDATX3(char *file,MDATX3 *data);

// Sorts into monotonically increasing order 
// sort energy list
// input parameter:
//    NMAX: size of E
//    E   : pointer to  list of energy
template <typename T>
void SORT(int NMAX, T *E);

// merges energy lists
template <typename T>
int MERGE(T* E1,int *K1,int *L1,int MMAX,T*  E2,int *K2,int *L2,int NMAX);

// Reverses the order of lists
template <typename T>
void REV(int NMAX,T *X);

//fits F as a function of X, and calculates cubic spline coefficients A,B,C and D
template <typename T>
int SCOF(T *X,T *F,const int NMAX,T *A,T *B,T *C,T *D);

// Evaluates cubic spline as function of S, to obtain fitted result G.
template <typename T>
int BSPOL(const T S,T *X,T *A,T *B,T *C,T *D,const int N,T &G);

// Linear interpolation routine
template <typename T>
int BLIN(const T S,T *X,T *Y,const int N,T &TY);

//////////////////////////////////////////////////////////////
// Read chemical formulas.
// input  :
//     W  :  chemical formulas
//     JZ :  list of atomic number
//     WT :  list of weight of element
int ParseFormulas(const char *W,int *JZ,float *WT);

//Read chemical formulas.
//  input   :
//  formula :  chemical formulas
//  Z       :  list of atomic number
//  Atoms   :  list of atoms of an element 
bool ParseChemicalFormula(string formula, vector<int>& Z, vector<double>& atoms);

//
bool GetFractionByWeight(vector<int>& Z, vector<double>& atoms, vector<double>& weight);

// Inititialize energy list 
// input parameter:
//   KMAX : number of elements
//   NZ   : list of atomic number of element
//   NEGO : flag of output status in energy list (1: default energy list(1 keV -- 100 GeV ) ; 2:default and add energy list; 3: only input energy list,)
//   JENG : number of add energy 
//   EAD  : list of add energy
//   EN   : list of energy wanted to output
//   KZ   : list of energy wanted to output flag
//   KM   : list of energy wanted to output flag
// return:  the number of output energy
int InitEnergyList(int KMAX,int *NZ,int NEGO,int JENG,float *EAD,float *EN,int *KZ,int *KM);

int GetEnergyList(int KMAX, int* NZ, int NEGO, int JENG, float* EAD, float* EN, int* KZ, int* KM);
// calculate the cross sections for the interactions of photons with element, compound or mixture
// input parameter:
//   KMAX  : number of element
//   NZ    : list of atomic number of element
//   WEIGHT: list of weight of element
//   NF    : Flag for select the unit of result (3:cm2/g) 2:b/atom in total cross
//   NEGO  : flag of output status in energy list (1: default energy list(1 keV -- 100 GeV ) ; 2:default and add energy list; 3: only input energy list,)
//   NENG  : number of output energy
//   EN    : list of energy
//   KZ    : list of shell flag of element
//   KM    : list
//   SCTCO : scattering with coherent
//   SCTIN : scattering with incoherent
//   PHT   : photo elctric
//   PRAT  : pair proudction in nuclear field
//   PREL  : pair proudction in electron field
//   PHDIF : shell jump energy
void Calculation(int KMAX,int *NZ,float *WEIGHT, int NF,int NEGO,int NENG,float *EN,int *KZ,int *KM,float *SCTCO,float *SCTIN,float *PHT,float *PRAT,float *PREL,float *PHDIF);

#endif
