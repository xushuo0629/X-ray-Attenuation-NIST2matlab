//**************************************************************************************************
//Calculate X-ray and Gamma ray Mass Attenuation coefficients on MATLAB platform with data from XCOM
//All rights reserved by ZhiZhen Zhang,Department of Engineering Physics,Tsinghua University
//Contact me: 
//  Email: 2245167804@qq.com
//  Phone: 13511055064
// 
//References: 
//[1] https://github.com/gsize/calMu
//[2] https://physics.nist.gov/cgi-bin/Xcom/xcom3_2
//**************************************************************************************************

#include "xcom.h"
#include "mex.h"

using namespace std;

template <typename T>
void CopyToStructMatrix(mxArray* dest, const char* fieldname, T* src, int length);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    //input and output energy units: keV
    //output attenuation: cm2/g
    if (nrhs == 0 || nrhs >= 3)
    {
        mexErrMsgTxt("Input Error: input parameters must be 1 or 2\n");
    }
    if (nlhs >= 2)
    {
        mexErrMsgTxt("Output Error: too much outputs\n");
    }

    //the first input must be an 1 by N char matrix or a 2 by N double matrix
    //if the first input is a 2 by N double matrix,the first row implies atomic number 
    //and the sencond row implies atoms in chemical formula
    std::string formula;
    size_t row = mxGetM(prhs[0]);
    size_t col = mxGetN(prhs[0]);
    vector<int> Z;
    vector<double> atoms;
    vector<double> weight;
    if (mxIsChar(prhs[0]))
    {
        //mexPrintf("Char Matrix row = %u, col = %u\n", row, col);
        if (row != 1u)
        {
            mexErrMsgTxt("Input Error : first input dimensions error and is not a 1 by N char matrix\n");
        }
        formula = mxArrayToString(prhs[0]);
        std::string temp = "formula = " + formula+"\n";
        mexPrintf(temp.c_str());
        if (!ParseChemicalFormula(formula, Z, atoms))
        {
            mexErrMsgTxt("Input Error: illegal chemical formula\n");
        }
        GetFractionByWeight(Z, atoms, weight);

    }
    else if (mxIsDouble(prhs[0]))
    {
        double* input0 = mxGetPr(prhs[0]);
        if (row != 2u)
        {
            mexErrMsgTxt("Input Error : first input dimension error and is not a 2 by N double matrix\n");
        }
        Z.resize(col);
        atoms.resize(col);
        weight.resize(col);
        for (int i = 0; i < col; ++i)
        {
            Z[i] = int(input0[i * row]);
            atoms[i] = input0[i * row + 1];
            if (Z[i] <= 0 || Z[i] > 100 || atoms[i] <= 0)
            {
                mexErrMsgTxt("Input Error : first input 2 by N double matrix has illegal number\n");
            }
            mexPrintf("Z[%d] = %f, atoms[%d] = %f\n", i, input0[i * row], i, input0[i * row + 1]);
        }
        GetFractionByWeight(Z, atoms, weight);
    }
    else
    {
        mexErrMsgTxt("Input Error : first input parameters must be a chemical formula that is 1 by N char matrix or a 2 by N double matrix\n");
    }
    int MMAX = static_cast<int>(Z.size());
    int* NZ = new int[MMAX];
    float* WT = new float[MMAX];

    for (int i = 0; i < MMAX; i++)
    {
        NZ[i] = Z[i];
        WT[i] = static_cast<float>(weight[i]);
    }

    for (int i = 0; i < MMAX; i++)
        mexPrintf("ele:%s MMAX: %d Z[%d]: %d   FractionByWeight[%d]: %f\n", formula.c_str(), MMAX, i, NZ[i], i, WT[i]);

    mexPrintf("-------------------------------------------------\n");

    int JENG = 1;
    int NEGO = 2; //3:only input energy list  2: show all energy list (default and adding) 
    int NF = 3;  //3: cm2/g
    int length = 0;
    int KMAX = MMAX;
    

    //输出结构体的总域名个数，包括CROSS_SECTION_TYPES个截面，1个能量列表和1个edge元素列表
    int nfields = CROSS_SECTION_TYPES + 2;
    const char** fieldnames = new const char* [nfields];
    for (int i = 0; i < nfields; ++i)
    {
        fieldnames[i] = OutputNames.at(i).c_str();
    }
    if (nlhs == 1)
    {
        plhs[0] = mxCreateStructMatrix(1, 1, nfields, fieldnames);
    }
    delete[] fieldnames;
    

    if (nrhs == 1)
    {
        JENG = 1;
        NEGO = 2; //3:only input energy list  2: show all energy list (default and adding) 
        NF = 3;  //3: cm2/g
        float EAD[600];
        float EN[1200];
        int KZ[1200], KM[1200];
        //float SCTCO[ME], SCTIN[ME], PHT[ME], PRAT[ME], PREL[ME], PHDIF[ME];
        EAD[0] = 511000.0;
        length = InitEnergyList(KMAX, NZ, NEGO, JENG, EAD, EN, KZ, KM);
        float* SCTCO = new float[length];
        float* SCTIN = new float[length];
        float* PHT = new float[length];
        float* PRAT = new float[length];
        float* PREL = new float[length];
        float* PHDIF = new float[length];
        Calculation(KMAX, NZ, WT, NF, NEGO, length, EN, KZ, KM, SCTCO, SCTIN, PHT, PRAT, PREL, PHDIF);
        
        if (nlhs == 0)
        {
            mexPrintf("          Photon      Scattering        Photo-     Pair production    Total attenuation\n");
            mexPrintf("          energy   coherent  incoher.  electric     in        in       with     without\n");
            mexPrintf("                                       absorption  nuclear electron  coherent  coherent\n");
            mexPrintf("                                                    field    field     scatt.    scatt.\n");
            mexPrintf("          (MeV)    (cm2/g)    (cm2/g)   (cm2/g)    (cm2/g)  (cm2/g)   (cm2/g)   (cm2/g)\n");
            for (int j = 0; j < length; j++)
            {
                float crossAll = SCTCO[j] + SCTIN[j] + PHT[j] + PRAT[j] + PREL[j];
                double ca_noco = SCTIN[j] + PHT[j] + PRAT[j] + PREL[j];
                if (KZ[j] > 0)
                {
                    mexPrintf("       %10.3E%10.3E%10.3E%10.3E%10.3E%10.3E%10.3E%10.3E\n", EN[j] * 1.0e-6, SCTCO[j], SCTIN[j], PHT[j] - PHDIF[j], PRAT[j], PREL[j], crossAll - PHDIF[j], ca_noco - PHDIF[j]);
                    mexPrintf("%7d%10.3E%10.3E%10.3E%10.3E%10.3E%10.3E%10.3E%10.3E\n", KZ[j], EN[j] * 1.0e-6, SCTCO[j], SCTIN[j], PHT[j], PRAT[j], PREL[j], crossAll, ca_noco);
                }
                else
                    mexPrintf("       %10.3E%10.3E%10.3E%10.3E%10.3E%10.3E%10.3E%10.3E\n", EN[j] * 1.0e-6, SCTCO[j], SCTIN[j], PHT[j] - PHDIF[j], PRAT[j], PREL[j], crossAll - PHDIF[j], ca_noco - PHDIF[j]);
            }
        }
        else if (nlhs == 1)
        {

            int LEN = length;
            for (int i = 0; i < length; i++)
            {
                if (KZ[i] > 0)
                {
                    LEN++;
                }
            }
            double* Energy = new double[LEN];
            double* Coh = new double[LEN];
            double* Incoh = new double[LEN];
            double* PhotonEA = new double[LEN];
            double* NuclearFieldPairProduction = new double[LEN];
            double* ElectronPairFieldProduction = new double[LEN];
            double* TotalWithCoh = new double[LEN];
            double* TotalWithoutCoh = new double[LEN];

            int j = 0;
            for (int i = 0; i < length; ++i)
            {
                if (KZ[i] > 0)
                {
                    Energy[j] = EN[i] * 1E-3;
                    Energy[j + 1] = EN[i] * 1E-3;
                    Coh[j] = SCTCO[i];
                    Coh[j + 1] = SCTCO[i];
                    Incoh[j] = SCTIN[i];
                    Incoh[j + 1] = SCTIN[i];
                    PhotonEA[j] = PHT[i] - PHDIF[i];
                    PhotonEA[j + 1] = PHT[i];
                    NuclearFieldPairProduction[j] = PRAT[i];
                    NuclearFieldPairProduction[j + 1] = PRAT[i];
                    ElectronPairFieldProduction[j] = PREL[i];
                    ElectronPairFieldProduction[j + 1] = PREL[i];
                    TotalWithCoh[j] = SCTCO[i] + SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                    TotalWithCoh[j + 1] = SCTCO[i] + SCTIN[i] + PHT[i] + PRAT[i] + PREL[i];
                    TotalWithoutCoh[j] = SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                    TotalWithoutCoh[j + 1] = SCTIN[i] + PHT[i] + PRAT[i] + PREL[i];

                    j = j + 2;
                }
                else
                {
                    Energy[j] = EN[i] * 1E-3;
                    Coh[j] = SCTCO[i];
                    Incoh[j] = SCTIN[i];
                    PhotonEA[j] = PHT[i] - PHDIF[i];
                    NuclearFieldPairProduction[j] = PRAT[i];
                    ElectronPairFieldProduction[j] = PREL[i];
                    TotalWithCoh[j] = SCTCO[i] + SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                    TotalWithoutCoh[j] = SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                    j = j + 1;
                }
                
            }
            CopyToStructMatrix(plhs[0], OutputNames.at(0).c_str(), KZ, length);
            CopyToStructMatrix(plhs[0], OutputNames.at(1).c_str(), Energy, LEN);
            CopyToStructMatrix(plhs[0], OutputNames.at(2).c_str(), Coh, LEN);
            CopyToStructMatrix(plhs[0], OutputNames.at(3).c_str(), Incoh, LEN);
            CopyToStructMatrix(plhs[0], OutputNames.at(4).c_str(), PhotonEA, LEN);
            CopyToStructMatrix(plhs[0], OutputNames.at(5).c_str(), NuclearFieldPairProduction, LEN);
            CopyToStructMatrix(plhs[0], OutputNames.at(6).c_str(), ElectronPairFieldProduction, LEN);
            CopyToStructMatrix(plhs[0], OutputNames.at(7).c_str(), TotalWithCoh, LEN);
            CopyToStructMatrix(plhs[0], OutputNames.at(8).c_str(), TotalWithoutCoh, LEN);

            delete[] Energy;
            delete[] Coh;
            delete[] Incoh;
            delete[] PhotonEA;
            delete[] NuclearFieldPairProduction;
            delete[] ElectronPairFieldProduction;
            delete[] TotalWithCoh;
            delete[] TotalWithoutCoh;
            
        }
        delete[] SCTCO;
        delete[] SCTIN;
        delete[] PHT;
        delete[] PRAT;
        delete[] PREL;
        delete[] PHDIF;
    }
    else
    {
        if (mxIsDouble(prhs[1]))
        {
            if (mxGetM(prhs[1]) != 1u)
            {
                mexErrMsgTxt("Input Error : second input dimensions error and is not a 1 by N double matrix\n");
            }
            else
            {
                size_t NumOfPoints = mxGetN(prhs[1]);
                //float SCTCO[ME], SCTIN[ME], PHT[ME], PRAT[ME], PREL[ME], PHDIF[ME];

                JENG = int(NumOfPoints);
                float* EN = new float[JENG];
                int* KZ = new int[JENG];
                int* KM = new int[JENG];
                mxDouble* input1 = mxGetPr(prhs[1]);
                for (int i = 0; i < JENG; ++i)
                {
                    KZ[i] = -1;
                    KM[i] = 0;
                    EN[i] = static_cast<float>(input1[i]*1E3);//keV->eV
                }
                NEGO = 3; //3:only input energy list  2: show all energy list (default and adding) 
                NF = 3;  //3: cm2/g
                length = JENG;
                float* SCTCO = new float[length];
                float* SCTIN = new float[length];
                float* PHT = new float[length];
                float* PRAT = new float[length];
                float* PREL = new float[length];
                float* PHDIF = new float[length];
                Calculation(KMAX, NZ, WT, NF, NEGO, length, EN, KZ, KM, SCTCO, SCTIN, PHT, PRAT, PREL, PHDIF);
                int LEN = length;
                for (int i = 0; i < length; i++)
                {
                    if (KZ[i] > 0)
                    {
                        LEN++;
                    }
                }
                double* Energy = new double[LEN];
                double* Coh = new double[LEN];
                double* Incoh = new double[LEN];
                double* PhotonEA = new double[LEN];
                double* NuclearFieldPairProduction = new double[LEN];
                double* ElectronPairFieldProduction = new double[LEN];
                double* TotalWithCoh = new double[LEN];
                double* TotalWithoutCoh = new double[LEN];

                int j = 0;
                for (int i = 0; i < length; ++i)
                {
                    if (KZ[i] > 0)
                    {
                        Energy[j] = EN[i] * 1E-3;
                        Energy[j + 1] = EN[i] * 1E-3;
                        Coh[j] = SCTCO[i];
                        Coh[j + 1] = SCTCO[i];
                        Incoh[j] = SCTIN[i];
                        Incoh[j + 1] = SCTIN[i];
                        PhotonEA[j] = PHT[i] - PHDIF[i];
                        PhotonEA[j + 1] = PHT[i];
                        NuclearFieldPairProduction[j] = PRAT[i];
                        NuclearFieldPairProduction[j + 1] = PRAT[i];
                        ElectronPairFieldProduction[j] = PREL[i];
                        ElectronPairFieldProduction[j + 1] = PREL[i];
                        TotalWithCoh[j] = SCTCO[i] + SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                        TotalWithCoh[j + 1] = SCTCO[i] + SCTIN[i] + PHT[i] + PRAT[i] + PREL[i];
                        TotalWithoutCoh[j] = SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                        TotalWithoutCoh[j + 1] = SCTIN[i] + PHT[i] + PRAT[i] + PREL[i];

                        j = j + 2;
                    }
                    else
                    {
                        Energy[j] = EN[i] * 1E-3;
                        Coh[j] = SCTCO[i];
                        Incoh[j] = SCTIN[i];
                        PhotonEA[j] = PHT[i] - PHDIF[i];
                        NuclearFieldPairProduction[j] = PRAT[i];
                        ElectronPairFieldProduction[j] = PREL[i];
                        TotalWithCoh[j] = SCTCO[i] + SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                        TotalWithoutCoh[j] = SCTIN[i] + PHT[i] + PRAT[i] + PREL[i] - PHDIF[i];
                        j = j + 1;
                    }

                }
                CopyToStructMatrix(plhs[0], OutputNames.at(0).c_str(), KZ, length);
                CopyToStructMatrix(plhs[0], OutputNames.at(1).c_str(), Energy, LEN);
                CopyToStructMatrix(plhs[0], OutputNames.at(2).c_str(), Coh, LEN);
                CopyToStructMatrix(plhs[0], OutputNames.at(3).c_str(), Incoh, LEN);
                CopyToStructMatrix(plhs[0], OutputNames.at(4).c_str(), PhotonEA, LEN);
                CopyToStructMatrix(plhs[0], OutputNames.at(5).c_str(), NuclearFieldPairProduction, LEN);
                CopyToStructMatrix(plhs[0], OutputNames.at(6).c_str(), ElectronPairFieldProduction, LEN);
                CopyToStructMatrix(plhs[0], OutputNames.at(7).c_str(), TotalWithCoh, LEN);
                CopyToStructMatrix(plhs[0], OutputNames.at(8).c_str(), TotalWithoutCoh, LEN);

                delete[] Energy;
                delete[] Coh;
                delete[] Incoh;
                delete[] PhotonEA;
                delete[] NuclearFieldPairProduction;
                delete[] ElectronPairFieldProduction;
                delete[] TotalWithCoh;
                delete[] TotalWithoutCoh;

                delete[] SCTCO;
                delete[] SCTIN;
                delete[] PHT;
                delete[] PRAT;
                delete[] PREL;
                delete[] PHDIF;

                delete[] EN;
                delete[] KZ;
                delete[] KM;
            }
        }
        else
        {   
            mexErrMsgTxt("Input Error : second input is not a double matrix\n");
        }
    }
    
    delete[] NZ;
    delete[] WT;
}

template <typename T>
void CopyToStructMatrix(mxArray* dest, const char* fieldname, T* src, int length)
{
    mxArray* temp;
    temp = mxCreateDoubleMatrix(length, 1, mxREAL);
    mxDouble* temp_ptr = mxGetPr(temp);
    for (int j = 0; j < length; j++)
    {
        temp_ptr[j] = mxDouble(src[j]);
    }
    mxSetField(dest, 0, fieldname, temp);

}