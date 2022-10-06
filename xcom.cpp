#include "xcom.h"
#include <map>
#include <stack>
#include <vector>
#include <regex>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"
// read MDATX3 file 
int ReadMDATX3(char *file,MDATX3 *p)
{
    FILE *fin;
    fin = fopen(file,"r");
    if(!fin)
    {
        //printf("read data %s failed\n",file);
        string filename(file);
        string msg = "read data" + filename + "failed\n";
        mexErrMsgTxt(msg.c_str());
        return -1;
    }
    fscanf(fin,"%d%f",&(p->IZ),&(p->ATWT));
    //ATWT =ATWTS[[K]];
    fscanf(fin,"%d%d",&(p->MAXEDG),&(p->MAXE));
    if(p->MAXEDG >0 )
    {
        for(int i=0;i<p->MAXEDG;i++)
            fscanf(fin,"%d",p->IDG + p->MAXEDG-i-1);
        for(int i=0;i<p->MAXEDG;i++)
        {
            char strBuf[20];
            fscanf(fin,"%s",strBuf);
            (p->ADG)[i] = strBuf;
        }
        for(int i=0;i<p->MAXEDG;i++)
            fscanf(fin,"%f",p->EDGEN + p->MAXEDG-i-1);
    }
    for(int M=0;M<p->MAXE;M++)
    {
        fscanf(fin,"%f",p->E+M);
    }
    for(int M=0;M<p->MAXE;M++)
    {
        fscanf(fin,"%f",p->SCATCO+M);
    }
    for(int M=0;M<p->MAXE;M++)
    {
        fscanf(fin,"%f",p->SCATIN+M);
    }
    for(int M=0;M<p->MAXE;M++)
    {
        fscanf(fin,"%f",p->PHOT+M);
    }
    for(int M=0;M<p->MAXE;M++)
    {
        fscanf(fin,"%f",p->PAIRAT+M);
    }
    for(int M=0;M<p->MAXE;M++)
    {
        fscanf(fin,"%f",p->PAIREL+M);
    }
    if(p->MAXEDG> 0)
    {
        fscanf(fin,"%d",&(p->LAX));
        for(int L=0;L<p->LAX;L++)
        {
            fscanf(fin,"%d",p->KMX+L);
        }
        for(int L=0;L<p->LAX;L++)
        {
            int IMAX = (p->KMX)[L];
            for(int I=0;I<IMAX;I++)
            {
                float fd;
                fscanf(fin,"%f",&fd);
                (p->ENG)[L][I] =fd;
            }
        }
        for(int L=0;L<p->LAX;L++)
        {
            int IMAX = (p->KMX)[L];
            for(int I=0;I<IMAX;I++)
            {
                float fd;
                fscanf(fin,"%f",&fd);
                (p->PHC)[L][I] =fd;
            }
        }
    }
    fclose(fin);
    return 0;
}

/////////////////////////////////////////////////
//  Initial the parameters about the searched material 
//  Input:
//  - KMAX: numbers of element in material
//  - NZ  : lists of Z of element in material
//  - NEGO: Flag how to print energy list(3:only print energy list inputed; 2:add input energy list in default)
//  - JENG: number of energy list inputed
//  - EAD : list of energy inputed
//  - EN  : list of energy outputed
//  - KZ  : flag of 
//  - KM  : flag of 
int InitEnergyList(int KMAX,int *NZ,int NEGO,int JENG,/*int *JZ,int *JM,*/float *EAD,float *EN,int *KZ,int *KM)
{
    int lenEN =0;
    if(NEGO ==3) //only input energy list
    {
        for(int i=0;i<JENG;i++)
        {
            EN[i] =EAD[i];
            KZ[i] = -1;
            KM[i] = 0;
        }
        // MMAX =1;
        //LEN[0] = JENG;
        lenEN =JENG;
    }
    else
    {
        int  NENG= 80;//number of energy list 
        float ENB[80] ={1.0E+03,1.5E+03,2.0E+03,3.0E+03,4.0E+03,5.0E+03,
            6.0E+03,8.0E+03,1.0E+04,1.5E+04,2.0E+04,3.0E+04,4.0E+04,
            5.0E+04,6.0E+04,8.0E+04,1.0E+05,1.5E+05,2.0E+05,3.0E+05,
            4.0E+05,5.0E+05,6.0E+05,8.0E+05,1.0E+06,1.022E+06,1.25E+06,
            1.5E+06,2.0E+06,2.044E+06,3.0E+06,4.0E+06,5.0E+06,6.0E+06,
            7.0E+06,8.0E+06,9.0E+06,1.0E+07,1.1E+07,1.2E+07,1.3E+07,
            1.4E+07,1.5E+07,1.6E+07,1.8E+07,2.0E+07,2.2E+07,2.4E+07,
            2.6E+07,2.8E+07,3.0E+07,4.0E+07,5.0E+07,6.0E+07,8.0E+07,
            1.0E+08,1.5E+08,2.0E+08,3.0E+08,4.0E+08,5.0E+08,6.0E+08,
            8.0E+08,1.0E+09,1.5E+09,2.0E+09,3.0E+09,4.0E+09,5.0E+09,
            6.0E+09,8.0E+09,1.0E+10,1.5E+10,2.0E+10,3.0E+10,4.0E+10,
            5.0E+10,6.0E+10,8.0E+10,1.0E+11};
        for(int i=0;i<NENG;i++) //default energy list
        {
            EN[i] = ENB[i];
            KZ[i] = 0;
            KM[i] = i+1;
        }
        lenEN = NENG;
        if(NEGO ==2) // add input energy list
        {
            int JZ[MEA],JM[MEA];
            for(int i=0;i<JENG;i++)
            {
                JZ[i] = -1;
                JM[i] = 0;
            }
            lenEN = MERGE(EN,KZ,KM,NENG,EAD,JZ,JM,JENG);
        }
        for(int i=0;i<KMAX;i++) // add nuclear shell energy
        {
            char file[30];
            int LZ[14],LM[14];
            sprintf(file,"./data/MDATX3.%03d",NZ[i]);
            printf("KMAX: %d  NZ[%d]:%d  file: %s\n",KMAX,i,NZ[i],file);
            MDATX3 data;
            if (ReadMDATX3(file, &data)!=0)
            {

            }
            for(int j=0;j<data.MAXEDG;j++)
            {
                LZ[j] = data.IZ;
                LM[j] = j+80+1;
            }
            lenEN = MERGE(EN,KZ,KM,lenEN,data.EDGEN,LZ,LM,data.MAXEDG);
        }

        int KOUNT =0;
        float EAS[MEA];
        int SZ[MEA],SM[MEA];
        if(KZ[1]>0)
        {
            EAS[KOUNT] = sqrt(EN[0]*EN[1]);
            KOUNT++;
        }
        for(int N=1;N<lenEN;N++)
        {
            if(KZ[N]>0&& KZ[N-1]>0)
                KOUNT++;
            else 
                continue;

            if(fabs(EN[N]-EN[N-1]) < 0.15)
            {	
                EAS[KOUNT-1] = EN[N];
                EN[N-1] = EN[N-1] *0.99995;
                EN[N] = EN[N] *1.00005;
            }else{
                EAS[KOUNT-1] = sqrt(EN[N]*EN[N-1]);
            }
        }
        if(KOUNT>0)
        {
            for(int i=0;i<KOUNT;i++)
            {
                SZ[i] = -2;
                SM[i] = 0;
            }
            lenEN = MERGE(EN,KZ,KM,lenEN,EAS,SZ,SM,KOUNT);
        }

        int MX=0;
        int JDG[MEA];
        for(int N=0;N<lenEN;N++)
        {
            if(KZ[N]>0)
            {
                JDG[MX]=N+1;
                MX++;
            }
        }
        int LEN[MEA];
        int MMAX =MX +1;
        if(MMAX > 1)
        {
            int MXED =MX;
            LEN[0] = JDG[0];
            for(int M=1;M<MX;M++)
                LEN[M] = JDG[M]- JDG[M-1]+1;
            LEN[MMAX-1] = lenEN -JDG[MXED-1]+1;
        }
        else
            LEN[0] = NENG;
    }
    return lenEN;
}



///////////////////////////////////////////////////////////////
//  Calculate X-ray and gamma Phnton cross sections for  material 
//  Input:
//  - KMAX: numbers of element in material
//  - NZ  : lists of Z of element in material
//  - WEIGHT : ratio of each element in material
//  - NF   : Flag for select the unit of result (3:cm2/g) 2:b/atom in total cross
//  - NEGO: Flag how to print energy list(3:only print energy list inputed;2:add input energy list in default)
//  - NENG  : number of energy list to calculation
//  -  EN  : lists of energy 
//  - KZ  :
//  - KM  :
//  - SCTCO  : scattering coherent cross sections
//  - SCTIN  : scattering incoherent cross sections
//  - PHT  : photo-electric absorption
//  - PRAT  : pair production in nuclear field
//  - PREL  : pair production in electron field
//  - PHDIF : shell jump energy
void Calculation(int KMAX, int* NZ, float* WEIGHT, int NF, int NEGO, int NENG, float* EN, int* KZ, int* KM, float* SCTCO, float* SCTIN, float* PHT, float* PRAT, float* PREL, float* PHDIF)
{
    /*float ENL[ME];
    float AT[ME], ATNC[ME];
    string ALAB[ME];*/
    float* ENL = new float[NENG];
    float* AT = new float[NENG];
    float* ATNC = new float[NENG];
    string* ALAB = new string[NENG];

    string EDGE[94];

    float AFIT[MEB], BFIT[MEB], CFIT[MEB], DFIT[MEB];
    float SATCO[94], SATIN[94], POT[94], PDIF[94], PAIRT[94], PAIRL[94];
    float PR1[55], PR2[51];

    for (int i = 0; i < NENG; i++)
    {
        ENL[i] = log(EN[i]);
        SCTCO[i] = 0.0f;
        SCTIN[i] = 0.0f;
        PHT[i] = 0.0f;
        PRAT[i] = 0.0f;
        PREL[i] = 0.0f;
        PHDIF[i] = 0.0f;
    }
    for (int K = 0; K < KMAX; K++)
    {
        MDATX3 data;
        char file[30];
        sprintf(file, "./data/MDATX3.%03d", NZ[K]);
        //	printf("%d %s\n",K,file);
        ReadMDATX3(file, &data);

        float ATWT = data.ATWT;
        int MAXK = 0;
        if (data.MAXEDG > 0) {
            MAXK = data.MAXE - data.IDG[data.MAXEDG - 1] + 1;
        }
        if (data.MAXEDG <= 0)
        {
            for (int i = 0; i < data.MAXE; i++)
            {
                SATCO[i] = data.SCATCO[i];
                SATIN[i] = data.SCATIN[i];
                POT[i] = data.PHOT[i];
                PDIF[i] = 0.0;
                PAIRT[i] = data.PAIRAT[i];
                PAIRL[i] = data.PAIREL[i];
            }
        }
        else {
            int IRV = data.MAXEDG - 1;
            for (int i = 0; i < data.MAXEDG; i++)
            {
                int IG = data.IDG[i] - 1;
                int IP = i + 80;
                SATCO[IP] = data.SCATCO[IG];
                SATIN[IP] = data.SCATIN[IG];
                POT[IP] = data.PHOT[IG];
                PDIF[IP] = data.PHOT[IG] - data.PHOT[IG - 1];
                PAIRT[IP] = data.PAIRAT[IG];
                PAIRL[IP] = data.PAIREL[IG];
                EDGE[IP] = data.ADG[IRV];
                IRV--;
            }
            int MB = 0;
            for (int M = 0; M < data.MAXE; M++)
            {
                int I = 0;
                for (I = 0; I < data.MAXEDG; I++)
                {
                    if (M - data.IDG[I] + 2 != 0)
                    {
                        if (M - data.IDG[I] + 1 != 0)
                            continue;
                        else
                            break;
                    }
                    else
                        break;
                }
                if (I < data.MAXEDG && ((M - data.IDG[I] + 2 == 0) || (M - data.IDG[I] + 1 == 0)))
                    continue;

                SATCO[MB] = data.SCATCO[M];
                SATIN[MB] = data.SCATIN[M];
                POT[MB] = data.PHOT[M];
                PDIF[MB] = 0.0;
                PAIRT[MB] = data.PAIRAT[M];
                PAIRL[MB] = data.PAIREL[M];

                MB++;
            }
            int MS = 0;
            for (int M = 0; M < data.MAXE; M++)
            {
                int I = 0;
                for (I = 0; I < data.MAXEDG; I++)
                {
                    if (M - data.IDG[I] + 2 != 0)
                        continue;
                    else
                        break;
                }
                if (I < data.MAXEDG && (M - data.IDG[I] + 2 == 0))
                    continue;

                data.E[MS] = data.E[M];
                data.SCATCO[MS] = data.SCATCO[M];
                data.SCATIN[MS] = data.SCATIN[M];
                data.PHOT[MS] = data.PHOT[M];
                data.PAIRAT[MS] = data.PAIRAT[M];
                data.PAIREL[MS] = data.PAIREL[M];
                MS++;
            }
            data.MAXE = MS;
        }

        REV(data.MAXE, data.E);
        REV(data.MAXE, data.SCATCO);
        REV(data.MAXE, data.SCATIN);
        REV(data.MAXE, data.PHOT);
        REV(data.MAXE, data.PAIRAT);
        REV(data.MAXE, data.PAIREL);

        data.E[50] = EPAIR2;
        data.E[54] = EPAIR1;
        for (int i = 0; i < 54; i++)
        {
            float TERM = (data.E[i] - EPAIR1) / data.E[i];
            PR1[i] = log(data.PAIRAT[i] / pow(TERM, 3.));
        }
        PR1[54] = 3.006275 * PR1[53] - 2.577757 * PR1[52] + 0.571482 * PR1[51];
        for (int i = 0; i < 50; i++)
        {
            float TERM = (data.E[i] - EPAIR2) / data.E[i];
            PR2[i] = log(data.PAIREL[i] / pow(TERM, 3.));
        }
        PR2[50] = 3.006275 * PR2[49] - 2.577757 * PR2[48] + 0.571482 * PR2[47];
        for (int i = 0; i < data.MAXE; i++)
        {
            data.E[i] = log(data.E[i]);
            data.SCATCO[i] = log(data.SCATCO[i]);
            data.SCATIN[i] = log(data.SCATIN[i]);
            data.PHOT[i] = log(data.PHOT[i]);
        }
        float FRAC = 1.0;
        if (NF == 3)
            FRAC = AVOG * WEIGHT[K] / data.ATWT;
        float ECUT = 0.0;
        if (data.MAXEDG > 0)
            ECUT = data.EDGEN[data.MAXEDG - 1];
        if (NENG == 80 && NEGO != 3)
        {
            for (int i = 0; i < NENG; i++)
            {
                SCTCO[i] = SCTCO[i] + FRAC * SATCO[i];
                SCTIN[i] = SCTIN[i] + FRAC * SATIN[i];
                PHT[i] = PHT[i] + FRAC * POT[i];
                PRAT[i] = PRAT[i] + FRAC * PAIRT[i];
                PREL[i] = PREL[i] + FRAC * PAIRL[i];
            }
            continue;
        }

        int IMP = 0;
        for (int N = 0; N < NENG; N++)
        {
            if (KZ[N] < 0 || (KZ[N] > 0 && KZ[N] - data.IZ != 0))
            {
                float T = EN[N];
                float TL = ENL[N];
                float RES;

                SCOF(data.E, data.SCATCO, data.MAXE, AFIT, BFIT, CFIT, DFIT);
                BSPOL(TL, data.E, AFIT, BFIT, CFIT, DFIT, data.MAXE, RES);
                SCTCO[N] = SCTCO[N] + FRAC * exp(RES);

                SCOF(data.E, data.SCATIN, data.MAXE, AFIT, BFIT, CFIT, DFIT);
                BSPOL(TL, data.E, AFIT, BFIT, CFIT, DFIT, data.MAXE, RES);
                SCTIN[N] = SCTIN[N] + FRAC * exp(RES);

                if ((T - EPAIR1) > 0)
                {
                    float TERM = (T - EPAIR1) / T;
                    SCOF(data.E, PR1, 55, AFIT, BFIT, CFIT, DFIT);
                    BSPOL(TL, data.E, AFIT, BFIT, CFIT, DFIT, 55, RES);
                    PRAT[N] = PRAT[N] + FRAC * (pow(TERM, 3.) * exp(RES));

                    if ((T - EPAIR2) > 0)
                    {
                        TERM = (T - EPAIR2) / T;
                        SCOF(data.E, PR2, 51, AFIT, BFIT, CFIT, DFIT);
                        BSPOL(TL, data.E, AFIT, BFIT, CFIT, DFIT, 51, RES);
                        PREL[N] = PREL[N] + FRAC * (pow(TERM, 3.) * exp(RES));
                    }
                }

                if (data.MAXEDG <= 0)
                {
                    SCOF(data.E, data.PHOT, data.MAXE, AFIT, BFIT, CFIT, DFIT);
                    BSPOL(TL, data.E, AFIT, BFIT, CFIT, DFIT, data.MAXE, RES);
                    PHT[N] = PHT[N] + FRAC * exp(RES);
                    continue;
                }
                else {
                    if ((T - ECUT) >= 0)
                    {
                        SCOF(data.E, data.PHOT, MAXK, AFIT, BFIT, CFIT, DFIT);
                        BSPOL(TL, data.E, AFIT, BFIT, CFIT, DFIT, MAXK, RES);
                        PHT[N] = PHT[N] + FRAC * exp(RES);
                        continue;
                    }
                    else {
                        while (T - data.EDGEN[IMP] >= 0)
                            IMP++;

                        int MAXX = data.KMX[IMP];
                        float ENGL[35], PHCL[35];

                        for (int M = 0; M < MAXX; M++)
                        {
                            ENGL[M] = log(1.0E+06) + log(data.ENG[IMP][M]);
                            PHCL[M] = log(data.PHC[IMP][M]);
                        }
                        BLIN(TL, ENGL, PHCL, MAXX, RES);
                        PHT[N] = PHT[N] + FRAC * exp(RES);
                    }
                }
            }
            if (KZ[N] == 0 || KZ[N] == data.IZ)
            {
                int NN = KM[N] - 1;
                if (KZ[N] == data.IZ)
                {
                    PHDIF[N] = FRAC * PDIF[NN];
                    ALAB[N] = EDGE[NN];
                }
                SCTCO[N] = SCTCO[N] + FRAC * SATCO[NN];
                SCTIN[N] = SCTIN[N] + FRAC * SATIN[NN];
                PHT[N] = PHT[N] + FRAC * POT[NN];
                PRAT[N] = PRAT[N] + FRAC * PAIRT[NN];
                PREL[N] = PREL[N] + FRAC * PAIRL[NN];

                continue;
            }
        }
    }
    delete[] ENL;
    delete[] AT; 
    delete[] ATNC;
    delete[] ALAB;
}



//  Reads element symbols or chemical formulas.
//  Input:
//  - W   : C风格字符串数组(以空字符结尾)
//          支持形如 {元素1 原子数1 元素2 原子数2 元素3 原子数3 ...} 的化学式输入
//          元素符号大小写必须正确，否则PbEu(铅和铕1:1合金)与 PBeU(磷、铍和铀1:1:1混合物)无法区分
//          原子数为正浮点数或整数,采用形如"x.y"表达,其中x表示小数点前的数，y表示小数点后的数
//          x为0时可以省略 即 "0.y" 等价于 ".y" ; 原子数为1时可以省略数字1.不可采用科学计数法.
//          例 H2O C2H4 C6H12O6  C27H22O4S  C0.6H1.2O.6
//          
//  - JZ  : 按化学式中元素出现顺序排列的原子序数的数组
//  - WT  : 与JZ中原子序数一一对应的元素质量比例
//  - return : numbers of element in a chemical formulas, i.e. KMAX
int ParseFormulas(const char *W,int *JZ,float *WT)
{
    int MASH1[26]={0,5,6,0,0,9,0,1,53,0,19,0,0,7,8,15,0,0,16,0,92,23,74,0,39,0};
    int MASH2[418] ={70,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
        ,54,0,0,0,73,0,0,0,0,0,0,0,0,65,0,0,0,0,0,0,0,0,43,51,88
            ,0,0,0,0,0,0,0,21,37,0,0,0,0,0,0,52,0,0,0,91,0,0,0,0,0,34,0,0,82,0,0,0,0,0,0,75,30,0,0,11,0,0,
            90,0,0,0,46,0,41
                ,0,0,22,0,0,0,0,0,0,0,57,0,14,45,0,0,0,60,0,0,0,0,0,40
                ,0,0,10,0,0,81,0,0,0,0,0,0,0,0,69,0,0,0,0,0,0,0,0,0,62,0,0,0,0,0,12,0,0,50,0,0,31,0,28,0,0,0,0
                ,86,0,0,0,0,0,0,0,0,0,0,61,0,0,0,3,0,0,0,2,64,0,0,0,0,0,38
                ,0,72,84,0,0,0,20,0,0,0,80,0,26,0,0,0,56,0,0,0,0,0,0
                ,25,0,0,0,0,0,59,0,93,42,48,0,0,44,0,0,0,0,0,58,0,89,0,0,78,76,0,0,98
                ,4,0,0,0,94,0,0,0,0,0,0,49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,47,0,67
                ,0,100,0,0,0,83,0,0,0,0,0,0,0,71,0,0,77,0,0,0,0,0,17,97
                ,0,0,0,0,0,0,0,96,0,0,0,0,0,0,0,0,0,0,13,0,0,0,87,0,0,27,0,95,0,0,0,0,68,0,0,0,0,0,0,0,0,99,0,0,0,0,0,0,0,0,0,0,24
                ,0,0,0,0,0,0,63,0,55,35,0,0,0,0,0,0,0,0,0,18,0,0,0,0,0,0,29,0,33,0,0,0,0,0,0,0,0,85,0,0,0,0,0,0,0,0,79,0,0,0,0,0,66
    };
    //		printf("MASH2[110]:%d ,MASH2[204]:%d ,MASH2[287]:%d ,MASH2[417]:%d\n",MASH2[110],MASH2[204],MASH2[287],MASH2[417]);

    //		printf("MASH1[24]:%d\n",MASH1[24]);


    char IC[72],K[72];
    int NZ[100];

    for(int i=0;i<72;i++)
    {
        IC[i] = W[i];
        if(IC[i]<32) 
            K[i] =1;//no alphabeta digit 
        else if(IC[i]==32)
            K[i]=2; // space
        else if(IC[i]<48)
            K[i] =1;
        else if(IC[i]<58)
            K[i] =3; // digit
        else if(IC[i]<65)
            K[i] =1;
        else if(IC[i]<91)
            K[i] =4; //upper case
        else if(IC[i]<97)
            K[i] =1;
        else if(IC[i]<123)
            K[i] =5;//lower case
        else 
            K[i] = 1;
    }
    /*
       printf("W : %s\n",W);
       for(int i=0;i<72;i++)
       if(K[i]>2)
       printf("w:%c  IC:%d K:%d\n",W[i],IC[i],K[i]);
       */
    int L=0; // index of formulas word
    int M=0; // index of element

    while(K[L]<=2)
        L++;

    int LMIN =L;

    while(1)
    {
        int KG = K[L];
        if(L <=LMIN &&(KG==1 ||KG==3 ||KG==5))
            return M;
        if(L>LMIN && KG==2)
        {
            break;
        }
        if(KG ==4)
        {
            int KG1 = K[L+1];
            if(KG1==1)
                return -M;
            if(KG1==2||KG1==3 ||KG1==4)
            {
                int ICC = IC[L]-64;
                int JT = MASH1[ICC-1];

                if(JT<=0)
                    return -M;

                JZ[M] = JT;
                M++;
                if(KG1 ==2)
                {
                    NZ[M-1]=1;
                    break;
                }
                if(KG1==3)
                {
                    int INN = L+1;
                    int IS =0;
                    int MS[100];

                    while(K[INN]==3)
                    {
                        MS[IS] = IC[INN] -48;
                        IS++;
                        INN++;
                    }
                    L = INN;

                    int ISN =IS;
                    int KFAC =1;

                    NZ[M-1]=0;
                    do{
                        NZ[M-1] += KFAC*MS[IS-1];
                        KFAC *=10;
                        IS--;
                        if(IS==0)
                        {
                            if(NZ[M-1]<=0)
                                return -(M-1);

                            continue;
                        }
                    }while(IS>0);
                }
                if(KG1==4)
                {
                    NZ[M-1]=1;
                    L++;
                    continue;
                }
            }
            if(KG1==5)
            {
                int ICC =9*IC[L+1]-10*IC[L]+9;

                if(ICC<1||ICC>418)
                    return -M;

                if(ICC == 208)
                {
                    if(IC[L]==71)
                    {
                        JZ[M]=32;
                    }else{
                        JZ[M] =84;
                    }
                }else{
                    int JT =MASH2[ICC-1];
                    if(JT<=0)
                        return -M;
                    JZ[M]=JT;
                }
                M++;
                int KG2=K[L+2];
                if(KG2 ==1 ||KG2==5)
                    return -(M-1);
                if(KG2 ==2)
                {
                    NZ[M-1] =1;
                    break;	
                }
                if(KG2==3)
                {
                    int INN = L+2;
                    int IS =0;
                    int MS[100];

                    while(K[INN]==3)
                    {
                        MS[IS] = IC[INN] -48;
                        IS++;
                        INN++;
                    }
                    L = INN;

                    int ISN =IS;
                    int KFAC =1;

                    NZ[M-1]=0;
                    do{
                        NZ[M-1] += KFAC*MS[IS-1];
                        KFAC *=10;
                        IS--;
                        if(IS==0)
                        {
                            if(NZ[M-1]<=0)
                                return -(M-1);

                            continue;
                        }
                    }while(IS>0);
                }
                if(KG2==4)
                {
                    NZ[M-1]=1;
                    L+= 2;
                    continue;
                }
            }
        }
    }

    if(M>0)
    {
        float ASUM =0.0;
        for(int i=0;i<M;i++)
            ASUM =ASUM+ATWTS[JZ[i]-1]*NZ[i];
        for(int i=0;i<M;i++)
            WT[i] = ATWTS[JZ[i]-1]*NZ[i]/ASUM;
    }
    return M;
}

//Read chemical formulas.
//  input   :
//  formula :  chemical formulas,for example , H2O , Fe(OH)3, CuSO4(H2O)5 , Fe1.5Al.6Cu0.2 (i.e. Fe : Al : Cu = 1.5 : 0.6 : 0.2)
//  Z       :  list of atomic number
//  Atoms   :  list of atoms of an element 
//  return  :  true or false
bool ParseChemicalFormula(string formula, vector<int>& Z, vector<double>& atoms)
{
    //Mul: 存储每个')'后的数字因子
    //s: 当前的元素名
    //d: 当前元素名后面的数字因子
    //cnt: 当前元素名需要的总乘积因子
    //ElementMap: 化学式中的元素和相应原子个数成对形成的表
    formula = '(' + formula + ')';//加一层括号为了处理方便. 如果不加括号，形如"aa(...)","b",这样的错误识别不出来，
    //必须在循环最后一步处理进行处理，比较麻烦
    string s;
    string d;
    double cnt = 1.0;
    stack<double> Mul;
    stack<char> Brackets;
    map<string, double> ElementMap;
    regex PositiveNumber("[1-9]\\d*|[1-9]\\d*\\.\\d*|0?\\.\\d*[1-9]\\d*");
    //从末尾到开始进行遍历，每个大写字母可能对应一种元素
    for (auto i = formula.rbegin(); i != formula.rend(); ++i)
    {
        if (isdigit(*i) || (*i) == '.')
        {
            d += *i;

        }
        else if (islower(*i))
        {
            s += *i;
        }
        else if (isupper(*i))
        {
            s += *i;
            reverse(s.begin(), s.end());
            reverse(d.begin(), d.end());
            if (PeriodicTable.count(s) == 0 || (!regex_match(d, PositiveNumber) && !d.empty()))
            {
                return false;
            }
            double temp_d = d.empty() ? 1.0 : stod(d);
            ElementMap[s] += temp_d * cnt;
            d.clear();
            s.clear();
        }
        else if (*i == ')')
        {
            reverse(d.begin(), d.end());
            if (!s.empty() || (!regex_match(d, PositiveNumber) && !d.empty()))
            {
                return false;
            }
            double temp_d = d.empty() ? 1.0 : stod(d);
            cnt *= temp_d;
            Mul.push(temp_d);
            Brackets.push(')');
            d.clear();
        }
        else if (*i == '(')
        {

            if (!s.empty() || !d.empty() || Brackets.empty() || Brackets.top() != ')')
            {
                return false;
            }

            double temp_d = Mul.top();
            cnt /= temp_d;
            Mul.pop();
            Brackets.pop();
        }
        else
        {
            return false;
        }
    }
    if (!Brackets.empty() || ElementMap.empty())
    {
        return false;
    }
    Z.resize(ElementMap.size());
    atoms.resize(ElementMap.size());
    int j = 0;
    for (auto i = ElementMap.begin(); i != ElementMap.end(); ++i)
    {
        Z[j] = PeriodicTable.at(i->first);
        atoms[j] = i->second;
        j++;
    }
    return true;
}

bool GetFractionByWeight(vector<int>& Z, vector<double>& atoms, vector<double>& weight)
{
    size_t n = Z.size();
    if (n == 0)
    {
        return false;
    }
    weight.resize(n);
    double ASUM = 0.0;
    for (int i = 0; i < n; ++i)
    {
        ASUM += ATWTS[Z[i] - 1] * atoms[i];
    }
    for (int i = 0; i < n; ++i)
    {
        weight[i] = ATWTS[Z[i] - 1] * atoms[i] / ASUM;
    }

    return true;
}
// Sorts into monotonically increasing order 
template <typename T>
void SORT(int NMAX, T* E)
{
    float EBIC = 1.0E+20;
    for (int i = 0; i < NMAX - 1; i++)
    {
        int NS;
        T EMIN = EBIC;
        for (int j = i; j < NMAX; j++)
        {
            if (E[j] > EMIN)
                continue;
            EMIN = E[j];
            NS = j;
        }
        E[NS] = E[i];
        E[i] = EMIN;
    }
}
// merges energy lists
template <typename T>
int MERGE(T* E1, int* K1, int* L1, int MMAX, T* E2, int* K2, int* L2, int NMAX)
{
    int MLIM = 200;
    for (int i = 0; i < NMAX; i++)
    {
        int M = 1;
        while (E2[i] - E1[M] > 0.)
            M++;

        int	MC = M;
        int	MF = M + 1;
        MMAX++;

        if (MMAX - MLIM <= 0)
        {
            for (int j = MMAX - 1; j >= MF; j--)
            {
                E1[j] = E1[j - 1];
                K1[j] = K1[j - 1];
                L1[j] = L1[j - 1];
            }
            E1[MC] = E2[i];
            K1[MC] = K2[i];
            L1[MC] = L2[i];
        }
        else
        {
            printf("MERGE stop! MMAX:%d  MLIM:%d\n", MMAX, MLIM);
            MMAX = -1;
            break;
        }
    }

    return MMAX;
}

// Reverses the order of lists
template <typename T>
void REV(int NMAX, T* X)
{
    int NH = NMAX / 2;
    for (int i = 0; i <= NH; i++)
    {
        int i1 = NMAX - i - 1;
        T DT = X[i1];
        X[i1] = X[i];
        X[i] = DT;
    }
}

//fits F as a function of X, and calculates cubic spline coefficients A,B,C and D
template <typename T>
int SCOF(T* X, T* F, const int NMAX, T* A, T* B, T* C, T* D)
{
    int M1 = 1;
    int M2 = NMAX - 2;
    T S = 0.0;
    T R = 0.0;

    for (int i = 0; i <= M2; i++)
    {
        D[i] = X[i + 1] - X[i];
        R = (F[i + 1] - F[i]) / D[i];
        C[i] = R - S;
        S = R;
    }

    S = 0.0;
    R = 0.0;
    C[0] = 0.0;
    C[NMAX - 1] = 0.0;
    for (int i = M1; i <= M2; i++)
    {
        C[i] += R * C[i - 1];
        B[i] = (X[i - 1] - X[i + 1]) * 2.0 - R * S;
        S = D[i];
        R = S / B[i];
    }

    for (int i = M2; i >= M1; i--)
        C[i] = (D[i] * C[i + 1] - C[i]) / B[i];

    for (int i = 0; i <= M2; i++)
    {
        S = D[i];
        R = C[i + 1] - C[i];
        D[i] = R / S;
        C[i] = 3.0 * C[i];
        B[i] = (F[i + 1] - F[i]) / S - (C[i] + R) * S;
        A[i] = F[i];
    }
    return 0;
}

// Evaluates cubic spline as function of S, to obtain fitted result G.
template <typename T>
int BSPOL(const T S, T* X, T* A, T* B, T* C, T* D, const int N, T& G)
{
    int IDIR, MLB, MUB, MU, ML;
    if (X[0] > X[N - 1])
    {
        IDIR = 1;
        MLB = N;
        MUB = 0;
    }
    else
    {
        IDIR = 0;
        MLB = 0;
        MUB = N;
    }

    if (S > X[MUB + IDIR - 1])
        MU = MUB + 2 * IDIR - 2;
    else if (S < X[MLB - IDIR])
        MU = MLB - 2 * IDIR;
    else {
        ML = MLB;
        MU = MUB;

        while (abs(MU - ML) > 1)
        {
            int MAV = (ML + MU) / 2;
            if (S < X[MAV - 1])
                MU = MAV;
            else
                ML = MAV;
        }
        MU = MU + IDIR - 2;
    }

    T Q = S - X[MU];
    G = ((D[MU] * Q + C[MU]) * Q + B[MU]) * Q + A[MU];

    return MU;
}

// Linear interpolation routine
template <typename T>
int BLIN(const T S, T* X, T* Y, const int N, T& TY)
{
    int IDIR, MLB, MUB, MU, ML;

    if (X[0] > X[N - 1])
    {
        IDIR = 1;
        MLB = N;
        MUB = 0;
    }
    else
    {
        IDIR = 0;
        MLB = 0;
        MUB = N;
    }

    if (S > X[MUB + IDIR - 1])
    {
        MU = MUB + 2 * IDIR - 2;
    }
    else if (S < X[MLB - IDIR])
    {
        MU = MLB - 2 * IDIR;
    }
    else {
        ML = MLB;
        MU = MUB;

        while (abs(MU - ML) > 1)
        {
            int MAV = (ML + MU) / 2;
            if (S < X[MAV - 1])
                MU = MAV;
            else
                ML = MAV;
        }
        MU = MU + IDIR - 2;
    }

    T Q = S - X[MU];
    TY = Y[MU] + Q * (Y[MU + 1] - Y[MU]) / (X[MU + 1] - X[MU]);

    return MU;
}

