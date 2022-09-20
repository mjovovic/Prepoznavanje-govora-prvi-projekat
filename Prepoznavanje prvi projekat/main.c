#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <Python.h>
#include <string.h>
#include "matplotlib-cpp/matplotlibcpp.h"
#include "kissfft/kiss_fft.h"

namespace plt = matplotlibcpp;

int noiseLimit(SF_INFO* fileInfo,  short* sampleVal){
    //10ms for num of samples
    int sampleNum = fileInfo->samplerate/10;
    int midVal = 0;

    for (int i = 0; i < sampleNum; i++) {
        midVal += *(sampleVal+i);
    }
    
    midVal = midVal / sampleNum;
 
    int standardDev = 0;
 
    for (int i = 0; i < sampleNum; i++){
        standardDev += (*(sampleVal+i)-midVal) * (*(sampleVal+i)-midVal);
    }
    
    standardDev = sqrt( standardDev / sampleNum );
    
    return (midVal + 2 * standardDev);
    
}
//plotter
void plottingWord(int limit, short* items, int size) {

    int start = 0, end = 0;
    int startFlag = 0;
    for (int i = 0; i < size; i++) {
        if (startFlag == 0 && limit < abs(*(items+i))) {
            start = i;
            startFlag = 1;
        }
        if (startFlag == 1 && limit < abs(*(items+i))) {
            end = i;
        }
        
    }
    printf("start %d\nfinish %d\nlimit %d", start, end, limit);
    
    if (start < 5 ) {
        printf("JUST NOISE");
    }
    
    std::vector<int> pltVec(items, items + size);

    plt::plot( pltVec );

    plt::axvline(start);
    plt::axvline(end);
    plt::axhline(limit);
    plt::axhline(-limit);

    
    plt::show();

    
}

//konverzija u mono
void conversionMono(short* items, SF_INFO* info, short* itemsMono){

    short holder;
    int j = 0;
    for (int i = 1; i < ( info->frames*info->channels ); i+=2) {        

        holder = *(items+i-1);
        *(itemsMono + j) = (holder + *(items + i)) / 2;
        j++;
        
    }
    free(items);     
}

void oneWindowSpec(SF_INFO* info, short* items){
     //samplerate/10 for 100ms
     
    kiss_fft_cfg cfg = kiss_fft_alloc( info->samplerate/10, 0, NULL, NULL );
    
    kiss_fft_cpx * pInput = (kiss_fft_cpx*)malloc(info->samplerate/10 * sizeof(kiss_fft_cpx));
    kiss_fft_cpx * pOutput = (kiss_fft_cpx*)malloc(info->samplerate/10 * sizeof(kiss_fft_cpx));

    for (int i = 0; i < info->samplerate/10; i++){
        pInput[i].r = *(items+i);
        pInput[i].i = 0;
    }
    

    kiss_fft(cfg, pInput, pOutput);

    std::vector<double> pltVec(items, items + info->samplerate/10);
    for (int i = 0; i < info->samplerate/10; i++){
        pltVec[i] = sqrt(pOutput[i].r*pOutput[i].r + pOutput[i].i*pOutput[i].i)/info->frames;//euklidska da bi se otarasili imaginarnog dela i delimo sa brojem framejvoa za normalizaciju
    }

    free(pInput);
    plt::plot( pltVec, 10 );
    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Magnitude");
    plt::show();
}


void sonogram(SF_INFO* info, short* items){
    kiss_fft_cfg cfg = kiss_fft_alloc( info->samplerate/10, 0, NULL, NULL );
    
    kiss_fft_cpx * pInput = (kiss_fft_cpx*)malloc(info->samplerate/10 * sizeof(kiss_fft_cpx));
    kiss_fft_cpx * pOutput = (kiss_fft_cpx*)malloc(info->samplerate/10 * sizeof(kiss_fft_cpx));
    std::vector<double> pltVec(items, items + info->frames);
    
  
   
    for (int i = 0; i < info->frames; i+=info->samplerate/10){

        for (int j = 0; j < info->samplerate/10; j++) {
            if(i+j < info->frames) {
                pInput[j].r = *(items+j+i);
                pInput[j].i = 0;
            }
        }
        

        kiss_fft(cfg, pInput, pOutput);
    
        for (int k = 0; k < info->samplerate/10; k++) {
            if(i+k < info->frames)
                pltVec[i+k] = sqrt(pOutput[k].r*pOutput[k].r + pOutput[k].i*pOutput[k].i)/info->frames;//euklidska da bi se otarasili imaginarnog dela i delimo sa brojem framejvoa za normalizaciju
        }
    }


    std::vector<int>a(1,0);
    std::vector<int>b(1,0);
    char output[50];

    for (int i = 0; i < info->frames/2; i+=info->samplerate/10) {
        a[0] = i/(info->samplerate/10);
        
        int max = 0, min = 0x7fffffff;
    
        for (int j = 0; j < (info->samplerate/20) && i+j < info->frames/2; j++) {
            
            if(pltVec[j+i] > max) 
                max = pltVec[j];
            
            if(pltVec[j+i] < min)
                min = pltVec[j];

        }
        

        for (int j = 0; j < (info->samplerate/20) && i+j < info->frames/2; j++) {
            b[0] = pltVec[i+j];
            if(b[0] > 0){

                
                double alpha = (double)(b[0] - min) / (max-min);
            
                sprintf(output, "%lf", alpha);
                std::string str(output);
                std::map<std::string, std::string> map = {{"alpha", str}};    

                b[0] = j*10;
                
                plt::scatter(a, b, 2.0, map);

            }

        }
        
    }
    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Magnitude");
    plt::show();

    free(pInput);
}

/*
    printf("%d - %d - %d - %d - %d - %d\n", pMaleInfo->channels, pMaleInfo->frames, pMaleInfo->format, pMaleInfo->samplerate, pMaleInfo->sections, pMaleInfo->seekable);
    printf("%d - %d - %d - %d - %d - %d\n", pFemaleInfo->channels, pFemaleInfo->frames, pFemaleInfo->format, pFemaleInfo->samplerate, pFemaleInfo->sections, pFemaleInfo->seekable);
    printf("%d - %d - %d - %d - %d - %d\n", pGuitarInfo->channels, pGuitarInfo->frames, pGuitarInfo->format, pGuitarInfo->samplerate, pGuitarInfo->sections, pGuitarInfo->seekable);
    printf("%d - %d - %d - %d - %d - %d\n", pWhitenoiseInfo->channels, pWhitenoiseInfo->frames, pWhitenoiseInfo->format, pWhitenoiseInfo->samplerate, pWhitenoiseInfo->sections, pWhitenoiseInfo->seekable);
    2 - 27958 - 65538 - 44100 - 1 - 1
    2 - 54455 - 65538 - 48000 - 1 - 1
    2 - 97374 - 65538 - 44100 - 1 - 1
    1 - 441001 - 65538 - 44100 - 1 - 1
  */
int main(int argc, char**argv){
    int a;
    
    printf("chose sound file:\n1.male\n2.female\n3.white noise\n4.gutar\n0.EXIT\n");
    while(1){
        
        scanf("%d", &a);

        if(a == 1) {

            SF_INFO * pMaleInfo = (SF_INFO*)malloc( sizeof(SF_INFO) );
            SNDFILE * pMale = sf_open("./samples/male.wav", SFM_READ, pMaleInfo );
            short * pMaleItemsStereo = (short*)malloc(pMaleInfo->frames*pMaleInfo->channels * sizeof(short));
            short * pMaleItems = (short*)malloc(pMaleInfo->frames * sizeof(short));
            sf_readf_short( pMale, pMaleItemsStereo, pMaleInfo->frames);
            conversionMono(pMaleItemsStereo, pMaleInfo, pMaleItems);
            int maleLimit = noiseLimit(pMaleInfo, pMaleItems);
            plottingWord(maleLimit, pMaleItems, pMaleInfo->frames);
            oneWindowSpec(pMaleInfo, pMaleItems);
          // ne radi zbog labela kad se skinu radi
           // sonogram(pMaleInfo, pMaleItems);

        } else if (a == 2) { 

            SF_INFO * pFemaleInfo = (SF_INFO*)malloc(sizeof(SF_INFO));
            SNDFILE * pFemale = sf_open("./samples/female.wav", SFM_READ, pFemaleInfo);
            short * pFemaleItemsStereo = (short*)malloc(pFemaleInfo->frames*pFemaleInfo->channels * sizeof(short));
            short * pFemaleItems = (short*)malloc(pFemaleInfo->frames * sizeof(short));
            sf_readf_short( pFemale, pFemaleItemsStereo, pFemaleInfo->frames);
            conversionMono(pFemaleItemsStereo, pFemaleInfo, pFemaleItems);
            int femaleLimit = noiseLimit(pFemaleInfo, pFemaleItems);
            plottingWord(femaleLimit, pFemaleItems, pFemaleInfo->frames);
            oneWindowSpec(pFemaleInfo, pFemaleItems);
            sonogram(pFemaleInfo, pFemaleItems);
        
        } else if  (a == 3) {
            
            SF_INFO * pWhitenoiseInfo = (SF_INFO*)malloc(sizeof(SF_INFO));
            SNDFILE * pWhitenoise = sf_open("./samples/whitenoise.wav", SFM_READ, pWhitenoiseInfo);
            short * pWhitenoiseItemsStereo = (short*)malloc(pWhitenoiseInfo->frames*pWhitenoiseInfo->channels * sizeof(short));
            short * pWhitenoiseItems = (short*)malloc(pWhitenoiseInfo->frames * sizeof(short));
            sf_readf_short( pWhitenoise, pWhitenoiseItemsStereo, pWhitenoiseInfo->frames);
            conversionMono(pWhitenoiseItemsStereo, pWhitenoiseInfo, pWhitenoiseItems);
            int whitenoiseLimit = noiseLimit(pWhitenoiseInfo, pWhitenoiseItems);
            plottingWord(whitenoiseLimit, pWhitenoiseItems, pWhitenoiseInfo->frames);
           // oneWindowSpec(pWhitenoiseInfo, pWhitenoiseItems);
           // sonogram(pWhitenoiseInfo, pWhitenoiseItems);

        } else if (a == 4) {

            SF_INFO * pGuitarInfo = (SF_INFO*)malloc(sizeof(SF_INFO));
            SNDFILE * pGuitar = sf_open("./samples/guitar.wav", SFM_READ, pGuitarInfo);
            short * pGuitarItemsStereo = (short*)malloc(pGuitarInfo->frames*pGuitarInfo->channels * sizeof(short));  
            short * pGuitarItems = (short*)malloc(pGuitarInfo->frames * sizeof(short));
            sf_readf_short( pGuitar, pGuitarItemsStereo, pGuitarInfo->frames);
            conversionMono(pGuitarItemsStereo, pGuitarInfo, pGuitarItems);
            int guitarLimit = noiseLimit(pGuitarInfo, pGuitarItems);
            plottingWord(guitarLimit, pGuitarItems, pGuitarInfo->frames);
            oneWindowSpec(pGuitarInfo, pGuitarItems);
            sonogram(pGuitarInfo, pGuitarItems);

        } else if(a == 0){

            return 0;

        } else {

            printf("sound not chosen try again");

        }

   } 
    
}
