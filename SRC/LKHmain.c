#include "LKH.h"
#include "Genetic.h"
#include "LKHmain.h"



/*
 * This file contains the main function of the program.
 */

 int add(int i, int j) {
    return i + j;
}

 int getCost(double** ppArr) {
    printf("cost is %f\n", ppArr[1][3]);
    return ppArr[2][6];
}


void aPythonInput(int* sequence, double** travelCost, int n_stops)
 {

     printf("cost is %f\n", travelCost[1][3]);

     int DIMENSION = n_stops;


     MaxMatrixDimension = 20000;
     MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
         MergeWithTourGPX2;

     customInput(n_stops, travelCost);
     printf("Read Good\n");
     arunLKH();
     printf("solved, LKH\n");

     for (int i = 0; i < n_stops; i++) {

         sequence[i] = BestTour[i];
     }
     //Do not commenents the following out. 
     FreeStructures();
     BestTour = 0;
     BetterTour = 0;
     FirstNode = 0;
     FirstActive = 0;
     LastActive = 0;
     return;


 }


 void arunLKH() {
     GainType Cost, OldOptimum;
     double Time, LastTime;

     StartTime = LastTime = GetTime();

     AllocateStructures();
     CreateCandidateSet();
     InitializeStatistics();

     if (Norm != 0)
         BestCost = PLUS_INFINITY;
     else {
         /* The ascent has solved the problem! */
         Optimum = BestCost = (GainType)LowerBound;
         UpdateStatistics(Optimum, GetTime() - LastTime);
         RecordBetterTour();
         RecordBestTour();
         WriteTour(OutputTourFileName, BestTour, BestCost);
         WriteTour(TourFileName, BestTour, BestCost);
         Runs = 0;
         printf("The ascent has solved the problem!");
     }
     /* Find a specified number (Runs) of local optima */
     for (Run = 1; Run <= Runs; Run++) {
         LastTime = GetTime();
         if (LastTime - StartTime >= TimeLimit) {
             if (TraceLevel >= 1)
                 printff("*** Time limit exceeded ***\n");
             break;
         }
         Cost = FindTour();      /* using the Lin-Kernighan heuristic */

         if (Cost < BestCost) {
             BestCost = Cost;
             RecordBetterTour();
             RecordBestTour();

         }

     }
     //PrintStatistics();
 }








 void PythonInput(int* sequence, double** travelCost, int n_stops)
{

    printf("cost is %f\n", travelCost[1][3]);

	int DIMENSION = n_stops;

    
    MaxMatrixDimension = 20000;
    MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
        MergeWithTourGPX2;

    customInput(n_stops, travelCost);
    printf("Read Good\n");
    runLKH();
    printf("solved, LKH\n");
	
    for (int i = 0; i < n_stops; i++){

        sequence[i] = BestTour[i];
    }
    //Do not commenents the following out. 
    FreeStructures();
    BestTour = 0;
    BetterTour = 0;
    FirstNode = 0;
   FirstActive = 0;
   LastActive = 0;
    return;


}

 void freeme(int* sequence, double** travelCost, int n_stops)
{
    for (int i = 0; i < n_stops; i++) {
        free(travelCost[i]);
    }
    free(travelCost);
    free(sequence);
    travelCost = 0;
    sequence = 0;
    printf("Free all pointers");
}

void runLKH() {
    GainType Cost, OldOptimum;
    double Time, LastTime;

    StartTime = LastTime = GetTime();
    //ReadProblem();

  /*  if (SubproblemSize > 0) {
        if (DelaunayPartitioning)
            SolveDelaunaySubproblems();
        else if (KarpPartitioning)
            SolveKarpSubproblems();
        else if (KCenterPartitioning)
            SolveKCenterSubproblems();
        else if (KMeansPartitioning)
            SolveKMeansSubproblems();
        else if (RohePartitioning)
            SolveRoheSubproblems();
        else if (MoorePartitioning || SierpinskiPartitioning)
            SolveSFCSubproblems();
        else
            SolveTourSegmentSubproblems();
        return EXIT_SUCCESS;
    }*/
    AllocateStructures();
    CreateCandidateSet();
    InitializeStatistics();

    if (Norm != 0)
        BestCost = PLUS_INFINITY;
    else {
        /* The ascent has solved the problem! */
        Optimum = BestCost = (GainType)LowerBound;
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        WriteTour(OutputTourFileName, BestTour, BestCost);
        WriteTour(TourFileName, BestTour, BestCost);
        Runs = 0;
        printf("The ascent has solved the problem!");
    }
    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        if (LastTime - StartTime >= TimeLimit) {
            if (TraceLevel >= 1)
                printff("*** Time limit exceeded ***\n");
            break;
        }
        Cost = FindTour();      /* using the Lin-Kernighan heuristic */
      //  if (MaxPopulationSize > 1) {
      //      ///* Genetic algorithm */
      //      //int i;
      //      //for (i = 0; i < PopulationSize; i++) {
      //      //    GainType OldCost = Cost;
      //      //    Cost = MergeTourWithIndividual(i);
      //      //    if (TraceLevel >= 1 && Cost < OldCost) {
      //      //        printff("  Merged with %d: Cost = " GainFormat, i + 1,
      //      //            Cost);
      //      //        if (Optimum != MINUS_INFINITY && Optimum != 0)
      //      //            printff(", Gap = %0.4f%%",
      //      //                100.0 * (Cost - Optimum) / Optimum);
      //      //        printff("\n");
      //      //    }
      //      //}
      //      //if (!HasFitness(Cost)) {
      //      //    if (PopulationSize < MaxPopulationSize) {
      //      //        AddToPopulation(Cost);
      //      //        if (TraceLevel >= 1)
      //      //            PrintPopulation();
      //      //    }
      //      //    else if (Cost < Fitness[PopulationSize - 1]) {
      //      //        i = ReplacementIndividual(Cost);
      //      //        ReplaceIndividualWithTour(i, Cost);
      //      //        if (TraceLevel >= 1)
      //      //            PrintPopulation();
      //      //    }
      //      //}
      //  }
      //  else if (Run > 1)
      //      Cost = MergeTourWithBestTour();
        if (Cost < BestCost) {
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
            /*WriteTour(OutputTourFileName, BestTour, BestCost);
            WriteTour(TourFileName, BestTour, BestCost);*/
        }
      //  OldOptimum = Optimum;
      //  if (Cost < Optimum) {
      //      if (FirstNode->InputSuc) {
      //          Node* N = FirstNode;
      //          while ((N = N->InputSuc = N->Suc) != FirstNode);
      //      }
      //      Optimum = Cost;
      //      printff("*** New optimum = " GainFormat " ***\n", Optimum);
      //  }
      //  Time = fabs(GetTime() - LastTime);
      //  UpdateStatistics(Cost, Time);
      //  if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
      //      printff("Run %d: Cost = " GainFormat, Run, Cost);
      //      if (Optimum != MINUS_INFINITY && Optimum != 0)
      //          printff(", Gap = %0.4f%%",
      //              100.0 * (Cost - Optimum) / Optimum);
      //      printff(", Time = %0.2f sec. %s\n\n", Time,
      //          Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
      //  }
      //  if (StopAtOptimum && Cost == OldOptimum && MaxPopulationSize >= 1) {
      //      Runs = Run;
      //      break;
      //  }
      ///*  if (PopulationSize >= 2 &&
      //      (PopulationSize == MaxPopulationSize ||
      //          Run >= 2 * MaxPopulationSize) && Run < Runs) {
      //      Node* N;
      //      int Parent1, Parent2;
      //      Parent1 = LinearSelection(PopulationSize, 1.25);
      //      do
      //          Parent2 = LinearSelection(PopulationSize, 1.25);
      //      while (Parent2 == Parent1);
      //      ApplyCrossover(Parent1, Parent2);
      //      N = FirstNode;
      //      do {
      //          if (ProblemType != HCP && ProblemType != HPP) {
      //              int d = C(N, N->Suc);
      //              AddCandidate(N, N->Suc, d, INT_MAX);
      //              AddCandidate(N->Suc, N, d, INT_MAX);
      //          }
      //          N = N->InitialSuc = N->Suc;
      //      } while (N != FirstNode);
      //  }*/
      //  SRandom(++Seed);
    }
    //PrintStatistics();
}



void test1() {
    int DIMENSION = 6;
    double tempMatrix[6][6] = {
        {0,200,300,200,150,600},
        {100,0.0,300,189,134,400},
        {140,240,0.0,120,150,300},
        {155,130,300,0.0,165,105.6},
        {100,200,300,200,0.0,190.8},
        {100,200,300,215.6,150,0.0},

    };
    double** weightMatrix;
    weightMatrix = (double**)malloc(DIMENSION * sizeof(double*));

    for (int row = 0; row < DIMENSION; row++) {
        weightMatrix[row] = (double*)malloc(DIMENSION * sizeof(double));
    }
    for (int i = 0; i < DIMENSION; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            weightMatrix[i][j] = tempMatrix[i][j];
        }
    }
    int* optimalSequence;
    optimalSequence = (int*)malloc(DIMENSION * sizeof(int));

    PythonInput(optimalSequence, weightMatrix, DIMENSION);
    free(weightMatrix);
    weightMatrix = 0;
    free(optimalSequence);
    optimalSequence = 0;
}
void test2() {
    int DIMENSION = 5;
    double** weightMatrix;
    weightMatrix = (double**)malloc(DIMENSION * sizeof(double*));

    for (int row = 0; row < DIMENSION; row++) {
        weightMatrix[row] = (double*)malloc(DIMENSION * sizeof(double));
    }
    double tempMatrix[5][5] = {
        {0.0,200,300,200,150},
        {100,0.0,300,189,134},
        {140,240,0.0,120,150},
        {155,130,300,0.0,165},
        {100,200,300,200,0.0},

    };
    for (int i = 0; i < DIMENSION; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            weightMatrix[i][j] = tempMatrix[i][j];
        }
    }
    int* optimalSequence;
    optimalSequence = (int*)malloc(DIMENSION * sizeof(int));

    PythonInput(optimalSequence, weightMatrix, DIMENSION);
    for (int i = 0; i < DIMENSION; i++) {
        free(weightMatrix[i]);
    }
    free(weightMatrix);
    weightMatrix = 0;
    free(optimalSequence);
    optimalSequence = 0;
}

int main(int argc, char *argv[])
{
   
    for (int i = 1; i <= 10000; i++) {
        //test1();
        test2();
    }

    system("pause");
    return EXIT_SUCCESS;
}
