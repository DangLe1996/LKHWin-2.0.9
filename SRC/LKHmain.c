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


 void PythonInput(int* sequence, double** travelCost, int n_stops)
{

	int DIMENSION = n_stops;
    int** weightMatrix;
    weightMatrix = (int**)malloc(DIMENSION * sizeof(int*));

    for (int row = 0; row < DIMENSION; row++) {
        weightMatrix[row] = (int*)malloc(DIMENSION * sizeof(int));
    }

	for(int i = 0; i < n_stops; i++){
		
		for(int j = 0; j < n_stops; j++){
			
			weightMatrix[i][j] = travelCost[i][j];
		}
	}
    customInput(n_stops, weightMatrix);
    printf("Read Good\n");

    runLKH();
    printf("solved, LKH\n");
	
    for (int i = 0; i < n_stops; i++){

        sequence[i] = BestTour[i];
    }
	FreeStructures();

}

void runLKH() {
    GainType Cost, OldOptimum;
    double Time, LastTime;

    StartTime = LastTime = GetTime();
    MaxMatrixDimension = 20000;
    MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
        MergeWithTourGPX2;
    //ReadProblem();

    if (SubproblemSize > 0) {
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
    }
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
        if (MaxPopulationSize > 1) {
            /* Genetic algorithm */
            int i;
            for (i = 0; i < PopulationSize; i++) {
                GainType OldCost = Cost;
                Cost = MergeTourWithIndividual(i);
                if (TraceLevel >= 1 && Cost < OldCost) {
                    printff("  Merged with %d: Cost = " GainFormat, i + 1,
                        Cost);
                    if (Optimum != MINUS_INFINITY && Optimum != 0)
                        printff(", Gap = %0.4f%%",
                            100.0 * (Cost - Optimum) / Optimum);
                    printff("\n");
                }
            }
            if (!HasFitness(Cost)) {
                if (PopulationSize < MaxPopulationSize) {
                    AddToPopulation(Cost);
                    if (TraceLevel >= 1)
                        PrintPopulation();
                }
                else if (Cost < Fitness[PopulationSize - 1]) {
                    i = ReplacementIndividual(Cost);
                    ReplaceIndividualWithTour(i, Cost);
                    if (TraceLevel >= 1)
                        PrintPopulation();
                }
            }
        }
        else if (Run > 1)
            Cost = MergeTourWithBestTour();
        if (Cost < BestCost) {
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
            WriteTour(OutputTourFileName, BestTour, BestCost);
            WriteTour(TourFileName, BestTour, BestCost);
        }
        OldOptimum = Optimum;
        if (Cost < Optimum) {
            if (FirstNode->InputSuc) {
                Node* N = FirstNode;
                while ((N = N->InputSuc = N->Suc) != FirstNode);
            }
            Optimum = Cost;
            printff("*** New optimum = " GainFormat " ***\n", Optimum);
        }
        Time = fabs(GetTime() - LastTime);
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
            printff("Run %d: Cost = " GainFormat, Run, Cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff(", Gap = %0.4f%%",
                    100.0 * (Cost - Optimum) / Optimum);
            printff(", Time = %0.2f sec. %s\n\n", Time,
                Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
        }
        if (StopAtOptimum && Cost == OldOptimum && MaxPopulationSize >= 1) {
            Runs = Run;
            break;
        }
        if (PopulationSize >= 2 &&
            (PopulationSize == MaxPopulationSize ||
                Run >= 2 * MaxPopulationSize) && Run < Runs) {
            Node* N;
            int Parent1, Parent2;
            Parent1 = LinearSelection(PopulationSize, 1.25);
            do
                Parent2 = LinearSelection(PopulationSize, 1.25);
            while (Parent2 == Parent1);
            ApplyCrossover(Parent1, Parent2);
            N = FirstNode;
            do {
                if (ProblemType != HCP && ProblemType != HPP) {
                    int d = C(N, N->Suc);
                    AddCandidate(N, N->Suc, d, INT_MAX);
                    AddCandidate(N->Suc, N, d, INT_MAX);
                }
                N = N->InitialSuc = N->Suc;
            } while (N != FirstNode);
        }
        SRandom(++Seed);
    }
    PrintStatistics();
}


int main(int argc, char *argv[])
{
    int DIMENSION = 5;
    int** weightMatrix;
    weightMatrix = (int**)malloc(DIMENSION * sizeof(int*));

    for (int row = 0; row < DIMENSION; row++) {
        weightMatrix[row] = (int*)malloc(DIMENSION * sizeof(int));
    }
    int tempMatrix[5][5] = {
        {100,200,300,200,150},
        {100,200,300,189,134},
        {140,240,300,120,150},
        {155,130,300,340,165},
        {100,200,300,200,150},

    };

    for (int i = 0; i < DIMENSION; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            weightMatrix[i][j] = tempMatrix[i][j];
        }
    }

    GainType Cost, OldOptimum;
    double Time, LastTime;

    /* Read the specification of the problem */
    if (argc >= 2)
        ParameterFileName = argv[1];
    int* optimalSequence;
    optimalSequence = (int*)malloc(DIMENSION * sizeof(int));

    PythonInput(optimalSequence, weightMatrix, DIMENSION);
    //customInput(weightMatrix, DIMENSION); //Custom input for Amazon Challenge 
    //////ReadParameters(); // This is initial function to read in problem parameters. 
    //runLKH();
    system("pause");
    return EXIT_SUCCESS;
}
