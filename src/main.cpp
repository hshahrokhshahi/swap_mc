#define NDEBUG
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <sys/stat.h>
#include "system.h"
#include "config.h"
#include "timer.h"
#include "export2D.h"
#include <random>

void MonteCarlo(Config config, bool usePreviousStates,
                std::string simulationID, std::string previousID,
                int sampleBegin, std::string simulation_num);

void MC_to_eq(Config config, bool usePreviousStates,
              std::string simulationID, std::string previousID,
              int sampleBegin, std::string simulation_num);

void PrintAcceptanceInfo(const System& system, int numIterations);
std::string GetCurrentTime(std::chrono::high_resolution_clock::time_point current,
                           double offset);
void CheckFileExistence(std::string fileName);

// Argument options:
// 2: sample delay & path
// 3: sample delay & path & previous path
int main(int argc, char* argv[])
{
    Timer timer;

    try
    {
        bool usePreviousStates = false;

        std::string previousPath;
        std::string currentPath;
        std::string simulation_num;
        long long int sampleBegin;
        if(argc==3)
        {
            currentPath = argv[2];
            std::cout<<"Current path: "<<currentPath<<std::endl;
            sampleBegin = std::stoll(argv[1]);
        }
        else if(argc==4)
        {
            currentPath = argv[2];
            std::cout<<"Current path: "<<currentPath<<std::endl;
            previousPath = argv[3];
            std::cout<<"Previous path: "<<previousPath<<std::endl;

            std::string previousConfigFile = previousPath + "/config.txt";
            CheckFileExistence(previousConfigFile);

            usePreviousStates = true;

            sampleBegin = std::stoll(argv[1]);
        }
        else if(argc==5)
        {
            currentPath = argv[2];
            std::cout<<"now we have 5 argument!"<<std::endl;
            std::cout<<"Current path: "<<currentPath<<std::endl;
            previousPath = argv[3];
            std::cout<<"Previous path: "<<previousPath<<std::endl;
            simulation_num = argv[4];

            std::string previousConfigFile = previousPath + "/config.txt";
            CheckFileExistence(previousConfigFile);

            usePreviousStates = true;

            sampleBegin = std::stoll(argv[1]);
        }
        else
        {
            std::cout<<argc<<std::endl;
            throw std::out_of_range("Wrong number of arguments.");
        }


        std::string configFile = currentPath + "/config.txt";
        CheckFileExistence(configFile);
        Config config(configFile);

        std::cout<<"\nPosition sampling starts at iteration: "<<sampleBegin<<std::endl;

        MC_to_eq(config, usePreviousStates, currentPath, previousPath, sampleBegin, simulation_num);

        MonteCarlo(config, usePreviousStates, currentPath, previousPath, sampleBegin, simulation_num);
    }
    catch(std::out_of_range& e)
    {
        std::cout<<e.what()<<std::endl;
    }
    catch(std::invalid_argument& e)
    {
        std::cout<<e.what()<<std::endl;
    }
}

void MC_to_eq(Config config, bool usePreviousStates,
                std::string currentPath, std::string previousPath,
                int sampleBegin, std::string simulation_num)
{
    System system(config, usePreviousStates, previousPath);
    const long long int numIterations = 100000;
    const double swapProbability = config.GetSwapProbability();

    long long int attemptedSwaps = 0;
    long long int attemptedTranslations = 0;

    for(long long int it=0; it<numIterations; ++it)
    {

        if(system.IsChosenWithProbability(swapProbability))
        {
            system.AttemptSwap();
            ++attemptedSwaps;
        }
        else
        {
            system.AttemptTranslation();
            ++attemptedTranslations;
        }
    }


}
void MonteCarlo(Config config, bool usePreviousStates,
                std::string currentPath, std::string previousPath,
                int sampleBegin, std::string simulation_num)
{
    System system(config, usePreviousStates, previousPath);

    const long long int numIterations = config.GetNumIterations();
    const double swapProbability = config.GetSwapProbability();

    std::string outputIterationsFile = currentPath + "/iterations.txt";
    std::string outputEnergyFile = currentPath + "/energy.txt";
    std::string outputSwapFile = currentPath + "/swapAcceptance.txt";
    std::string outputTranslationFile = currentPath + "/translationAcceptance.txt";
    ClearContents(outputIterationsFile);
    ClearContents(outputEnergyFile);
    ClearContents(outputSwapFile);
    ClearContents(outputTranslationFile);
    std::vector<std::vector<double>> exportedStates;

    // Record start time.
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    std::cout<<"Simulation started at "<<GetCurrentTime(start, 0)<<std::endl<<std::endl;
///// THIS PART I HAVE ADDED THE OUTPUT FILES CREATOR
    const long long int n_data_per_graph_points = 2000; // this should be either 1000 or 10000
    int n_graph_points = (int) (5 * log10(numIterations));
    long long int max_iteration = numIterations;
    long int delta_time_arr[n_graph_points];
    int reminder, quitinet; // these two number set the delta_time for each j by using dt = log10(10^(x/5)*2*(x%5))
    std::mt19937 generator ((unsigned int)time(NULL));
    long int random_position_arr[n_graph_points][n_data_per_graph_points];
    for(int i=0; i < n_graph_points; ++i){
        reminder = i % 5;
        quitinet = (int) (i / 5);
        if(reminder == 0)
            delta_time_arr[i] = (long int) round(pow(10, quitinet));
        if(reminder == 1)
            delta_time_arr[i] = (long int) round(pow(10, quitinet + 0.2));
        if(reminder == 2)
            delta_time_arr[i] = (long int) round(pow(10, quitinet+ 0.4));
        if(reminder == 3)
            delta_time_arr[i] = (long int) round(pow(10, quitinet+ 0.6));
        if(reminder == 4)
            delta_time_arr[i] = (long int) round(pow(10, quitinet+ 0.8));
        long int max_rand = (long int) (max_iteration - delta_time_arr[i])/ n_data_per_graph_points;
        std::cout<< i << " = i , max rand = " << max_rand << " ,dt = "<<delta_time_arr[i]<<std::endl;
        int sum_rand = 0;
        std::uniform_int_distribution<long int> distribution(0, (int) max_rand);
        std::string stringpath = "./" + std::to_string(i);
        int status = mkdir(stringpath.c_str(), 0777);
        for(int j=0; j<n_data_per_graph_points; ++j){
            sum_rand += distribution(generator);
            random_position_arr[i][j] = sum_rand;
            std::string stringpath2 = "./" + std::to_string(i) + "/" + std::to_string(j);
            int status2 = mkdir(stringpath2.c_str(), 0777);
            std::string outputStatesFile = stringpath2 + "/requested_configs.txt";
            ClearContents(outputStatesFile);
//            std::cout << "i, j = " << i << j << " , rand = " << sum_rand << std::endl;
        }
    }
    const long int rand_arr_size = n_graph_points * n_data_per_graph_points;
    long int start_iteration[rand_arr_size];
    long int finish_iteration[rand_arr_size];
    long int i_start_index[rand_arr_size];
    long int j_start_index[rand_arr_size];
    long int i_finish_index[rand_arr_size];
    long int j_finish_index[rand_arr_size];
    long long int index;
    for(int i = 0; i < n_graph_points; ++i) {
        for (int j = 0; j < n_data_per_graph_points; j++) {
            index = i * n_data_per_graph_points + j;
            start_iteration[index] = random_position_arr[i][j];
            finish_iteration[index] = random_position_arr[i][j] + delta_time_arr[i];
            i_start_index[index] = i;
            j_start_index[index] = j;
            i_finish_index[index] = i;
            j_finish_index[index] = j;
//            std::cout << "i, j, index = " << i<<j << " " <<index<< " , min_arr[i][j] = " << start_iteration[index] << " | "<< finish_iteration[index] << std::endl;
        }
    }
    /// this part we sort the start & finish arr so we have 2 increasing arr in order to have more ef alg to write xyz
    long int temp, temp_i, temp_j;
    for (int i = 0; i < rand_arr_size - 1; ++i) {
        int min = i;
        for (int j = i + 1; j < rand_arr_size; ++j)
            if (start_iteration[j] < start_iteration[min])
                min = j;
        temp = start_iteration[i];
        temp_i = i_start_index[i];
        temp_j = j_start_index[i];
        start_iteration[i] = start_iteration[min];
        i_start_index[i] = i_start_index[min];
        j_start_index[i] = j_start_index[min];
        start_iteration[min] = temp;
        i_start_index[min] = temp_i;
        j_start_index[min] = temp_j;
    }
    for (int i = 0; i < rand_arr_size - 1; ++i) {
        int min = i;
        for (int j = i + 1; j < rand_arr_size; ++j)
            if (finish_iteration[j] < finish_iteration[min])
                min = j;
        temp = finish_iteration[i];
        temp_i = i_finish_index[i];
        temp_j = j_finish_index[i];
        finish_iteration[i] = finish_iteration[min];
        i_finish_index[i] = i_finish_index[min];
        j_finish_index[i] = j_finish_index[min];
        finish_iteration[min] = temp;
        i_finish_index[min] = temp_i;
        j_finish_index[min] = temp_j;
    }
//    for(int i = 0; i < rand_arr_size; ++i)
//        std::cout << "index = " << i << " , start_arr[i]= " << start_iteration[i] << " ,with i, j =" << i_start_index[i] << j_start_index[i] << std::endl;
//    for(int i = 0; i < rand_arr_size; ++i)
//        std::cout << "index = " << i << " , finish_arr[i]= " << finish_iteration[i] << " ,with i, j =" << i_finish_index[i] << j_finish_index[i] << std::endl;
    long long int attemptedSwaps = 0;
    long long int attemptedTranslations = 0;
    int whichPrint = 0;
    long long int logScaler = 2;
    long int start_index = 0;
    long int finish_index = 0;
    long int i, j;
    for(long long int it=0; it<numIterations; ++it)
    {
//        std::cout << "it = " << it << std::endl;
        //// this part is dedicated to writing the output
        if(it == start_iteration[start_index]){
            while(it == start_iteration[start_index]){
                i = i_start_index[start_index];
                j = j_start_index[start_index];
                std::string stringpath2 = "./" + std::to_string(i) + "/" + std::to_string(j);
                std::string outputStatesFile = stringpath2 + "/requested_configs.txt";
                exportedStates = system.GetStates();
                Export2D(exportedStates, outputStatesFile, it);
                start_index += 1;
            }
        }

        if(it == finish_iteration[finish_index]){
            while(it == finish_iteration[finish_index]){
                i = i_finish_index[finish_index];
                j = j_finish_index[finish_index];
                std::string stringpath2 = "./" + std::to_string(i) + "/" + std::to_string(j);
                std::string outputStatesFile = stringpath2 + "/requested_configs.txt";
                exportedStates = system.GetStates();
                Export2D(exportedStates, outputStatesFile, it);
                finish_index += 1;
            }
        }
        ////
        if(it%((long long int) logScaler/2)==0)
        {
            ExportItem(it, outputIterationsFile);
            ExportItem(system.GetTotalEnergy(), outputEnergyFile);

            // Export MC acceptance ratios.
            double swapRatio = (double) system.GetAcceptedSwaps()/attemptedSwaps;
            double translationRatio = (double) system.GetAcceptedTranslations()/attemptedTranslations;

            ExportItem(swapRatio, outputSwapFile);
            ExportItem(translationRatio, outputTranslationFile);
        }
        if(it==logScaler)
        {
            logScaler *= 2;
        }

        if(system.IsChosenWithProbability(swapProbability))
        {
            system.AttemptSwap();
            ++attemptedSwaps;
        }
        else
        {
            system.AttemptTranslation();
            ++attemptedTranslations;
        }

        // Printing progress and ETA.
        int numProgressUpdates = 10;
        if(it%((long long int) (numIterations-1)/numProgressUpdates)==0)
        {
            if(it!=0)
            {
                whichPrint++;
                auto current = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> timeSinceStart = current - start;
                int numProgressUpdatesToDo = numProgressUpdates - whichPrint;
                auto estimatedTimeOfCompletion = timeSinceStart/whichPrint*numProgressUpdatesToDo;

                long long int progress = (long long int) 100*it/(numIterations-1);

                auto dateTime = GetCurrentTime(current, 0);

                double secondsLeft = estimatedTimeOfCompletion.count();

                std::cout<<std::setw(3)<<progress<<"%"
                         << std::setw(15)<<secondsLeft/3600<<"h"
                         <<" at "<<dateTime<<std::endl;
                std::cout<<std::setw(24)<<"ETA: "<<GetCurrentTime(current, secondsLeft)
                         <<std::endl;
            }
        }
    }
}



std::string GetCurrentTime(std::chrono::high_resolution_clock::time_point current, double offset)
{
    std::time_t currentTime = std::chrono::system_clock::to_time_t(current) + offset;
    auto dateTimeInfo = std::ctime(&currentTime);
    std::string dateTimeInfoString(dateTimeInfo);
    auto dateTime = dateTimeInfoString.substr(4,15);
    return dateTime;
}

void CheckFileExistence(std::string fileName)
{
    std::ifstream file(fileName);
    if(!file)
    {
        std::string message = fileName + " does not exist.";
        throw std::invalid_argument(message);
    }
}
