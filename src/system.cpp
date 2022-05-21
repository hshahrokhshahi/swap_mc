#define NDEBUG

#include <iostream>
#include <cmath>
#include <random>
#include <cassert>
#include <algorithm>
#include <string>
#include "timer.h"
#include "system.h"
#include "config.h"
#include "states.h"
#include "distribution.h"

System::System(const Config& config, const bool usePreviousStates, std::string previousPath)
        :   numSpheres(config.GetNumSpheres()),
            spheres(numSpheres),
            temperatureFixed(config.GetTemperatureFixed()),
            numDensity(config.GetNumDensity()),
            volumeBox(numSpheres/numDensity),
            lengthBox(cbrt(volumeBox)),
            mersenneTwister(config.GetRandomSeed()),
            randomDouble(0,1),
            randomPosNegDouble(-1,1),
            randomParticle(0,numSpheres-1),
            maxdV(config.GetmaxdV())
{
    Timer timer;

    // Get the number of spheres along one dimension.
    int numSpheres1D;
    for(int i=2; i<=numSpheres; ++i)
    {
        if(i*i*i == numSpheres)
        {
            numSpheres1D = i;
            break;
        }
    }
    double lengthUnit = lengthBox/numSpheres1D;

    maxTranslationDistance = config.GetMaxTranslationDistanceInLengthUnits() * lengthUnit;

    if(usePreviousStates)
    {
        std::cout<<"Using previous states"<<std::endl;
        std::string previousOutputFile = previousPath + "/lastState.txt";
        States states(previousOutputFile, numSpheres);

        std::vector<std::vector<double>> spherePositions(numSpheres);
        spherePositions = states.GetSample();

        for(int i=0; i<numSpheres; ++i)
        {
            spheres[i].position.x = spherePositions[i][0];
            spheres[i].position.y = spherePositions[i][1];
            spheres[i].position.z = spherePositions[i][2];
            spheres[i].radius = spherePositions[i][3];
        }
    }
    else
    {
        std::cout<<"Not using previous states"<<std::endl;
        int numSpecies = config.GetNumSpecies();
        int n = numSpecies - 1;

        Distribution distribution(numSpecies);

        double sigmaMax = distribution.GetSigmaMax();
        double sigmaMin = distribution.GetSigmaMin();

        // Function constants.
        double a = distribution.GetConstA();
        double b = distribution.GetConstB();
        double c = distribution.GetConstC();

        double interval = (sigmaMax - sigmaMin)/n;

        std::vector<double> sigmas(numSpecies);
        std::vector<double> weights(numSpecies);

        for(int i=0; i<numSpecies; ++i)
        {
            sigmas[i] = i*interval + sigmaMin;
            weights[i] = a/(pow((b+c*i),3));
        }

        std::discrete_distribution<> randomSigma(weights.begin(), weights.end());

        double cumulativeRadius=0;
        for(int i=0; i<numSpheres; ++i)
        {
            int index = randomSigma(mersenneTwister);
            double radius = 0.5*sigmas[index];
            spheres[i].radius = radius;

            cumulativeRadius += radius;
        }
        std::cout<<"Average radius: "<<cumulativeRadius/numSpheres<<std::endl;

        // Initialize positions.
        int numPlacedSpheres = 0;
        for(int i=0; i<numSpheres1D; ++i)
        {
            for(int j=0; j<numSpheres1D; ++j)
            {
                for(int k=0; k<numSpheres1D; ++k)
                {
                    if(numPlacedSpheres<numSpheres)
                    {
                        spheres[numPlacedSpheres].position.x = (i + 0.5)*lengthUnit;
                        spheres[numPlacedSpheres].position.y = (j + 0.5)*lengthUnit;
                        spheres[numPlacedSpheres].position.z = (k + 0.5)*lengthUnit;
                    }
                    else
                    {
                        break;
                    }
                    ++numPlacedSpheres;
                }
            }
        }
    }
}

std::vector<std::vector<double>> System::GetStates() const
{
    std::vector<std::vector<double>> sphereStates(numSpheres);
    std::vector<double> sphereState(4);

    for(int i=0; i<numSpheres; ++i)
    {
        sphereState[0] = spheres[i].position.x;
        sphereState[1] = spheres[i].position.y;
        sphereState[2] = spheres[i].position.z;
        sphereState[3] = spheres[i].radius;
        sphereStates[i] = sphereState;
    }
    return sphereStates;
}

void System::AttemptTranslation()
{
    int randomParticleIndex = ChooseRandomParticle();
    double energy = CalculateEnergy(randomParticleIndex, spheres[randomParticleIndex]);
//
    Sphere newSphere = spheres[randomParticleIndex];
    newSphere.position.x += maxTranslationDistance*randomPosNegDouble(mersenneTwister);
    newSphere.position.y += maxTranslationDistance*randomPosNegDouble(mersenneTwister);
    newSphere.position.z += maxTranslationDistance*randomPosNegDouble(mersenneTwister);

    double energyNew = CalculateEnergy(randomParticleIndex, newSphere);
    double energyDifference = energyNew - energy;

    double acceptProbability = std::min(1.0,
                                        std::exp(-energyDifference/(boltzmannConstant * temperatureFixed)));

    if(IsChosenWithProbability(acceptProbability))
    {
        CorrectForPeriodicSphere(newSphere);
        spheres[randomParticleIndex] = newSphere;
        ++acceptedTranslations;
    }
}

void System::AttemptSwap()
{
    double p = 7.23;
    double energyOld = GetTotalEnergy();
    double v_old = volumeBox;
    double l_old = lengthBox;
    double dv = maxdV * randomPosNegDouble(mersenneTwister);
    volumeBox += dv;
    lengthBox = pow(volumeBox, 0.333);
    double l_coef = pow(dv, 0.333);
//    std::cout << dv << "-" << l_coef << std::endl;
    for(int i = 0; i < numSpheres; ++i){
        spheres[i].position.x *= lengthBox / l_old;
        spheres[i].position.y *= lengthBox / l_old;
        spheres[i].position.z *= lengthBox / l_old;
    }
    double energyNew = GetTotalEnergy();
    double dU = energyNew - energyOld;
    double weight = std::exp(-(dU + p * dv - numSpheres*temperatureFixed*std::log(volumeBox/v_old))/temperatureFixed);
    if(weight < 1) {
        double x = (randomPosNegDouble(mersenneTwister) + 1) / 2;
        if (weight < x) {
            volumeBox -= dv;
            l_old = lengthBox;
            lengthBox = pow(volumeBox, 0.333);
            for (int i = 0; i < numSpheres; ++i) {
                spheres[i].position.x *= lengthBox / l_old;
                spheres[i].position.y *= lengthBox / l_old;
                spheres[i].position.z *= lengthBox / l_old;
            }
        } else {
            acceptedSwaps++;
        }
    }
    else{
        acceptedSwaps++;
    }


}

int System::ChooseRandomParticle()
{
    int randomParticleIndex = randomParticle(mersenneTwister);
    return randomParticleIndex;
}

double System::CalculateEnergy(const int index, const Sphere& sphere)
{
    double energy = 0;
    for(int i=0; i<numSpheres; ++i)
    {
        if(i!=index)
        {
//            if(sphere.radius + spheres[i].radius == 0.907268){
//                type = 1;
//            }
//            else if(sphere.radius + spheres[i].radius == 1.460254){
//                type = 2;
//            }else if(sphere.radius + spheres[i].radius == 2.01324){
//                type = 3;
//            }else{
//                std::cout << "hahahaha" << std::endl;
//            }
            energy += PotentialSRPP(RadiusSumOf(sphere,spheres[i]),
                                    DistanceBetween(sphere,spheres[i]));
        }
    }
    return energy;
}

double System::PotentialSRPP(const double sigmaSummedRadius, const double distanceBetweenSpheres) const
{
    double potential;
//    std::cout << sigmaSummedRadius << std::endl;
    const double cutOffDistance = 1.25 * sigmaSummedRadius;
    if(distanceBetweenSpheres > cutOffDistance)
    {
        potential = 0;
        return potential;
    }

    else {
        if (sigmaSummedRadius < 2.95) {
            const double sigma = 1.414;
            const double epsilon = 0.23852;;
            potential = 4 * epsilon * (pow(sigma / distanceBetweenSpheres, 12) -
                                       pow(sigma / distanceBetweenSpheres, 6));  // we have added this potential
            return potential;
        } else if (sigmaSummedRadius < 2.99) {
            const double sigma = 1.0;
            const double epsilon = 0.23852;;
            potential = 4 * epsilon * (pow(sigma / distanceBetweenSpheres, 12) -
                                       pow(sigma / distanceBetweenSpheres, 6));  // we have added this potential
            return potential;

        } else {
            const double sigma = 1.414;
            const double epsilon = 0.23852;
            potential = 4 * epsilon * (pow(sigma / distanceBetweenSpheres, 12) -
                                       pow(sigma / distanceBetweenSpheres, 6));  // we have added this potential
            return potential;
        }
    }
        const int n = 12;
        const double c0 = -1.92415;
        const double c2 = 2.11106;
        const double c4 = -0.591097;
//        const double sigma = 1.0;
//        const double epsilon = 1.0;
//        potential = 4 * epsilon * (pow(sigma/distanceBetweenSpheres, 12) - pow(sigma/distanceBetweenSpheres, 6));  // we have added this potential
////        potential = pow((sigmaSummedRadius/distanceBetweenSpheres),n)
////                    + c4 * pow((distanceBetweenSpheres/sigmaSummedRadius),4)
////                    + c2 * pow((distanceBetweenSpheres/sigmaSummedRadius),2)
////                    + c0;
//    }
//    return potential;
}

double System::RadiusSumOf(const Sphere& sphere1, const Sphere& sphere2) const
{
    double nonAdditivityConstant = 0.2;

    return (sphere1.radius + sphere2.radius)*
           (1 - 2*nonAdditivityConstant*abs(sphere1.radius - sphere2.radius));
}

double System::DistanceBetween(const Sphere& sphere1, const Sphere& sphere2)
{
    double diffX = sphere1.position.x - sphere2.position.x;
    double diffY = sphere1.position.y - sphere2.position.y;
    double diffZ = sphere1.position.z - sphere2.position.z;

    CorrectForPeriodicDistance(diffX);
    CorrectForPeriodicDistance(diffY);
    CorrectForPeriodicDistance(diffZ);

    return sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
}

void System::CorrectForPeriodicDistance(double& length) const
{
    bool negativeOutsideBoundary = (length < -0.5*lengthBox);
    bool positiveOutsideBoundary = (length > 0.5*lengthBox);
    if(positiveOutsideBoundary)
    {
        length -= lengthBox;
    }
    else if(negativeOutsideBoundary)
    {
        length += lengthBox;
    }
    assert(!(length < -0.5*lengthBox || length > 0.5*lengthBox));
}

void System::CorrectForPeriodicSphere(Sphere& sphere)
{
    CorrectForPeriodicCoordinate(sphere.position.x);
    CorrectForPeriodicCoordinate(sphere.position.y);
    CorrectForPeriodicCoordinate(sphere.position.z);
}

void System::CorrectForPeriodicCoordinate(double& coordinate) const
{
    while(coordinate < 0.0)
    {
        coordinate += lengthBox;
    }
    while(coordinate >= lengthBox)
    {
        coordinate -= lengthBox;
    }
}

bool System::IsChosenWithProbability(const double probabilityReference)
{
    double probabilityRandom = randomDouble(mersenneTwister);

    if(probabilityRandom<=probabilityReference)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int System::GetAcceptedTranslations() const
{
    return acceptedTranslations;
}
int System::GetAcceptedSwaps() const
{
    return acceptedSwaps;
}

double System::GetTotalEnergy()
{
    double energy = 0;
    for(int i=0; i < (numSpheres-1); ++i)
    {
        for(int j=i+1; j < numSpheres; ++j)
        {
            energy += PotentialSRPP(RadiusSumOf(spheres[i],spheres[j]),
                                    DistanceBetween(spheres[i],spheres[j]));
        }
    }
    return energy;
}
