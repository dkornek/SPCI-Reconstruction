// ROOT Macro to calculate the intricacy given a system matrix

// author: Dominik Kornek <dominik.kornek@gmail.com>
// last modified: 19-06-20

void getNumberOfDetectorsAndBins(TFile* f, Int_t& nDet, Int_t& nBin){
    // extract index of both the detectors from the name of the spectrum

    Int_t detectorD;
    Int_t detectorC;

    TString nameOfLastVoxel = f->GetListOfKeys()->Last()->GetName();
    TDirectory* dirOfVoxel = (TDirectory*)f->Get(nameOfLastVoxel);
    TString nameOfLastSpectrum = dirOfVoxel->GetListOfKeys()->Last()->GetName();
    nBin = ((TH1F*)dirOfVoxel->Get(nameOfLastSpectrum))->GetNbinsX();

    detectorD = ((TString)nameOfLastSpectrum(3, 2)).Atoi();
    detectorC = ((TString)nameOfLastSpectrum(5, 2)).Atoi();
    nDet = std::max(detectorD, detectorC) + 1;

    Int_t nElements = 0.5 * (nDet * (nDet - 1)) * nBin;
    cout << "Number of Detectors:\t\t\t" << nDet <<
    "\nNumber of Energy Bins:\t\t\t" << nBin << 
    "\nNumber of Measurement Elements:\t\t" << nElements << "\n";
}

Int_t getNumberOfVoxels(TFile* f){

    Int_t nVoxel = 0;

    TIter nextVoxel(f->GetListOfKeys());
    TKey* keyVoxel;
    while ((keyVoxel = (TKey*)nextVoxel())){
        ++nVoxel;
    }
    delete keyVoxel;

    cout << "Number of Voxels:\t\t\t" << nVoxel << "\n";

    return nVoxel;
}

void calculateWorstIntricacy(Int_t nVx, Double_t zS){
    // calculate the worst intricacy

    Double_t worstZeta = TMath::Log2(nVx);

    cout << "\nFor the best imaging system:\t0\n";
    cout << "For the worst imaging system:\t" << 1 << "\n";
    cout << "System Intricacy:\t\t" << zS / worstZeta << "\n";
}

std::vector<std::array<Int_t, 3> > getDetectorElements(Int_t nDet, Int_t nBin){

    std::vector<std::array<Int_t, 3> > elements;
    for (Int_t d = 0; d < nDet; ++d){
        for (Int_t c = d + 1; c < nDet; ++c){
            for (Int_t b = 1; b <= nBin; ++b){
                std::array<Int_t, 3> e = { d, c, b };
                elements.push_back(e);
            }
        }
    }

    return elements;
}

std::vector<std::vector<Double_t> > getSMCoefficients(TFile* f, Int_t nVx, Int_t nDet, Int_t nBin){

    std::vector<std::vector<Double_t> > SM;

    for (Int_t d = 0; d < nDet; ++d){
        for (Int_t c = d + 1; c < nDet; ++c){
            for (Int_t b = 1; b <= nBin; ++b){
                std::vector<Double_t> p_dcbv;
                SM.push_back(p_dcbv);
            }
        }
    }

    for (Int_t v = 0; v < nVx; ++v){

	Int_t indexOfMeasurementElement = 0;

        TString nameOfVoxel;
        nameOfVoxel.Form("%03i", v);
        TDirectory* directoryOfVoxel = (TDirectory*)f->Get(nameOfVoxel);
	
        for (Int_t d = 0; d < nDet; ++d){
            for (Int_t c = d + 1; c < nDet; ++c){
		
                TString nameOfSpectrum;
                nameOfSpectrum.Form("%03i%02i%02i", v, d, c);
                TH1F* spectrum = (TH1F*)directoryOfVoxel->Get(nameOfSpectrum);
		for (Int_t b = 1; b <= nBin; ++b){

                    Double_t binContent = spectrum->GetBinContent(b);
		    SM.at(indexOfMeasurementElement).push_back(binContent);
                    ++indexOfMeasurementElement;
                }
		
		delete spectrum;
            }
        }
    } 

    return SM;
}

std::vector<Double_t> calculateSumOfSMCoef(std::vector<std::vector<Double_t> > SMCoef){

    std::vector<Double_t> sumVector;
    std::for_each(SMCoef.begin(), SMCoef.end(), [&] (std::vector<Double_t> c){
    
        Double_t sumOverVoxels = 0;
	std::for_each(c.begin(), c.end(), [&] (Double_t n){

	    sumOverVoxels += n;

	});

	sumVector.push_back(sumOverVoxels);
    });

    return sumVector;
}

std::vector<Double_t> calculateIntricacies(std::vector<std::vector<Double_t> > SMCoef, 						   std::vector<Double_t> sumOfSMElements){

    std::vector<Double_t> zetas;
    Int_t nEl = 0;
    std::for_each(SMCoef.begin(), SMCoef.end(), [&] (std::vector<Double_t> c){
    
        Double_t zeta = 0;
	
	Int_t voxel = 0;
	std::for_each(c.begin(), c.end(), [&] (Double_t n){

	    if (n != 0){

		Double_t s = sumOfSMElements.at(nEl);
		if (s != 0){

		    Double_t z = n / s;
		    z *= TMath::Log2(z);
		    zeta += z;
		}
	    }

	    ++voxel;
	});

	if (zeta != 0){
	    zeta *= -1;
	    zetas.push_back(zeta);
	}
	++nEl;
    });

    return zetas;
}

Double_t calculateMeanZeta(std::vector<Double_t> RORzetas){

    Int_t K = RORzetas.size();
    Double_t meanZeta = 0;

    std::for_each(RORzetas.begin(), RORzetas.end(), [&] (Double_t z){

	meanZeta += z;
    });

    meanZeta = meanZeta / K;
    return meanZeta;
}


void Intricacy(){

    TString pathToSystemMatrix = "../../data/SystemMatrix/Original/SPCIBase441.root";
    TFile* input = new TFile(pathToSystemMatrix, "READ");

    if (!input->IsOpen()){
        cout << "File not found!\nShut down macro.\n";
    }

    // *******************************************************************************
    // get detector info
    Int_t numberOfDetectors;
    Int_t numberOfBins;
    getNumberOfDetectorsAndBins(input, numberOfDetectors, numberOfBins);

    // get image space info
    Int_t numberOfVoxels = getNumberOfVoxels(input);
    

    // get system matrix info
    Int_t numberOfSMCoefficients = 0.5 * (numberOfDetectors *
        (numberOfDetectors - 1)) * numberOfBins * numberOfVoxels;
    cout << "Number of SystemMatrix Elements:\t" << numberOfSMCoefficients << "\n";

    // *******************************************************************************
    // get measurement elements
    cout << "\nGet measurement elements ...\n";

    std::vector<std::array<Int_t, 3> > detectorElements =
        getDetectorElements(numberOfDetectors, numberOfBins);

    // get system matrix coefficients for measurement elements
    cout << "\nGet coefficients for measurement elements ...\n";

    std::vector<std::vector<Double_t> > SMCoefInMeasurementElements =
        getSMCoefficients(input, numberOfVoxels, numberOfDetectors, numberOfBins);

    // calculate sum for every measurement element over all voxels
    cout << "\nGet sum over voxels for measurement elements ...\n";

    std::vector<Double_t> sumOfME = calculateSumOfSMCoef(SMCoefInMeasurementElements);

    // calculate ROR intricacies
    cout << "\nCalculate ROR intricacies ...\n";

    std::vector<Double_t> zetaOfMeasurementElements = 
	calculateIntricacies(SMCoefInMeasurementElements, sumOfME);

    // calculate system intricacy
    cout << "\nCalculate system intricacy ...\n";

    Double_t zetaOfSystem = calculateMeanZeta(zetaOfMeasurementElements);
    calculateWorstIntricacy(numberOfVoxels, zetaOfSystem);

    // *******************************************************************************
    // shut down macro
    input->Close();
}















