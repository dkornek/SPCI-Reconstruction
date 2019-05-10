// this macro reads a tree of histograms and rebins them. The new file is saved.

TString newFileName(TString str, int bins){
    // create new file name based on the original name and the new number of bins

    // get rid of ".root"
    TString newName = str;
    Int_t extension = newName.Last('.');
    Int_t stringLength = newName.Length();
    newName.Remove(extension, stringLength - extension);

    // append number of bins + ".root"
    TString suffix;
    suffix.Form("_%i_bins.root", bins);
    newName.Append(suffix.Data());

    return newName;
}

Bool_t isBinningConsistent(int nSpectra, int nBins, double& intPart){

    double ratio;
    ratio = (double)nBins / (double)nSpectra;
    double doublePart = std::modf(ratio, &intPart);

    if (doublePart != 0.0){
        cout << "Binning pattern not consistent.\nMacro aborted.\n";
        return kFALSE;
    } else{
        return kTRUE;
    }
}

int getFactor(double intPart){
    int rebinningFactor;  // how many bins will be combined into one

    cout << "There are " << intPart << " bins in each histogram.\n";
    cout << "Type rebinning factor: ";
    cin >> rebinningFactor;
    return rebinningFactor;
}

int getNewBinNumber(double intPart, double rebinningFactor){
    int newBins;
    newBins = intPart / rebinningFactor;
    cout << "The new histograms will contain " << newBins << " bins each.\n";

    return newBins;
}

void RebinningMeasurements(TFile& file, TString path){

    // 1) check if all histograms follow the same binning pattern
    TList* contents = file.GetListOfKeys();
    TH1F* spectrum;

    int nSpectra = 0, nBins = 0;
    TIter next(contents);
    TKey* key;
    while ((key = (TKey*)next())){
        spectrum = (TH1F*)file.Get(key->GetName());
        nBins += spectrum->GetNbinsX();
        ++nSpectra;
    }

    double intPart;
    if (!isBinningConsistent(nSpectra, nBins, intPart)){
        return;
    }

    // 2) do the actual rebinning and save the file
    int rebinningFactor = getFactor(intPart);
    int newBins = getNewBinNumber(intPart, rebinningFactor);

    TString newName = newFileName(path, newBins);
    TFile* newFile = new TFile(newName, "RECREATE");

    next.Reset();
    while ((key = (TKey*)next())){
        spectrum = (TH1F*)file.Get(key->GetName());
        spectrum->Rebin(rebinningFactor);
        spectrum->Write(key->GetName());

        delete spectrum;
    }

    newFile->Write();
    cout << "File saved.\n";
    delete newFile;
}

void RebinningProjections(TFile& file, TString path){

    // 1) check if all histograms follow the same binning pattern
    TList* voxels = file.GetListOfKeys();
    TDirectory* spectra;
    TH1F* spectrum;

    int nSpectra = 0, nBins = 0;
    TIter next(voxels);
    TKey* key;
    TKey* keyVoxel;
    while ((key = (TKey*)next())){
        spectra = (TDirectory*)file.Get(key->GetName());

        TIter nextVoxel(spectra->GetListOfKeys());
        while ((keyVoxel = (TKey*)nextVoxel())){
            spectrum = (TH1F*)spectra->Get(keyVoxel->GetName());
            nBins += spectrum->GetNbinsX();
            ++nSpectra;

            delete spectrum;
        }
    }

    double intPart;
    if (!isBinningConsistent(nSpectra, nBins, intPart)){
        return;
    }

    // 2) do the actual rebinning and save the file
    int rebinningFactor = getFactor(intPart);
    int newBins = getNewBinNumber(intPart, rebinningFactor);

    TString newName = newFileName(path, newBins);
    TFile* newFile = new TFile(newName, "RECREATE");

    next.Reset();
    while ((key = (TKey*)next())){
        spectra = (TDirectory*)file.Get(key->GetName());
        TDirectory* subDir = newFile->mkdir(key->GetName(), key->GetTitle());
        subDir->cd();

        TIter nextVoxel(spectra->GetListOfKeys());
        while ((keyVoxel = (TKey*)nextVoxel())){
            spectrum = (TH1F*)spectra->Get(keyVoxel->GetName());
            spectrum->Rebin(rebinningFactor);
            spectrum->Write(keyVoxel->GetName());

            delete spectrum;
        }
    }

    newFile->Write();
    cout << "File saved.\n";
    delete newFile;
}

void RebinningSpectra(){
    int choice;
    cout << "Do you want to rebin the measurements (1) or the base file (2)?: ";
    cin >> choice;

    TString pathToFile;
    cout << "Type path to file: ";
    cin >> pathToFile;

    // 0) Open file
    TFile* pFile = new TFile(pathToFile, "READ");

    if (!pFile->IsOpen()){
        cout << "File not found.\n";
        return;
    }
    if (choice == 1){
        RebinningMeasurements(*pFile, pathToFile);
    }

    if (choice == 2){
        RebinningProjections(*pFile, pathToFile);
    }

    delete pFile;
    return;
}
