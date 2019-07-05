// ROOT Macro to delete specified voxels of the image space in System Matrix data

template <typename T>
std::pair<bool, int> findInVector(const std::vector<T>& vecOfElements, const T& element)
{
        std::pair<bool, int> result;

        // Find given element in vector
        auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

        if (it != vecOfElements.end())
        {
                result.second = distance(vecOfElements.begin(), it);
                result.first = true;
        }
        else
        {
                result.first = false;
                result.second = -1;
        }

        return result;
}

void ChooseVoxelsInSM(){

    TString pathToSM = "../../data/SystemMatrix/Original/SPCIBase441_above.root";
    TFile* input = new TFile(pathToSM, "UPDATE");

    // choose voxels to be left in the SM
    std::vector<Int_t> vxVector = { 154, 155, 156, 157, 158, 159, 160,
                                    175, 176, 177, 178, 179, 180, 181,
                                    196, 197, 198, 199, 200, 201, 202,
                                    217, 218, 219, 220, 221, 222, 223,
                                    238, 239, 240, 241, 242, 243, 244,
                                    259, 260, 261, 262, 263, 264, 265,
                                    280, 281, 282, 283, 284, 285, 286 };
    Int_t vxSize = vxVector.size();

    TIter nextVoxel(input->GetListOfKeys());
    TKey* keyVoxel;
    while ((keyVoxel = (TKey*)nextVoxel())){
        // iterate through every voxel
        TString nameOfVoxel = keyVoxel->GetName();
        Int_t voxel = nameOfVoxel.Atoi();

        std::pair<bool, int> result = findInVector<int>(vxVector, voxel);

        if (!result.first){

            TDirectory* dir = (TDirectory*)input->Get(nameOfVoxel);
            dir->cd();
            dir->Delete("T*;*");
            input->cd();
            TString n;
            n.Form("%s;*", nameOfVoxel.Data());
            input->Delete(n);
            cout << "Deleted voxel " << nameOfVoxel << "...\n";
        }
    }

    input->Close();
    cout << "File closed!\n";
}
