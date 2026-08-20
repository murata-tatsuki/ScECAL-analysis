#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using namespace std;

namespace {

struct SampleInput {
  string label;
  int nEnergy = 0;
  vector<double> energies;
  vector<int> filesPerEnergy;
  vector<string> files;
};

double parseEnergyToken(const string& token) {
  char* endPtr = nullptr;
  const double value = strtod(token.c_str(), &endPtr);
  if (endPtr == token.c_str()) {
    return 0.0;
  }
  if (value < 1.0) {
    return 0.5;
  }
  return value;
}

bool parseSample(int argc, char* argv[], int& argIndex, const string& label, SampleInput& sample) {
  if (argIndex >= argc) {
    cerr << "[ERROR] " << label << ": missing Nenergy." << endl;
    return false;
  }

  sample.label = label;
  sample.nEnergy = atoi(argv[argIndex++]);
  if (sample.nEnergy <= 0) {
    cerr << "[ERROR] " << label << ": Nenergy must be > 0." << endl;
    return false;
  }

  sample.energies.reserve(sample.nEnergy);
  sample.filesPerEnergy.reserve(sample.nEnergy);

  for (int i = 0; i < sample.nEnergy; ++i) {
    if (argIndex >= argc) {
      cerr << "[ERROR] " << label << ": missing energy token." << endl;
      return false;
    }
    sample.energies.push_back(parseEnergyToken(argv[argIndex++]));
  }

  int totalFiles = 0;
  for (int i = 0; i < sample.nEnergy; ++i) {
    if (argIndex >= argc) {
      cerr << "[ERROR] " << label << ": missing file-count token." << endl;
      return false;
    }
    const int nFiles = atoi(argv[argIndex++]);
    if (nFiles < 0) {
      cerr << "[ERROR] " << label << ": file count must be >= 0." << endl;
      return false;
    }
    sample.filesPerEnergy.push_back(nFiles);
    totalFiles += nFiles;
  }

  sample.files.reserve(totalFiles);
  for (int i = 0; i < totalFiles; ++i) {
    if (argIndex >= argc) {
      cerr << "[ERROR] " << label << ": missing input ROOT file path." << endl;
      return false;
    }
    sample.files.push_back(argv[argIndex++]);
  }

  return true;
}

int validateRootFiles(const SampleInput& sample, TTree* summaryTree, string& outSample, string& outPath,
                      double& outEnergy, int& outEntries, int& outStatus) {
  int readOk = 0;
  int fileOffset = 0;

  for (int iEnergy = 0; iEnergy < sample.nEnergy; ++iEnergy) {
    const double beamEnergy = sample.energies[iEnergy];
    const int nFiles = sample.filesPerEnergy[iEnergy];

    for (int iFile = 0; iFile < nFiles; ++iFile) {
      const string& filePath = sample.files[fileOffset + iFile];
      outSample = sample.label;
      outPath = filePath;
      outEnergy = beamEnergy;
      outEntries = 0;
      outStatus = 0;

      TFile* fileIn = TFile::Open(filePath.c_str(), "READ");
      if (!fileIn || fileIn->IsZombie()) {
        cerr << "[WARN] failed to open: " << filePath << endl;
        outStatus = 1;
        summaryTree->Fill();
        if (fileIn) {
          fileIn->Close();
          delete fileIn;
        }
        continue;
      }

      TTree* tree = dynamic_cast<TTree*>(fileIn->Get("Calib_Hit"));
      if (!tree) {
        cerr << "[WARN] missing TTree Calib_Hit: " << filePath << endl;
        outStatus = 2;
        summaryTree->Fill();
        fileIn->Close();
        delete fileIn;
        continue;
      }

      outEntries = static_cast<int>(tree->GetEntries());
      double totalEnergyDep = 0.0;
      tree->SetBranchAddress("TotalEnergyDep", &totalEnergyDep);
      if (outEntries > 0) {
        tree->GetEntry(0);
      }

      outStatus = 10;
      summaryTree->Fill();
      ++readOk;

      fileIn->Close();
      delete fileIn;
    }

    fileOffset += nFiles;
  }

  return readOk;
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "usage: ./between_files_discrepancy output.root "
         << "Nenergy_data E1...En_data E1Nfiles...EnNfiles_data data_input1.root ... "
         << "Nenergy_sim E1...En_sim E1Nfiles...EnNfiles_sim sim_input1.root ..." << endl;
    return 1;
  }

  cout << "=====> " << argv[1] << endl;

  int argIndex = 2;
  SampleInput dataSample;
  SampleInput simSample;

  if (!parseSample(argc, argv, argIndex, "data", dataSample)) {
    return 1;
  }
  if (!parseSample(argc, argv, argIndex, "simulation", simSample)) {
    return 1;
  }

  if (argIndex != argc) {
    cerr << "[ERROR] unexpected extra arguments detected from index " << argIndex << endl;
    return 1;
  }

  TFile outFile(argv[1], "RECREATE");
  if (outFile.IsZombie()) {
    cerr << "[ERROR] failed to create output file: " << argv[1] << endl;
    return 1;
  }

  string sampleName;
  string filePath;
  double energy = 0.0;
  int nEntries = 0;
  int status = 0;

  TTree summary("read_summary", "ROOT file read summary");
  summary.Branch("sample", &sampleName);
  summary.Branch("file_path", &filePath);
  summary.Branch("energy_GeV", &energy);
  summary.Branch("entries", &nEntries);
  summary.Branch("status", &status);

  const int dataOk = validateRootFiles(dataSample, &summary, sampleName, filePath, energy, nEntries, status);
  const int simOk = validateRootFiles(simSample, &summary, sampleName, filePath, energy, nEntries, status);

  outFile.cd();
  summary.Write();
  outFile.Close();

  cout << "[INFO] data files readable: " << dataOk << "/" << dataSample.files.size() << endl;
  cout << "[INFO] simulation files readable: " << simOk << "/" << simSample.files.size() << endl;
  cout << "[INFO] wrote read summary to " << argv[1] << endl;

  return 0;
}
