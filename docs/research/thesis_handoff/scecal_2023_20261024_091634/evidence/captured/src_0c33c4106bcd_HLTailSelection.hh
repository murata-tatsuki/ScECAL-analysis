#ifndef RESOLUTION_HL_TAIL_SELECTION_HH
#define RESOLUTION_HL_TAIL_SELECTION_HH

#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"

namespace hl_tail {
namespace fs = std::filesystem;
inline void require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error("HL tail: " + message);
}
inline std::unique_ptr<TFile> open(const fs::path& path) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  require(file && !file->IsZombie(), "cannot read " + path.string());
  return file;
}
template<class T> inline T* object(TFile& file, const char* name) {
  auto* value = dynamic_cast<T*>(file.Get(name));
  hl_tail::require(value, std::string("missing ") + name + " in " + file.GetName());
  return value;
}
template<class T> inline void bind(TTree& tree, const char* name, T* value) {
  require(tree.GetBranch(name), std::string("missing event branch ") + name);
  tree.SetBranchStatus(name, 1);
  require(tree.SetBranchAddress(name, value) >= 0, std::string("invalid event branch ") + name);
}

// Match the exact calibrated input, never a basename or a TriggerID alone.
class Selection {
  struct Event { int run, time, trigger; bool rich; };
  using FileKey = std::pair<int, int>;
  std::map<std::string, FileKey> inputs_;
  std::map<FileKey, std::map<Long64_t, Event>> events_;
  std::string definition_;
  fs::path source_;
public:
  explicit Selection(const fs::path& directory) : source_(fs::canonical(directory)) {
    for (int sample = 0; sample < 2; ++sample) {
      const auto manifest = source_ / "inputs" / (sample ? "mc.tsv" : "data.tsv");
      std::ifstream input(manifest);
      require(bool(input), "cannot read " + manifest.string());
      std::string line;
      int index = 0;
      while (std::getline(input, line)) {
        const auto tab = line.find('\t');
        require(tab != std::string::npos && tab + 1 < line.size(), "invalid manifest " + manifest.string());
        const auto path = fs::canonical(line.substr(tab + 1)).string();
        require(inputs_.emplace(path, FileKey{sample, index++}).second, "duplicate manifest input " + path);
      }
    }
    auto file = open(source_ / "tail_study.root");
    require(object<TParameter<int>>(*file, "min_tail_channels")->GetVal() == 3,
            "expected tail-rich threshold of 3 channels");
    require(object<TParameter<double>>(*file, "nsigma")->GetVal() == 5,
            "expected low-tail threshold of -5 sigma");
    definition_ = object<TNamed>(*file, "tail_definition")->GetTitle();
    auto* tree = object<TTree>(*file, "events");
    tree->SetBranchStatus("*", 0);
    int sample = 0, index = 0, run = 0, time = 0, trigger = 0, channels = 0;
    Long64_t entry = 0;
    bool rich = false;
    bind(*tree, "sample", &sample); bind(*tree, "file_index", &index);
    bind(*tree, "entry", &entry); bind(*tree, "Run_Num", &run);
    bind(*tree, "Event_Time", &time); bind(*tree, "TriggerID", &trigger);
    bind(*tree, "tail_rich", &rich); bind(*tree, "low_tail_channels", &channels);
    // Deliberately do not read 'selected': resolution retains its own CoG/event cuts.
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
      require(tree->GetEntry(i) > 0, "cannot read tail event");
      require((sample == 0 || sample == 1) && index >= 0 && entry >= 0 && channels >= 0,
              "invalid tail event key");
      require(rich == (channels >= 3), "inconsistent tail-rich flag");
      require(events_[{sample, index}].emplace(entry, Event{run, time, trigger, rich}).second,
              "duplicate tail event entry");
    }
    tree->ResetBranchAddresses();
  }

  std::vector<unsigned char> mask(const std::string& calibrated) const {
    const auto path = fs::canonical(calibrated).string();
    const auto found = inputs_.find(path);
    require(found != inputs_.end(), "input is absent from tail manifest: " + path +
            "; use matching --tail-results or --keep-tail-events");
    const auto foundEvents = events_.find(found->second);
    auto file = open(path);
    auto* tree = object<TTree>(*file, "Calib_Hit");
    const Long64_t count = tree->GetEntries();
    require(foundEvents != events_.end() && Long64_t(foundEvents->second.size()) == count,
            "incomplete tail coverage for " + path + "; run hg_lg_tail on all events first");
    tree->SetBranchStatus("*", 0);
    int run = 0, time = 0, trigger = 0;
    bind(*tree, "Run_Num", &run); bind(*tree, "Event_Time", &time); bind(*tree, "Event_Num", &trigger);
    std::vector<unsigned char> result(count);
    for (Long64_t entry = 0; entry < count; ++entry) {
      require(tree->GetEntry(entry) > 0, "cannot read calibrated event");
      const auto event = foundEvents->second.find(entry);
      require(event != foundEvents->second.end(), "missing tail entry in " + path);
      const auto& expected = event->second;
      require(run == expected.run && time == expected.time && trigger == expected.trigger,
              "event identity mismatch in " + path + " at entry " + std::to_string(entry));
      result[entry] = expected.rich;
    }
    tree->ResetBranchAddresses();
    return result;
  }
  const std::string& definition() const { return definition_; }
};

// MultiEnergy consumes fitted histograms. It must never pretend to cut events there.
inline void checkAnalysis(TFile& file, bool exclude) {
  auto* flag = dynamic_cast<TParameter<int>*>(file.Get("hl_tail_excluded"));
  if (!flag) {
    require(!exclude, std::string("no tail-selection metadata in ") + file.GetName() +
            "; rerun SingleEnergyAnalysis, or use --keep-tail-events for legacy results");
    return;
  }
  require(flag->GetVal() == int(exclude), std::string("tail selection mismatch in ") + file.GetName());
  auto* status = object<TNamed>(file, "hl_tail_status");
  require(std::string(status->GetTitle()) == "complete", std::string("incomplete analysis ") + file.GetName());
}
} // namespace hl_tail
#endif
