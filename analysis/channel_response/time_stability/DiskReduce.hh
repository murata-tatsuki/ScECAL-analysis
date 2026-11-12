#pragma once
#include <filesystem>
#include <fstream>
#include <map>
#include <queue>
#include <functional>
#include <type_traits>
#include <unistd.h>
#include <zlib.h>

namespace ts {
// Large sequential TSV writes avoid tiny synchronous writes on the shared filesystem.
class BufferedOutput:public std::ofstream {
 std::unique_ptr<char[]> buffer;
 public:
 explicit BufferedOutput(const std::filesystem::path&p,std::ios::openmode mode=std::ios::out):buffer(new char[1024*1024]){rdbuf()->pubsetbuf(buffer.get(),1024*1024);open(p,mode|std::ios::out);require(bool(*this),"Cannot open output "+p.string());}
 ~BufferedOutput(){if(is_open())close();}
};

// Bounded ordered reduction. Records are private, same-process scratch data, not a cache format.
// Flush sorted aggregates, then merge at most 16 streams at a time. No event-sized vector.
class WorkDirectory {
 public:
 std::filesystem::path path;
 explicit WorkDirectory(const std::filesystem::path& parent) {
  std::string pattern=(parent/".work-XXXXXX").string();
  std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);
  require(mkdtemp(name.data())!=nullptr,"Cannot create analysis scratch directory");path=name.data();
 }
 ~WorkDirectory(){std::error_code e;std::filesystem::remove_all(path,e);}
};
template<class K,class V> class DiskReduce {
 static_assert(std::is_trivially_copyable<K>::value&&std::is_trivially_copyable<V>::value,"Scratch records must be POD");
 struct Record {K key;V value;};
 std::filesystem::path dir;std::map<K,V> cache;size_t limit,serial=0;
 std::vector<std::filesystem::path> runs;std::map<std::filesystem::path,size_t> runCounts;
 struct BinaryFile {
  gzFile f=nullptr;
  BinaryFile()=default;BinaryFile(const BinaryFile&)=delete;BinaryFile&operator=(const BinaryFile&)=delete;
  ~BinaryFile(){if(f)gzclose(f);}
  void open(const std::filesystem::path&p,bool writing=false){f=gzopen(p.c_str(),writing?"wb1":"rb");require(f!=nullptr,"Cannot open scratch run "+p.string());require(gzbuffer(f,256*1024)==0,"Cannot buffer scratch run");}
  void close(){int result=gzclose(f);f=nullptr;require(result==Z_OK,"Scratch close failed (check disk space)");}
 };
 static bool read(BinaryFile&f,Record&r){int n=gzread(f.f,&r,sizeof r);if(n==0&&gzeof(f.f))return false;require(n==sizeof r,"Truncated/unreadable scratch record");return true;}
 static void write(BinaryFile&f,const Record&r){require(gzwrite(f.f,&r,sizeof r)==sizeof r,"Scratch write failed (check disk space)");}
 std::filesystem::path fresh(){return dir/(std::to_string(serial++)+".bin");}
 std::filesystem::path merge(const std::vector<std::filesystem::path>& paths){
  struct Reader{BinaryFile f;Record r;};std::vector<Reader> readers(paths.size());
  auto greater=[&](size_t a,size_t b){return readers[b].r.key<readers[a].r.key;};
  std::priority_queue<size_t,std::vector<size_t>,decltype(greater)> q(greater);
  for(size_t i=0;i<paths.size();++i){readers[i].f.open(paths[i]);if(read(readers[i].f,readers[i].r))q.push(i);}
  auto p=fresh();BinaryFile f;f.open(p,true);Record sum{};bool active=false;size_t count=0;
  while(!q.empty()){auto i=q.top();q.pop();auto&r=readers[i].r;if(active&&sum.key==r.key)sum.value.merge(r.value);else{if(active){write(f,sum);++count;}sum=r;active=true;}if(read(readers[i].f,r))q.push(i);}
  if(active){write(f,sum);++count;}
  f.close();runCounts[p]=count;
  for(auto&p:paths){std::filesystem::remove(p);runCounts.erase(p);}
  return p;
 }
 void flush(){if(cache.empty())return;auto p=fresh();BinaryFile f;f.open(p,true);for(auto&[k,v]:cache)write(f,Record{k,v});f.close();runCounts[p]=cache.size();cache.clear();runs.push_back(p);}
 public:
 DiskReduce(const std::filesystem::path&parent,const std::string&name,size_t bytes):dir(parent/name),limit(std::max<size_t>(1,bytes/(sizeof(Record)+128))){require(std::filesystem::create_directory(dir),"Scratch reducer already exists");}
 V& get(const K&k){auto it=cache.find(k);if(it!=cache.end())return it->second;if(cache.size()>=limit)flush();return cache[k];}
 void finish(){flush();while(runs.size()>1){std::vector<std::filesystem::path>next;for(size_t i=0;i<runs.size();i+=16){size_t n=std::min(runs.size(),i+16);if(n-i==1)next.push_back(runs[i]);else next.push_back(merge({runs.begin()+i,runs.begin()+n}));}runs.swap(next);}}
 size_t size(){finish();return runs.empty()?0:runCounts.at(runs[0]);}
 template<class F>void each(F fn){finish();if(runs.empty())return;BinaryFile f;f.open(runs[0]);Record r;while(read(f,r))fn(r.key,r.value);}
 class Cursor {
  BinaryFile f;
  public:explicit Cursor(const std::filesystem::path&p){if(!p.empty())f.open(p);}
  bool next(K&k,V&v){if(!f.f)return false;Record r;if(!read(f,r))return false;k=r.key;v=r.value;return true;}
 };
 Cursor cursor(){finish();return Cursor(runs.empty()?std::filesystem::path{}:runs[0]);}
 void clear(){cache.clear();for(auto&p:runs)std::filesystem::remove(p);runs.clear();runCounts.clear();}
 ~DiskReduce(){std::error_code e;std::filesystem::remove_all(dir,e);}
};
inline void mergeMoment(Moment&a,const Moment&b){a.n+=b.n;a.s+=b.s;a.q+=b.q;}
inline size_t memoryBudget(const char*value){long long mb=std::stoll(value);require(mb>=16&&mb<=8192,"Memory budget must be 16..8192 MiB");return size_t(mb)*1024*1024;}
}
