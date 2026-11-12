#pragma once
#include <poll.h>
#include <sys/wait.h>
#include <signal.h>
#include <cerrno>
// Process workers receive aggregate distributions only. They do not read raw
// input and never write ROOT/TSV files. The parent writes all results in key order.
class FitPool {
 struct Worker{pid_t pid=-1;int request=-1,response=-1;};
 struct Wire{double values[7];long long nfit;int code,cov;char status[64];};
 std::vector<Worker>workers;const ts::Config&config;
 static bool readAll(int fd,void*buffer,size_t bytes){auto*p=static_cast<char*>(buffer);size_t done=0;while(done<bytes){ssize_t n=::read(fd,p+done,bytes-done);if(n<0&&errno==EINTR)continue;if(!n&&done==0)return false;ts::require(n>0,"Truncated fit-worker response");done+=n;}return true;}
 static void writeAll(int fd,const void*buffer,size_t bytes){auto*p=static_cast<const char*>(buffer);size_t done=0;while(done<bytes){ssize_t n=::write(fd,p+done,bytes-done);if(n<0&&errno==EINTR)continue;ts::require(n>0,"Fit worker pipe closed");done+=n;}}
 [[noreturn]]static void child(int request,int response,const ts::Config&c){try{Histogram<540>b;while(readAll(request,&b,sizeof b)){TH1D h("worker_adc","",540,-100,2600);h.SetDirectory(nullptr);b.restore(h);auto r=fit(h,b,c);Wire w{};double v[7]={r.mpv,r.me,r.w,r.we,r.sigma,r.se,r.chi};std::copy_n(v,7,w.values);w.nfit=r.nfit;w.code=r.code;w.cov=r.cov;ts::require(r.status.size()<sizeof w.status,"Worker status too long");std::copy(r.status.begin(),r.status.end(),w.status);writeAll(response,&w,sizeof w);}::close(request);::close(response);_exit(0);}catch(const std::exception&e){std::cerr<<"Fit worker: "<<e.what()<<std::endl;_exit(1);}}
 public:explicit FitPool(const ts::Config&c,int count):config(c){ts::require(count>=1&&count<=32,"fit-workers must be 1..32");if(count==1)return;signal(SIGPIPE,SIG_IGN);std::cout.flush();std::cerr.flush();workers.reserve(count);for(int i=0;i<count;++i){int req[2],res[2];ts::require(pipe(req)==0&&pipe(res)==0,"Cannot create fit-worker pipes");pid_t pid=fork();ts::require(pid>=0,"Cannot fork fit worker");if(!pid){for(auto&w:workers){::close(w.request);::close(w.response);}::close(req[1]);::close(res[0]);child(req[0],res[1],c);}::close(req[0]);::close(res[1]);workers.push_back({pid,req[1],res[0]});}}
 struct Task{size_t position;const Histogram<540>*value;};
 void run(const std::vector<Task>&tasks,std::vector<Fitted>&results){if(workers.empty()){for(auto&t:tasks){TH1D h("worker_adc","",540,-100,2600);h.SetDirectory(nullptr);t.value->restore(h);results[t.position]=fit(h,*t.value,config);}return;}
  std::vector<long long>active(workers.size(),-1);std::vector<pollfd>fds(workers.size());size_t sent=0,completed=0;
  auto send=[&](size_t i){if(sent==tasks.size()){active[i]=-1;fds[i]={-1,0,0};return;}active[i]=sent;writeAll(workers[i].request,tasks[sent++].value,sizeof(Histogram<540>));fds[i]={workers[i].response,POLLIN,0};};
  for(size_t i=0;i<workers.size();++i)send(i);
  while(completed<tasks.size()){int ready=poll(fds.data(),fds.size(),-1);if(ready<0&&errno==EINTR)continue;ts::require(ready>0,"Fit-worker poll failed");for(size_t i=0;i<workers.size();++i){if(active[i]<0||!fds[i].revents)continue;Wire w;ts::require(readAll(workers[i].response,&w,sizeof w),"Fit worker exited before returning result");auto&r=results[tasks[active[i]].position];r.mpv=w.values[0];r.me=w.values[1];r.w=w.values[2];r.we=w.values[3];r.sigma=w.values[4];r.se=w.values[5];r.chi=w.values[6];r.nfit=w.nfit;r.code=w.code;r.cov=w.cov;ts::require(w.status[sizeof w.status-1]==0,"Invalid fit-worker status");r.status=w.status;++completed;send(i);}}
 }
 ~FitPool(){for(auto&w:workers){if(w.request>=0)::close(w.request);if(w.response>=0)::close(w.response);}for(auto&w:workers)if(w.pid>0){int status;while(waitpid(w.pid,&status,0)<0&&errno==EINTR){}}}
};
// Same early conditions as fit(), without invoking the minimizer.
inline bool requiresFit(TH1D&h,const Summary&b,const ts::Config&c){if(!b.threshold)return false;double low=std::max(5.,b.low),high=std::min(2495.,b.high);if(high-low<100)return false;long long n=0;for(int i=h.FindFixBin(low);i<=h.FindFixBin(high);++i)n+=std::llround(h.GetBinContent(i));return n>=c.integer("min-fit-hits",100);}
