// Phase 1: sparse circle kNN layers and incremental bridge contraction.
// Rebuild this portable C++17 helper; no binary is included in the release.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;

struct BridgeForest {
    vector<int> ecc, cc, size, parent, last;
    int components, bridges=0, iteration=0;
    explicit BridgeForest(int n):ecc(n),cc(n),size(n,1),parent(n,-1),last(n,0),components(n) {
        iota(ecc.begin(),ecc.end(),0); iota(cc.begin(),cc.end(),0);
    }
    int find_ecc(int v) {
        if(v<0)return -1;
        return ecc[v]==v?v:ecc[v]=find_ecc(ecc[v]);
    }
    int find_cc(int v) {
        v=find_ecc(v);
        return cc[v]==v?v:cc[v]=find_cc(cc[v]);
    }
    void root_at(int v) {
        v=find_ecc(v); int root=v, child=-1;
        while(v!=-1) {
            int p=find_ecc(parent[v]);
            parent[v]=child; cc[v]=root; child=v; v=p;
        }
        size[root]=size[child];
    }
    void close_cycle(int a,int b) {
        ++iteration;vector<int> pa,pb;int lca=-1;
        while(lca==-1) {
            if(a!=-1) {
                a=find_ecc(a);pa.push_back(a);
                if(last[a]==iteration){lca=a;break;}
                last[a]=iteration;a=parent[a];
            }
            if(b!=-1) {
                b=find_ecc(b);pb.push_back(b);
                if(last[b]==iteration){lca=b;break;}
                last[b]=iteration;b=parent[b];
            }
        }
        for(auto path:{pa,pb})for(int v:path) {
            ecc[v]=lca;if(v==lca)break;--bridges;
        }
    }
    void add(int a,int b) {
        a=find_ecc(a);b=find_ecc(b);if(a==b)return;
        int ca=find_cc(a),cb=find_cc(b);
        if(ca!=cb) {
            ++bridges;--components;
            if(size[ca]>size[cb]){swap(a,b);swap(ca,cb);}
            root_at(a);parent[a]=b;cc[a]=b;size[cb]+=size[a];
        } else close_cycle(a,b);
    }
};

struct Directed { uint64_t pair;int rank; };
struct Edge { int a,b,key_or,key_mu; };
template<class T> T read_one(istream& f) {
    T x;f.read(reinterpret_cast<char*>(&x),sizeof(x));
    if(!f)throw runtime_error("truncated input");
    return x;
}

int main(int argc,char**argv) {
    try {
        if(argc<2)throw runtime_error("need input file");
        ifstream input(argv[argc-1],ios::binary);
        int n=read_one<int32_t>(input),cap=read_one<int32_t>(input);
        vector<Edge> edges;
        auto begin=chrono::steady_clock::now();
        if(argc==3 && string(argv[1])=="--generic") {
            int m=read_one<int32_t>(input);
            for(int i=0;i<m;++i) {
                int a=read_one<int32_t>(input),b=read_one<int32_t>(input);
                int x=read_one<int32_t>(input),y=read_one<int32_t>(input);
                edges.push_back({a,b,x,y});
            }
        } else {
            if(n<3||cap<1||cap>=n)throw runtime_error("invalid circle dimensions");
            vector<double>s(n);input.read(reinterpret_cast<char*>(s.data()),n*sizeof(double));
            if(!input||!is_sorted(s.begin(),s.end()))throw runtime_error("invalid sorted sample");
            vector<Directed> directed;directed.reserve(size_t(n)*cap);
            for(int i=0;i<n;++i) {
                int left=i-1,right=i+1;
                for(int rank=1;rank<=cap;++rank) {
                    int jl=(left%n+n)%n,jr=right%n;
                    double dl=s[i]-(s[jl]+(left<0?-1:0));
                    double dr=(s[jr]+(right>=n?1:0))-s[i];
                    int j;
                    if(dl<dr||(dl==dr&&jl<jr)){j=jl;--left;}
                    else {j=jr;++right;}
                    int a=min(i,j),b=max(i,j);
                    directed.push_back({uint64_t(a)*n+b,rank});
                }
            }
            sort(directed.begin(),directed.end(),[](auto&a,auto&b){return a.pair<b.pair;});
            edges.reserve(directed.size());
            for(size_t i=0;i<directed.size();) {
                size_t j=i+1;
                while(j<directed.size()&&directed[j].pair==directed[i].pair)++j;
                if(j-i>2)throw runtime_error("duplicate neighbor");
                int x=directed[i].rank,y=j-i==2?directed[i+1].rank:cap+1;
                edges.push_back({int(directed[i].pair/n),int(directed[i].pair%n),min(x,y),max(x,y)});
                i=j;
            }
        }
        auto built=chrono::steady_clock::now();
        for(int method=0;method<2;++method) {
            vector<vector<pair<int,int>>> layers(cap+1);
            int close=0;
            vector<bool> ring(n,false);
            for(auto e:edges) {
                int k=method==0?e.key_or:e.key_mu;
                if(k>cap)continue;
                layers[k].push_back({e.a,e.b});
                if(e.b==e.a+1||(e.a==0&&e.b==n-1)) {
                    ring[e.b==e.a+1?e.a:n-1]=true;close=max(close,k);
                }
            }
            if(!all_of(ring.begin(),ring.end(),[](bool x){return x;}))close=-1;
            BridgeForest forest(n);
            int conn=-1,br=-1,at_conn=-1;size_t inserted=0,at_br=0;
            for(int k=1;k<=cap;++k) {
                for(auto e:layers[k]){forest.add(e.first,e.second);++inserted;}
                // Test only after the complete tied integer-k layer.
                if(conn<0&&forest.components==1){conn=k;at_conn=forest.bridges;}
                if(br<0&&forest.components==1&&forest.bridges==0){br=k;at_br=inserted;break;}
            }
            if(method)cout<<" ";
            cout<<conn<<" "<<br<<" "<<close<<" "<<at_conn<<" "<<at_br;
        }
        auto end=chrono::steady_clock::now();
        cout<<" "<<edges.size()<<" "
            <<chrono::duration<double>(built-begin).count()<<" "
            <<chrono::duration<double>(end-built).count()<<"\n";
    } catch(exception&e){cerr<<e.what()<<"\n";return 1;}
}
