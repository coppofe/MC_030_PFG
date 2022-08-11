#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "mygraphlib.h"
using namespace std;
int main(int argc, char *argv[]) 
{
  int n;
  srand48(clock());
  if (argc!=2) {cout<<"Usage: "<< argv[0]<<" <number_of_nodes>"<<endl; exit(0);} 
  n = atoi(argv[1]);
  Graph g;
  NodeStringMap vname(g);
  NodePosMap px(g);
  NodePosMap py(g);
  EdgeValueMap weight(g);
  GenerateRandomEuclideanGraph(g,vname,px,py,weight,n,1000,1000);
  cout << "nnodes narcs type" << endl;
  cout << n << " " << n*(n-1)/2 << " digraph" << endl;
  cout << "nodename" << endl;
  for (NodeIt v(g);v!=INVALID;++v) 
    cout<<vname[v]<<endl;
  cout << "tail head weight" << endl;
  for (EdgeIt e(g);e!=INVALID;++e){
    cout<<vname[g.u(e)]<<" "<<vname[g.v(e)]<<" "<<int(weight[e])<<endl;
    cout<<vname[g.v(e)]<<" "<<vname[g.u(e)]<<" "<<int(weight[e])<<endl;
  }
}
