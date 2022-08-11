#include <lemon/list_graph.h>
#include <chrono>
#include "mygraphlib.h"
#include "myutils.h"
#include "solver.h"


int main(int argc, char *argv[]){
  char name[1000];
  int time_limit = 3600;
  Digraph g;
  ArcValueMap weight(g);
  DNodeStringMap vname(g);
  DNodePosMap   posx(g),posy(g);
  string filename,graphtitle;
  DNode sourcenode,targetnode;
  int NNodes = countNodes(g);
  int seed=1;

  srand48(seed);

//   /***********************************************************************/
//   /***************************Input Checking******************************/
  if (argc!=4) {cout<<endl<<
      "Program to obtain a shortest path from node s to a node t," << endl <<
      "Usage: "<< argv[0]<<"  <graph>  <source_node_name>  <target_node_name>"<<
      endl << endl;
    cout << "Example:" << endl <<
      "\t"<<argv[0]<<" "<<getpath(argv[0])+"../instances/t_100.dig 12 50"<<endl<<endl;
    exit(0);}

  filename = argv[1];
  
  if (!ReadDigraph(filename,g,vname,posx,posy,weight))   // Read the graph
    {cout<<"Error reading digraph file "<<argv[1]<<"."<<endl;exit(0);}
//   /**************************End Input Check***********************************/
//   /****************************************************************************/

  Digraph::ArcMap<GRBVar> x(g);
  Digraph::ArcMap<GRBVar> x_0(g);
  Digraph::NodeMap<GRBVar> r(g);
  Digraph::NodeMap<GRBVar> r_0(g);
  DNodeIntMap demand(g);
  DNodeIntMap demand_type(g);
  const vector<int> a = {100, 50};
  const DNode depot = g.nodeFromId(0);
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
  model.getEnv().set(GRB_IntParam_Seed, seed);
  model.set(GRB_StringAttr_ModelName, "HVRP Formulation Golden et al"); // name to the problem
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem
  vector<int> max_vehicles;
  vector<DNode> BestCircuit;
  max_vehicles.push_back(std::atoi(argv[2]));
  max_vehicles.push_back(std::atoi(argv[3]));

  int q = a[1];
  // Randomly assign demand values and demand types to the nodes
  for(DNodeIt node(g); node!=INVALID; ++node){
    // Define the demand of each node so that it never exceeds maximum demand
    q = q - rand() % q;
    demand_type[node] = rand() % 2;
    demand[node] = q;
  }

  demand[depot] = 0;

for (int vehicle_type = 0; vehicle_type < 2; vehicle_type++) {  
    for(DNodeIt e(g); e!=INVALID; ++e){
      sprintf(name,"r_%s",vname[e].c_str());
      r[e] = model.addVar(0, a[vehicle_type], 0, GRB_INTEGER, name);
    }

    srand(seed);
    for (ArcIt e(g); e!=INVALID; ++e) {
      //Name structure x_i_j_k
      // i-> origin node
      // j-> destination node
      // k-> vehicle type

      if(g.source(e) != g.target(e)){
        sprintf(name,"x_%s_%s_%i",vname[g.source(e)].c_str(),vname[g.target(e)].c_str(), vehicle_type);
        x[e] = model.addVar(0.0, 1.0, weight[e],GRB_BINARY,name);  
      }
    }
    model.update(); // run update to use model inserted variables

  // Constraint 3
  for (DNodeIt p(g); p != INVALID; ++p){
    if(p!=depot){
      GRBLinExpr expr3;
      for(InArcIt in(g, p); in!=INVALID; ++in){
        expr3 += x[in];
      }
      for(OutArcIt out(g, p); out!=INVALID; ++out){
        expr3 -= x[out];  
      }
      model.addConstr(expr3 == 0);
    }
  }

  // Constraint 4
  GRBLinExpr expr4;
  for(OutArcIt out(g, depot); out!=INVALID; ++out){
    expr4 += x[out];
  }
  model.addConstr(expr4 <= max_vehicles[vehicle_type]);
  model.addConstr(1 <= expr4);
  cout << max_vehicles[vehicle_type] <<endl;

  // Create an auxiliary copy of the VarMaps
  if(vehicle_type == 0){
    for(ArcIt arc(g); arc!=INVALID; ++arc){
      x_0[arc] = x[arc];
    }
    for(DNodeIt node(g); node!=INVALID; ++node){
      r_0[node] = r[node];
    }
  }
}

// Constraint 2
for(DNodeIt j(g); j!=INVALID; ++j) { 
  GRBLinExpr expr2;
  for (int vehicle_type = 0; vehicle_type < 2; vehicle_type++){
    for (InArcIt e(g,j); e!=INVALID; ++e){
      if(g.target(e) != depot){
        if(vehicle_type == 0){
          expr2 += x_0[e];
        }else{
          expr2 += x[e];
        }
      }
    }
  }
  if(expr2.size() > 0){
    model.addConstr(expr2 == 1 );  
  }
}

// Constraint 5
for(int vehicle_type = 0; vehicle_type < 2; vehicle_type++){
  for(ArcIt ij(g); ij!=INVALID; ++ij){
    GRBLinExpr expr7; 
    if(g.target(ij)!=depot){
      for(int k = 0; k < 2; k++){
        if(k == 0){
          expr7 += (demand[g.target(ij)] + a[vehicle_type]) * x_0[ij];
        }else{
          expr7 += (demand[g.target(ij)] + a[vehicle_type]) * x[ij];
        }
      }
      expr7 -= a[vehicle_type];
    }
    if(expr7.size() > 0){
      if(vehicle_type == 0){
        model.addConstr(expr7 <= (r_0[g.target(ij)] - r_0[g.source(ij)]));
      }else{
        model.addConstr(expr7 <= (r[g.target(ij)] - r[g.source(ij)]));
      }
    } 
  } 
}

// Constraint 6
for(int k = 0; k < 2; k++){
  for (DNodeIt j(g); j!=INVALID; ++j) { 
    GRBLinExpr expr8;
    for (int vehicle_type = 0; vehicle_type < 2; vehicle_type++){
      for (InArcIt e(g,j); e!=INVALID; ++e){
        if(g.target(e) != depot){
          if(vehicle_type == 0){
            expr8 += a[k] * x_0[e];
            cout << "Varname: " << x_0[e].get(GRB_StringAttr_VarName) << endl;
          }else{
            expr8 += a[k] * x[e];
            cout << "Varname: " << x_0[e].get(GRB_StringAttr_VarName) << endl;
          }
        }
      }
    }
    if(expr8.size() > 0){
      if(k == 0){
        model.addConstr(r_0[j] <= expr8);  
      }else{
        model.addConstr(r[j] <= expr8);  
      }
    }
  }
}


// Make it so only one type of vehicle can pass through certain nodes
for(DNodeIt i(g); i!=INVALID; ++i){
  if(demand_type[i] == 1 && i != depot){
    for(InArcIt in(g, i); in !=INVALID; ++in){
      if (g.source(in) != depot){
        model.addConstr(x[in] == 0);
      }
    }
    for(OutArcIt out(g, i); out != INVALID; ++out){
      if(g.target(out) != depot){
        model.addConstr(x[out] == 0);
      }
    }
  }
  
}

// Show the percentage distribution of nodes that can only be accessed by ont type of vehicle
cout << "Node types" <<endl;
cout << "Type 0" <<endl;
int count_0 = 0;
int count_1 = 0;
for(DNodeIt node(g); node!=INVALID; ++node){
  count_0++;
  if(demand_type[node] == 0){
    cout << "Node: " << vname[node].c_str() << " Type:" << demand_type[node] << endl;
  }
}
cout << "Type 1" << endl;
for(DNodeIt node(g); node!=INVALID; ++node){
  count_1++;
  if(demand_type[node] == 1){
    cout << "Node: " << vname[node].c_str() << " Type: " << demand_type[node] << endl;
  }
}
cout << "End node types" <<endl;

cout << "Type 0 percentage: " << (count_0 / (count_0+ count_1)) <<endl;
cout << "Type 1 percentage: " << (count_1 / (count_0+ count_1)) <<endl;

try{
  if (time_limit >= 0){ 
    model.getEnv().set(GRB_DoubleParam_TimeLimit,time_limit);
  }

  // Time the algorithm execution
  auto start_time = chrono::high_resolution_clock::now();
  model.optimize();
  auto end_time = chrono::high_resolution_clock::now();
  auto time = end_time - start_time;
  cout << "Running Time in microseconds: " << chrono::duration_cast<chrono::microseconds>(time).count() << endl;

  // If the model has obtained at least one solution, print its members
  if (model.get(GRB_IntAttr_SolCount) > 0){
    cout << "Begin Solution" << endl;
    for(ArcIt e(g); e!=INVALID; ++e){
      if (BinaryIsOne(x_0[e].get(GRB_DoubleAttr_X))){
        cout << " -" << x_0[e].get(GRB_StringAttr_VarName) << endl;
      }
    }
    for(ArcIt e(g); e!=INVALID; ++e){
      if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))){
        cout << " -" << x[e].get(GRB_StringAttr_VarName) << endl;
      }
    }
    cout << "End of solution" << endl;
  }
}
catch(...)
{

} 
}