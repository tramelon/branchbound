B&B1 B&B2: templates for branch and bound

defined_matrix.h: head file for defining cplex matrix, not used

parinit.h: designed to initiate the parameters, not used, might be changed to struct

parinit.cpp: for above

callb.h: not used

callb.cpp: not used

definedmatrix.h, definedmatrix.cpp: used to define cplex matrix

calub.h + calub.cpp: used to declare and define function to solve the switched order 

B&B: template from Yu AN

branbond1.cpp: branch and bound for sub problem only, no robust uc master problem

branbond3.cpp: same to branbond1.cpp

branbond2.cpp: added something

ruc.cpp: robust uc with master problem and switched order subproblem solved by branch and bound

subduality.cpp
RUC separation(sub problem)


ruc_laz_dec.cpp

ruc_checkbound.cpp  
robust uc with all options( strong duality, ccg ans switched order ccg)


ruc_laz.cpp
master problem using lazy constraints