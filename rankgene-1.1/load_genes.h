// Copyright (c) 2002 T. M. Murali
//
// File: load_genes.h
// Author: T. M. Murali <murali@bu.edu>
// Created: 10/24/02
//
// Purpose: function declarations for load_genes.C
//

// $Id: load_genes.h,v 1.2 2002/12/05 18:42:16 murali Exp $
//


#ifndef _LOAD_GENES_H
#define _LOAD_GENES_H

#include <map>
#include <string>
#include <vector>

using namespace std;

// forward declarations. struct unidim is defined in oc1.h
struct unidim;


//void load_classes(char *class_data, int **classes, int &no_of_classes);
void load_classes(char *class_data, map< string, int > &classes, int &no_of_classes);
void load_genes(char *train_data, const map< string, int > &classes,
                vector< struct unidim * > &genes,
                vector< string > &gene_names, vector< string > &gene_ids,
                int &no_of_genes, int &no_of_points);



#endif // _LOAD_GENES_H 
