// Copyright (c) 2002 T. M. Murali
//
// File: load_genes.C
// Author: T. M. Murali <murali@bu.edu>
// Created: 10/24/02
//
// Purpose: load class and genes information for RankGene.
//


// static char load_genes_id[] = "$Id: load_genes.C,v 1.9 2002/12/09 22:50:59 murali Exp $";

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <stdio.h>

extern "C"
{
  
#include "oc1.h"
  
}

#include "load_genes.h"

const int UNKNOWN_CLASS = -1;

// print vectors.
template < class item >
ostream& operator<<(ostream& ostr, const vector< item >& v)
{
  for (unsigned int i = 0; i < v.size(); i++)
    ostr << v[i] << " ";
  ostr << endl;
  return(ostr);
}

// does c equal any char in chars?
inline bool match_char(char c,
                       string chars)
//                       const vector< char >& chars)
{
  return(string::npos != chars.find(c));
}

// simple function to keep reading into name from ostr until one of
// the chars in limits is read. return the read char. i use this
// function mainly to read strings e.g., gene names with spaces in
// them in a tab-delimited file. the function *READS* the char its
// returns (the char is no longer in istr).
char read_till_char(istream &istr, string &name, string limits)
{
  char c;
  string local_name;
  while (true)
    {
      // peek at the next character. i am peeking because i want to
      // get a character only if it is whitespace. otherwise, i can
      // get a on-character word in a gene name (e.g. 'E' in
      // "Semaphorin E", add that to the name and then also add the
      // next word into the name, which is incorrect).
      c = istr.peek();
      if (match_char(c, limits))
        // i should stop reading now.
        {
          c = istr.get();
          return(c);
        }
      if (' ' == c)
        // sometimes a gene name ends in a space. what i should do
        // then is to read the next character and check if it is a tab
        // rather than reading the next string.
        {
          // add c to name, otherwise i lose spaces.
          c = istr.get();
          name += c;
          continue;
        }
      // at this stage, i am sure that the next character is not
      // white-space. i can add the next word to the name.
      istr >> local_name;
//      cout << "local name is " << local_name << endl;
      name += local_name;
    }
}

  
/* read class information from class_file. the file lists the classes
   that the samples belong to. each line contains a sample name and a
   class name. the function assigns class names to unique integers and
   returns a map from sample names to classes (the integers). */
void load_classes(char *class_data, map< string, int > &classes,
                  int &no_of_classes)
{
  ifstream classfile;
  string class_name, sample_name, rest_of_line;
  int no_of_classes_allocated= 0;
  vector< string > class_names;
  /* class indices start at 1 in rankgene. */
  int class_index = 1;
  int i;
  char c;
  
  map< string, int > class_map;
  classfile.open(class_data);
  if (!classfile)
    error("Load_Classes: Class file can not be opened");
  
  while (-1 != classfile.peek())
    {
      // reset.
      sample_name = class_name = rest_of_line = "";
      read_till_char(classfile, sample_name, "\t");
      // sometimes the class name is followed by a \n. in that case,
      // don't call the next read_till_char.
      c = read_till_char(classfile, class_name, "\t\n");
      // ignore rest of line.
      if ('\n' != c)
        read_till_char(classfile, rest_of_line, "\n");
      
      class_names.push_back(class_name);
      if (class_map.end() == class_map.find(class_name))
        {
          // i map a class name of "NA" to UNKNOWN_CLASS.
          if ("NA" == class_name)
            class_map[class_name] = UNKNOWN_CLASS;
          else
            {
              class_map[class_name] = class_index;
              class_index++;
              no_of_classes++;
            }
        }
      classes[sample_name] = class_map[class_name];
//       cout << "Sample " << sample_name << " has class " << class_name
//            << " (mapped to " << class_map[class_name] << ")" << endl;
    }
  classfile.close();
  cout << "Read class names for " << classes.size() << " samples." << endl;
}

// replaces indices with missing values with the average of the other
// indices. if all indices are missing, return false, otherwise,
// returns true.
bool fix_missing_values(int num_values, struct unidim *gene,
                        const vector< int > &missing_values)
{
  float sum = 0, mean;
  int i, current_missing;

//  cout << "Fixing missing values." << endl;
  
  if ((num_values - 1) == missing_values.size())
    // all values are missing. this gene has no hope.
    return(false);
  current_missing = 0;
  for (i = 1; i <= num_values; i++)
    {
      if (i == missing_values[current_missing])
        // this index is missing.
        current_missing++;
      else
        // i can add this value to sum.
        sum += gene[i].value;
    }
  // the number of values present is num_values - missing_values.size() - 1.
  mean = sum/(num_values - missing_values.size() - 1);
  for (i = 0; i < missing_values.size(); i++)
    // fill in missing vale.
    gene[missing_values[i]].value = mean;
  return(true);
}

// simple function to return the class name of a sample.
int get_class_name(string sample, const map< string, int > &classes)
{
  return((*classes.find(sample)).second);
}

// does the sample have a class?
bool has_class_name(string sample, const map< string, int > &classes)
{
  return(classes.end() != classes.find(sample));
}

// does the sample have a known class? answer is NO if the class is
// UNKNOWN_CLASS.
bool is_class_known(string sample, const map< string, int > &classes)
{
  return(UNKNOWN_CLASS != get_class_name(sample, classes));
}


// load_genes reads a file containing gene expression data. the first
// line contains the column names: first the name of the column with
// gene ids, then the name of the column with the gene names, and then
// the names of the samples. each succeeding line gives the gene id,
// the gene name and the expression values of the samples.
void load_genes(char *gene_data,// int *classes,
                const map< string, int > &classes,
                vector < struct unidim * > &genes,
                vector< string > &gene_names,
                vector< string > &gene_ids,
                int &no_of_genes, int &no_of_points)
{
  ifstream genefile;
  genefile.open(gene_data);
  if (!genefile)
    error("Load_Genes: Gene file can not be opened.");

  string gene_id_column_name, gene_name_column_name;
  vector< string > sample_names;
  string gene_id, gene_name, sample_name;
  string limits("\t\n");
  
  unsigned int no_of_samples;
  no_of_genes = 0;
  bool end_of_line = false;
  float value;
  char c;
  int num_values;
  struct unidim *gene = NULL;
  bool reading_first_line = true;

  // HACKHACKHACKHACK. store a list of possible names of gene id and
  // gene name column names so that i know which column contains
  // which. these strings will have to be constantly updated for newer
  // formats.
  string possible_gene_id_column_names("Gene Accession Number GID Image Id. ");
  string possible_gene_name_column_names("Gene Description NAME ");
  string first_name, second_name, missing_value;
  int gene_name_column, gene_id_column, num_missing, num_genes_with_missing;
  vector< int > missing_indices;
  bool fixed;
  
  // push in a trivial gene since the rest of the oc1 code assumes
  // that indices start at 1.
  genes.push_back(gene);
  gene_names.push_back("");
  gene_ids.push_back("");
  sample_names.push_back("");

  num_genes_with_missing = 0;
  
  while (-1 != genefile.peek())
    // peek() returns -1 when it doesn't see anything.
    {
      // reset gene_id and gene_name.
      gene_id = "";
      gene_name = "";
      first_name = "";
      second_name = "";
      // assume that the gene's missing values are fixed (unless set
      // to false by fix_missing_values().
      fixed = true;

      // read the first two strings in the column.
      c = read_till_char(genefile, first_name, limits);
      c = read_till_char(genefile, second_name, limits);      

      
      // reset information on missing values.
      missing_indices.clear();
      num_missing = 0;

      if (reading_first_line)
        {
          if ((string::npos != possible_gene_id_column_names.find(first_name)) &&
              (string::npos != possible_gene_name_column_names.find(second_name)))
            // the first column corresponds to the gene id and the
            // second column to the gene name.
            {
              gene_id_column = 1;
              gene_name_column = 2;
            }
          else if ((string::npos != possible_gene_id_column_names.find(second_name)) &&
                   (string::npos != possible_gene_name_column_names.find(first_name)))
            // the first column corresponds to the gene name and the
            // second column to the gene id.
            {
              gene_id_column = 2;
              gene_name_column = 1;
            }
          else
            {
              cerr << "Name of first column (in the first row) is \""
                   << first_name
                   << "\" and name of second column (in the first row) is \""
                   << second_name << "\"" << endl;
              error("Unknown column names for gene names and gene ids");
            }
        }
      else
        // assign the gene name and gene id.
        {
          // +1 because OC1 arrays start at 1.
          gene = (struct unidim *)malloc((1 + no_of_samples)*sizeof(struct unidim));
          if (1 == gene_id_column)
            {
              gene_id = first_name;
              gene_name = second_name;
            }
          else
            {
              gene_id = second_name;
              gene_name = first_name;
            }
//          cout << "Gene id = " << gene_id << " and gene name = " << gene_name << endl;
          gene_names.push_back(gene_name);
          gene_ids.push_back(gene_id);
        }
      // since OC1's array start at index 1, start at 1.
      num_values = 1;
      // now read the rest of the line.
      while (true)
        {
          if (reading_first_line)
            // rest of the line should be made up of strings. read the
            // next sample name.
            {
              c = read_till_char(genefile, sample_name, limits);
//                  cout << "Read sample name " << sample_name << endl;

              // check that sample_name has an associated class.
              if (!get_class_name(sample_name, classes))
                {
                  cerr << "Sample \"" << sample_name << "\" has no associated class." << endl;
                  error("All samples do not have class labels. Please remove samples without class labels from the input file or provide class labels for them in the class file.");
                }
              // remember all sample names. if i only remember samples
              // whose classes are known, then when i read a value for
              // that sample, i don't have a corresponding sample name
              // and i get a nice seg fault.
              sample_names.push_back(sample_name);
              // reset sample_name.
              sample_name = "";
            }
          else
            // i am not reading the first line. the rest of the line
            // should be made up of numbers and strings indicating
            // missing values. to check if the next string is a
            // missing value, i will peek at it and guess.
            {
              c = genefile.peek();
              if (('?' == c) || ('N' == c) || ('\t' == c) || ('\n' == c))
                {
                  if (('\t' != c) && ('\n' != c))
                    // read the missing value only if it not the empty
                    // string.
                    {
                      // i have to read the missing value.
                      missing_value = "";
                      c = read_till_char(genefile, missing_value, limits);
                    }
                  num_missing++;
                  missing_indices.push_back(num_values);
                }
              else
                genefile >> value;
              // check if i managed to read a good value.
              if (!genefile)
                cerr << "Did not read a valid value." << endl;
//              cout << gene[num_values].value << " ";
              if (is_class_known(sample_names[num_values], classes))
                // i know this sample's class name, so i can store it.
                {
                  gene[num_values].value = value;
                  gene[num_values].cat =
                    get_class_name(sample_names[num_values], classes);
                }
            }
          num_values++;
          // i have read a sample name or value. check if the line has
          // ended by examining the next character.
          if (!reading_first_line)
            //  if i am reading the first line, i am getting the value
            //  of c from the previous call to read_till_char.
            c = genefile.get();
          if ('\t' == c)
            // check if the next character is a newline. if it is, get
            // it and break.
            {
              char oldc = c;
              c = genefile.get();
              if ('\n' == c)
                {
                  // there are two cases here. the last \t might just
                  // be an extra one tagged on at the end of the line.
                  // in that case, i just break. otherwise, the empty
                  // string between the \t and \n might be a missing
                  // value. in that case, i increment do the whole
                  // missing value thing.
                  if (!reading_first_line &&
                      ((num_values - 1) != no_of_samples))
                    {
                      missing_indices.push_back(num_values);
                      num_missing++;
                      gene[num_values].cat =
                        (*classes.find(sample_names[num_values])).second;
                      num_values++;
                    }
                  break;
                }
              
              // next character is a regular character.
              genefile.putback(c);
              c = oldc;
              continue;
            }
          
          if ('\n' == c)
            break;
          // c seems to be a regular ascii character. put it back in
          // the stream.
          genefile.putback(c);
        } // finished reading a line.
      if (reading_first_line)
        {
          reading_first_line = false;
          // -1 because num_values started at 1.
          no_of_samples  = num_values - 1;
//          cout << "Sample names are " << sample_names << endl;
          cout << "Read first line. There are " << no_of_samples << " samples." << endl;
        }
      else
        {
//          if ((num_values + num_missing - 1) != no_of_samples) 
          if ((num_values - 1) != no_of_samples)
           {
              cerr << "Read a gene with " << num_values - 1
                   << " sample values." << endl;
              error("Gene with wrong number of sample values.");
            }
          if (0 != missing_indices.size())
            {
              // check.
              if (num_missing != missing_indices.size())
                cerr << "#missing values = " << num_missing << " but there are" << missing_indices.size() << " missing indices." << endl;
              fixed = fix_missing_values(num_values, gene, missing_indices);
              num_genes_with_missing++;
            }
          if (fixed)
            genes.push_back(gene);
          else
            // all the values for this gene are missing. forget about it. 
            {
              gene_id = gene_ids.back();
              gene_ids.pop_back();
              gene_name = gene_names.back();
              gene_names.pop_back();
              cerr << "Discarding gene with id \"" << gene_id << "\" and name \""
                   << gene_name << "\"" << endl;
            }
          // print progress info.
          if (0 == genes.size()%1000)
            cout << "Read " << genes.size() << " genes.\n";
        }
    }
  // -1 because of the first push_back() of an empty gene.
  no_of_genes = genes.size() - 1;
  // -1 because of the first push_back() of an empty sample name.
  no_of_points = sample_names.size() - 1;
  cout << "Read " << no_of_genes << " genes and " << no_of_points << " samples." << endl;
  cout << num_genes_with_missing << " genes had some missing values." << endl;
  genefile.close();
}




