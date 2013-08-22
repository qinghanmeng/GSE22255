/****************************************************************/
/* Copyright 2002                                               */
/* Bioinformatics Program			                */
/* Boston University		                                */
/****************************************************************/
/* Contact : yangsu@bu.edu					*/
/****************************************************************/
/* File Name : rankgene.c                                       */
/* Author : Yang Su					        */
/* Last modified : Aug 2002					*/
/* The main module in Rankgene program                          */
/****************************************************************/

#include <unistd.h>

extern "C"
{
  
#include "oc1.h"
  
}

#include "load_genes.h"

char *pname;
char train_data[LINESIZE],test_data[LINESIZE], class_data[LINESIZE];
char test_out[LINESIZE],gene_rank[LINESIZE];

int no_of_dimensions=0,no_of_genes=0,no_of_categories=0;
int *left_count=NULL, *right_count=NULL;
int split_option=1;
int no_of_train_points=0,no_of_test_points=0;
int no_of_missing_values=0,unlabeled=FALSE;;
int no_of_wanted_genes=100;
int sv1_class,sv2_class;
float onedsv1,onedsv2;
float svm_c=1;
float t_test;
float *t_test_array;

vector < struct unidim * > genes;
vector< string > gene_names;
vector< string > gene_ids;
map< string, int > classes;

unsigned int use_oc1_format = 0;

float compute_impurity();



void srand48();

struct unidim *candidates;
struct unidim *gene_list;
struct support_vector *svectors;
struct test_outcome estimate_accuracy();

POINT **train_points=NULL,**test_points=NULL;


/************************************************************************/
/************************************************************************/

void read_data(char *input_file, int no_of_points)
{
  
 FILE *infile;
 int i,j,k,count;
 //POINT **points,**allocate_point_array();
 POINT **points,**allocate_point_array(POINT **array_name, int size, int prev_size); 
 if (strlen(input_file) == 0 )
   error("Read_Data : No data filename specified.");
 if (no_of_points < -1)
   error("Read_Data : Invalid number of points to be loaded.");
 if ((infile = fopen(input_file,"r")) == NULL)
   error("Read_Data : Data file can not be opened.");
 
 count = load_points(infile,&points);
 if (no_of_points != -1)  shuffle_points (points,count);
 fclose(infile);
 
 if (no_of_points > count)
  error("Read_Data : Insufficient data in input file.");
 
 if (no_of_points == 0 || no_of_points == count)
  {
    no_of_test_points = 0;
    no_of_train_points = count;
  }
 else if (no_of_points == -1)
  {
    no_of_test_points = count;
    no_of_train_points = 0;
  }
 else
  {
   no_of_test_points = count - no_of_points;
   no_of_train_points = no_of_points;
  }
 
 if (no_of_train_points)
  {
   train_points = allocate_point_array(train_points,no_of_train_points,0);
   for (i=1;i<=no_of_train_points;i++)
    {
     for (j=1;j<=no_of_dimensions;j++)
      train_points[i]->dimension[j] = points[i]->dimension[j];
     train_points[i]->category = points[i]->category;
     train_points[i]->val = points[i]->val;
    }
  }
 
 if (no_of_test_points)
  {
   test_points = allocate_point_array(test_points,no_of_test_points,0);
   for (i=no_of_train_points+1;i<=count;i++)
    {
     k = i - no_of_train_points;
 
     for (j=1;j<=no_of_dimensions;j++)
      test_points[k]->dimension[j] = points[i]->dimension[j];
     test_points[k]->category = points[i]->category;
     test_points[k]->val = points[i]->val;
    }
  }
 
 for (i=1;i<=count;i++)
  {
     free_vector(points[i]->dimension,1,no_of_dimensions);
     free((char *)points[i]);
  }
 free((char *)(points+1));
 
}

void allocate_structures(int no_of_points)
{
  int i;
  
  no_of_genes = no_of_dimensions;
  left_count = allocate_ivector(1,no_of_categories);
  right_count = allocate_ivector(1,no_of_categories);
  t_test_array = allocate_vector(1,no_of_genes);
  if (use_oc1_format)
    {
      candidates = (struct unidim *)malloc((unsigned)no_of_points *
                                           sizeof(struct unidim));
      candidates -= 1;
    }
  gene_list = (struct unidim *)malloc((unsigned)no_of_genes * sizeof(struct unidim));
  gene_list -= 1;
  svectors = (struct support_vector *)malloc((unsigned)no_of_genes * sizeof(struct support_vector));
  svectors -= 1;
}


void deallocate_structures(int no_of_points)
{
  free_ivector(left_count,1,no_of_categories);
  free_ivector(right_count,1,no_of_categories);
  free_vector(t_test_array,1,no_of_genes);
  if (use_oc1_format)
    free((candidates+1));
  else
    // free the genes. start at 1 because genes[0] is NULL (see load_genes()).
//    for (unsigned int i = 1; i <= no_of_genes; i++)
    for (unsigned int i = no_of_genes; i > 0; i--)
      free(genes[i]);
  free((gene_list+1));
  free((svectors+1));
}

/************************************************************************/
/************************************************************************/

void axis_parallel_split(POINT **cur_points, int cur_no_of_points)
{
  int i,j,cur_gene;
  float linear_split(int);
  FILE *testfile;
  
   
  for (cur_gene=1;cur_gene<=no_of_genes;cur_gene++)
    {
      if (use_oc1_format)
        {
          for (j=1;j<=cur_no_of_points;j++) 
            { 
              candidates[j].value = cur_points[j]->dimension[cur_gene]; 
              candidates[j].cat   = cur_points[j]->category; 
          
            }
        }
      else
        candidates = genes[cur_gene];
      gene_list[cur_gene].value = (float)linear_split(cur_no_of_points);      
      gene_list[cur_gene].cat = cur_gene;
      t_test_array[cur_gene] = t_test;
      svectors[cur_gene].vector1 = onedsv1;
      svectors[cur_gene].vector2 = onedsv2;
      svectors[cur_gene].class1 = sv1_class;
      svectors[cur_gene].class2 = sv2_class;
     }
   
}


//int compare(struct unidim *ptr1, struct unidim *ptr2) 
int compare(const void *ptr1, const void *ptr2) 
{
  float x;
  
  x = (*(struct unidim *)ptr1).value - (*(struct unidim *)ptr2).value;
  
  if (x > 0) return(1); 
  else if (x) return(-1);
  else return(0);
}

void reset_counts()
{
  int i;
  
  for (i=1;i<=no_of_categories;i++)
    left_count[i] = right_count[i] = 0;
}


float linear_split(int no_of_eff_points) 
{
  int i,j,from,to,count1,count2,position1,position2;
  int correct = 0;
  int position = 0;
  float temp,impurity_1d;
  float smallest,largest;
  float mean1,mean2,var1,var2,temp1,temp2,t;
  float *array1,*array2;
  struct oned_svm *candidates1,*candidates2;
  
  float compute_impurity(int);
  float compute_lp(int position, struct oned_svm *class1, struct oned_svm *class2,
                   int c1, int c2);
  float unnormalize(float value, float small, float large);
  int check_sv(int position, struct oned_svm *class1, struct oned_svm *class2,
               int c1, int c2);

  if (split_option <= 6) {
   candidates += 1;
   qsort(candidates,no_of_eff_points,sizeof(struct unidim),compare);
   candidates -= 1;
  
   reset_counts();
   for (i=1;i<=no_of_eff_points;i++)
     right_count[candidates[i].cat]++;
  
   impurity_1d = compute_impurity(no_of_eff_points);
  

   for (i=1;i<=no_of_eff_points;i++)
    {
      from = i;
      for (to=from+1;to<=no_of_eff_points && candidates[to].value ==
	   candidates[from].value;to++);
      to -= 1;
      
      for (j=from;j<=to;j++)
	{
	  left_count[candidates[j].cat]++;
	  right_count[candidates[j].cat]--;
	}

      i = to;
      temp = compute_impurity(no_of_eff_points);
      
      if (temp < impurity_1d ||
	  (temp == impurity_1d && myrandom(0.0,1.0) < 0.5))
	{
	  impurity_1d = temp;
          onedsv1 = candidates[from].value;
          onedsv2 = candidates[to+1].value;           
	  if (impurity_1d == 0) break;
	}
    } 
  }
  
  else if (split_option == 7) {
   candidates += 1;
   qsort(candidates,no_of_eff_points,sizeof(struct unidim),compare);
   candidates -= 1;
   
   count1 = 0;
   count2 = 0;
   for(i=1;i<=no_of_eff_points;i++) {
       if(candidates[i].cat == 1)
          count1++;
       else
          count2++;
   }
   array1 = allocate_vector(1,count1);
   array2 = allocate_vector(1,count2);
   
   from = 1;
   to = 1;
   for(i=1;i<=no_of_eff_points;i++) {
       if(candidates[i].cat == 1)  {
       	   array1[from] = candidates[i].value;
       	   from++;
       }
       else {
           array2[to] = candidates[i].value;
       	   to++;	  
       }
   }
   
   mean1 = average(array1,count1);
   mean2 = average(array2,count2);
   temp = 0;
   for(i=1;i<=count1;i++) 
      temp += (array1[i]-mean1)*(array1[i]-mean1);
   var1 = temp/((count1-1)*count1);
   temp = 0;
   for(i=1;i<=count2;i++)
      temp += (array2[i]-mean2)*(array2[i]-mean2);
   var2 = temp/((count2-1)*count2);
   temp = sqrt(var1+var2);
   t = (mean1 - mean2)/temp;
   t_test = t;
   if(t < 0)
      t = -1*t;
   impurity_1d = 1/t;
   
   free_vector(array1,1,count1);
   free_vector(array2,1,count2);
  }  
  
  
  else {
   candidates += 1;
   qsort(candidates,no_of_eff_points,sizeof(struct unidim),compare);
   candidates -= 1;
   smallest = candidates[1].value;
   largest = candidates[no_of_eff_points].value;  
 
   count1 = 0;
   count2 = 0;
   for(i=1;i<=no_of_eff_points;i++) {
       if(candidates[i].cat == 1)
          count1++;
       else
          count2++;
   }
   candidates1 = (struct oned_svm *)malloc(count1*sizeof(struct oned_svm));
   candidates2 = (struct oned_svm *)malloc(count2*sizeof(struct oned_svm));
   
   from = 0;
   to = 0;
   for(i=0;i<no_of_eff_points;i++) {
       if(candidates[i+1].cat == 1) {
       	  candidates1[from].cat = 1;
       	  if(candidates[i+1].value == smallest) 
       	     candidates1[from].value = 0;
       	  else if(candidates[i+1].value == largest)
       	     candidates1[from].value = 1;
       	  else {
       	     temp = candidates[i+1].value;
       	     candidates1[from].value = (temp - smallest)/(largest - smallest);
       	  }
       	  from++;
       }
       else {
       	  candidates2[to].cat = -1;
       	  if(candidates[i+1].value == smallest) 
       	     candidates2[to].value = 0;
       	  else if(candidates[i+1].value == largest)
       	     candidates2[to].value = 1;
       	  else {
       	     temp = candidates[i+1].value;
       	     candidates2[to].value = (temp - smallest)/(largest - smallest);
       	  }
       	  to++;
       }
   }
   for(i=0;i<count1;i++) {
       if(i == 0) {
       	  candidates1[i].left_total = 0;
       	  candidates1[count1-1].right_total = 0;
       }
       else {
       	  candidates1[i].left_total = candidates1[i-1].left_total + candidates1[i-1].value;
       	  candidates1[count1-i-1].right_total = candidates1[count1-i].right_total + candidates1[count1-i].value;
       }
   }
   for(i=0;i<count2;i++) {
       if(i == 0) {
       	  candidates2[i].left_total = 0;
       	  candidates2[count2-1].right_total = 0;
       }
       else {
       	  candidates2[i].left_total = candidates2[i-1].left_total + candidates2[i-1].value;
       	  candidates2[count2-i-1].right_total = candidates2[count2-i].right_total + candidates2[count2-i].value;
       }
   }
   
   
   if(count1 <= count2) {
      impurity_1d = compute_lp(0,candidates1,candidates2,count1,count2);  
      for(i=1;i<count1;i++) {
          temp = compute_lp(i,candidates1,candidates2,count1,count2);
          if(impurity_1d == 0) {
             impurity_1d = temp;
             position = i;
          }
          else {
	    if(temp > 0 && temp < impurity_1d) { 
                 impurity_1d = temp;
	         position = i;
            }
          }
      }
      position1 = position;
      onedsv1 = unnormalize(candidates1[position1].value,smallest,largest);
      sv1_class = 1;
      position2 = check_sv(position1,candidates1,candidates2,count1,count2);
      onedsv2 = unnormalize(candidates2[position2].value,smallest,largest);
      sv2_class = -1;
   }
   else {
      impurity_1d = compute_lp(0,candidates2,candidates1,count2,count1);
      for(i=1;i<count2;i++) {
          temp = compute_lp(i,candidates2,candidates1,count2,count1);
          if(impurity_1d == 0) {
	     impurity_1d = temp;
             position = i;
          }
          else {
	    if(temp > 0 && temp < impurity_1d) {
	         impurity_1d = temp;
                 position = i;
            }
	  }
      }
      position2 = position;
      onedsv2 = unnormalize(candidates2[position2].value,smallest,largest);
      sv2_class = -1;
      position1 = check_sv(position2,candidates2,candidates1,count2,count1);
      onedsv1 = unnormalize(candidates1[position1].value,smallest,largest);
      sv1_class = 1;
   }
   free(candidates1);
   free(candidates2);
  }
   	
  return(impurity_1d);
}


float unnormalize(float value, float small, float large)
{
      return(small + value*(large - small));
}

int check_sv(int position, struct oned_svm *class1, struct oned_svm *class2,
             int c1, int c2)
{
  int i;
  float w,temp;
  float temp1 = 0;
  float temp2 = 0;
  
  i = position;
  if(class1[i].value < class2[c1-i-1].value) {
      w = 2/(class1[i].value - class2[c1-i-1].value);
      temp1 = w*w/2 + svm_c*(2*(c1-i-1) - w*(class1[i].right_total - class2[c1-i-1].left_total));
  }
  if(class1[i].value > class2[c2-i-1].value) {
      w = 2/(class1[i].value - class2[c2-i-1].value);
      temp2 = w*w/2 + svm_c*(2*i - w*(class1[i].left_total - class2[c2-i-1].right_total));
  }
  
  if(temp1 == 0)
     return(c2-i-1);
  else if(temp2 == 0)
     return(c1-i-1);
  else {
     if(temp1 < temp2)
        return(c1-i-1);
     else
        return(c2-i-1);
  }
}

float compute_lp(int position, struct oned_svm *class1, struct oned_svm *class2,
                 int c1, int c2)
{
  int i;
  float w,temp;
  float temp1 = 0;
  float temp2 = 0;
  
  i = position;
  if(class1[i].value < class2[c1-i-1].value) {
      w = 2/(class1[i].value - class2[c1-i-1].value);
      temp1 = w*w/2 + svm_c*(2*(c1-i-1) - w*(class1[i].right_total - class2[c1-i-1].left_total));
  }
  if(class1[i].value > class2[c2-i-1].value) {
      w = 2/(class1[i].value - class2[c2-i-1].value);
      temp2 = w*w/2 + svm_c*(2*i - w*(class1[i].left_total - class2[c2-i-1].right_total));
  }
  
  if(temp1 == 0) {
     if(temp2 > 0) return(temp2);
     else return(0);
  }
  else if(temp2 == 0) {
     if(temp1 > 0) return(temp1);
     else return(0);
  }
  else {
     if(temp1 < temp2) return(temp1);
     else return(temp2);
  }
}

int largest_element(int *array, int count)
{
  int i,major;
  
  major = 1;
  for (i=2;i<=count;i++)
    if (array[i] > array[major]) major = i;
  return(major);
}


float compute_impurity(int cur_no_of_points)
{
  float temp;
   
  switch(split_option)
    {
       case 1:
         temp = info_gain();
         break;
       case 2:
         temp = twoing();
         break;
       case 3:
         temp = summinority();
         break;
       case 4:
         temp = maxminority();
         break;
       case 5:
         temp = gini_index();
         break;
       case 6:
         temp = variance();
         break;
      }
  return(temp);
}
  

/************************************************************************/
/************************************************************************/



void print_list()
{
  int i=0;
  int rank = 1;
  int j;
  FILE *outfile;
  

  if ((outfile = fopen(gene_rank,"w")) == NULL)
    {
      fprintf(stderr,"rankgene: Output file cannot be written to.\n");
      exit(0);
    }
  
  if (strlen(train_data)) fprintf(outfile,"Gene expression data: %s\n",train_data);
  fprintf(outfile,"Gene ranking measures: 1-information gain; 2-twoing rule; 3-summinority; 4-maxminority; 5-gini index; 6-sum of variances; 7-t_test; 8-1D support vector machine.\n");
  if (split_option) fprintf(outfile,"Using gene ranking measure: %d\n",split_option);
  fprintf(outfile,"Data has %d genes with %d classes.\n",
	  no_of_dimensions,no_of_categories);
  
  gene_list += 1;
  qsort(gene_list,no_of_genes,sizeof(struct unidim),compare);
  gene_list -= 1;
  if (use_oc1_format)
    fprintf(outfile, "Rank\tIndex\tMeasure value\n");
  else
    fprintf(outfile, "Rank\tIndex\tGene Id\tGene Name\tMeasure value\n");
  for(i=1;i <= no_of_wanted_genes;i++)
    {
       if(i > 1) {
         if(gene_list[i].value != gene_list[i-1].value)
	    rank++;
       }

      if (use_oc1_format)
        // don't print gene name and id since i don't know what they
        // are in this format.
        fprintf(outfile,"%i\t%i\t%f", rank, gene_list[i].cat,
               gene_list[i].value);
      else
       fprintf(outfile,"%i\t%i\t%s\t%s\t",rank,gene_list[i].cat,
               gene_ids[gene_list[i].cat].c_str(),
               gene_names[gene_list[i].cat].c_str());
       
       j = gene_list[i].cat;                
       if(split_option < 3)
         fprintf(outfile,"%f\t%f\t%f\n",-gene_list[i].value,svectors[j].vector1,
                   svectors[j].vector2);
       else if(split_option < 7)
         fprintf(outfile,"%f\t%f\t%f\n",gene_list[i].value,svectors[j].vector1,
                   svectors[j].vector2);
       else if(split_option == 7) {
         fprintf(outfile,"%f\n",t_test_array[j]);
       }
       else {
         fprintf(outfile,"%f\t%i: %f\t%i: %f\n",gene_list[i].value,svectors[j].class1,svectors[j].vector1,svectors[j].class2,svectors[j].vector2);
       }
    }

  fclose(outfile);
}



/************************************************************************/
/************************************************************************/ 
void exit_with_help()
{
    printf( "\nUsage: rankgene -m8 -n<number of genes> -i<input file> -c<class file> -o<output file> -w<weight parameter> \n"
            "\t\tOR\n"
            "       rankgene -m<1 to 7> -n<number of genes> -i<input file> -c<class file> -o<output file>\n"
            "Options:\n"
            "\tc: class filename.\n"
            "\ti: input filename.\n"
            "\tm: information measures selection: 1-6 OC1 measures; 7 t test; 8 1D SVM;\n"
            "\tn: number of genes that will be listed in the output file;\n"
            "\to: output filename;\n"
            "\tR: input file is in old RankGene (oc1) format;\n"
            "\tw: weight parameter for SVM, can be any positive number;\n"
            "Please read README for more details.\n\n"
           );
    exit(1);
 }

void check_parameter()
{
  if(svm_c <= 0) {
    printf("option -w error: weight parameter cannot be negative.\n");
    exit_with_help();
  }
  if(split_option < 1 || split_option > 8) {
    printf("option -m error: please specify options from 1 to 8.\n");
    exit_with_help();
  }
  if(no_of_wanted_genes < 1) {
    printf("option -n error: number of listed genes should be at least 1.\n");
    exit_with_help();
  }
  if(!strlen(train_data)) {
    printf("option -i error: no input data file name available.\n");
    exit_with_help();
  }
  if ((!strlen(class_data)) && (!use_oc1_format)) {
    printf("option -c error: no class file name available.\n");
    exit_with_help();
  }
  if(!strlen(gene_rank)) {
    printf("option -o error: please specify output file name.\n");
    exit_with_help();
  }
}

int main (int argc, char *argv[])
{
  extern char *optarg;
  extern int optind;
  int c1,i,j;
  
  strcpy(train_data,"\0");
  strcpy(gene_rank,"\0");
  
  pname = argv[0];
  if (argc == 1) 
     exit_with_help();
  while ((c1 = 
	  getopt (argc, argv, "c:w:m:n:o:i:R")) 
         != EOF)
  
  switch (c1)
     {
      case 'c':                      /*Class file name. */
        strcpy(class_data,optarg);
	break;
      case 'i':                      /*Input file name. */
        strcpy(train_data,optarg);
	break;
      case 'm':                      /*i=1-6 oc1 impurity measure,7 t test,8 1-d algorithm*/
        split_option = atoi(optarg);
        break;
      case 'n':
        no_of_wanted_genes = atoi(optarg);  /*The number of listed genes*/
        break;
      case 'o':
        strcpy(gene_rank,optarg);    /*Output file name*/
        break;
       case 'R':
         use_oc1_format = 1;         /* expect files in old RankGene (oc1) format. */
         break;
       case 'w':                      /*w is user-defined cost parameter*/
        svm_c = atof(optarg);
        break;
      default:
        exit_with_help();
     }

  check_parameter();

  if (!use_oc1_format)
    {
      
      load_classes(class_data, classes, no_of_categories);
      
      load_genes(train_data, classes, genes, gene_names, gene_ids,
                 no_of_dimensions, no_of_train_points);
    }

  else
    read_data(train_data,0);
  // correct value of no_of_wanted_genes if the data has fewer genes
  // than the number wanted.
  if (no_of_wanted_genes > no_of_dimensions)
    no_of_wanted_genes = no_of_dimensions;
  // commenting this code out now since we are checking if
  // no_of_wanted_genes is >= 1 in check_parameter(). will figure out
  // a way later to ask for all the genes to ranked.
  //
//   else if (-1 == no_of_wanted_genes)
//     // if the command line specifies -1, rank all the genes.
//     no_of_wanted_genes = no_of_dimensions;
    
  if (split_option >=7 && no_of_categories > 2) {
    printf("Error: Measures 7 (t-test) and 8 (one-dimensional support vector machines) only accept data sets with two classes.\n");
    exit(1);
  }
  
  allocate_structures(no_of_train_points);
  
  axis_parallel_split(train_points,no_of_train_points);
  
  print_list();
  
  deallocate_structures(no_of_train_points);
   
}
