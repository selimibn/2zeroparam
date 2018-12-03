//#include <gsl/gsl_multimin.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "parameters.h" //Parametros Experimentales
#include "parametrizations/par88.h" //Definición de funciones  
#include "chi2alt.h" //Definición del chi cuadrado

#define mu1 mqu/mqt
#define mu2 mqc/mqt
#define md1 mqd/mqb
#define md2 mqs/mqb
#define chimin 6E6
#define uparam 8
#define dparam 8


int main(){

  int grid_size=5;

  FILE *seed_file;
  seed_file = fopen("dat88.dat", "w");

  /* int grid_size = 25; */
  /* int uparam = 8; */
  /* int dparam = 1; */
  int iter;
  int points;
  double pu[2], pd[2], ph[2];
  points = pow(grid_size, 6);
  double eps = 0;
  double chi2;


  double uparrange[8][2][2] = {
    { {mu1,mu2},  {mu2,1} },		/* parametrization 1 */
    { {-mu1,mu2}, {mu2,1}},  		/* parametrization 2 */
    { {-mu2,mu1}, {mu1,1}},  		/* parametrization 3 */
    { {-mu2,-mu1}, {-mu1,1}},  		/* parametrization 4 */
    { {-1,mu1}, {mu1,mu2}},  		/* parametrization 5 */
    { {-1,-mu1}, {-mu1,mu2}},  		/* parametrization 6 */
    { {-1,-mu2}, {-mu2,mu1}},  		/* parametrization 7 */
    { {-1,-mu2}, {-mu2,-mu1}}  		/* parametrization 8 */
  };

  double dparrange[8][2][2] = {
    { {md1,md2},  {md2,1} },		/* parametrization 1 */
    { {-md1,md2}, {md2,1}},  		/* parametrization 2 */
    { {-md2,md1}, {md1,1}},  		/* parametrization 3 */
    { {-md2,-md1}, {-md1,1}},  		/* parametrization 4 */
    { {-1,md1}, {md1,md2}},  		/* parametrization 5 */
    { {-1,-md1}, {-md1,md2}},  		/* parametrization 6 */
    { {-1,-md2}, {-md2,md1}},  		/* parametrization 7 */
    { {-1,-md2}, {-md2,-md1}}  		/* parametrization 8 */
  };

  double Aumin = uparrange[uparam-1][1][0];
  double Aumax = uparrange[uparam-1][1][1];
  double Admin = dparrange[dparam-1][1][0];
  double Admax = dparrange[dparam-1][1][1];
  double fumin = uparrange[uparam-1][0][0];
  double fumax = uparrange[uparam-1][0][1];
  double fdmin = dparrange[dparam-1][0][0];
  double fdmax = dparrange[dparam-1][0][1];
  double ph1min = 0;
  double ph1max = 2*Pi;
  double ph2min = 0;
  double ph2max = 2*Pi;

  double VCKM[3][3]; //Arreglo para contener los valores de CKM para una semilla en Particular

#pragma omp parallel for
  for(int i = 1; i < points ; ++i)
    {
      int grid_index[6];
      double texture_parameter[6];
      grid_index[0] = i%grid_size ;
      grid_index[1] = (i/grid_size)%grid_size ;
      grid_index[2] = (i/(grid_size*grid_size))%grid_size ;
      grid_index[3] = (i/(grid_size*grid_size*grid_size))%grid_size ;
      grid_index[4] = (i/(grid_size*grid_size*grid_size*grid_size))%grid_size ;
      grid_index[5] = (i/(grid_size*grid_size*grid_size*grid_size*grid_size))%grid_size ;

      texture_parameter[2] = Admin*(1+eps) + (1-eps)*(grid_index[0]+1) * (Admax-Admin)/grid_size; //Interval for Ad
      texture_parameter[0] = Aumin*(1+eps) + (1-eps)*(grid_index[1]+1) * (Aumax-Aumin)/grid_size; //Interval for Au
      texture_parameter[3] = fdmin*(1+eps) + (1-eps)*(grid_index[2]+1) * (fdmax-fdmin)/grid_size; //Interval for fd
      texture_parameter[1] = fumin*(1+eps) + (1-eps)*(grid_index[3]+1) * (fumax-fumin)/grid_size; //Interval for fu
      texture_parameter[4] = ph1min*(1+eps) + (1-eps)*(grid_index[4]+1) * (ph1max-ph1min)/grid_size; //Interval for ph1
      texture_parameter[5] = ph2min*(1+eps) + (1-eps)*(grid_index[5] +1)* (ph2max-ph2min)/grid_size; //Interval for ph2
      /* printf("Linear index %i is grid position (%i,%i,%i,%i,%i,%i) with parameter values (%lf,%lf, %lf, %lf, %lf, %lf)\n", i, grid_index[0], grid_index[1], grid_index[2], grid_index[3], grid_index[4], grid_index[5], texture_parameter[0], texture_parameter[1], texture_parameter[2], texture_parameter[3], texture_parameter[4], texture_parameter[5]); */
#pragma omp critical(dataupdate)
      {         write_VCKM(VCKM, texture_parameter);
	chi2 = chi_sq(VCKM,JarlskogTH(VCKM,texture_parameter));
	if(chi2<chimin){
	  fprintf(seed_file,"%.12lf %.12lf  %.12lf  %.12lf  %.12lf  %.12lf  %.12lf\n", texture_parameter[0],  texture_parameter[1],  texture_parameter[2] , texture_parameter[3],  texture_parameter[4], texture_parameter[5], chi2);
	  printf("%.12lf %.12lf  %.12lf %.12lf %.12lf %.12lf %.12lf \n", texture_parameter[0],  texture_parameter[1],  texture_parameter[2] , texture_parameter[3],  texture_parameter[4],  texture_parameter[5], chi2);
	}
      }
    }

  printf("orden de los parámetros: Au, fu,  Ad, fd, phi1, ph2\n");
  printf("El parámetro Au va de %.12lf a %.12lf\n", Aumin, Aumax);
  printf("El parámetro fu va de %.12lf a %.12lf\n", fumin, fumax);
  printf("El parámetro Ad va de %.12lf a %.12lf\n", Admin, Admax);
  printf("El parámetro fd va de %.12lf a %.12lf\n", fdmin, fdmax);
  return 0 ;
  fclose(seed_file);
}       
