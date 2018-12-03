//Función que escribe VCKM a partir de la semilla
void write_VCKM(double V[3][3], double seed[6])
{
  int i = 0;
  int j = 0;
  double mu[2], md[2];
  memset(mu, 0, 2*sizeof(mu[0]));
  memset(md, 0, 2*sizeof(mu[0]));
  double pu[2], pd[2], ph[2];
  memset(pu, 0, 2*sizeof(mu[0]));
  memset(pd, 0, 2*sizeof(mu[0]));
  memset(ph, 0, 2*sizeof(mu[0]));
  mu[0] = mqu/mqt;
  mu[1] = mqc/mqt;
  md[0] = mqd/mqb;
  md[1] = mqs/mqb; //Masas escaladas a m3
  double Ou[3][3],Od[3][3];
  memset(Ou, 0, 3*3*sizeof(Ou[0][0]));
  memset(Od, 0, 3*3*sizeof(Ou[0][0]));
  
  pu[0] = seed[0]; //Au
  pu[1] = seed[1]; //fu 
  pd[0] = seed[2]; //Ad
  pd[1] = seed[3]; //fd
  ph[0] = seed[4];//ph1
  ph[1] = seed[5];//ph2
  //Definción de la Matriz Ortogonal
  Ou[0][0] = O11u(pu,mu);
  Ou[0][1] = O12u(pu,mu);
  Ou[0][2] = O13u(pu,mu);
  Ou[1][0] = O21u(pu,mu);
  Ou[1][1] = O22u(pu,mu);
  Ou[1][2] = O23u(pu,mu);
  Ou[2][0] = O31u(pu,mu);
  Ou[2][1] = O32u(pu,mu);
  Ou[2][2] = O33u(pu,mu);
  Od[0][0] = O11d(pd,md);
  Od[0][1] = O12d(pd,md);
  Od[0][2] = O13d(pd,md);
  Od[1][0] = O21d(pd,md);
  Od[1][1] = O22d(pd,md);
  Od[1][2] = O23d(pd,md);
  Od[2][0] = O31d(pd,md);
  Od[2][1] = O32d(pd,md);
  Od[2][2] = O33d(pd,md);

  //Cálculo de la matriz CKM
  for (i = 0; i <=2; ++i)
    {
      for( j = 0; j <= 2; ++j)
	{
	  V[i][j] = sqrt(pow(Ou[0][i]*Od[0][j] + Ou[1][i]*Od[1][j] * cos(ph[0]) + Ou[2][i]*Od[2][j] * cos(ph[0] + ph[1]),2) + pow(Ou[0][i]*Od[0][j] + Ou[1][i]*Od[1][j] * sin(ph[0]) + Ou[2][i]*Od[2][j] * sin(ph[0] + ph[1]) ,2));
	}   
    }
}

double JarlskogTH(double V[3][3], double seed[6])
{
  double mu[2], md[2];
  double pu[2], pd[2], ph[2];
  mu[0] = mqu/mqt;
  mu[1] = mqc/mqt;
  md[0] = mqd/mqb;
  md[1] = mqs/mqb; //Masas escaladas a m3
  double Ou[3][3],Od[3][3];

  pu[0] = seed[0]; //Au
  pu[1] = seed[1]; //fu
  pd[0] = seed[2]; //Ad
  pd[1] = seed[3]; //fd
  ph[0] = seed[4];//ph1
  ph[1] = seed[5];//ph2
  //Definción de la Matriz Ortogonal
  Ou[0][0] = O11u(pu,mu);
  Ou[0][1] = O12u(pu,mu);
  Ou[0][2] = O13u(pu,mu);
  Ou[1][0] = O21u(pu,mu);
  Ou[1][1] = O22u(pu,mu);
  Ou[1][2] = O23u(pu,mu);
  Ou[2][0] = O31u(pu,mu);
  Ou[2][1] = O32u(pu,mu);
  Ou[2][2] = O33u(pu,mu);
  Od[0][0] = O11d(pd,md);
  Od[0][1] = O12d(pd,md);
  Od[0][2] = O13d(pd,md);
  Od[1][0] = O21d(pd,md);
  Od[1][1] = O22d(pd,md);
  Od[1][2] = O23d(pd,md);
  Od[2][0] = O31d(pd,md);
  Od[2][1] = O32d(pd,md);
  Od[2][2] = O33d(pd,md);

  double A22, A21, Jarls;
  A22 = atan((Ou[1][2] * Od[1][2] + Ou[2][2] * Od[2][2] * sin(ph[0]) + Ou[3][2] * Od[3][2] * sin(ph[0] + ph[2]))/(Ou[1][2] * Od[1][2] + Ou[2][2] * Od[2][2] * cos(ph[0]) + Ou[3][2] * Od[3][2] * cos(ph[0] + ph[2])));
  A21 = atan((Ou[1][2] * Od[1][1] + Ou[2][2] * Od[2][1] * sin(ph[0]) + Ou[3][2] * Od[3][1] * sin(ph[0] + ph[2]))/(Ou[1][2] * Od[1][1] + Ou[2][2] * Od[2][1] * cos(ph[0]) + Ou[3][2] * Od[3][1] * cos(ph[0] + ph[2])));
  Jarls = V[1][1] * V[2][2] * V[1][2] * V[2][1] * sin(A22- A21);
  
  return Jarls;
}

double chi_sq(double V[3][3], double Jarlskog)
{
  int j;
  double chi_squared = 0.0;
  double VckmExp[9][2] = {{Vud , ErrVud },
			  {Vus , ErrVus },
			  {Vub , ErrVub },
			  {Vcd , ErrVcd },
			  {Vcs , ErrVcs },
			  {Vcb , ErrVcb },
			  {Vtd , ErrVtd },
			  {Vts , ErrVts },
			  {Vtb , ErrVtb }};
  //Cálculo de Chi cuadrada
  for (j=0;j<3;++j)
    {
      chi_squared = chi_squared + pow(VckmExp[j][0] - V[0][j],2)/pow(VckmExp[j][1],2);
    }
  chi_squared = chi_squared + pow(JS - Jarlskog,2)/pow(ErrJS,2);
  return chi_squared;
}

